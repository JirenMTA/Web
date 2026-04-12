
function [assocIdx, unassocDet, predCache] = associateDetections(tracks, detections, peakList, dt, Q, gate, rgpoCfg)
% Association with tracker-side RGPO support.
%
% Input:
%   tracks    : tracker struct
%   detections: Nx2 [range, velocity]
%   peakList  : Nx6 [ir, iv, pow, th, marginDb, clusterSize]
%   dt, Q     : motion model params
%   gate      : struct with sigR, sigV, threshold
%   rgpoCfg   : struct, example:
%       rgpoCfg.enable = true;
%       rgpoCfg.victimTrackId = 1;           % track to be pulled off
%       rgpoCfg.velLane = 1.5;               % m/s, same Doppler lane
%       rgpoCfg.minPullOff = 2.0;            % m, farther than predicted by at least this much
%       rgpoCfg.maxBackStep = 1.0;           % do not accept much nearer detections
%       rgpoCfg.minVictimMarginDb = 0;       % optional extra filter
%       rgpoCfg.confirmedOnly = true;
%
% Output:
%   assocIdx(k) = detection index for track k, 0 = miss
%   unassocDet  = detections still unused
%   predCache(k).state, predCache(k).P = predicted state/cov for update step

    nTracks = numel(tracks.list);
    nDets   = size(detections, 1);

    assocIdx = zeros(1, nTracks);
    usedDet  = false(1, nDets);
    predCache = repmat(struct('state', [], 'P', []), 1, max(nTracks,1));

    if nTracks == 0 || nDets == 0
        unassocDet = 1:nDets;
        return;
    end

    % ------------------------------------------------------------
    % 1) Predict ALL tracks first
    % ------------------------------------------------------------
    for k = 1:nTracks
        if strcmp(tracks.list(k).status, 'deleted')
            predCache(k).state = tracks.list(k).state;
            predCache(k).P     = tracks.list(k).P;
            continue;
        end

        [sp, Pp] = kalmanPredict(tracks.list(k).state, tracks.list(k).P, dt, Q);
        predCache(k).state = sp;
        predCache(k).P     = Pp;
    end

    % ------------------------------------------------------------
    % 2) Optional RGPO association for victim track ONLY
    % ------------------------------------------------------------
    if nargin >= 7 && isfield(rgpoCfg, 'enable') && rgpoCfg.enable
        for k = 1:nTracks
            if strcmp(tracks.list(k).status, 'deleted')
                continue;
            end

            if ~isVictimTrack(tracks.list(k), rgpoCfg)
                continue;
            end

            if isfield(rgpoCfg, 'confirmedOnly') && rgpoCfg.confirmedOnly && ...
               ~strcmp(tracks.list(k).status, 'confirmed')
                continue;
            end

            jVictim = chooseRgpoDetection(predCache(k).state, predCache(k).P, ...
                                          detections, peakList, usedDet, gate, rgpoCfg);
            if jVictim > 0
                assocIdx(k) = jVictim;
                usedDet(jVictim) = true;
            end
        end
    end

    % ------------------------------------------------------------
    % 3) Global greedy NN on REMAINING pairs
    % Better than per-track greedy because it is order-independent enough
    % for small scenes, even though it is not full Hungarian assignment.
    % ------------------------------------------------------------
    pairList = [];
    for k = 1:nTracks
        if strcmp(tracks.list(k).status, 'deleted') || assocIdx(k) > 0
            continue;
        end

        sp = predCache(k).state;
        Pp = predCache(k).P;

        for j = 1:nDets
            if usedDet(j)
                continue;
            end

            [dist, insideGate] = gatedMahalanobis(sp, Pp, detections(j,:)', gate);
            if insideGate
                pairList(end+1,:) = [k, j, dist]; %#ok<AGROW>
            end
        end
    end

    if ~isempty(pairList)
        [~, order] = sort(pairList(:,3), 'ascend');
        pairList = pairList(order,:);

        assignedTrack = false(1, nTracks);
        for ii = 1:size(pairList,1)
            k = pairList(ii,1);
            j = pairList(ii,2);

            if assignedTrack(k) || usedDet(j) || assocIdx(k) > 0
                continue;
            end

            assocIdx(k) = j;
            usedDet(j) = true;
            assignedTrack(k) = true;
        end
    end

    unassocDet = find(~usedDet);
end


function [state_pred, P_pred] = kalmanPredict(state, P, dt, Q)
    F = [1, dt; 0, 1];
    state_pred = F * state;
    P_pred     = F * P * F' + Q;
end


function tf = isVictimTrack(track, rgpoCfg)
    tf = false;

    if isfield(rgpoCfg, 'victimTrackId') && track.id == rgpoCfg.victimTrackId
        tf = true;
        return;
    end

    if isfield(track, 'isVictim') && track.isVictim
        tf = true;
    end
end

function jBest = chooseRgpoDetection(statePred, Ppred, detections, peakList, usedDet, gate, rgpoCfg)
% For a victim track, prefer the farther-range detection in the same Doppler lane
% if it is still inside gate and not too weak.

    jBest = 0;
    cand = [];

    rPred = statePred(1);
    vPred = statePred(2);

    if ~isfield(rgpoCfg, 'velLane'), rgpoCfg.velLane = 1.5; end
    if ~isfield(rgpoCfg, 'minPullOff'), rgpoCfg.minPullOff = 2.0; end
    if ~isfield(rgpoCfg, 'maxBackStep'), rgpoCfg.maxBackStep = 1.0; end
    if ~isfield(rgpoCfg, 'minVictimMarginDb'), rgpoCfg.minVictimMarginDb = -inf; end

    for j = 1:size(detections,1)
        if usedDet(j)
            continue;
        end

        meas = detections(j,:)';
        [dist, insideGate] = gatedMahalanobis(statePred, Ppred, meas, gate);
        if ~insideGate
            continue;
        end

        dv = abs(meas(2) - vPred);
        if dv > rgpoCfg.velLane
            continue;
        end

        dr = meas(1) - rPred;
        if dr < -rgpoCfg.maxBackStep
            continue;
        end

        marginDb = -inf;
        if ~isempty(peakList) && size(peakList,1) >= j
            marginDb = peakList(j,5);
        end
        if marginDb < rgpoCfg.minVictimMarginDb
            continue;
        end

        cand(end+1,:) = [j, dr, dist, marginDb, meas(1), meas(2)]; %#ok<AGROW>
    end

    if isempty(cand)
        return;
    end

    % First preference: farther detections that represent pull-off
    farMask = cand(:,2) >= rgpoCfg.minPullOff;
    if any(farMask)
        farCand = cand(farMask,:);
        % choose farthest range; if tie, smaller distance wins
        [~, idx] = sortrows([ -farCand(:,5), farCand(:,3) ], [1 2]);
        jBest = farCand(idx(1),1);
        return;
    end

    % Fallback: normal nearest detection inside same Doppler lane
    [~, idx] = min(cand(:,3));
    jBest = cand(idx,1);
end

function [dist, insideGate] = gatedMahalanobis(statePred, Ppred, meas, gate)
    innov = meas - statePred;
    S = Ppred + diag([gate.sigR^2, gate.sigV^2]);
    dist = innov' / S * innov;
    insideGate = (dist < gate.threshold);
end
