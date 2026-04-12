function tracks = updateTracks(tracks, detections, peakList, assocIdx, ...
                               unassocDet, predCache, R, tmCfg)
% Update existing tracks + filtered track initiation.
%
% New tmCfg fields you should set:
%   tmCfg.minInitMarginDb = 15;
%   tmCfg.minInitRangeSep = 4.0;   % m
%   tmCfg.minInitVelSep   = 1.5;   % m/s

    if nargin < 8
        error('updateTracks now expects predCache before R and tmCfg.');
    end

    % ------------------------------------------------------------
    % 1) Update / miss existing tracks using PRECOMPUTED predictions
    % ------------------------------------------------------------
    for k = 1:numel(tracks.list)
        if strcmp(tracks.list(k).status, 'deleted')
            continue;
        end

        sp = predCache(k).state;
        Pp = predCache(k).P;

        if assocIdx(k) > 0
            meas = detections(assocIdx(k), 1:2)';
            [su, Pu, ~] = kalmanUpdate(sp, Pp, meas, R);

            tracks.list(k).state  = su;
            tracks.list(k).P      = Pu;
            tracks.list(k).hits   = tracks.list(k).hits + 1;
            tracks.list(k).misses = 0;

            if strcmp(tracks.list(k).status, 'tentative') && ...
               tracks.list(k).hits >= tmCfg.confirmHits
                tracks.list(k).status = 'confirmed';
                fprintf('[TM] Track %d CONFIRMED\n', tracks.list(k).id);
            end
        else
            tracks.list(k).state  = sp;
            tracks.list(k).P      = Pp;
            tracks.list(k).misses = tracks.list(k).misses + 1;

            if tracks.list(k).misses >= tmCfg.maxMisses
                tracks.list(k).status = 'deleted';
                fprintf('[TM] Track %d DELETED (miss=%d)\n', ...
                    tracks.list(k).id, tracks.list(k).misses);
            end
        end

        if ~strcmp(tracks.list(k).status, 'deleted')
            tracks.list(k).history(end+1, :) = tracks.list(k).state';
        end
    end

    % ------------------------------------------------------------
    % 2) Initiate new tracks more carefully
    % Do not start a new track from every leftover detection.
    % ------------------------------------------------------------
    if ~isfield(tmCfg, 'minInitMarginDb'), tmCfg.minInitMarginDb = 15; end
    if ~isfield(tmCfg, 'minInitRangeSep'), tmCfg.minInitRangeSep = 4.0; end
    if ~isfield(tmCfg, 'minInitVelSep'),   tmCfg.minInitVelSep   = 1.5; end

    for j = unassocDet
        if ~isempty(peakList)
            marginDb = peakList(j,5);
            if marginDb < tmCfg.minInitMarginDb
                continue;
            end
        end

        meas = detections(j,1:2);
        if isTooCloseToActiveTrack(tracks, meas, tmCfg.minInitRangeSep, tmCfg.minInitVelSep)
            continue;
        end

        newTrack.id      = tracks.nextID;
        tracks.nextID    = tracks.nextID + 1;
        newTrack.state   = meas(:);
        newTrack.P       = tmCfg.initCov;
        newTrack.hits    = 1;
        newTrack.misses  = 0;
        newTrack.status  = 'tentative';
        newTrack.history = meas;

        % Optional victim flag support
        newTrack.isVictim = false;

        if numel(tracks.list) == 0
            tracks.list = newTrack;
        else
            newTrack = orderfields(newTrack, tracks.list(1));
            tracks.list(end+1) = newTrack;
        end

        fprintf('[TM] Track %d INITIATED r=%.1f m v=%.1f m/s\n', ...
            newTrack.id, newTrack.state(1), newTrack.state(2));
    end
end


function [state_upd, P_upd, innov] = kalmanUpdate(state_pred, P_pred, meas, R)
    H = eye(2);
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    innov      = meas - H * state_pred;
    state_upd  = state_pred + K * innov;
    P_upd      = (eye(2) - K * H) * P_pred;
end


function tf = isTooCloseToActiveTrack(tracks, meas, minRangeSep, minVelSep)
    tf = false;
    for k = 1:numel(tracks.list)
        if strcmp(tracks.list(k).status, 'deleted')
            continue;
        end

        st = tracks.list(k).state(:);
        if abs(st(1) - meas(1)) <= minRangeSep && abs(st(2) - meas(2)) <= minVelSep
            tf = true;
            return;
        end
    end
end
