function tracks = initTrackManager()
    tracks.list = struct( ...
        'id',           {}, ...   % Track ID
        'state',        {}, ...   % [r; v] — range & velocity
        'P',            {}, ...   % Covariance matrix 2x2
        'hits',         {}, ...   % số lần detect liên tiếp
        'misses',       {}, ...   % số lần miss liên tiếp
        'status',       {}, ...   % 'tentative' | 'confirmed' | 'deleted'
        'history',      {}  ...   % lịch sử state để plot
    );
    tracks(1).nextID = 1;
end

function [state_pred, P_pred] = kalmanPredict(state, P, dt, Q)
% Predict bước tiếp theo
% State: [range; velocity]
% Model: range += velocity * dt (constant velocity)

    F = [1, dt; ...   % transition matrix
         0,  1];

    state_pred = F * state;
    P_pred     = F * P * F' + Q;
end


function [state_upd, P_upd, innov] = kalmanUpdate(state_pred, P_pred, meas, R)
% Update với measurement mới
% meas = [r_meas; v_meas]

    H = eye(2);                          % observation matrix (đo trực tiếp r,v)
    S = H * P_pred * H' + R;            % innovation covariance
    K = P_pred * H' / S;                % Kalman gain

    innov      = meas - H * state_pred; % innovation
    state_upd  = state_pred + K * innov;
    P_upd      = (eye(2) - K * H) * P_pred;
end

function [assocIdx, unassocDet] = associateDetections(tracks, detections, gate)
% GNN — Global Nearest Neighbor association
% gate: ngưỡng khoảng cách Mahalanobis để chấp nhận cặp (track, detect)
%
% Output:
%   assocIdx(k) = chỉ số detection gán cho track k (0 = miss)
%   unassocDet  = detection chưa được gán (candidates for new track)

    nTracks = numel(tracks);
    nDets   = size(detections, 1);   % [r, v] mỗi hàng

    assocIdx  = zeros(1, nTracks);
    usedDet   = false(1, nDets);

    for k = 1:nTracks
        if strcmp(tracks(k).status, 'deleted'), continue; end

        bestDist = inf;
        bestJ    = 0;

        for j = 1:nDets
            if usedDet(j), continue; end

            meas  = detections(j, :)';
            innov = meas - tracks(k).state;

            % Mahalanobis distance
            S    = tracks(k).P + diag([gate.sigR^2, gate.sigV^2]);
            dist = innov' / S * innov;

            if dist < gate.threshold && dist < bestDist
                bestDist = dist;
                bestJ    = j;
            end
        end

        if bestJ > 0
            assocIdx(k) = bestJ;
            usedDet(bestJ) = true;
        end
    end

    unassocDet = find(~usedDet);
end

function tracks = updateTracks(tracks, detections, assocIdx, ...
                               unassocDet, dt, Q, R, tmCfg)
% tmCfg: track management config
%   .confirmHits   — số hit để confirm (vd: 3)
%   .maxMisses     — số miss để xóa   (vd: 3)
%   .initCov       — P khởi tạo

    %% --- Update các track đã associated ---
    for k = 1:numel(tracks)
        if strcmp(tracks(k).status, 'deleted'), continue; end

        % Predict
        [sp, Pp] = kalmanPredict(tracks(k).state, tracks(k).P, dt, Q);

        if assocIdx(k) > 0
            % Hit: update với measurement
            meas = detections(assocIdx(k), :)';
            [su, Pu, ~] = kalmanUpdate(sp, Pp, meas, R);

            tracks(k).state  = su;
            tracks(k).P      = Pu;
            tracks(k).hits   = tracks(k).hits + 1;
            tracks(k).misses = 0;

            % Promote tentative → confirmed
            if strcmp(tracks(k).status, 'tentative') && ...
               tracks(k).hits >= tmCfg.confirmHits
                tracks(k).status = 'confirmed';
                fprintf('[TM] Track %d CONFIRMED\n', tracks(k).id);
            end
        else
            % Miss: chỉ predict, không update
            tracks(k).state  = sp;
            tracks(k).P      = Pp;
            tracks(k).misses = tracks(k).misses + 1;

            % Xóa nếu miss quá nhiều
            if tracks(k).misses >= tmCfg.maxMisses
                tracks(k).status = 'deleted';
                fprintf('[TM] Track %d DELETED (miss=%d)\n', ...
                    tracks(k).id, tracks(k).misses);
            end
        end

        % Lưu lịch sử
        if ~strcmp(tracks(k).status, 'deleted')
            tracks(k).history(end+1, :) = tracks(k).state';
        end
    end

    %% --- Khởi tạo track mới từ unassociated detections ---
    persistent nextID;
    if isempty(nextID), nextID = 1; end

    for j = unassocDet
        newTrack.id      = nextID;   nextID = nextID + 1;
        newTrack.state   = detections(j, :)';
        newTrack.P       = tmCfg.initCov;
        newTrack.hits    = 1;
        newTrack.misses  = 0;
        newTrack.status  = 'tentative';
        newTrack.history = detections(j, :);

        tracks(end+1) = newTrack;
        fprintf('[TM] Track %d INITIATED at r=%.1f m, v=%.1f m/s\n', ...
            newTrack.id, newTrack.state(1), newTrack.state(2));
    end
end

