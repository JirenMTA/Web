function tracks = updateTracks(tracks, detections, assocIdx, ...
                               unassocDet, dt, Q, R, tmCfg)
    %% Update existing tracks
    for k = 1:numel(tracks.list)
        if strcmp(tracks.list(k).status, 'deleted'), continue; end

        [sp, Pp] = kalmanPredict(tracks.list(k).state, ...
                                 tracks.list(k).P, dt, Q);

        if assocIdx(k) > 0
            meas = detections(assocIdx(k), 1:2)';
            [su, Pu, ~] = kalmanUpdate(sp, Pp, meas, R);

            tracks.list(k).state  = su;
            tracks.list(k).P      = Pu;
            tracks.list(k).hits   = tracks.list(k).hits + 1;   % ✓ sửa bug
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

    %% Khởi tạo track mới
    for j = unassocDet
        newTrack.id      = tracks.nextID;
        tracks.nextID    = tracks.nextID + 1;
        newTrack.state   = detections(j, 1:2)';
        newTrack.P       = tmCfg.initCov;
        newTrack.hits    = 1;
        newTrack.misses  = 0;
        newTrack.status  = 'tentative';
        newTrack.history = detections(j, 1:2);

        if numel(tracks.list) == 0
            tracks.list = newTrack;                          % list rỗng
        else
            newTrack = orderfields(newTrack, tracks.list(1)); % khớp field order
            tracks.list(end+1) = newTrack;
        end

        fprintf('[TM] Track %d INITIATED r=%.1f m v=%.1f m/s\n', ...
            newTrack.id, newTrack.state(1), newTrack.state(2));
    end
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