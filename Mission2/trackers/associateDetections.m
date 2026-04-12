function [assocIdx, unassocDet] = associateDetections(tracks, detections, gate)
% GNN — Global Nearest Neighbor association
% gate: ngưỡng khoảng cách Mahalanobis để chấp nhận cặp (track, detect)
%
% Output:
%   assocIdx(k) = chỉ số detection gán cho track k (0 = miss)
%   unassocDet  = detection chưa được gán (candidates for new track)

    nTracks = numel(tracks.list);
    nDets   = size(detections, 1);   % [r, v] mỗi hàng

    assocIdx  = zeros(1, nTracks);
    usedDet   = false(1, nDets);

    for k = 1:nTracks
        if strcmp(tracks.list(k).status, 'deleted'), continue; end

        bestDist = inf;
        bestJ    = 0;

        for j = 1:nDets
            if usedDet(j), continue; end

            meas  = detections(j, :)';
            innov = meas - tracks.list(k).state;

            % Mahalanobis distance
            S    = tracks.list(k).P + diag([gate.sigR^2, gate.sigV^2]);
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