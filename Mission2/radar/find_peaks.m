function [detCellsFinal, peakList, detections] = find_peaks(RDpow, r_axis, v_axis)
    numTrain = [10 4];
    numGuard = [6 2];

    Nr = size(RDpow, 1);
    Nd = size(RDpow, 2);

    rMargin = numTrain(1) + numGuard(1);
    dMargin = numTrain(2) + numGuard(2);

    if Nr <= 2*rMargin || Nd <= 2*dMargin
        error('RD map qua nho cho cau hinh CFAR hien tai.');
    end

    rValid = (rMargin + 1) : (Nr - rMargin);
    dValid = (dMargin + 1) : (Nd - dMargin);
    [rr, dd] = ndgrid(rValid, dValid);
    cutidx = [rr(:).'; dd(:).'];

    trainCellsTotal = (2*(numTrain(1)+numGuard(1))+1) * ...
                      (2*(numTrain(2)+numGuard(2))+1) - ...
                      (2*numGuard(1)+1) * (2*numGuard(2)+1);
    rankOS = round(0.75 * trainCellsTotal);

    cfar2d = phased.CFARDetector2D( ...
        'Method', 'OS', ...
        'Rank', rankOS, ...
        'GuardBandSize', numGuard, ...
        'TrainingBandSize', numTrain, ...
        'ProbabilityFalseAlarm', 1e-4, ...
        'ThresholdOutputPort', true);

    [detFlags, th] = cfar2d(RDpow, cutidx);
    detCells = cutidx(:, detFlags);

    if isempty(detCells)
        fprintf('No detection.\n');
        detCellsFinal = [];
        peakList = [];
        fprintf('So detections sau OS-CFAR + clustering + filter = 0\n');
        return;
    end

    linIdx = sub2ind(size(RDpow), detCells(1,:), detCells(2,:));
    detMap = false(size(RDpow));
    detMap(linIdx) = true;

    CC = bwconncomp(detMap, 8);
    peakList = [];
    thDet = th(detFlags);

    for ic = 1:CC.NumObjects
        pix = CC.PixelIdxList{ic};
        [~, kmax] = max(RDpow(pix));
        bestLin = pix(kmax);
        [ir, iv] = ind2sub(size(RDpow), bestLin);
        bestPow = RDpow(bestLin);

        match = find(detCells(1,:) == ir & detCells(2,:) == iv, 1, 'first');
        if isempty(match)
            bestTh = min(thDet);
        else
            bestTh = thDet(match);
        end

        marginDb = 10 * log10(bestPow / (bestTh + eps));
        clusterSize = numel(pix);
        peakList(end+1, :) = [ir, iv, bestPow, bestTh, marginDb, clusterSize]; %#ok<AGROW>
    end

    marginDbMin = 10;
    clusterSizeMin = 5;

    % Luu y:
    % RGPO/VGPO thuong lam nang luong peak bi trai rong hon false target tinh.
    % Neu can debug jammer, hay thu tam thoi giam:
    %   marginDbMin = 6;
    %   clusterSizeMin = 1;
    keep = peakList(:,5) >= marginDbMin & peakList(:,6) >= clusterSizeMin;
    peakList = peakList(keep, :);

    if isempty(peakList)
        detCellsFinal = [];
        fprintf('So detections sau OS-CFAR + clustering + filter = 0\n');
        return;
    end

    [~, order] = sort(peakList(:,3), 'descend');
    peakList = peakList(order, :);
    detCellsFinal = peakList(:,1:2).';

    fprintf('So detections sau OS-CFAR + clustering + filter = %d\n', size(detCellsFinal,2));
    
    detections = zeros(size(detCellsFinal,2), 2);

    for k = 1:size(detCellsFinal,2)
        ir = detCellsFinal(1,k);
        iv = detCellsFinal(2,k);
        detections(k,:) = [r_axis(ir), v_axis(iv)];
        fprintf(['Detection %d: R = %.2f m, V = %.2f m/s, P = %.3e, ' ...
                 'Margin = %.2f dB, ClusterSize = %d\n'], ...
                 k, r_axis(ir), v_axis(iv), peakList(k,3), peakList(k,5), peakList(k,6));
    end
end