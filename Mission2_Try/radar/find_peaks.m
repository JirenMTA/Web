function [detCellsFinal, peakList, detections] = find_peaks(RDpow, r_axis, v_axis)
    % ---------------------------------------------------------------------
    % CLEAN / NEUTRAL peak finder
    %
    % Muc tieu:
    % - tra ve cac peak sach tu RD map
    % - khong bias theo false target / RGPO
    % - phu hop hon khi co nhieu target that
    %
    % Dau ra:
    %   detCellsFinal : 2xK [rangeIdx; dopplerIdx]
    %   peakList      : Kx6 [ir, iv, bestPow, bestTh, marginDb, clusterSize]
    %   detections    : Kx2 [range_m, vel_mps]
    % ---------------------------------------------------------------------

    detCellsFinal = [];
    peakList      = [];
    detections    = [];

    % -----------------------------
    % CFAR parameters
    % -----------------------------
    numTrain = [10 4];
    numGuard = [6 2];

    Nr = size(RDpow, 1);
    Nd = size(RDpow, 2);

    rMargin = numTrain(1) + numGuard(1);
    dMargin = numTrain(2) + numGuard(2);

    if Nr <= 2*rMargin || Nd <= 2*dMargin
        error('RD map qua nho cho cau hinh CFAR hien tai.');
    end

    % -----------------------------
    % Peak clean-up tuning
    % -----------------------------
    marginDbMin = 10;      % bo peak qua sat nguong
    suppressRangeBins = 2; % NMS theo range
    suppressDoppBins  = 1; % NMS theo Doppler

    % -----------------------------
    % CUT list for CFAR
    % -----------------------------
    rValid = (rMargin + 1) : (Nr - rMargin);
    dValid = (dMargin + 1) : (Nd - dMargin);
    [rr, dd] = ndgrid(rValid, dValid);
    cutidx = [rr(:).'; dd(:).'];

    % -----------------------------
    % OS-CFAR
    % -----------------------------
    trainCellsTotal = (2*(numTrain(1)+numGuard(1))+1) * ...
                      (2*(numTrain(2)+numGuard(2))+1) - ...
                      (2*numGuard(1)+1) * (2*numGuard(2)+1);
    rankOS = max(1, min(trainCellsTotal, round(0.75 * trainCellsTotal)));

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
        fprintf('So detections sau OS-CFAR + clustering + filter = 0\n');
        return;
    end

    % -----------------------------
    % Detection map va threshold map
    % -----------------------------
    detMap = false(size(RDpow));
    thMap  = nan(size(RDpow));

    linIdx = sub2ind(size(RDpow), detCells(1,:), detCells(2,:));
    detMap(linIdx) = true;
    thMap(linIdx)  = th(detFlags);

    % -----------------------------
    % Chi giu local maxima trong cac cell da qua CFAR
    % -----------------------------
    localMaxMap = imregionalmax(RDpow, 8) & detMap;
    if ~any(localMaxMap(:))
        localMaxMap = detMap;
    end

    [candR, candD] = find(localMaxMap);
    if isempty(candR)
        fprintf('So detections sau OS-CFAR + clustering + filter = 0\n');
        return;
    end

    % -----------------------------
    % Lap danh sach candidate peak
    % -----------------------------
    numCand = numel(candR);
    candPow    = zeros(numCand,1);
    candTh     = zeros(numCand,1);
    candMargin = zeros(numCand,1);

    for i = 1:numCand
        ir = candR(i);
        iv = candD(i);
        candPow(i)    = RDpow(ir, iv);
        candTh(i)     = thMap(ir, iv);
        candMargin(i) = 10 * log10(candPow(i) / (candTh(i) + eps));
    end

    % -----------------------------
    % Bo candidate yeu gan nguong
    % -----------------------------
    keep0 = candMargin >= marginDbMin;
    candR      = candR(keep0);
    candD      = candD(keep0);
    candPow    = candPow(keep0);
    candTh     = candTh(keep0);
    candMargin = candMargin(keep0);

    if isempty(candR)
        fprintf('So detections sau OS-CFAR + clustering + filter = 0\n');
        return;
    end

    % -----------------------------
    % NMS trung tinh: giu peak manh nhat trong vung lan can
    % -----------------------------
    [~, order] = sort(candPow, 'descend');

    selR = [];
    selD = [];
    selPow = [];
    selTh = [];
    selMargin = [];

    for kk = 1:numel(order)
        idx = order(kk);
        ir = candR(idx);
        iv = candD(idx);

        if isempty(selR)
            conflict = false;
        else
            conflict = any(abs(selR - ir) <= suppressRangeBins & ...
                           abs(selD - iv) <= suppressDoppBins);
        end

        if ~conflict
            selR(end+1,1)      = ir; %#ok<AGROW>
            selD(end+1,1)      = iv; %#ok<AGROW>
            selPow(end+1,1)    = candPow(idx); %#ok<AGROW>
            selTh(end+1,1)     = candTh(idx); %#ok<AGROW>
            selMargin(end+1,1) = candMargin(idx); %#ok<AGROW>
        end
    end

    if isempty(selR)
        fprintf('So detections sau OS-CFAR + clustering + filter = 0\n');
        return;
    end

    % -----------------------------
    % ClusterSize = 1 vi moi peak cuoi cung la 1 local maxima sau NMS
    % -----------------------------
    selClust = ones(size(selR));

    % Sort lai theo power giam dan
    [~, finalOrder] = sort(selPow, 'descend');
    selR      = selR(finalOrder);
    selD      = selD(finalOrder);
    selPow    = selPow(finalOrder);
    selTh     = selTh(finalOrder);
    selMargin = selMargin(finalOrder);
    selClust  = selClust(finalOrder);

    peakList = [selR, selD, selPow, selTh, selMargin, selClust];
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
