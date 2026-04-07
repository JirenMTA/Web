function [RDdB, r_axis, v_axis, detCellsFinal, peakList] = radar(waveform, targets, fc, Nchirp, Nsample, jammerCfg)
    % =========================
    % 1) Sinh data cube từ mô phỏng radar
    % =========================
    sig_all = waveform();
    sig_all = reshape(sig_all, Nsample, Nchirp);

    fs = waveform.SampleRate;
    Tchirp = waveform.SweepTime;
    B = waveform.SweepBandwidth;

    tx = phased.Transmitter('PeakPower', 1, 'Gain', 30);
    rx = phased.ReceiverPreamp('Gain', 30, 'NoiseFigure', 4, 'SampleRate', fs);
    ant = phased.IsotropicAntennaElement;
    radiator = phased.Radiator('Sensor', ant, 'OperatingFrequency', fc);
    collector = phased.Collector('Sensor', ant, 'OperatingFrequency', fc);
    channel = phased.FreeSpace( ...
        'SampleRate', fs, ...
        'OperatingFrequency', fc, ...
        'TwoWayPropagation', true);

    tx_pos = [0;0;0];
    tx_vel = [0;0;0];

    numTargets = numel(targets);
    tgt_platform = cell(1, numTargets);

    for i = 1:numTargets
        tgt_platform{i} = phased.Platform( ...
            'InitialPosition', targets(i).Range, ...
            'Velocity', targets(i).Velocity);
    end

    data_cube = zeros(Nsample, Nchirp);

    for i = 1:Nchirp
        txsig = tx(sig_all(:, i));
        txsig = radiator(txsig, [0;0]);
        total_rxsig = complex(zeros(size(txsig)));

        for k = 1:numTargets
            [tgt_pos, tgt_vel] = tgt_platform{k}(Tchirp);
            sig_k = channel(txsig, tx_pos, tgt_pos, tx_vel, tgt_vel);
            sig_k = targets(k).Object(sig_k);
            total_rxsig = total_rxsig + sig_k;
        end
        
        if jammerCfg.enable
            raw_jam_sig = drfm_jammer_module(txsig, fs, fc, i, Tchirp, jammerCfg);
            jam_sig = drfm_apply_jsr(raw_jam_sig, total_rxsig, jammerCfg.JSRdB);
            total_rxsig = total_rxsig + jam_sig;
        end
        

        total_rxsig = collector(total_rxsig, [0;0]);
        total_rxsig = rx(total_rxsig);
        data_cube(:, i) = dechirp(total_rxsig, sig_all(:, i));
    end

    % =========================
    % 2) Build RD map
    % =========================
    [RDdB, r_axis, v_axis, RDpow] = build_rd_map(data_cube, fc, B, Tchirp);

    % =========================
    % 3) Find peaks
    % =========================
    [detCellsFinal, peakList] = find_peaks(RDpow, r_axis, v_axis);

    % =========================
    % 4) Plot detections
    % =========================
    plot_peaks_to_rdmap(RDdB, r_axis, v_axis, detCellsFinal);
end

function [RDdB, r_axis, v_axis, RDpow] = build_rd_map(data_cube, fc, B, Tchirp)
    c = physconst('LightSpeed');
    lambda = c / fc;

    [Nsample, Nchirp] = size(data_cube);

    % =========================
    % FFT Range
    % =========================
    winK = hann(Nsample);
    data_cube = data_cube .* winK;

    sig_range = fft(data_cube, Nsample, 1);
    sig_range = sig_range / sum(winK);
    sig_range = sig_range(1:Nsample/2, :);

    % =========================
    % FFT Doppler
    % =========================
    winL = hann(Nchirp)';
    sig_range = sig_range .* winL;

    sig_rd = fft(sig_range, Nchirp, 2);
    sig_rd = sig_rd / sum(winL);
    sig_rd = fftshift(sig_rd, 2);

    % =========================
    % Axes
    % =========================
    dR = c / (2 * B);
    r_axis = (0:Nsample/2 - 1) * dR;

    PRF = 1 / Tchirp;
    fd_axis = (-Nchirp/2 : Nchirp/2 - 1) * (PRF / Nchirp);
    v_axis = fd_axis * lambda / 2;

    % =========================
    % RD map
    % =========================
    RDpow = abs(sig_rd).^2;
    RDdB  = 10 * log10(RDpow + 1e-12);
end

function [detCellsFinal, peakList] = find_peaks(RDpow, r_axis, v_axis)
    % =========================
    % OS-CFAR parameters
    % =========================
    numTrain = [10 4];
    numGuard = [6 2];

    Nr = size(RDpow, 1);
    Nd = size(RDpow, 2);

    rMargin = numTrain(1) + numGuard(1);
    dMargin = numTrain(2) + numGuard(2);

    if Nr <= 2*rMargin || Nd <= 2*dMargin
        error('RD map quá nhỏ cho cấu hình CFAR hiện tại. Hãy giảm numTrain/numGuard.');
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

    % =========================
    % Binary detection map
    % =========================
    linIdx = sub2ind(size(RDpow), detCells(1,:), detCells(2,:));
    detMap = false(size(RDpow));
    detMap(linIdx) = true;

    % =========================
    % Clustering detections
    % =========================
    CC = bwconncomp(detMap, 8);

    peakList = [];   % [ir iv pow th marginDb clusterSize]
    thDet = th(detFlags);

    for ic = 1:CC.NumObjects
        pix = CC.PixelIdxList{ic};

        % local max mạnh nhất trong cụm
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

    % =========================
    % Filter peaks
    % =========================
    marginDbMin = 10;
    clusterSizeMin = 5;

    keep = peakList(:,5) >= marginDbMin & peakList(:,6) >= clusterSizeMin;
    peakList = peakList(keep, :);

    if isempty(peakList)
        detCellsFinal = [];
        fprintf('So detections sau OS-CFAR + clustering + filter = 0\n');
        return;
    end

    % sort theo power giảm dần
    [~, order] = sort(peakList(:,3), 'descend');
    peakList = peakList(order, :);
    detCellsFinal = peakList(:,1:2).';

    % =========================
    % Print detection info
    % =========================
    fprintf('So detections sau OS-CFAR + clustering + filter = %d\n', size(detCellsFinal,2));
    for k = 1:size(detCellsFinal,2)
        ir = detCellsFinal(1,k);
        iv = detCellsFinal(2,k);
        fprintf(['Detection %d: R = %.2f m, V = %.2f m/s, P = %.3e, ' ...
                 'Margin = %.2f dB, ClusterSize = %d\n'], ...
                 k, r_axis(ir), v_axis(iv), ...
                 peakList(k,3), peakList(k,5), peakList(k,6));
    end
end

function plot_peaks_to_rdmap(RDdB, r_axis, v_axis, detCellsFinal)
    figure;
    imagesc(v_axis, r_axis, RDdB);
    axis xy;
    colorbar;
    colormap jet;
    xlabel('Velocity (m/s)');
    ylabel('Range (m)');
    title('Range-Doppler Map với OS-CFAR + clustering');
    hold on;

    if ~isempty(detCellsFinal)
        plot(v_axis(detCellsFinal(2,:)), r_axis(detCellsFinal(1,:)), ...
            'ws', 'MarkerSize', 8, 'LineWidth', 1.5);

        for k = 1:size(detCellsFinal, 2)
            text(v_axis(detCellsFinal(2,k)), r_axis(detCellsFinal(1,k)), ...
                sprintf(' %d', k), ...
                'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
        end
    end
end