function [RDdB, detections, peakList] = radar( ...
    waveform, targets, fc, Nchirp, Nsample, jammerCfg)

    % =========================
    % 1) Generate transmit burst
    % =========================
    sig_all = waveform();
    sig_all = reshape(sig_all, Nsample, Nchirp);

    fs = waveform.SampleRate;
    Tchirp = waveform.SweepTime;
    B = waveform.SweepBandwidth;

    % Luu y: function nay co 6 input, nen phai check nargin < 6
    if nargin < 6 || isempty(jammerCfg)
        jammerCfg = struct('enable', false);
    end
    if ~isfield(jammerCfg, 'enable')
        jammerCfg.enable = false;
    end
    if ~isfield(jammerCfg, 'TchirpForPhase')
        jammerCfg.TchirpForPhase = Tchirp;
    end

    % =========================
    % 2) Radar front-end
    % =========================
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

    % =========================
    % 3) Real target states at current scan time
    % =========================
    numTargets = numel(targets);
    tgt_platform = cell(1, numTargets);
    for k = 1:numTargets
        tgt_platform{k} = phased.Platform( ...
            'InitialPosition', targets(k).Range, ...
            'Velocity', targets(k).Velocity);
    end

    % =========================
    % 4) Simulate data cube
    % =========================
    data_cube = complex(zeros(Nsample, Nchirp));

    for i = 1:Nchirp
        txsig = tx(sig_all(:, i));
        txsig = radiator(txsig, [0;0]);

        total_rxsig = complex(zeros(size(txsig)));
        victim_echo = complex(zeros(size(txsig)));

        % ---- real echoes ----
        for k = 1:numTargets
            [tgt_pos, tgt_vel] = tgt_platform{k}(Tchirp);
            sig_k = channel(txsig, tx_pos, tgt_pos, tx_vel, tgt_vel);
            sig_k = targets(k).Object(sig_k);

            total_rxsig = total_rxsig + sig_k;

            % target 1 la victim mac dinh
            if k == 1
                victim_echo = sig_k;
            end
        end

        % ---- jammer ----
        if jammerCfg.enable
            raw_jam_sig = drfm_jammer_module(txsig, fs, fc, i, jammerCfg);
            jam_sig = drfm_apply_jsr(raw_jam_sig, victim_echo, jammerCfg.JSRdB);
            total_rxsig = total_rxsig + jam_sig;
        end

        total_rxsig = collector(total_rxsig, [0;0]);
        total_rxsig = rx(total_rxsig);
        data_cube(:, i) = dechirp(total_rxsig, sig_all(:, i));
    end

    % =========================
    % 5) Build RD map
    % =========================
    [RDdB, r_axis, v_axis, RDpow] = build_rd_map(data_cube, fc, B, Tchirp, fs);

    % =========================
    % 6) Detect peaks
    % =========================
    [~, peakList, detections] = find_peaks(RDpow, r_axis, v_axis);
end
