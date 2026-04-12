function main()
    clc;
    close all;

    addpath("drfm\");
    addpath("signal\");
    addpath("radar\");
    addpath("drfm\");
    addpath("trackers\");

    % =========================
    % 1) Radar parameters
    % =========================
    c = physconst('LightSpeed');
    fc = 24e9;
    B = 40e6;
    fs = 80e6;
    Tchirp = 30e-6;
    Nchirp = 128;
    Nsample = round(Tchirp * fs);

    % Chu kỳ giữa 2 lần scan
    numScans = 200;
    scanPeriod = 0.3;   % giây

    % =========================
    % 2) Waveform
    % =========================
    waveform = simulate_signal(Tchirp, B, fs, Nchirp);

    % =========================
    % 3) Real targets at t = 0
    %    Range là vector vị trí ban đầu [x;y;z]
    % =========================
    targets(1).Object = phased.RadarTarget('MeanRCS', 10, 'OperatingFrequency', fc);
    targets(1).Velocity = [14; 0; 0];
    targets(1).Range    = [253; 0; 0];

    targets(2).Object = phased.RadarTarget('MeanRCS', 20, 'OperatingFrequency', fc);
    targets(2).Velocity = [38; 0; 0];
    targets(2).Range    = [213; 0; 0];
    % 
    % targets(3).Object = phased.RadarTarget('MeanRCS', 13, 'OperatingFrequency', fc);
    % targets(3).Velocity = [32; 0; 0];
    % targets(3).Range    = [321; 0; 0];

    % =========================
    % 4) Jammer config + state
    % =========================
    dt = 0.2;
    sigma_a = 2;
    Q = sigma_a^2 * [dt^4/4, dt^3/2; dt^3/2, dt^2];
    R = diag([2.0, 0.5]);

    tmCfg.confirmHits = 3;
    tmCfg.maxMisses   = 3;
    tmCfg.initCov     = diag([100, 10]);

    gate.sigR      = 5;
    gate.sigV      = 1;
    gate.threshold = 9.21;

    [jammerCfg, jammerState] = init_config_and_state(fc, Tchirp);

    % Chọn policy update theo thời gian
    tracks = initTrackManager();

    % [tmCfg,  gate, tracks] = init_atracker();
    
    for scanIdx = 1:numScans
        currentTime = (scanIdx - 1) * dt;

        [jammerCfg, jammerState] = update_jammer_strategy( ...
            jammerCfg, jammerState, scanIdx, currentTime, targets);

        [~, detections] = radar( ...
            waveform, targets, fc, Nchirp, Nsample, currentTime, jammerCfg);
        
        if isempty(detections)
            detections = zeros(0, 2);
        end
        
        [assocIdx, unassocDet] = associateDetections(tracks, detections, gate);
        tracks = updateTracks(tracks, detections, assocIdx, ...
                              unassocDet, dt, Q, R, tmCfg);
        show_tracks(tracks);
        fprintf("------------------------------------------\n");
    end
end