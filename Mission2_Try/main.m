function main()
    clc;
    close all;

    addpath("drfm\");
    addpath("signal\");
    addpath("radar\");
    addpath("drfm\");
    addpath("target\");
    addpath("trackers\");

    fc = 24e9;
    B = 30e6;
    fs = 60e6;
    Tchirp = 30e-6;
    Nchirp = 100;
    Nsample = round(Tchirp * fs);
    waveform = simulate_signal(Tchirp, B, fs, Nchirp);

    targets(1).Object = phased.RadarTarget('MeanRCS', 10, 'OperatingFrequency', fc);
    targets(1).Velocity = [14; 0; 0];
    targets(1).Range    = [253; 0; 0];

    % targets(2).Object = phased.RadarTarget('MeanRCS', 20, 'OperatingFrequency', fc);
    % targets(2).Velocity = [25; 0; 0];
    % targets(2).Range    = [213; 0; 0];

    dt = 0.2;
    sigma_a = 10;
    Q = sigma_a^2 * [dt^4/4, dt^3/2; dt^3/2, dt^2];
    R = diag([1.0, 0.5]);

    tmCfg.confirmHits    = 3;
    tmCfg.maxMisses      = 3;
    tmCfg.initCov        = diag([100, 10]);
    tmCfg.minInitMarginDb = 15;
    tmCfg.minInitRangeSep = 4.0;
    tmCfg.minInitVelSep   = 1.5;

    gate.sigR      = 8;
    gate.sigV      = 3;
    gate.threshold = 9.21;

    rgpoCfg.enable            = true;
    rgpoCfg.victimTrackId     = 1;
    rgpoCfg.velLane           = 1.5;
    rgpoCfg.minPullOff        = 2.0;
    rgpoCfg.maxBackStep       = 1.0;
    rgpoCfg.minVictimMarginDb = 0;
    rgpoCfg.confirmedOnly     = true;

    jammerCfg = drfm_default_config(fc, 'rgpo', Tchirp);
    tracks = initTrackManager();

    for scanIdx = 1:200
        currentTime = (scanIdx - 1) * dt;

        [~, detections, peakList] = radar( ...
            waveform, targets, fc, Nchirp, Nsample, jammerCfg);

        jammerCfg = update_jammer_strategy(jammerCfg, scanIdx, currentTime);

        if isempty(detections)
            detections = zeros(0, 2);
        end
        if isempty(peakList)
            peakList = zeros(0, 6);
        end

        [assocIdx, unassocDet, predCache] = associateDetections( ...
            tracks, detections, peakList, dt, Q, gate, rgpoCfg);

        tracks = updateTracks( ...
            tracks, detections, peakList, assocIdx, unassocDet, predCache, R, tmCfg);

        for i = 1:numel(targets)
            fprintf("Current [Target]: r=%0.2f; v=%0.2f \n", targets(i).Range(1), targets(i).Velocity(1));
        end
        fprintf("Current [Jammer]: r=%0.2f; v=%0.2f \n", jammerCfg.falseRange, jammerCfg.falseVelocity);
        show_tracks(tracks);

        for i = 1:numel(targets)
            targets(i) = change_target_position(targets(i), dt);
        end

        fprintf("------------------------------------------\n");
    end
end
