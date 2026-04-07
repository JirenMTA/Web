function main()
    clc
    close all
    addpath("radar\");
    addpath("target\");
    addpath("signal\");
    addpath("drfm\");

    B = 40e6;
    fs = 80e6;
    Tchirp = 30e-6;
    Nchirp = 128;
    fc = 24e9;
    Nsample = Tchirp*fs;

    waveform = simulate_signal(Tchirp, B, fs, Nchirp);    
    targets(1).Object = phased.RadarTarget('MeanRCS', 10, 'OperatingFrequency', fc);
    targets(1).Velocity = [14; 0; 0];
    targets(1).Range = [253; 0; 0];

    targets(2).Object = phased.RadarTarget('MeanRCS', 20, 'OperatingFrequency', fc);
    targets(2).Velocity = [38; 0; 0];
    targets(2).Range = [213; 0; 0];

    targets(3).Object = phased.RadarTarget('MeanRCS', 13, 'OperatingFrequency', fc);
    targets(3).Velocity = [32; 0; 0];
    targets(3).Range = [321; 0; 0];
    
    jammerCfg.enable = true;
    jammerCfg.mode = 'single_false';
    jammerCfg.JSRdB = -20;          % bắt đầu âm, đừng dùng dương ngay
    jammerCfg.falseRanges = [280];
    jammerCfg.falseVelocities = [5];
    jammerCfg.TchirpForPhase = Tchirp;
    jammerCfg.initialPhase = 0;

    dt = 0.2;
    for i=1:1000
        radar(waveform, targets, fc, Nchirp, Nsample, jammerCfg);
        % for k=1:numel(targets)
        %     targets(k) = change_target_position(targets(k), dt);
        % end
        pause(dt);
    end
end