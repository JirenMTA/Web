function cfg = drfm_default_config(fc, mode, Tchirp)
    c = physconst('LightSpeed');
    cfg.enable = false;
    cfg.mode = mode;
    cfg.lambda = c / fc;
    cfg.JSRdB = -18;
    cfg.TchirpForPhase = Tchirp;
    cfg.initialPhase = 0;

    % single/multi false
    cfg.falseRanges = 120;
    cfg.falseVelocities = 0;

    % rgpo/vgpo seeded values at scan-level
    cfg.baseRange = 200;
    cfg.baseVelocity = 20;
    cfg.pullOffRate = 5;       % m/s  (default vua phai de false target van gom peak tot)
    cfg.pullOffAccel = 0.02;    % m/s^2
    cfg.velocityWalk = 0.08;   % m/s/chirp
    cfg.velocityAccel = 0.002; % m/s/chirp^2

    % range smear
    cfg.centerRange = 140;
    cfg.spanRange = 20;
    cfg.numCopies = 6;
    cfg.smearVelocity = 0;
end