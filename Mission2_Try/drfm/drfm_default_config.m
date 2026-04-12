function cfg = drfm_default_config(fc, mode, Tchirp)
    c = physconst('LightSpeed');

    cfg.enable = true;
    cfg.mode   = lower(mode);
    cfg.lambda = c / fc;

    cfg.JSRdB        = 14;
    cfg.initialJSRdB = cfg.JSRdB;

    cfg.TchirpForPhase = Tchirp;
    cfg.initialPhase   = 0;
    cfg.scanStartTime  = 0;

    cfg.victimIndex   = 1;
    cfg.victimTrackId = [];

    cfg.baseRange     = 253;
    cfg.baseVelocity  = 14;
    cfg.falseRange    = cfg.baseRange;
    cfg.falseVelocity = cfg.baseVelocity;

    cfg.rgpoStartTime = 0.0;
    cfg.maxPullOff    = 150;

    cfg.pullOffRate   = 8.0;    % m/s
    cfg.pullOffAccel  = 0.02;   % m/s^2

    cfg.velocityWalk  = 0.0;    % m/s per scan-step
    cfg.velocityAccel = 0.0;    % m/s^2

    cfg.lastUpdateTime   = 0;
    cfg.isLockedToVictim = false;

    switch cfg.mode
        case 'rgpo'
            % keep defaults

        case 'vgpo'
            cfg.pullOffRate   = 0.0;
            cfg.pullOffAccel  = 0.0;
            cfg.velocityWalk  = 0.5;
            cfg.velocityAccel = 0.02;

        case 'single_false'
            cfg.pullOffRate   = 0.0;
            cfg.pullOffAccel  = 0.0;
            cfg.velocityWalk  = 0.0;
            cfg.velocityAccel = 0.0;
    end
end