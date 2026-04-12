function [tmCfg,  gate, tracks] = init_tracker(dt)
    tmCfg.dt = dt;          % thời gian giữa 2 scan
    tmCfg.Q  = diag([0.1, 0.01]);        % process noise [range, velocity]
    tmCfg.R  = diag([2.0, 0.5]);         % measurement noise
    
    % --- Cấu hình Track Management ---
    tmCfg.confirmHits = 3;
    tmCfg.maxMisses   = 3;
    tmCfg.initCov     = diag([100, 10]);
    
    % --- Cấu hình Association Gate ---
    gate.sigR      = 10;    % m
    gate.sigV      = 2;    % m/s
    gate.threshold = 9.21; % chi2 95% với dof=2
    
    % --- Khởi tạo ---
    tracks = initTrackManager();
end 