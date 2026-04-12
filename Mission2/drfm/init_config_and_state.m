function [jammerCfg, jammerState] = init_config_and_state(fc, Tchirp)
    jammerCfg = drfm_default_config(fc, 'rgpo', Tchirp);
    jammerState.policy = 'rgpo_time_seeded';

    % Ví dụ update theo thời gian mô phỏng:
    jammerState.initialFalseRange = 280;
    jammerState.falseRangeRate = 20;          % m/s
    jammerState.initialFalseVelocity = 10;
    jammerState.falseVelocityRate = 1.5;      % m/s^2 nếu hiểu v(t)=v0+a*t
    jammerState.falseVelocityConst = 6;       % dùng cho policy cố định vận tốc
    jammerState.initialJSRdB = -18;
    jammerState.jsrRatePerSecond = 2;

    % Nếu muốn RGPO/VGPO seeded theo thời gian:
    jammerState.pullOffRate = 0.8;            % m/chirp
    jammerState.pullOffAccel = 0.01;          % m/chirp^2
    jammerState.velocityWalk = 0.08;          % m/s/chirp
    jammerState.velocityAccel = 0.002;        % m/s/chirp``````````^2
end