function y = drfm_false_target(x, fs, lambda, Rfalse, Vfalse, amp, chirpIdx, cfg)
    c = physconst('LightSpeed');
    N = numel(x);
    t_fast = (0:N-1).' / fs;

    tau = 2 * Rfalse / c;
    fd = 2 * Vfalse / lambda;

    x_delayed = fractional_delay(x, fs, tau);
    phase_fast = exp(1j * 2*pi * fd * t_fast);
    t_slow = (chirpIdx - 1) * cfg.TchirpForPhase;
    phase_slow = exp(1j * 2*pi * fd * t_slow);

    if isfield(cfg, 'initialPhase')
        phi0 = cfg.initialPhase;
    else
        phi0 = 0;
    end

    y = amp * x_delayed .* phase_fast .* phase_slow * exp(1j * phi0);
end
