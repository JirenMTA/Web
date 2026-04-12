function y = drfm_false_target(x, fs, lambda, Rfalse, Vfalse, amp, chirpIdx, cfg)
    c = physconst('LightSpeed');
    N = numel(x);
    t_fast = (0:N-1).' / fs;

    if isfield(cfg, 'TchirpForPhase')
        Tslow = cfg.TchirpForPhase;
    else
        Tslow = 0;
    end

    if isfield(cfg, 'scanStartTime')
        t0 = cfg.scanStartTime;
    else
        t0 = 0;
    end

    if isfield(cfg, 'initialPhase')
        phi0 = cfg.initialPhase;
    else
        phi0 = 0;
    end

    tau = 2 * Rfalse / c;
    fd  = 2 * Vfalse / lambda;

    % delay for false range
    x_delayed = fractional_delay(x, fs, tau);

    % total time for Doppler coherence
    t = t0 + (chirpIdx - 1) * Tslow + t_fast;

    % carrier propagation phase for false range
    phase_prop = exp(-1j * 4*pi * Rfalse / lambda);

    % Doppler phase
    phase_dopp = exp(-1j * 2*pi * fd * t);

    y = amp * x_delayed .* phase_prop .* phase_dopp * exp(1j * phi0);
end