function [cfg, state] = update_jammer_strategy(cfg, state, scanIdx, currentTime, targets)
    if nargin < 2 || isempty(state)
        state = init_jammer_state();
    end

    state.scanIdx = scanIdx;
    state.currentTime = currentTime;

    switch lower(state.policy)
        case 'static'
            % Giữ nguyên cfg như hiện tại

        case 'single_false_linear_range'
            cfg.mode = 'single_false';
            cfg.falseRanges(1) = state.initialFalseRange + state.falseRangeRate * currentTime;
            cfg.falseVelocities(1) = state.falseVelocityConst;
            cfg.JSRdB = state.initialJSRdB;

        case 'single_false_linear_range_doppler'
            cfg.mode = 'single_false';
            cfg.falseRanges(1) = state.initialFalseRange + state.falseRangeRate * currentTime;
            cfg.falseVelocities(1) = state.initialFalseVelocity + state.falseVelocityRate * currentTime;
            cfg.JSRdB = state.initialJSRdB;

        case 'jsr_ramp'
            cfg.JSRdB = state.initialJSRdB + state.jsrRatePerSecond * currentTime;

        case 'mode_schedule_time'
            if currentTime < 0.9
                cfg.mov = 'single_false';
                cfg.falseRanges = 120;
                cfg.falseVelocities = 0;
                cfg.JSRdB = -18;
            elseif currentTime < 1.8
                cfg.mode = 'multi_false';
                cfg.falseRanges = [110 150 190];
                cfg.falseVelocities = [0 3 -2];
                cfg.JSRdB = -12;
            else
                cfg.mode = 'rgpo';
                cfg.baseRange = 120;
                cfg.baseVelocity = 0;
                cfg.pullOffRate = 0.8;
                cfg.pullOffAccel = 0.00;
                cfg.JSRdB = -10;
            end

        case 'targeted_follow_time'
            targetIdx = min(state.targetIdx, numel(targets));
            Rtrue = norm(targets(targetIdx).Range + targets(targetIdx).Velocity * currentTime);
            Vtrue = norm(targets(targetIdx).Velocity);

            cfg.mode = 'single_false';
            cfg.falseRanges(1) = Rtrue + state.rangeBias0 + state.rangeBiasRate * currentTime;
            cfg.falseVelocities(1) = Vtrue + state.velBias0 + state.velBiasRate * currentTime;
            cfg.JSRdB = state.initialJSRdB + state.jsrRatePerSecond * currentTime;

        case 'rgpo_time_seeded'
            cfg.mode = 'rgpo';
            cfg.baseRange = state.initialFalseRange + state.falseRangeRate * currentTime;
            cfg.baseVelocity = state.falseVelocityConst;
            cfg.pullOffRate = state.pullOffRate;      % m/s
            cfg.pullOffAccel = state.pullOffAccel;    % m/s^2
            cfg.JSRdB = state.initialJSRdB;

        case 'vgpo_time_seeded'
            cfg.mode = 'vgpo';
            cfg.baseRange = state.initialFalseRange;
            cfg.baseVelocity = state.initialFalseVelocity + state.falseVelocityRate * currentTime;
            cfg.velocityWalk = state.velocityWalk;
            cfg.velocityAccel = state.velocityAccel;
            cfg.JSRdB = state.initialJSRdB;

        case 'rgpo_vgpo_time_seeded'
            cfg.mode = 'rgpo_vgpo';
            cfg.baseRange = state.initialFalseRange + state.falseRangeRate * currentTime;
            cfg.baseVelocity = state.initialFalseVelocity + state.falseVelocityRate * currentTime;
            cfg.pullOffRate = state.pullOffRate;
            cfg.pullOffAccel = state.pullOffAccel;
            cfg.velocityWalk = state.velocityWalk;
            cfg.velocityAccel = state.velocityAccel;
            cfg.JSRdB = state.initialJSRdB + state.jsrRatePerSecond * currentTime;

        otherwise
            error('Khong ho tro state.policy = %s', state.policy);
    end

    if isfield(cfg, 'falseRanges')
        cfg.falseRanges = max(cfg.falseRanges, 0);
    end
    if isfield(cfg, 'baseRange')
        cfg.baseRange = max(cfg.baseRange, 0);
    end
end