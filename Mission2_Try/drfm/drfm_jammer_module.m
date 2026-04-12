function jam_sig = drfm_jammer_module(txsig, fs, fc, chirpIdx, cfg)
    c = physconst('LightSpeed');
    lambda = c / fc;

    % Bien doi JSR (power ratio, dB) -> he so bien do
    amp = 10^(cfg.JSRdB / 20);

    switch lower(cfg.mode)
        case 'single_false'
            Rf = cfg.falseRanges(1);
            if isfield(cfg, 'falseVelocities') && ~isempty(cfg.falseVelocities)
                Vf = cfg.falseVelocities(1);
            else
                Vf = 0;
            end
            jam_sig = drfm_false_target(txsig, fs, lambda, Rf, Vf, amp, chirpIdx, cfg);

        case 'multi_false'
            jam_sig = complex(zeros(size(txsig)));
            nT = numel(cfg.falseRanges);
            for k = 1:nT
                Rf = cfg.falseRanges(k);
                if isfield(cfg, 'falseVelocities') && numel(cfg.falseVelocities) >= k
                    Vf = cfg.falseVelocities(k);
                else
                    Vf = 0;
                end

                % Chia nang luong cho nhieu false target
                amp_k = amp / max(1, sqrt(nT));
                jam_sig = jam_sig + drfm_false_target(txsig, fs, lambda, Rf, Vf, amp_k, chirpIdx, cfg);
            end

        case 'rgpo'
            % ===========================================================
            % RGPO:
            % - Khong tu tinh lai range gia theo chirpIdx o day
            % - Range gia phai duoc cap nhat truoc do boi update_jammer_strategy
            %   va luu trong cfg.falseRange
            % - Velocity nen giu gan muc tieu that de deception on dinh
            % ===========================================================
            if ~isfield(cfg, 'falseRange') || isempty(cfg.falseRange)
                Rg = cfg.baseRange;
            else
                Rg = cfg.falseRange;
            end

            if ~isfield(cfg, 'falseVelocity') || isempty(cfg.falseVelocity)
                Vg = cfg.baseVelocity;
            else
                Vg = cfg.falseVelocity;
            end

            jam_sig = drfm_false_target(txsig, fs, lambda, Rg, Vg, amp, chirpIdx, cfg);

        case 'vgpo'
            % ===========================================================
            % VGPO:
            % - Range giu nguyen
            % - Velocity gia co the duoc cap nhat san trong cfg.falseVelocity
            % - Neu chua co, moi fallback sang cong thuc cu
            % ===========================================================
            if isfield(cfg, 'falseRange') && ~isempty(cfg.falseRange)
                Rg = cfg.falseRange;
            else
                Rg = cfg.baseRange;
            end

            if isfield(cfg, 'falseVelocity') && ~isempty(cfg.falseVelocity)
                Vg = cfg.falseVelocity;
            else
                Vg = cfg.baseVelocity + ...
                     cfg.velocityWalk  * (chirpIdx - 1) + ...
                     cfg.velocityAccel * (chirpIdx - 1)^2;
            end

            jam_sig = drfm_false_target(txsig, fs, lambda, Rg, Vg, amp, chirpIdx, cfg);

        case 'rgpo_vgpo'
            % ===========================================================
            % RGPO + VGPO:
            % - Neu cfg da co falseRange / falseVelocity thi uu tien dung
            % - Khong tinh lai song song mot logic RGPO khac tai day
            % ===========================================================
            if isfield(cfg, 'falseRange') && ~isempty(cfg.falseRange)
                Rg = cfg.falseRange;
            else
                Rg = cfg.baseRange;
            end

            if isfield(cfg, 'falseVelocity') && ~isempty(cfg.falseVelocity)
                Vg = cfg.falseVelocity;
            else
                Vg = cfg.baseVelocity + ...
                     cfg.velocityWalk  * (chirpIdx - 1) + ...
                     cfg.velocityAccel * (chirpIdx - 1)^2;
            end

            jam_sig = drfm_false_target(txsig, fs, lambda, Rg, Vg, amp, chirpIdx, cfg);

        case 'range_smear'
            jam_sig = drfm_range_smear(txsig, fs, fc, cfg.centerRange, cfg.spanRange, ...
                cfg.numCopies, cfg.smearVelocity, chirpIdx, cfg);
            jam_sig = amp * jam_sig;

        otherwise
            error('Che do DRFM khong hop le: %s', cfg.mode);
    end
end
