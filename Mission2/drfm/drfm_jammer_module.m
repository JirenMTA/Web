function jam_sig = drfm_jammer_module(txsig, fs, fc, chirpIdx, cfg)
    c = physconst('LightSpeed');
    lambda = c / fc;

    switch lower(cfg.mode)
        case 'single_false'
            Rf = cfg.falseRanges(1);
            if isfield(cfg, 'falseVelocities') && ~isempty(cfg.falseVelocities)
                Vf = cfg.falseVelocities(1);
            else
                Vf = 0;
            end
            jam_sig = drfm_false_target(txsig, fs, lambda, Rf, Vf, 1, chirpIdx, cfg);

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
                amp_k = 1 / max(1, sqrt(nT));
                jam_sig = jam_sig + drfm_false_target(txsig, fs, lambda, Rf, Vf, amp_k, chirpIdx, cfg);
            end

        case 'rgpo'
            % RGPO: range gia thay doi theo thoi gian cham.
            % DUNG MO HINH VAT LY HON:
            %   Rg(t) = R0 + v_pull*t + 0.5*a_pull*t^2
            %   Vg(t) = Vbase + v_pull + a_pull*t
            % Neu chi doi R ma khong doi V tuong ung, peak jammer se mat coherence.
            t_slow = (chirpIdx - 1) * cfg.TchirpForPhase;
            Rg = cfg.baseRange + cfg.pullOffRate * t_slow + 0.5 * cfg.pullOffAccel * t_slow^2;
            Vg = cfg.baseVelocity + cfg.pullOffRate + cfg.pullOffAccel * t_slow;
            jam_sig = drfm_false_target(txsig, fs, lambda, Rg, Vg, 1, chirpIdx, cfg);

        case 'vgpo'
            Rg = cfg.baseRange;
            Vg = cfg.baseVelocity + cfg.velocityWalk*(chirpIdx - 1) + cfg.velocityAccel*(chirpIdx - 1)^2;
            jam_sig = drfm_false_target(txsig, fs, lambda, Rg, Vg, 1, chirpIdx, cfg);

        case 'rgpo_vgpo'
            % RGPO + VGPO: range va velocity deu duoc keo.
            % Phan RGPO phai coherent ve mat dao ham range -> velocity.
            t_slow = (chirpIdx - 1) * cfg.TchirpForPhase;
            Rg = cfg.baseRange + cfg.pullOffRate * t_slow + 0.5 * cfg.pullOffAccel * t_slow^2;
            Vg_rgpo = cfg.baseVelocity + cfg.pullOffRate + cfg.pullOffAccel * t_slow;
            Vg_vgpo = cfg.velocityWalk*(chirpIdx - 1) + cfg.velocityAccel*(chirpIdx - 1)^2;
            Vg = Vg_rgpo + Vg_vgpo;
            jam_sig = drfm_false_target(txsig, fs, lambda, Rg, Vg, 1, chirpIdx, cfg);

        case 'range_smear'
            jam_sig = drfm_range_smear(txsig, fs, fc, cfg.centerRange, cfg.spanRange, ...
                cfg.numCopies, cfg.smearVelocity, chirpIdx, cfg);

        otherwise
            error('Che do DRFM khong hop le: %s', cfg.mode);
    end
end
