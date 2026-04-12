function cfg = update_jammer_strategy(cfg, scanIdx, currentTime)
    cfg.scanIdx     = scanIdx;
    cfg.currentTime = currentTime;

    switch lower(cfg.mode)
        case 'rgpo'
            if currentTime < cfg.rgpoStartTime
                tPull = 0;
                deltaR = 0;
            else
                tPull = currentTime - cfg.rgpoStartTime;
                deltaR = cfg.pullOffRate * tPull + 0.5 * cfg.pullOffAccel * tPull^2;
            end

            % Gioi han muc pull-off
            if isfield(cfg, 'maxPullOff') && ~isempty(cfg.maxPullOff)
                deltaR = min(deltaR, cfg.maxPullOff);
            end

            % -----------------------------------------------------------
            % 3) Cap nhat false target
            %    RGPO: doi range la chinh, velocity nen giu gan target that
            % -----------------------------------------------------------
            cfg.falseRange       = cfg.baseRange + cfg.falseVelocity*currentTime + deltaR;


        otherwise
            error('Khong ho tro cfg.mode = %s', cfg.mode);
    end
end