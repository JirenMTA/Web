% ========================================================================
% KIẾN TRÚC ĐỀ XUẤT
% ========================================================================
% 1) radar() giữ nguyên nhiệm vụ mô phỏng radar thật.
% 2) DRFM là module ngoài, ví dụ file: drfm_jammer.m
% 3) Trong radar(), tại mỗi chirp chỉ cần:
%
%    if jammerCfg.enable
%        jam_sig = drfm_jammer_module(txsig, fs, fc, i, Tchirp, jammerCfg);
%        total_rxsig = total_rxsig + jam_sig;
%    end
%
% 4) Như vậy:
%    - radar() không chứa logic chi tiết của DRFM
%    - dễ thay đổi thuật toán jammer
%    - dễ bật/tắt jammer để so sánh
%    - dễ mở rộng sau này thành nhiều loại ECM khác

% ========================================================================
% VÍ DỤ CÁCH GỌI TRONG radar()
% ========================================================================
% Chèn vào bên trong vòng lặp chirp của radar():
%
%   for i = 1:Nchirp
%       txsig = tx(sig_all(:, i));
%       txsig = radiator(txsig, [0;0]);
%       total_rxsig = complex(zeros(size(txsig)));
%
%       % Echo thật
%       for k = 1:numTargets
%           [tgt_pos, tgt_vel] = tgt_platform{k}(Tchirp);
%           sig_k = channel(txsig, tx_pos, tgt_pos, tx_vel, tgt_vel);
%           sig_k = targets(k).Object(sig_k);
%           total_rxsig = total_rxsig + sig_k;
%       end
%
%       % Gọi module DRFM bên ngoài
%       if jammerCfg.enable
%           raw_jam_sig = drfm_jammer_module(txsig, fs, fc, i, Tchirp, jammerCfg);
%           jam_sig = drfm_apply_jsr(raw_jam_sig, total_rxsig, jammerCfg.JSRdB);
%           total_rxsig = total_rxsig + jam_sig;
%       end
%
%       total_rxsig = collector(total_rxsig, [0;0]);
%       total_rxsig = rx(total_rxsig);
%       data_cube(:, i) = dechirp(total_rxsig, sig_all(:, i));
%   end

% ========================================================================
% FILE 1: drfm_jammer_module.m
% ========================================================================
function jam_sig = drfm_jammer_module(txsig, fs, fc, chirpIdx, Tchirp, cfg)
%DRFM_JAMMER_MODULE Module DRFM độc lập, được radar() gọi khi cần.
%
% Lưu ý quan trọng:
% - Hàm này chỉ sinh jammer thô (raw jammer).
% - Không nên dùng JSRdB ở đây để nhân biên độ trực tiếp với txsig,
%   vì txsig là tín hiệu phát rất mạnh, còn echo thật đã bị suy hao lớn qua kênh.
% - JSR nên được áp sau đó, bằng cách scale jammer theo echo thật hiện có.

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
            Rg = cfg.baseRange + ...
                 cfg.pullOffRate  * (chirpIdx - 1) + ...
                 cfg.pullOffAccel * (chirpIdx - 1)^2;
            Vg = cfg.baseVelocity;
            jam_sig = drfm_false_target(txsig, fs, lambda, Rg, Vg, 1, chirpIdx, cfg);

        case 'vgpo'
            Rg = cfg.baseRange;
            Vg = cfg.baseVelocity + ...
                 cfg.velocityWalk  * (chirpIdx - 1) + ...
                 cfg.velocityAccel * (chirpIdx - 1)^2;
            jam_sig = drfm_false_target(txsig, fs, lambda, Rg, Vg, 1, chirpIdx, cfg);

        case 'rgpo_vgpo'
            Rg = cfg.baseRange + ...
                 cfg.pullOffRate  * (chirpIdx - 1) + ...
                 cfg.pullOffAccel * (chirpIdx - 1)^2;
            Vg = cfg.baseVelocity + ...
                 cfg.velocityWalk  * (chirpIdx - 1) + ...
                 cfg.velocityAccel * (chirpIdx - 1)^2;
            jam_sig = drfm_false_target(txsig, fs, lambda, Rg, Vg, 1, chirpIdx, cfg);

        case 'range_smear'
            jam_sig = drfm_range_smear( ...
                txsig, fs, fc, ...
                cfg.centerRange, cfg.spanRange, cfg.numCopies, ...
                cfg.smearVelocity, chirpIdx, cfg);

        otherwise
            error('Che do DRFM khong hop le.');
    end
end


% ========================================================================
% FILE 2: drfm_false_target.m
% ========================================================================
function y = drfm_false_target(x, fs, lambda, Rfalse, Vfalse, amp, chirpIdx, cfg)
%DRFM_FALSE_TARGET Tạo 1 mục tiêu giả tại range/velocity chỉ định.
%
% Sửa quan trọng:
% - Không random pha mới ở mỗi chirp.
% - Pha Doppler phải liên tục theo thời gian chậm để jammer không bị bè ra toàn trục Doppler.

    c = physconst('LightSpeed');
    N = numel(x);
    t_fast = (0:N-1).' / fs;

    tau = 2 * Rfalse / c;
    fd  = 2 * Vfalse / lambda;

    x_delayed = fractional_delay(x, fs, tau);

    % Pha nhanh trong một chirp
    phase_fast = exp(1j * 2*pi * fd * t_fast);

    % Pha chậm giữa các chirp để giữ tính coherent
    t_slow = (chirpIdx - 1) * cfg.TchirpForPhase;
    phase_slow = exp(1j * 2*pi * fd * t_slow);

    % Pha cố định ban đầu nếu chưa có thì lấy 0 hoặc cfg.initialPhase
    if isfield(cfg, 'initialPhase')
        phi0 = cfg.initialPhase;
    else
        phi0 = 0;
    end

    y = amp * x_delayed .* phase_fast .* phase_slow * exp(1j*phi0);
end

% ========================================================================
% FILE 3: fractional_delay.m
% ========================================================================
function y = fractional_delay(x, fs, tau)
%FRACTIONAL_DELAY Tạo trễ phân số cho tín hiệu replay của DRFM.

    N = numel(x);
    t  = (0:N-1).' / fs;
    ts = t - tau;

    y_re = interp1(t, real(x), ts, 'linear', 0);
    y_im = interp1(t, imag(x), ts, 'linear', 0);
    y = complex(y_re, y_im);
end

% ========================================================================
% FILE 4: drfm_range_smear.m
% ========================================================================
function jam_sig = drfm_range_smear(txsig, fs, fc, centerRange, spanRange, numCopies, vel, chirpIdx, cfg)
%DRFM_RANGE_SMEAR Tạo nhiều bản sao gần nhau để làm bè mục tiêu theo range.

    c = physconst('LightSpeed');
    lambda = c / fc;

    ranges = linspace(centerRange - spanRange/2, centerRange + spanRange/2, numCopies);
    jam_sig = complex(zeros(size(txsig)));

    for k = 1:numCopies
        amp_k = 1 / sqrt(numCopies);
        jam_sig = jam_sig + drfm_false_target(txsig, fs, lambda, ranges(k), vel, amp_k, chirpIdx, cfg);
    end
end

% ========================================================================
% FILE 5: jsr_to_amplitude.m
% ========================================================================
function a = jsr_to_amplitude(JSRdB)
%JSR_TO_AMPLITUDE Đổi JSR từ dB sang biên độ tương đối.
    a = sqrt(10^(JSRdB/10));
end

% ========================================================================
% FILE 6: drfm_default_config.m
% ========================================================================
function cfg = drfm_default_config(fc, mode)
%DRFM_DEFAULT_CONFIG Tạo cấu hình mặc định cho module DRFM.

    c = physconst('LightSpeed');
    cfg.enable = true;
    cfg.mode = mode;
    cfg.lambda = c / fc;
    cfg.JSRdB = -10;           % bắt đầu rất nhỏ để test
    cfg.TchirpForPhase = 1e-4; % phải set đúng bằng Tchirp thực tế trong main/radar
    cfg.initialPhase = 0;

    % false target
    cfg.falseRanges = 80;
    cfg.falseVelocities = 0;

    % RGPO / VGPO
    cfg.baseRange = 35;
    cfg.baseVelocity = 0;
    cfg.pullOffRate = 4;
    cfg.pullOffAccel = 0.08;
    cfg.velocityWalk = 0.3;
    cfg.velocityAccel = 0.015;

    % range smear
    cfg.centerRange = 100;
    cfg.spanRange = 20;
    cfg.numCopies = 6;
    cfg.smearVelocity = 0;
end

% ========================================================================
% CÁCH DÙNG ĐỀ XUẤT
% ========================================================================
% Trong main hoặc script điều khiển:
%
%   jammerCfg = drfm_default_config(fc, 'rgpo');
%   [RDdB, r_axis, v_axis, detCellsFinal, peakList] = radar(waveform, targets, fc, Nchirp, Nsample, jammerCfg);
%
% Hoặc nếu muốn tắt jammer:
%
%   jammerCfg.enable = false;
%
% Nếu bạn chưa muốn đổi signature của radar(), cũng có thể cho jammerCfg là biến global
% hoặc biến struct khai báo trước trong workspace, nhưng cách tốt nhất vẫn là truyền qua tham số.

% ========================================================================
% GỢI Ý SỬA radar() CHO SẠCH
% ========================================================================
% Nên đổi đầu hàm thành:
%
% function [RDdB, r_axis, v_axis, detCellsFinal, peakList] = radar(waveform, targets, fc, Nchirp, Nsample, jammerCfg)
%
% Và ở đầu hàm:
%
% if nargin < 6 || isempty(jammerCfg)
%     jammerCfg.enable = false;
% end
%
% Như vậy radar() luôn chạy được ở cả 2 chế độ:
% - không jammer
% - có jammer

% ========================================================================
% Ý NGHĨA THIẾT KẾ NÀY
% ========================================================================
% Đây là cách đúng hơn về mặt kiến trúc:
% - radar = khối mô phỏng cảm biến
% - drfm_jammer_module = khối ECM / jammer
% - build_rd_map, find_peaks = khối xử lý tín hiệu và phát hiện
%
% Tức là mỗi module có nhiệm vụ riêng, dễ kiểm thử và dễ mở rộng.
