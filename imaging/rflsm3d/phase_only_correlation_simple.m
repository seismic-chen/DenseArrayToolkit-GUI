function [delay_time, poc_result, lag] = phase_only_correlation_simple(sig1, sig2, dt)
    sig1 = sig1(:); sig2 = sig2(:);
    N = min(length(sig1), length(sig2));
    sig1 = sig1(1:N); sig2 = sig2(1:N);
    sig1 = sig1 - mean(sig1); sig2 = sig2 - mean(sig2);
    win = hanning(N); 
    sig1 = sig1 .* win; sig2 = sig2 .* win;
    
    % 零填充到2N
    N_fft = 2 * N;
    F1 = fft(sig1, N_fft); 
    F2 = fft(sig2, N_fft);
    
    % 修正归一化公式
    cross_power = F2 .* conj(F1);
    R = cross_power ./ (abs(F1) .* abs(F2) + 1e-12);
    
    poc = real(ifft(R));
    poc = fftshift(poc);
    lag = (-floor(N_fft/2)) : (N_fft - floor(N_fft/2) - 1); 
    lag = lag(:);
    
    % 抛物线插值
    [~, max_idx] = max(poc);
    idx_range = max(1, max_idx-1) : min(N_fft, max_idx+1);
    p = polyfit(lag(idx_range), poc(idx_range), 2);
    peak_offset = -p(2) / (2 * p(1)); % 顶点位置
    
    delay_time = peak_offset * dt;
    poc_result = poc;
end