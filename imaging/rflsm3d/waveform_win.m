function [win] = waveform_win(data,t,t1,t2)
    % data:
    % t: time
    % tstart:
    % tend: 
    % td: begin to decrase (sec)

    alpha = 0.3; 
    t_rise = [t(1), t1];
    t_fall = [t2, t(end)];

    dt = t(2) - t(1);
    win = zeros(size(data(:), 1), 1);

    toIdx = @(tau) round((tau - t(1))/dt) + 1;

    i_rise_start = toIdx(t_rise(1));
    i_rise_end = toIdx(t_rise(2));
    i_fall_start = toIdx(t_fall(1));
    i_fall_end = toIdx(t_fall(2));

    L_rise = i_rise_end - i_rise_start + 1;
    L_fall = i_fall_end - i_fall_start + 1;

    w_up_full = tukeywin(2*L_rise - 1, alpha);
    w_up = w_up_full(1:L_rise);

    w_down_full = tukeywin(2*L_fall - 1, alpha);
    w_down = w_down_full(L_fall:end);

    i_plateau_start = i_rise_end + 1;
    i_plateau_end = i_fall_start - 1;

    N = length(win);
    i_rise_start = max(1, i_rise_start);
    i_rise_end = min(N, i_rise_end);
    i_plateau_start = max(1, i_plateau_start);
    i_plateau_end = min(N, i_plateau_end);
    i_fall_start = max(1, i_fall_start);
    i_fall_end = min(N, i_fall_end);

    if i_rise_start <= i_rise_end
        win(i_rise_start:i_rise_end) = w_up;
    end
    if i_plateau_start <= i_plateau_end
        win(i_plateau_start:i_plateau_end) = 1;
    end
    if i_fall_start <= i_fall_end
        win(i_fall_start:i_fall_end) = w_down;
    end
end