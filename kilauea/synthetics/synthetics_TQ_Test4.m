close all
clear
clc

% ---- specify parameters ----

window_multiplier = 10; % this determines where the time series gets trimmed

T_ac = 2.0; % acoustic mode T (sec)
omega1 = 2*pi/T_ac;
Q_ac = 10; % acoustic mode quality factor
alpha1 = omega1/(2*Q_ac);

T_cr = 40.0; % conduit reservoir T (sec)
omega2 = 2*pi/T_cr;
Q_cr = 5; % condiuit reservoir mode quality factor
alpha2 = omega2/(2*Q_cr);

dt = 0.01;
Ttotal = -log(0.001)/min(alpha1, alpha2);  % amount of decay
event_t = 0:dt:Ttotal;

noise_level = 35; % noise level for awgn function

%% ---- build time series with two pulses ----

% build pulses
F_event1 = sin(omega1.*event_t).*exp(-alpha1*event_t);
F_event2 = sin(omega2.*event_t).*exp(-alpha2*event_t);

% pad the signal
total_duration = 1000; % set to 9 minutes to give time for cr decay (sec)
times = round(total_duration / dt);

F_ac= zeros(1, times);
start_idx = floor((times - length(F_event1)) / 2) + 1;
F_ac(start_idx:start_idx + length(F_event1) - 1) = F_event1;
F_ac = awgn(F_ac,noise_level); % add gaussian white noise to ac timeseries

F_cr= zeros(1, times);
start_idx = floor((times - length(F_event2)) / 2) + 1;
F_cr(start_idx:start_idx + length(F_event2) - 1) = F_event2;

% combine signals
F = F_ac + F_cr;

% time vector
t = (0:times - 1) * dt;


% plot the signal

figure(1)
subplot(3,1,1)
plot(t, F)
hold on
xlabel('time (s)')
ylabel('amplitude')
title(sprintf('signal with two modes of resonance \n AC: T = %.1f s, Q = %.1f; CR: T = %.1f s, Q = %.1f', T_ac, Q_ac, T_cr, Q_cr))




%% --- find optimal window length for each pulse ---

% step 1: make two filtered copies of the signal, one for > 10 sec 
% and one for < 10 sec

    Fs = 1/dt;                % sampling frequency
    cutoff_period = 10;       % sec
    cutoff_freq = 1 / cutoff_period;  % hz

    % design butterworth filters (2nd order)
    % one low-pass one high-pass
    % we are in terms of freq now so the low-pass will be for our conduit
    % reservoir mode and high-pass for acoustic mode
    [b_low, a_low] = butter(2, cutoff_freq/(Fs/2), 'low');
    [b_high, a_high] = butter(2, cutoff_freq/(Fs/2), 'high');

    % make filtered copy for each mode
    F_low = filtfilt(b_low, a_low, F); % conduit reservoir mode
    F_high = filtfilt(b_high, a_high, F); % acoustic mode

    % plot
    subplot(3,1,2);
    plot(t, F_low);
    title('lowpass filtered signal (CR mode)');
    xlabel('time (s)'); ylabel('amplitude');

    subplot(3,1,3);
    plot(t, F_high);
    title('fighpass filtered signal (AC mode)');
    xlabel('time (s)'); ylabel('amplitude');
    hold on;


% STEP 2: calculate lta of unfiltered signalto set threshold
    
    lta = mean(F(F >= 0));  % long-term average (positive only)
    threshold = lta;
    
% STEP 3: time window onset for both signals
    % we will use the filtered signals to find the point just before the peak
    % amp of both signals
    % NOTE: we are making an assumption that each mode will be the dominant
    % amplitude within it's bandpassed timeseries. 

    % acoustic mode
    [~, max_idx_ac] = max(F_high);  % peak index in the highpass signal
    % conduit-reservoir mode
    [~, max_idx_cr] = max(F_low);  % peak index in the lowpass signal

% STEP 4: trim timeseries for each mode
   % now that we have the starting index, we are now back to using the unfiltered signal

    % acoustic mode onset/offset

    % find onset
    onset_idx = 1;
    for i = max_idx_ac:-1:1
        if F(i) < threshold
            onset_idx = i;
            break;
        end
    end

    % find offset
    window_duration = T_ac;  % based on original pulse
    window_samples = round(window_duration / dt);
    offset_idx = length(F);
    for i = max_idx_ac:(length(F) - window_samples)
        window = F(i:i + window_samples - 1);
        if mean(window < threshold) > 0.90
            offset_idx = (i + window_samples - 1) + window_multiplier*(T_ac/dt);
            break;
        end
    end

    % trim signal
    F_ac = F(onset_idx:offset_idx);
    t_ac = t(onset_idx:offset_idx);

    figure(2);
    subplot(4,1,1);
    plot(t_ac, F_ac);
    hold on;
    yline(lta, Color="red");
    title('trimmed signal (AC)');
    xlabel('time (s)'); ylabel('amplitude');
    legend("","LTA")
    hold on;


    % conduit reservoir mode onset/offset
    
    % find onset of eac
    onset_idx = 1;
    for i = max_idx_cr:-1:1
        if F(i) < threshold
            onset_idx = i;
            break;
        end
    end
    
    % find offset
    window_duration = T_cr;  % based on original pulse
    window_samples = round(window_duration / dt);
    offset_idx = length(F);
    for i = max_idx_cr:(length(F) - window_samples)
        window = F(i:i + window_samples - 1);
        if mean(window < threshold) > 0.90
            offset_idx = (i + window_samples - 1) + window_multiplier*(T_cr/dt);
            break;
        end
    end

    % trim signal
    F_cr = F(onset_idx:offset_idx);
    t_cr = t(onset_idx:offset_idx);

    subplot(4,1,2);
    plot(t_cr, F_cr);
    hold on;
    yline(lta, Color="red");
    title('trimmed signal (CR)');
    xlabel('time (s)'); ylabel('amplitude');
    legend("","LTA")
    hold on;



%% ---- compute spectrum and estimate T and Q for each mode ----

% STEP 5 (acoustic mode): find spectral peak of interest
    % fft
    [Fv_ac,FTs,Iv,spectrum_ac] = compute_fft(F_ac,dt);
    
    % find peaks
    [peaks, locs, widths, prominence] = findpeaks(spectrum_ac,'SortStr','descend');
    periods = 1 ./ Fv_ac(locs); % convert to period (for sanity)
    
    % filter for peaks with period < 10 sec
    valid_idx = find(periods < 10);
    if isempty(valid_idx)
        error('No spectral peaks in AC range.');
    end
    
    % choose largest peak among those in range
    [~, max_idx] = max(peaks(valid_idx));
    final_idx = valid_idx(max_idx);
    % peak = peaks(final_idx);
    loc = locs(final_idx);

% STEP 6 (acoustic mode): estimate T and Q in spectral domain

    width = widths(final_idx);
    dF = (Fv_ac(2) - Fv_ac(1)) / 2; % delta frequency
    T_pred_ac = 1 / Fv_ac(loc);
    Q_pred_ac = Fv_ac(loc) / (width * dF);

% STEP 5 (conduit reservoir mode): find spectral peak of interest
    % fft
    [Fv_cr,FTs,Iv,spectrum_cr] = compute_fft(F_cr,dt);
    
    % find peaks
    [peaks, locs, widths, prominence] = findpeaks(spectrum_cr,'SortStr','descend');
    periods = 1 ./ Fv_cr(locs); % convert to period (for sanity)
    
    % filter for peaks with period > 10 sec
    valid_idx = find(periods > 10);
    if isempty(valid_idx)
        error('No spectral peaks in CR range.');
    end
    
    % choose largest peak among those in range
    [~, max_idx] = max(peaks(valid_idx));
    final_idx = valid_idx(max_idx);
    % peak = peaks(final_idx);
    loc = locs(final_idx);

% STEP 6 (conduit reservoir mode): estimate T and Q in spectral domain

    width = widths(final_idx);
    dF = (Fv_cr(2) - Fv_cr(1)) / 2; % delta frequency
    T_pred_cr = 1 / Fv_cr(loc);
    Q_pred_cr = Fv_cr(loc) / (width * dF);


        
    subplot(4,1,3)
    findpeaks(spectrum_ac, Fv_ac, 'Annotate', 'extents')
    xlabel('frequency (hz)')
    ylabel('amplitude')
    title(sprintf('spectrum trimmed for acoustic mode \n Pred T: %.1f, Pred Q: %.1f', T_pred_ac, Q_pred_ac))
    set(gca, 'Xscale', 'log')
    hold on;
    
    subplot(4,1,4)
    findpeaks(spectrum_cr, Fv_cr, 'Annotate', 'extents')
    xlabel('frequency (hz)')
    ylabel('amplitude')
    title(sprintf('spectrum trimmed for cr mode \n Pred T: %.1f, Pred Q: %.1f', T_pred_cr, Q_pred_cr))
    set(gca, 'Xscale', 'log')

% --- Display results ---
disp(['Predicted period is: ', num2str(T_pred_ac)])
disp(['Predicted quality factor is: ', num2str(Q_pred_ac)])

disp(['Predicted period is: ', num2str(T_pred_cr)])
disp(['Predicted quality factor is: ', num2str(Q_pred_cr)])
