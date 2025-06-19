close all
clear
clc

%build synthetic pulse
window_multiplier = 10;
period = 2.0; %sec
omega = 2*pi/period; %angular frequency
Q = 10; %desired Quality factor
alpha = omega/(2*Q); %damping rate (1/s)

%make a signal
dt = 0.01;
%total time is where amplitude drops to 1% of original
Ttotal = -log(0.001)/alpha;  
event_t = 0:dt:Ttotal;
F_event = sin(omega.*event_t).*exp(-alpha*event_t);

% pad zeros to 2.5 minutes (150 seconds)
total_duration = 150; % seconds
times = round(total_duration / dt);
F = zeros(1, times);

% Insert event into the middle of the zero-padded array
start_idx = floor((times - length(F_event)) / 2) + 1;
F(start_idx:start_idx + length(F_event) - 1) = F_event;

% time vector for padded signal
t = (0:times - 1) * dt;

% add Gaussian white noise
F = awgn(F,25);


% ######### Amplitude-Based Detection (positive amplitudes only) #########
lta = mean(F(F >= 0));  % lta (positive only)
threshold = lta;

[~, max_idx] = max(F);  % peak index

% find onset (positive only)
onset_idx = 1;
for i = max_idx:-1:1
    if F(i) < threshold
        onset_idx = i;
        break;
    end
end

% ######### Step 3: Define how long the signal must stay below the threshold #########
window_duration = period;  % seconds
window_samples = round(window_duration / dt);

offset_idx = length(F);
for i = max_idx:(length(F) - window_samples)
    window = F(i:i + window_samples - 1);
    % Check if most of the window is below threshold and index the end of
    % that window if so
    if mean(window < threshold) > 0.90
        offset_idx = (i + window_samples - 1) + window_multiplier*(period/dt);
        break;
    end
end

% Convert to time
onset_time = t(onset_idx);
peak_time = t(max_idx);
offset_time = t(offset_idx);

% Step 3: Trim trace
event_data = F(onset_idx:offset_idx);
event_times = t(onset_idx:offset_idx);


%compute spectrum
[Fv,FTs,Iv,spectrum] = compute_fft(event_data,dt);

[peaks, locs, widths, prominence]=findpeaks(spectrum,'SortStr','descend');
peak = peaks(1);
loc = locs(1);
width = widths(1);

%recognize that the frequency spacing from compute_fft is not quite 1/dt.
%In fact it is close (but not exactly) (1/dt)/2 due to Nyquist.
%lets use it exactly

dF=(Fv(2)-Fv(1))/2;

% %predict T and Q based on spectrum
T_pred = 1/Fv(loc);
Q_pred = Fv(loc)/(width*dF);


%plot it
subplot(3,1,1)
plot(t, F)
hold on
yline(lta)
xlabel('Time (s)', 'FontSize', 16)
ylabel('Amplitude (m)', 'FontSize', 16)
title(sprintf('Synthetic Signal with T = %.3f sec, Q = %.3f',period,Q),'FontSize', 15)
subplot(3,1,2)
plot(event_times,event_data)
hold on
yline(lta)
xlabel('Time (s)', 'FontSize', 16)
ylabel('Amplitude (m)', 'FontSize', 16)
title(sprintf('Predicted T = %.3f sec, Q = %.3f',T_pred,Q_pred),'FontSize', 15)
subplot(3,1,3)
findpeaks(spectrum,Fv,'Annotate','extents')
xlabel('Time (s)', 'FontSize', 16)
ylabel('Amplitude (m)', 'FontSize', 16)
title(sprintf('Window extension: + %.0f * %.1f sec',window_multiplier, period),'FontSize', 15)
%semilogx(Fv, abs(spectrum),Fv, detrend(abs(spectrum)))
set(gca,'Xscale','log')
hold on
% xline(omega/(2*pi))



disp(['Predicted period is: ', num2str(T_pred)])
disp(['Predicted quality factor is: ', num2str(Q_pred)])
disp(['Signal length: ', num2str(length(event_times)*dt)])