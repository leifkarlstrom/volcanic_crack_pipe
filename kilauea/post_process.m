function [T_condres, Q_condres, T_ac, Q_ac] = post_process(t,amplitude,dt)


% break the resulting out.p_c into two parts, bandpass for 20-50s and then 1-20s, 
% do the trim thing for each, and then calc T and Q. make this all into one function 
% called post_process

% process both modes separately for T and Q estimation

%% conduit-reservoir mode %%

% isolate conduit-reservoir mode [20, 50] seconds
amp_condres = bandpass(amplitude, [1/50 1/20], 1/dt);

% trim timeseries to end at 10 percent max amp of the frequency of interest
[t_condres, amp_condres] = trim_timeseries(t, amp_condres);

%compute spectrum, note that there is a tapering and detrending in here!
[Fv,FTs,Iv,spectrum] = compute_fft(amp_condres,dt);

% the peaks will be listed in order of spectral amp
[peaks, locs, widths, prominence]=findpeaks(spectrum,'SortStr','descend');

%recognize that the frequency spacing from compute_fft is not quite 1/dt.
%In fact it is close (but not exactly) (1/dt)/2 due to Nyquist.
%lets use it exactly
dF=(Fv(2)-Fv(1))/2;

%predict T and Q based on spectrum
% listed in order of spectral amp
% Only keep the first one
T_condres = 1./Fv(locs);
T_condres = T_condres(1);
Q_condres = Fv(locs)./(widths*dF);
Q_condres = Q_condres(1);

% %plot it
% figure()
% subplot(4,1,1)
% subtitle('conduit reservoir mode')
% plot(t_condres,amp_condres)
% subplot(4,1,2)
% findpeaks(spectrum,Fv,'Annotate','extents')
% %semilogx(Fv, abs(spectrum),Fv, detrend(abs(spectrum)))
% set(gca,'Xscale','log')
% hold on

%% acoustic mode %%

% isolate acoustic mode [20, 50] seconds
amp_ac = bandpass(amplitude, [1/20 1], 1/dt);

% trim timeseries to end at 10 percent max amp of the frequency of interest
[t_ac, amp_ac] = trim_timeseries(t, amp_ac);

%compute spectrum, note that there is a tapering and detrending in here!
[Fv,FTs,Iv,spectrum] = compute_fft(amp_ac,dt);

% the peaks will be listed in order of spectral amp
[peaks, locs, widths, prominence]=findpeaks(spectrum,'SortStr','descend');

%recognize that the frequency spacing from compute_fft is not quite 1/dt.
%In fact it is close (but not exactly) (1/dt)/2 due to Nyquist.
%lets use it exactly
dF=(Fv(2)-Fv(1))/2;

%predict T and Q based on spectrum
% listed in order of spectral amp
T_ac = 1./Fv(locs);
T_ac = T_ac(1);
Q_ac = Fv(locs)./(widths*dF);
Q_ac = Q_ac(1);

% %plot it
% subplot(4,1,3)
% subtitle('acoustic mode')
% plot(t_ac,amp_ac)
% subplot(4,1,4)
% findpeaks(spectrum,Fv,'Annotate','extents')
% %semilogx(Fv, abs(spectrum),Fv, detrend(abs(spectrum)))
% set(gca,'Xscale','log')


% disp('---------------------------------------------------------------------------')
% disp(['predicted period for conduit-reservoir mode is: ' num2str(T_condres)]) 
% disp(['predicted quality factor for conduit-reservoir mode is: ' num2str(Q_condres)]) 
% disp('---------------------------------------------------------------------------')
% disp(['predicted period for acoustic mode is: ' num2str(T_ac)]) 
% disp(['predicted quality factor for acoustic mode is: ' num2str(Q_ac)])
% disp('---------------------------------------------------------------------------')
% 
%plot it
subplot(2,1,1)
plot(t_ac,amp_ac)
xlabel('Time (s)', 'FontSize', 16)
ylabel('Amplitude (m)', 'FontSize', 16)
title('Signal', 'FontSize', 18)
subplot(2,1,2)
findpeaks(spectrum,Fv,'Annotate','extents')
xlabel('Time (s)', 'FontSize', 16)
ylabel('Amplitude (m)', 'FontSize', 16)
%semilogx(Fv, abs(spectrum),Fv, detrend(abs(spectrum)))
set(gca,'Xscale','log')
hold on
% xline(omega/(2*pi))



end

%% functions %%

function [t_trimmed, y_trimmed] = trim_timeseries(t, y)

    % Find the maximum amplitude
    max_amp = max(y);
    
    % Determine the threshold (10% of max amplitude)
    threshold = 0.1 * max_amp;
    
    % Find the last index where y is above the threshold
    last_idx = find(y >= threshold, 1, 'last');
    
    % Trim the time and amplitude arrays
    t_trimmed = t(1:last_idx);
    y_trimmed = y(1:last_idx);
end


