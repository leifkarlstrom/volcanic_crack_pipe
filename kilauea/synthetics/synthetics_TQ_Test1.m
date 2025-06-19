close all
clear

%build synthetic pulse

period = 2.5; %sec
omega = 2*pi/period; %angular frequency
Q = 12; %desired Quality factor

alpha = omega/(2*Q); %damping rate (1/s)

%make a signal
dt = 0.01;
%total time is where amplitude drops to 1% of original
Ttotal = -log(0.01)/alpha;  

t = 0:dt:Ttotal;

F = sin(omega.*t).*exp(-alpha*t);

%compute spectrum
[Fv,FTs,Iv,spectrum] = compute_fft(F,dt);

[peaks, locs, widths, prominence]=findpeaks(spectrum);


%recognize that the frequency spacing from compute_fft is not quite 1/dt.
%In fact it is close (but not exactly) (1/dt)/2 due to Nyquist.
%lets use it exactly

dF=(Fv(2)-Fv(1))/2;

%predict T and Q based on spectrum
T_pred = 1/Fv(locs);
Q_pred = Fv(locs)/(widths*dF);

%plot it
subplot(2,1,1)
plot(t,F)
xlabel('Time (s)', 'FontSize', 16)
ylabel('Amplitude (m)', 'FontSize', 16)
title(sprintf('Synthetic Signal with T = %.1f sec, Q = %.1f',period,Q),'FontSize', 18)
subplot(2,1,2)
findpeaks(spectrum,Fv,'Annotate','extents')
xlabel('Time (s)', 'FontSize', 16)
ylabel('Amplitude (m)', 'FontSize', 16)
%semilogx(Fv, abs(spectrum),Fv, detrend(abs(spectrum)))
set(gca,'Xscale','log')
hold on
% xline(omega/(2*pi))



disp(['Predicted period is: ', num2str(T_pred)])
disp(['Predicted quality factor is: ', num2str(Q_pred)])
% 
% 
% % noisy
% 
% F_n = awgn(F,15);
% 
% %compute spectrum
% [Fv,FTs,Iv,spectrum] = compute_fft(F_n,dt);
% 
% [peaks, locs, widths, prominence]=findpeaks(spectrum,'SortStr','descend');
% 
% %recognize that the frequency spacing from compute_fft is not quite 1/dt.
% %In fact it is close (but not exactly) (1/dt)/2 due to Nyquist.
% %lets use it exactly
% 
% dF=(Fv(2)-Fv(1))/2;
% 
% %predict T and Q based on spectrum
% T_pred_noise = 1./Fv(locs);
% T_pred_noise = T_pred_noise(1);
% Q_pred_noise = Fv(locs)./(widths*dF);
% Q_pred_noise = Q_pred_noise(1);
% 
% disp(['Predicted period is: ', num2str(T_pred_noise)])
% disp(['Predicted quality factor is: ', num2str(Q_pred_noise)])
% 
% subplot(4,1,3)
% plot(t,F_n)
% xlabel('Time (s)', 'FontSize', 16)
% ylabel('Amplitude (m)', 'FontSize', 16)
% title(sprintf('Synthetic Signal (noisy) with T = %.1f sec, Q = %.1f',period,Q),'FontSize', 18)
% subplot(4,1,4)
% findpeaks(spectrum,Fv,'Annotate','extents')
% xlabel('Time (s)', 'FontSize', 16)
% ylabel('Amplitude (m)', 'FontSize', 16)
% %semilogx(Fv, abs(spectrum),Fv, detrend(abs(spectrum)))
% set(gca,'Xscale','log')