function [Fv,FTs,Iv,spectrum] = compute_fft(amplitude,dt)

% this function is used for plotting, not predicting T and Q (we taper here)
% use post_process.m for spectral domain T and Q predictions

%compute fft and spectrum
L = length(amplitude);                                           % length of signal

%taper end with cosine for 8% of L
len = round(.05*L);
%amplitude(L-len+1:L)=cos([1:len] * pi/2 /len).*amplitude(L-len+1:L);

Fs = 1/dt;                                              % Make Up Sampling Frequency & Units (Hz)
Fn = Fs/2;                                              % Nyquist Frequency
FTs = fft(amplitude - mean(amplitude))/L;              % Subtract Mean
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                     % Frequency Vector
Iv = 1:numel(Fv);                                       % Index Vector

% Compute the two-sided spectrum
spectrum_full = abs(FTs);%abs(FTs / L);

% Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
spectrum = spectrum_full(1:L/2+1);
spectrum(2:end-1) = 2*spectrum(2:end-1);

FTs = FTs(1:L/2+1); %single sided


