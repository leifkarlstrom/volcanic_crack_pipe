function [t_filtered, amplitude_filtered] = filter_timeseries(t, amplitude, fs)
    % Filters the pressure time series between resonant periods of 1 and 10 seconds.
    %
    % Parameters:
    %    t (array): Time array.
    %    amplitude (array): Amplitude array.
    %    fs (scalar): Sampling frequency in Hz.
    %
    % Returns:
    %    t_filtered (array): Filtered time array.
    %    amplitude_filtered (array): Filtered amplitude array.

    % Define frequency bounds
    f_low = 1 / 10; % 0.1 Hz
    f_high = 1 / 1; % 1 Hz
    
    % Design a bandpass filter
    [b, a] = butter(4, [f_low, f_high] / (fs / 2), 'bandpass');
    
    % Apply the filter to the amplitude signal
    amplitude_filtered = filtfilt(b, a, amplitude);
    
    % Return the filtered time and amplitude
    t_filtered = t;
end
