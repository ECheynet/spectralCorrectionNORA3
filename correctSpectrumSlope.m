function signal_corr = correctSpectrumSlope(signal, fs, fmin,fmax)
%FILTER_SLOPE_5_3_THRESHOLD Filters a time series to enforce a -5/3 spectral slope 
% above a specific frequency threshold, while preserving the low-frequency range and phase.
%
%   filtered_signal = FILTER_SLOPE_5_3_THRESHOLD(signal, fs, freq_thresh)
%
%   Inputs:
%       - signal: Original time series (vector)
%       - fs: Sampling frequency (scalar)
%       - freq_thresh: Frequency threshold above which the -5/3 slope is applied (Hz)
%
%   Outputs:
%       - filtered_signal: Time series filtered to have a -5/3 spectral slope above freq_thresh
%
%   Example:
%       fs = 1.0; % Sampling frequency
%       signal = randn(1, 1024); % White noise
%       freq_thresh = 0.1; % Threshold frequency in Hz
%       filtered_signal = correctSpectrumSlope(signal, fs, freq_thresh);

    % Ensure signal is a row vector
    if iscolumn(signal)
        signal = signal';
    end

    N = length(signal);  % Number of points

    % Step 1: FFT of the signal
    fft_signal = fft(signal);
    freqs = (0:N-1) * (fs / N);  % Frequency vector

    % Adjust frequencies for two-sided FFT representation
    freqs_shifted = fftshift(freqs - fs * (freqs >= fs/2)); % Shift negative frequencies

    % Step 2: Separate amplitude and phase
    amplitude = abs(fft_signal); % Amplitude spectrum
    phase = angle(fft_signal);   % Phase spectrum

 
    idx = abs(freqs_shifted) > fmin & abs(freqs_shifted) < fmax;
    if nnz(idx) > 1
        [S,f] = pwelch(detrend(signal),[],[],[],fs);

        M = 60;
        newF0 = logspace(log10(f(2)*0.8), log10(f(end)*1.1), M);
        [S1,f1] = binAveraging(S,f,'newX',newF0);

        indF = f1>fmin & f1<fmax;
        x = log(f1(indF));
        y = log(S1(indF));

        % Robust linear regression: y = b0 + b1*x
        [b,stats] = robustfit(x,y,'bisquare');   % default is good; try 'huber' too
        current_slope = -b(2);                   % because S ~ f^{-a2} ⇒ log S = log a1 - a2 log f

        % (Optional) 95% CI:
        ci = -[b(2)+1.96*stats.se(2), b(2)-1.96*stats.se(2)];
    else
        warning('no frequency found in the selected interval')
        current_slope = 0;
    end



    
    
    idx = abs(freqs_shifted) > fmin;
    % Step 4: Create the amplitude filter
    filter_response = ones(size(freqs_shifted)); 
    % Flatten the slope first by compensating for the current slope
    filter_response(idx) = (abs(freqs_shifted(idx)) / fmin).^(current_slope / 2);
    
    % Apply -5/3 slope after flattening
    slope = -5/3;
    filter_response(idx) = filter_response(idx) .* (abs(freqs_shifted(idx)) / fmin).^(slope / 2);
    
    % Ensure DC component (freq = 0) is not modified
    filter_response(freqs_shifted == 0) = 1.0; 

    % Ensure symmetry in the filter
    filter_response = ifftshift(filter_response); % Undo frequency shift for symmetry

    % Step 5: Apply the filter to the amplitude spectrum symmetrically
    filtered_amplitude = amplitude .* filter_response;

    % Step 6: Reconstruct the filtered FFT signal using original phase
    filtered_fft_signal = filtered_amplitude .* exp(1i * phase);

    % Step 7: Inverse FFT to return to time domain
    signal_corr = real(ifft(filtered_fft_signal));

    signal_corr(signal_corr<0) = 0;
    % Optional Diagnostic Plots
    % Uncomment these for debugging and visualization

    % % Plot the filter response
    % figure;
    % plot(freqs_shifted, filter_response, 'k', 'LineWidth', 1.5);
    % title('Frequency Filter Response');
    % xlabel('Frequency (Hz)');
    % ylabel('Filter Response');
    % grid on;
    % xlim([-2*freq_thresh, 2*freq_thresh]);
    % 
    % % Plot the PSD
    % [so, fo] = pwelch(signal, [], [], [], fs);
    % [sf, ff] = pwelch(filtered_signal, [], [], [], fs);
    % figure;
    % loglog(fo, so, 'b', 'DisplayName', 'Original PSD');
    % hold on;
    % loglog(ff, sf, 'r', 'DisplayName', 'Filtered PSD (-5/3 Above Threshold)');
    % legend('show');
    % title('Power Spectral Density');
    % xlabel('Frequency (Hz)');
    % ylabel('Power');
    % grid on;
    % 
    % % Plot the amplitude spectra
    % figure;
    % plot(freqs_shifted, amplitude, 'b', 'DisplayName', 'Original Amplitude Spectrum');
    % hold on;
    % plot(freqs_shifted, filtered_amplitude, 'r', 'DisplayName', 'Filtered Amplitude Spectrum');
    % legend('show');
    % title('Amplitude Spectrum Comparison');
    % xlabel('Frequency (Hz)');
    % ylabel('Amplitude');
    % grid on;
    % 
    % % Plot the time series
    % figure;
    % plot((0:N-1)/fs, signal, 'b', 'DisplayName', 'Original Signal');
    % hold on;
    % plot((0:N-1)/fs, filtered_signal, 'r', 'DisplayName', 'Filtered Signal (-5/3 Above Threshold)');
    % legend('show');
    % title('Time Series');
    % xlabel('Time (s)');
    % ylabel('Amplitude');
    % grid on;
end
