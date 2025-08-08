function signal_corr = correctSpectrumSlope(signal, fs, fmin,fmax)
% correctSpectrumSlope  Adjust PSD slope within a specified frequency band.
%   signal_corr = correctSpectrumSlope(signal, fs, fmin, fmax) modifies the
%   amplitude spectrum of a uniformly sampled time series so that its power
%   spectral density (PSD) follows a −5/3 power law between the given lower
%   and upper frequency limits (Hz). Phases are preserved, the low-frequency
%   range below fmin is unchanged, and the correction tapers out smoothly
%   above fmax to avoid artifacts.
%
%   This is intended for mesoscale wind speed time series (hourly resolution)
%   to compensate for model underestimation of variance and extremes by
%   restoring the expected slope in a selected frequency band. The method is
%   inspired by Bastine et al. (2018) and is a home-made implementation of 
%   the method.
%
%   Inputs
%   ------
%   signal : [N×1] or [1×N] double/single
%       Time series (uniform sampling, NaNs already handled).
%   fs     : scalar double
%       Sampling frequency in Hz.
%   fmin   : scalar double
%       Lower bound of frequency band to be slope-corrected (Hz).
%   fmax   : scalar double
%       Upper bound of frequency band to be slope-corrected (Hz).
%
%   Output
%   ------
%   signal_corr : same size and class as signal
%       Corrected time series with PSD ~ f^(−5/3) on [fmin, fmax].
%
%   Method (overview)
%   -----------------
%   1) FFT of the input signal.
%   2) Estimate current spectral slope on [fmin,fmax] via log–log fit.
%   3) Build a gain function to adjust that slope to −5/3 over the band,
%      with smooth tapers at fmin and fmax.
%   4) Apply the gain to the amplitude spectrum (phase preserved) and
%      inverse FFT to obtain the corrected series.
%
%   Notes
%   -----
%   • The function assumes evenly spaced samples with fs > 0.
%   • 0 < fmin < fmax < fs/2 must hold.
%   • Detrend before calling if needed to avoid bias at low frequencies.
%
%   Example
%   -------
%   % Synthetic hourly series (10 years)
%   fs = 1/3600;                           % Hz (hourly data)
%   t  = (0:(10*365*24-1))'/fs;            % time vector in seconds
%   x  = randn(size(t));                   % synthetic signal
%   fmin = 1/(12*3600); fmax = 1/(6*3600); % 12 h to 6 h band
%   xcorr = correctSpectrumSlope(x, fs, fmin, fmax);
%
%   Reference
%   ---------
%   Bastine, D., Larsén, X., Witha, B., Dörenkämper, M., & Gottschall, J. (2018).
%   Extreme winds in the new European wind atlas. J. Phys.: Conf. Ser., 1102, 012006.
%
%  See also: PWELCH, DETREND, INPAINT_NANS, BINAVERAGING
%
% Author: E. Cheynet - UiB -  last modified: 08/08/2025

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
