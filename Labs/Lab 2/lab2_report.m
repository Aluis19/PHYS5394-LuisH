%% Lab 2: Sampling, FIR Filtering, and Time-Frequency Analysis
% Luis A Hernandez
% PHYS 5394 - Statistical Methods

%% 1. Sampling parameters
% For phase phi(t) in radians, instantaneous frequency is phi'(t)/(2*pi).
% The lab comparison uses Fs = 2*fmax, 4*fmax, and 20*fmax for a 1 s record.
% These rates are based on maximum instantaneous frequency. For modulated
% and finite-duration signals, this need not bound the full spectrum.
lab2_signalParameters;

%% 1.1 Sinusoid
% The signal has a constant frequency of 10 Hz. At Fs=20 Hz, a zero-phase
% sine is sampled at its zero crossings and disappears. Sampling exactly
% at the Nyquist boundary can therefore lose phase-dependent information.
reportNyquist(1);

%% 1.2 Linear chirp
% The instantaneous frequency rises from 10 to 30 Hz. Higher sampling rates
% describe the waveform more faithfully. The spectrum spans a frequency
% interval because the signal sweeps through multiple frequencies.
reportNyquist(2);

%% 1.3 Quadratic chirp
% The frequency is 10+6*t+9*t^2 Hz and reaches 25 Hz at t=1 s. The sweep
% accelerates, so the signal spends unequal amounts of time at each frequency.
reportNyquist(3);

%% 1.4 Transient linear chirp
% During 0.2 <= t <= 0.7 s, frequency rises from 10 to 30 Hz. The finite
% gate reduces record-averaged power and introduces spectral spreading.
reportNyquist(4);

%% 1.5 Sine-Gaussian
% A Gaussian envelope localizes a 30 Hz carrier around t=0.5 s. At Fs=60 Hz,
% the zero-phase carrier samples vanish. The envelope broadens the spectrum,
% so twice the carrier frequency is not a strict bandwidth-based Nyquist rate.
reportNyquist(5);

%% 1.6 Amplitude modulation
% cos(2*pi*5*t)*sin(2*pi*40*t) produces equal-amplitude sidebands at 35 and
% 45 Hz. The carrier-based rate of 80 Hz misses the upper sideband and gives
% zero samples for this phase. The actual two-tone Nyquist boundary is 90 Hz.
reportNyquist(6);

%% 1.7 Frequency modulation
% The instantaneous frequency is 40-20*sin(2*pi*5*t) Hz, ranging from 20 to
% 60 Hz. FM has sidebands beyond that range; twice the maximum instantaneous
% frequency is not an exact alias-free sampling guarantee.
reportNyquist(7);

%% 1.8 Combined AM-FM
% The envelope and instantaneous frequency both change periodically. The
% amplitude modulation shifts the FM sidebands by +/-5 Hz. Higher sampling
% rates improve the representation of the significant spectral components.
reportNyquist(8);

%% Sampling and PSD interpretation
% All records last 1 s, so their DFT spacing remains approximately 1 Hz.
% Increasing Fs increases sample density and Nyquist frequency, not frequency
% resolution at fixed duration. The one-sided PSD uses |FFT|^2/(Fs*N), with
% doubled interior positive-frequency bins. Summing PSD*df equals mean(x.^2).
% This Parseval check validates normalization, not absence of aliasing.

%% 2. Three-component signal and FIR filters
% The input contains 100, 200, and 300 Hz sinusoids with amplitudes 10, 5,
% and 2.5, and phases 0, pi/6, and pi/4. With N=2048 and Fs=1024 Hz,
% duration is 2 s, df=0.5 Hz, and Nyquist frequency is 512 Hz.
% Component mean-square powers are 50, 12.5, and 3.125: total 65.625.
% Peak PSD values are 100, 25, and 6.25 because each occupies a 0.5 Hz bin.
%
% Order-100 Hamming-window FIR filters use cutoffs of 150 Hz (lowpass),
% 150-250 Hz (bandpass), and 250 Hz (highpass). These lie between the tones.
% The 101 symmetric coefficients give a delay of 50 samples, or 0.048828 s.
lab2_filtering;

%% Filtering results and physical interpretation
% Lowpass retains the 100 Hz component with amplitude about 10.005.
% Bandpass retains the 200 Hz component with amplitude 5.
% Highpass retains the 300 Hz component with amplitude about 2.505.
% Unwanted tones are attenuated by roughly 58-67 dB. Small passband gain
% deviations are expected for these finite-order filters.
% Full-record output PSDs include startup transients, producing broad leakage.
% After the first 100 samples, predicted and computed waveforms agree to
% approximately 1e-12. Numerical checks verify selectivity and PSD scaling.

%% 3.1 Sinusoid spectrogram
% A stationary horizontal ridge is centered at 10 Hz. Its finite width
% comes from the analysis window; the original signal is a single tone.
reportSpectrogram('testcrcbgensinsig');

%% 3.2 Linear chirp spectrogram
% The rising ridge follows f(t)=10+20*t Hz and directly displays the sweep.
reportSpectrogram('testcrcbgenlcsig');

%% 3.3 Quadratic chirp spectrogram
% The curved ridge follows f(t)=10+6*t+9*t^2 Hz, showing an accelerating sweep.
reportSpectrogram('testcrcbgenqcsig');

%% 3.4 Transient chirp spectrogram
% Power is concentrated during the gate interval. Windows overlapping its
% boundaries smear timing, and abrupt gating introduces broadband content.
reportSpectrogram('testcrcbgenltcsig');

%% 3.5 Sine-Gaussian spectrogram
% The time-frequency patch is centered near 0.5 s and 30 Hz. Localization
% comes from the Gaussian pulse, with extra broadening from the STFT window.
reportSpectrogram('testcrcbgensgsig');

%% 3.6 AM spectrogram
% Horizontal ridges at 35 and 45 Hz verify the suppressed-carrier sidebands.
% The 0.4 s window emphasizes frequency separation and averages the envelope.
reportSpectrogram('testcrcbgenamsig');

%% 3.7 FM spectrogram
% A periodic ridge follows the 20-60 Hz instantaneous-frequency variation.
% The 0.08 s window reveals rapid changes but produces a broad frequency band.
reportSpectrogram('testcrcbgenfmsig');

%% 3.8 AM-FM spectrogram
% Frequency varies periodically while the envelope changes power. Short
% windows reveal time variation, although frequency detail is less resolved.
reportSpectrogram('testcrcbgenamfmsig');

%% Conclusions
% Sampling at the instantaneous-frequency-based boundary can distort or even
% erase a signal. Bandwidth and phase must also be considered. FIR filters
% isolate the three tones with predictable delay and strong rejection.
% Spectrograms distinguish stationary tones, sweeps, localized pulses, and
% modulation. Window duration controls the time-frequency resolution tradeoff;
% overlap increases time-grid density, while zero padding refines the frequency
% grid without improving physical frequency resolution. Color is relative to
% each plot's maximum and cannot be used to compare absolute power across plots.

%% Local report helpers
function reportNyquist(selectedSignal)
    sourcePath = which('lab2_nyquistComparison');
    assert(~isempty(sourcePath),'Place lab2_nyquistComparison.m on the path.');
    sourceText = fileread(sourcePath);
    selector = 'signalIdx\s*=\s*\d+\s*;';
    assert(numel(regexp(sourceText,selector))==1, ...
        'Expected one signal selector in lab2_nyquistComparison.m.');
    sourceText = regexprep(sourceText,selector, ...
        sprintf('signalIdx = %d;',selectedSignal));
    eval(sourceText);
end

function reportSpectrogram(testName)
    sourcePath = which(testName);
    assert(~isempty(sourcePath),'Missing test script: %s',testName);
    sourceText = fileread(sourcePath);
    firstPlot = strfind(sourceText,'%% Plot the sampled signal');
    specSection = strfind(sourceText,'%% Spectrogram');
    assert(~isempty(firstPlot) && ~isempty(specSection), ...
        'Use the updated Lab 2 test script for %s.',testName);
    % Generate the original signal, then execute only its spectrogram block.
    % Local function workspace prevents variables from leaking between tests.
    eval(sourceText(1:firstPlot(1)-1));
    eval(sourceText(specSection(1):end));
end
