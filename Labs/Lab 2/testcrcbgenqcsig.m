%% Test the quadratic chirp signal
% This script generates a quadratic chirp, plots the sampled signal,
% computes its FFT, and plots the magnitude of its positive-frequency
% components.

%% Signal parameters

% Quadratic phase coefficients
a1 = 10;
a2 = 3;
a3 = 3;

% Requested matched-filtering SNR
snr = 10;

% Signal duration in seconds
signalDuration = 1.0;

%% Determine the sampling frequency

% The instantaneous frequency of the quadratic chirp is:
%
%       f(t) = a1 + 2*a2*t + 3*a3*t^2.
%
% For positive coefficients, the maximum instantaneous frequency occurs
% at the end of the signal.

maxFreq = a1 ...
        + 2*a2*signalDuration ...
        + 3*a3*signalDuration^2;

% Nyquist requires the sampling frequency to be at least twice the
% maximum frequency contained in the signal.
nyqFreq = 2*maxFreq;

% Use a sampling frequency five times larger than the Nyquist requirement.
samplFreq = 5*nyqFreq;

% Sampling interval
samplIntrvl = 1/samplFreq;

%% Generate the time samples

% Sample the interval from 0 to signalDuration.
timeVec = 0:samplIntrvl:signalDuration;

% Number of time samples
nSamples = length(timeVec);

%% Generate the quadratic chirp

% The parameter vector contains the three quadratic phase coefficients:

qcCoefs = [a1,a2,a3];

sigVec = crcbgenqcsig(timeVec,snr,qcCoefs);

%% Plot the sampled signal

figure;

plot(timeVec,sigVec, ...
    'Marker','.', ...
    'MarkerSize',12);

xlabel('Time (sec)');
ylabel('Signal amplitude');
title('Quadratic Chirp Signal');
grid on;

%% Compute the periodogram

% DFT index corresponding to the highest nonnegative frequency
kNyq = floor(nSamples/2)+1;

% Exact frequency resolution of the DFT
freqRes = samplFreq/nSamples;

% Nonnegative Fourier frequencies
posFreq = (0:(kNyq-1))*freqRes;

% Compute the FFT of the signal
fftSig = fft(sigVec);

% Retain only the nonnegative-frequency components
fftSig = fftSig(1:kNyq);

% Magnitude of the positive-frequency FFT
periodogram = abs(fftSig);

%% Plot the periodogram

figure;

plot(posFreq,periodogram);

xlabel('Frequency (Hz)');
ylabel('|FFT|');
title('Periodogram of the Quadratic Chirp');
grid on;

%% Spectrogram: time-frequency representation (Lab 2)
% Requires Signal Processing Toolbox. A longer window separates nearby
% frequencies better, but averages the chirp over a longer time interval.
windowDuration = 0.40;                 % seconds
windowLength = round(windowDuration*samplFreq);
windowVec = hamming(windowLength,'periodic');
overlapLength = floor(0.90*windowLength);
nfft = 2^nextpow2(4*windowLength);       % finer grid, not finer resolution

[stft,fSpec,tSpec,pSpec] = spectrogram(sigVec,windowVec, ...
    overlapLength,nfft,samplFreq);

% Display power relative to the strongest time-frequency bin.
peakPower = max(pSpec(:));
relativePowerDB = 10*log10(max(pSpec/peakPower,1e-8));

figure;
imagesc(tSpec,fSpec,relativePowerDB);
axis xy;
ylim([0,1.5*maxFreq]);
xlim([timeVec(1),timeVec(end)]);
caxis([-50,0]);
colormap parula;
cb = colorbar;
cb.Label.String = 'Relative PSD (dB)';
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title('Quadratic Chirp Spectrogram');
hold on;
expectedFrequency = a1 + 2*a2*timeVec + 3*a3*timeVec.^2;
plot(timeVec,expectedFrequency,'w--','LineWidth',1.5);
hold off;

% The dashed curve is the analytic instantaneous frequency. STFT columns
% are centered on complete windows, so they do not extend to both endpoints.
fprintf('Spectrogram: window = %d samples (%.3f s), overlap = %d, nfft = %d\n', ...
    windowLength,windowLength/samplFreq,overlapLength,nfft);
fprintf('Time step = %.3f s; frequency grid spacing = %.3f Hz\n', ...
    (windowLength-overlapLength)/samplFreq,samplFreq/nfft);
assert(all(isfinite(pSpec(:))) && all(pSpec(:)>=0), ...
    'Spectrogram power must be finite and nonnegative.');
