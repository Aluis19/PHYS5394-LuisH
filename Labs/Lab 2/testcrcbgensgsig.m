%% Test the sine-Gaussian signal
% This script generates a sine-Gaussian signal, plots the sampled signal,
% computes its FFT, and plots the magnitude of its positive-frequency
% components.

%% Signal parameters

% Center time of the Gaussian envelope in seconds
t0 = 0.5;

% Standard deviation of the Gaussian envelope in seconds
sigma = 0.1;

% Frequency of the sinusoidal carrier in Hz
f0 = 30;

% Initial phase in radians
phi0 = 0;

% Requested matched-filtering SNR
snr = 10;

% Signal duration in seconds
signalDuration = 1.0;

%% Determine the sampling frequency

%       sigma_f = 1/(2*pi*sigma).

freqStd = 1/(2*pi*sigma);
maxFreq = f0 + 3*freqStd;
nyqFreq = 2*maxFreq;
samplFreq = 5*nyqFreq;

% Sampling interval
samplIntrvl = 1/samplFreq;

%% Generate the time samples

timeVec = 0:samplIntrvl:(signalDuration-samplIntrvl);

% Number of time samples
nSamples = length(timeVec);

%% Generate the sine-Gaussian signal

sgParams = [t0,sigma,f0,phi0];
sigVec = crcbgensgsig(timeVec,snr,sgParams);

%% Plot the sampled signal

figure;

plot(timeVec,sigVec, ...
    'Marker','.', ...
    'MarkerSize',8);

xlabel('Time (sec)');
ylabel('Signal amplitude');
title('Sine-Gaussian Signal');
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
title('Periodogram of the Sine-Gaussian Signal');
grid on;

%% Spectrogram (Lab 2): requires Signal Processing Toolbox
% Window duration controls the time/frequency resolution tradeoff.
% Zero padding refines the frequency grid, not the physical resolution.
windowDuration = 0.15;
windowLength = round(windowDuration*samplFreq);
windowVec = hamming(windowLength,'periodic');
overlapLength = floor(0.90*windowLength);
nfft = 2^nextpow2(4*windowLength);
[stft,fSpec,tSpec,pSpec] = spectrogram(sigVec,windowVec, ...
    overlapLength,nfft,samplFreq);
peakPower = max(pSpec(:));
relativePowerDB = 10*log10(max(pSpec/peakPower,1e-8));
figure;
imagesc(tSpec,fSpec,relativePowerDB);
axis xy;
xlim([timeVec(1),timeVec(end)]);
ylim([0,min(samplFreq/2,1.5*maxFreq)]);
caxis([-50,0]);
colormap parula;
cb = colorbar;
cb.Label.String = 'Relative PSD (dB)';
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title('Sine-Gaussian Spectrogram');
hold on;
expectedFrequency = f0*ones(size(timeVec));

plot(timeVec,expectedFrequency,'w--','LineWidth',1.5);
xline(t0,'w:');
% The Gaussian envelope localizes power near t0; carrier stays at f0.
hold off;
% Complete windows leave uncovered margins at both ends of the record.
fprintf('Spectrogram: window = %d samples (%.3f s), overlap = %d, nfft = %d\n', ...
    windowLength,windowLength/samplFreq,overlapLength,nfft);
fprintf('Time step = %.3f s; frequency grid spacing = %.3f Hz\n', ...
    (windowLength-overlapLength)/samplFreq,samplFreq/nfft);
assert(all(isfinite(pSpec(:))) && all(pSpec(:)>=0), ...
    'Spectrogram power must be finite and nonnegative.');
