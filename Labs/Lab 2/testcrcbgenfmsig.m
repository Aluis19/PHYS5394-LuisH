%% Test the frequency-modulated sinusoidal signal
% This script generates an FM sinusoidal signal, plots the sampled signal,
% computes its FFT, and plots the magnitude of its positive-frequency
% components.

%% Signal parameters

% FM modulation index in radians
b = 4;

% Carrier frequency in Hz
f0 = 40;

% Modulation frequency in Hz
f1 = 5;

% Requested matched-filtering SNR
snr = 10;

% Signal duration in seconds
signalDuration = 1.0;

%% Determine the sampling frequency

% The instantaneous frequency of the FM signal is:
%       f_inst(t) = f0 - b*f1*sin(2*pi*f1*t).

% Therefore, the peak frequency deviation is:
%       deltaF = |b|*f1.

deltaF = abs(b)*f1;

%       f0 - (deltaF + f1) <= f <= f0 + (deltaF + f1).

maxFreq = f0 + deltaF + f1;
nyqFreq = 2*maxFreq;
samplFreq = 5*nyqFreq;
samplIntrvl = 1/samplFreq;

%% Generate the time samples

timeVec = 0:samplIntrvl:(signalDuration-samplIntrvl);
nSamples = length(timeVec);

%% Generate the FM sinusoidal signal

fmParams = [b,f0,f1];
sigVec = crcbgenfmsig(timeVec,snr,fmParams);

%% Plot the sampled signal

figure;

plot(timeVec,sigVec);

xlabel('Time (sec)');
ylabel('Signal amplitude');
title('Frequency-Modulated Sinusoidal Signal');
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

plot(posFreq,periodogram, ...
    'Marker','.', ...
    'MarkerSize',10);

xlabel('Frequency (Hz)');
ylabel('|FFT|');
title('Periodogram of the FM Sinusoid');
grid on;

%% Spectrogram (Lab 2): requires Signal Processing Toolbox
% Window duration controls the time/frequency resolution tradeoff.
% Zero padding refines the frequency grid, not the physical resolution.
windowDuration = 0.08;
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
title('Frequency Modulation Spectrogram');
hold on;
expectedFrequency = f0 - b*f1*sin(2*pi*f1*timeVec);

plot(timeVec,expectedFrequency,'w--','LineWidth',1.5);
% Short windows resolve the modulation in time but broaden frequency.
% The dashed curve shows the analytic carrier instantaneous frequency;
% the finite-window spectral peak need not coincide with it exactly.
hold off;
% Complete windows leave uncovered margins at both ends of the record.
fprintf('Spectrogram: window = %d samples (%.3f s), overlap = %d, nfft = %d\n', ...
    windowLength,windowLength/samplFreq,overlapLength,nfft);
fprintf('Time step = %.3f s; frequency grid spacing = %.3f Hz\n', ...
    (windowLength-overlapLength)/samplFreq,samplFreq/nfft);
assert(all(isfinite(pSpec(:))) && all(pSpec(:)>=0), ...
    'Spectrogram power must be finite and nonnegative.');
