%% Test the sinusoidal signal
% This script generates a sinusoidal signal, plots the sampled signal,
% computes its FFT, and plots the magnitude of its positive-frequency
% components.

%% Signal parameters

% Sinusoidal frequency in Hz
f0 = 10;

% Initial phase in radians
phi0 = 0;

% Requested matched-filtering SNR
snr = 10;

% Signal duration in seconds
signalDuration = 1.0;

%% Determine the sampling frequency

maxFreq = f0;
nyqFreq = 2*maxFreq;

% Use a sampling frequency five times larger than the Nyquist requirement.
samplFreq = 5*nyqFreq;

% Sampling interval
samplIntrvl = 1/samplFreq;

%% Generate the time samples

timeVec = 0:samplIntrvl:(signalDuration-samplIntrvl);

% Number of time samples
nSamples = length(timeVec);

%% Generate the sinusoidal signal

sinParams = [f0,phi0];

sigVec = crcbgensinsig(timeVec,snr,sinParams);

%% Plot the sampled signal

figure;

plot(timeVec,sigVec, ...
    'Marker','.', ...
    'MarkerSize',12);

xlabel('Time (sec)');
ylabel('Signal amplitude');
title('Sinusoidal Signal');
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
    'MarkerSize',12);

xlabel('Frequency (Hz)');
ylabel('|FFT|');
title('Periodogram of the Sinusoidal Signal');
grid on;