%% Test the AM-FM sinusoidal signal
% This script generates an amplitude- and frequency-modulated sinusoidal
% signal, plots the sampled signal, computes its FFT, and plots the
% magnitude of its positive-frequency components.

%% Signal parameters

% FM modulation index in radians
b = 4;

% Carrier frequency in Hz
f0 = 40;

% Common AM and FM modulation frequency in Hz
f1 = 5;

% Requested matched-filtering SNR
snr = 10;

% Signal duration in seconds
signalDuration = 1.0;

%% Determine the sampling frequency

%       f_inst(t) = f0 - b*f1*sin(2*pi*f1*t).
%       deltaF = |b|*f1.

deltaF = abs(b)*f1;

%       maxFreq = f0 + deltaF + 2*f1.

maxFreq = f0 + deltaF + 2*f1;
nyqFreq = 2*maxFreq;
samplFreq = 5*nyqFreq;
samplIntrvl = 1/samplFreq;

%% Generate the time samples

timeVec = 0:samplIntrvl:(signalDuration-samplIntrvl);
nSamples = length(timeVec);

%% Generate the AM-FM sinusoidal signal

amfmParams = [b,f0,f1];
sigVec = crcbgenamfmsig(timeVec,snr,amfmParams);

%% Plot the sampled signal

figure;

plot(timeVec,sigVec);

xlabel('Time (sec)');
ylabel('Signal amplitude');
title('AM-FM Sinusoidal Signal');
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
title('Periodogram of the AM-FM Sinusoid');
grid on;