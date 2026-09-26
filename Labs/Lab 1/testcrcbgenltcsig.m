%% Test the linear transient chirp signal
% This script generates a linear transient chirp, plots the sampled signal,
% computes its FFT, and plots the magnitude of its positive-frequency
% components.

%% Signal parameters

% Starting time of the transient in seconds
ta = 0.2;

% Initial frequency of the transient chirp in Hz
f0 = 10;

% Quadratic phase coefficient in Hz/sec
f1 = 20;

% Initial phase in radians
phi0 = 0;

% Duration of the transient in seconds
L = 0.5;

% Requested matched-filtering SNR
snr = 10;

% Total observation duration in seconds
signalDuration = 1.0;

%% Check that the transient is inside the observation interval

% The transient must end before or at the end of the complete observation.
if (ta + L) > signalDuration
    error('The transient extends beyond the observation interval.');
end

%% Determine the sampling frequency

% Inside the transient interval:
%       tau = t - ta.
% The phase in cycles is:
%       Phi(tau) = f0*tau + f1*tau^2.
% Therefore, the instantaneous frequency is:
%       f(tau) = f0 + 2*f1*tau.

initialFreq = f0;
finalFreq = f0 + 2*f1*L;

maxFreq = max(abs([initialFreq,finalFreq]));
nyqFreq = 2*maxFreq;
samplFreq = 5*nyqFreq;
samplIntrvl = 1/samplFreq;

%% Generate the time samples

timeVec = 0:samplIntrvl:(signalDuration-samplIntrvl);
nSamples = length(timeVec);

%% Generate the linear transient chirp

ltcParams = [ta,f0,f1,phi0,L];
sigVec = crcbgenltcsig(timeVec,snr,ltcParams);

%% Plot the sampled signal

figure;

plot(timeVec,sigVec, ...
    'Marker','.', ...
    'MarkerSize',8);

xlabel('Time (sec)');
ylabel('Signal amplitude');
title('Linear Transient Chirp Signal');
grid on;

% Mark the beginning and end of the transient.
xline(ta,'--r','Start');
xline(ta+L,'--r','End');

%% Compute the periodogram

% DFT index corresponding to the highest nonnegative frequency
kNyq = floor(nSamples/2)+1;

% Exact frequency resolution of the DFT
freqRes = samplFreq/nSamples;

% Nonnegative Fourier frequencies
posFreq = (0:(kNyq-1))*freqRes;

% Compute the FFT of the complete signal
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
title('Periodogram of the Linear Transient Chirp');
grid on;