%% Test the amplitude-modulated sinusoidal signal
% This script generates an AM sinusoidal signal, plots the sampled signal,
% computes its FFT, and plots the magnitude of its positive-frequency
% components.

%% Signal parameters

% Carrier frequency in Hz
f0 = 40;

% Amplitude-modulation frequency in Hz
f1 = 5;

% Initial carrier phase in radians
phi0 = 0;

% Requested matched-filtering SNR
snr = 10;

% Signal duration in seconds
signalDuration = 1.0;

%% Determine the sampling frequency

% The AM signal is:
%       s(t) = cos(2*pi*f1*t)*sin(2*pi*f0*t + phi0).

maxFreq = f0 + f1;
nyqFreq = 2*maxFreq;
samplFreq = 5*nyqFreq;
samplIntrvl = 1/samplFreq;

%% Generate the time samples

%       t_n = n*Delta,  n = 0,1,...,N-1.

timeVec = 0:samplIntrvl:(signalDuration-samplIntrvl);
nSamples = length(timeVec);

%% Generate the AM sinusoidal signal

amParams = [f0,f1,phi0];
sigVec = crcbgenamsig(timeVec,snr,amParams);

%% Plot the sampled signal

figure;

plot(timeVec,sigVec);

xlabel('Time (sec)');
ylabel('Signal amplitude');
title('Amplitude-Modulated Sinusoidal Signal');
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
title('Periodogram of the AM Sinusoid');
grid on;

%% Spectrogram (Lab 2): requires Signal Processing Toolbox
% Window duration controls the time/frequency resolution tradeoff.
% Zero padding refines the frequency grid, not the physical resolution.
windowDuration = 0.40;
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
title('Amplitude Modulation Spectrogram');
hold on;
yline(f0-f1,'w--','LineWidth',1.5);
yline(f0+f1,'w--','LineWidth',1.5);
% Product modulation has two spectral lines at f0-f1 and f0+f1.
hold off;
% Complete windows leave uncovered margins at both ends of the record.
fprintf('Spectrogram: window = %d samples (%.3f s), overlap = %d, nfft = %d\n', ...
    windowLength,windowLength/samplFreq,overlapLength,nfft);
fprintf('Time step = %.3f s; frequency grid spacing = %.3f Hz\n', ...
    (windowLength-overlapLength)/samplFreq,samplFreq/nfft);
assert(all(isfinite(pSpec(:))) && all(pSpec(:)>=0), ...
    'Spectrogram power must be finite and nonnegative.');
