%% Lab Topic 2_2 - Part 2: Whitening noisy data
clear;
close all;
clc;

%% 1. Load data
% Column 1: sampling times (s)
% Column 2: measured data
data = load('testData.txt');

timeVec = data(:,1);
rawData = data(:,2);

assert(all(isfinite(data(:))), 'Data contains nonfinite values.');

% Check timestamp spacing
timeSteps = diff(timeVec);
assert(all(timeSteps > 0), ...
    'Sampling times must be strictly increasing.');

% Average spacing is less sensitive to text-file rounding.
dt = (timeVec(end)-timeVec(1))/(numel(timeVec)-1);

relativeVariation = max(abs(timeSteps-dt))/dt;

fprintf('Mean sampling interval      : %.10g s\n',dt);
fprintf('Smallest sampling interval  : %.10g s\n',min(timeSteps));
fprintf('Largest sampling interval   : %.10g s\n',max(timeSteps));
fprintf('Maximum relative variation  : %.6f %%\n', ...
    100*relativeVariation);

% Practical screening threshold: 0.1% of a sampling interval.
assert(relativeVariation < 1e-3, ...
    ['Timestamp variation exceeds 0.1%%. Inspect the printed ', ...
    'intervals before treating the data as uniformly sampled.']);

sampFreq = 1/dt;
nSamples = numel(rawData);

signalStart = 5;                 % Signal begins after 5 s
noiseMask = timeVec < signalStart;
nNoise = nnz(noiseMask);

assert(nNoise >= 256, 'Too few noise-only samples.');
assert(any(timeVec > signalStart), 'No data after 5 s.');

% Estimate the offset using ONLY the noise-only interval.
noiseMean = mean(rawData(noiseMask));
inputData = rawData - noiseMean;
noiseData = inputData(noiseMask);

fprintf('Sampling frequency : %.4f Hz\n',sampFreq);
fprintf('Total samples     : %d\n',nSamples);
fprintf('Noise-only samples: %d\n',nNoise);

%% 2. Estimate the noise PSD using Welch's method
% Aim for segments of about 0.5 s, with at least several
% segments available within the noise-only record.
windowLength = min(round(0.5*sampFreq),floor(nNoise/4));
windowLength = max(32,windowLength);

welchWindow = hamming(windowLength);
nOverlap = floor(windowLength/2);
nFFT = max(256,2^nextpow2(windowLength));

[noisePSD,freqVec] = pwelch( ...
    noiseData,welchWindow,nOverlap,nFFT,sampFreq);

%% 3. Design the FIR whitening filter
% A PSD floor prevents excessive gain where the estimated
% noise power is very small. This is regularization.
psdFloor = 1e-6*max(noisePSD);
regularizedPSD = max(noisePSD,psdFloor);

% Target one-sided PSD for discrete noise with variance 1
whitePSDLevel = 2/sampFreq;

whiteningMagnitude = sqrt(whitePSDLevel./regularizedPSD);

% At DC and Nyquist, one-sided PSD bins are not doubled.
whiteningMagnitude([1,end]) = sqrt( ...
    (1/sampFreq)./regularizedPSD([1,end]));

% Even filter order gives an integer group delay.
% Limit filter duration relative to available noise-only data.
fltrOrdr = min(500,2*floor(nNoise/8));
filterDelay = fltrOrdr/2;

b = fir2(fltrOrdr, ...
         freqVec/(sampFreq/2),whiteningMagnitude);

%% 4. Apply the filter and compensate its group delay
% Append zeros to retain the delayed output through the end.
paddedData = [inputData; zeros(fltrOrdr,1)];
filteredData = fftfilt(b,paddedData);

whiteData = filteredData(filterDelay+(1:nSamples));

% Delay compensation uses samples on both sides of each time.
% Exclude record edges affected by missing neighboring data.
validMask = false(nSamples,1);
validMask(filterDelay+1:nSamples-filterDelay) = true;

% Keep PSD validation away from the signal onset as well.
cleanNoiseMask = validMask ...
              & (timeVec < signalStart-filterDelay/sampFreq);

assert(nnz(cleanNoiseMask) >= windowLength, ...
       'Insufficient interior noise-only data for PSD validation.');

%% 5. Compare noise PSDs before and after whitening
% Use the same time interval for both estimates.
[beforePSD,fPSD] = pwelch( ...
    inputData(cleanNoiseMask),welchWindow,nOverlap,nFFT,sampFreq);

[afterPSD,~] = pwelch( ...
    whiteData(cleanNoiseMask),welchWindow,nOverlap,nFFT,sampFreq);

figure('Name','Noise PSD before and after whitening','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
semilogy(fPSD,max(beforePSD,realmin),'b');
xlabel('Frequency (Hz)');
ylabel('PSD (data units^2/Hz)');
title('Noise-only PSD before whitening');
grid on;
xlim([0 sampFreq/2]);

nexttile;
semilogy(fPSD,max(afterPSD,realmin),'b');
hold on;
yline(whitePSDLevel,'k--','Target level: 2/f_s');
xlabel('Frequency (Hz)');
ylabel('PSD (1/Hz)');
title('Noise-only PSD after whitening');
grid on;
xlim([0 sampFreq/2]);

%% 6. Compare time series
figure('Name','Time series before and after whitening','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
plot(timeVec(validMask),inputData(validMask),'b');
xline(signalStart,'r--','Signal onset');
xlabel('Time (s)');
ylabel('Amplitude (data units)');
title('Before whitening: noise offset removed');
grid on;

ax2 = nexttile;
plot(timeVec(validMask),whiteData(validMask),'b');
xline(signalStart,'r--','Signal onset');
xlabel('Time (s)');
ylabel('Whitened amplitude');
title('After whitening');
grid on;

linkaxes([ax1 ax2],'x');
xlim([timeVec(find(validMask,1,'first')), ...
      timeVec(find(validMask,1,'last'))]);

%% 7. Spectrograms before and after whitening
% Identical STFT settings allow comparison of time-frequency structure.
specLength = windowLength;
specWindow = hamming(specLength);
specOverlap = floor(0.9*specLength);
specNFFT = max(512,2^nextpow2(specLength));

validTime = timeVec(validMask);

[~,specFreq,specTime,pBefore] = spectrogram( ...
    inputData(validMask),specWindow,specOverlap,specNFFT,sampFreq);

[~,~,~,pAfter] = spectrogram( ...
    whiteData(validMask),specWindow,specOverlap,specNFFT,sampFreq);

% Spectrogram times are relative to the first retained sample.
specTime = specTime + validTime(1);

beforeDB = 10*log10(max(pBefore,realmin));
afterDB = 10*log10(max(pAfter,realmin));

figure('Name','Spectrograms before and after whitening','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(specTime,specFreq,beforeDB);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Before whitening: PSD in dB re 1 data unit^2/Hz');
xline(signalStart,'w--','LineWidth',1.5);
colorbar;

nexttile;
imagesc(specTime,specFreq,afterDB);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('After whitening: PSD in dB re 1/Hz');
xline(signalStart,'w--','LineWidth',1.5);
colorbar;

colormap(parula);

%% 8. Numerical diagnostics
fprintf('\nWhitening diagnostics:\n');
fprintf('Welch segment length        : %d samples\n',windowLength);
fprintf('Welch frequency spacing     : %.4f Hz\n',sampFreq/nFFT);
fprintf('FIR filter order            : %d\n',fltrOrdr);
fprintf('Compensated group delay     : %.4f s\n', ...
        filterDelay/sampFreq);
fprintf('Noise variance before       : %.4f\n', ...
        var(inputData(cleanNoiseMask),1));
fprintf('Noise variance after        : %.4f\n', ...
        var(whiteData(cleanNoiseMask),1));
fprintf('Approximate target variance : 1\n');
%% Improve spectrogram contrast
% Express each frequency bin relative to its mean noise-only power.
% This display normalization does not modify the filtered time series.

% Only use STFT windows entirely before the signal-affected region.
halfWindowTime = specLength/(2*sampFreq);

noiseColumns = specTime + halfWindowTime ...
    < signalStart - filterDelay/sampFreq;

assert(any(noiseColumns), ...
    'No uncontaminated spectrogram windows available.');

referencePower = mean(pAfter(:,noiseColumns),2);

relativePowerDB = 10*log10( ...
    max(pAfter,realmin)./max(referencePower,realmin));

figure('Name','Whitened spectrogram: noise-relative power','Color','w');

imagesc(specTime,specFreq,relativePowerDB);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Whitened data: power relative to noise-only mean');
xline(signalStart,'w--','Signal onset','LineWidth',1.5);

cb = colorbar;
cb.Label.String = 'Power / mean noise power (dB)';

clim([-10 10]);
colormap(parula);
