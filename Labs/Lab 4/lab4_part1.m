%% Lab Topic 2_2 - Part 1: Colored Gaussian Noise
clear;
close all;
clc;

% Reproducible random realizations
rng(5394);

%% 1. Parameters and target PSD
sampFreq = 1024;                 % Sampling frequency (Hz)
sampleCounts = 16384*[1, 2, 4]; % N, 2N, and 4N
fltrOrdr = 500;                  % FIR filter order

% Positive-frequency branch of the target TWO-SIDED PSD
targetPSD = @(f) (f >= 100 & f <= 300) ...
                    .* (f-100) .* (300-f)/10000;

% Frequency grid includes DC and Nyquist
freqVec = 0:0.1:sampFreq/2;
psdVec = targetPSD(freqVec);

% Convert target to ONE-SIDED PSD for comparison with pwelch
targetOneSided = 2*psdVec;
targetOneSided([1,end]) = psdVec([1,end]);

% Fixed Welch settings: longer records provide more averages
windowLength = 256;
welchWindow = hamming(windowLength);
nOverlap = windowLength/2;
nFFT = windowLength;

% Store generated records
noiseRecords = cell(size(sampleCounts));

%% 2. Generate noise and compare PSD estimates
figure('Name','PSD comparison','Color','w');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

for k = 1:numel(sampleCounts)

    nSamples = sampleCounts(k);

    % Generate extra samples so the startup transient can be discarded.
    % The retained record contains exactly nSamples samples.
    noiseWithTransient = statgaussnoisegen( ...
        nSamples + fltrOrdr, ...
        [freqVec(:), psdVec(:)], ...
        fltrOrdr, sampFreq);

    outNoise = noiseWithTransient(fltrOrdr+1:end);
    noiseRecords{k} = outNoise;

    % One-sided PSD estimate, in data units^2/Hz
    [pxx,f] = pwelch(outNoise,welchWindow, ...
                    nOverlap,nFFT,sampFreq);

    nexttile;
    plot(f,pxx,'b','LineWidth',1);
    hold on;
    plot(freqVec,targetOneSided,'k--','LineWidth',1.5);

    xlabel('Frequency (Hz)');
    ylabel('PSD (data units^2/Hz)');
    title(sprintf('N = %d; duration = %.0f s', ...
                  nSamples,nSamples/sampFreq));
    legend('Welch estimate','One-sided target','Location','northeast');
    xlim([0 sampFreq/2]);
    ylim([0 3.5]);
    grid on;
end

sgtitle('Colored Gaussian noise: effect of record length');

%% 3. Inspect the baseline time series
outNoise = noiseRecords{1};
nSamples = sampleCounts(1);
timeVec = (0:nSamples-1)/sampFreq;

figure('Name','Noise time series','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(timeVec,outNoise,'b');
xlabel('Time (s)');
ylabel('Amplitude (data units)');
title('Colored Gaussian noise: full baseline record');
xlim([0 timeVec(end)]);
grid on;

% A short interval reveals the temporal structure
nexttile;
plot(timeVec,outNoise,'b.-');
xlabel('Time (s)');
ylabel('Amplitude (data units)');
title('Colored Gaussian noise: first 0.1 s');
xlim([0 0.1]);
grid on;

%% 4. Compare with white Gaussian noise
% Match the sample standard deviation to compare temporal structure.
whiteNoise = randn(size(outNoise));
whiteNoise = whiteNoise - mean(whiteNoise);
whiteNoise = whiteNoise/std(whiteNoise)*std(outNoise);

commonLimit = 1.1*max(abs([outNoise(:); whiteNoise(:)]));

figure('Name','White versus colored noise','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(timeVec,whiteNoise,'r.-');
xlabel('Time (s)');
ylabel('Amplitude (data units)');
title('White Gaussian noise: first 0.1 s');
xlim([0 0.1]);
ylim([-commonLimit commonLimit]);
grid on;

nexttile;
plot(timeVec,outNoise,'b.-');
xlabel('Time (s)');
ylabel('Amplitude (data units)');
title('Colored Gaussian noise: first 0.1 s');
xlim([0 0.1]);
ylim([-commonLimit commonLimit]);
grid on;

%% 5. Histogram and Gaussian PDF
sampleMean = mean(outNoise);
sampleStd = std(outNoise);

amplitudeGrid = linspace(sampleMean-4*sampleStd, ...
                         sampleMean+4*sampleStd,500);

% Gaussian PDF using the measured mean and standard deviation
gaussianPDF = exp(-0.5*((amplitudeGrid-sampleMean)/sampleStd).^2) ...
              /(sqrt(2*pi)*sampleStd);

figure('Name','Amplitude distribution','Color','w');

histogram(outNoise,60,'Normalization','pdf', ...
          'FaceColor',[0.3 0.6 0.9],'EdgeColor','none');
hold on;
plot(amplitudeGrid,gaussianPDF,'k','LineWidth',1.8);

xlabel('Noise amplitude (data units)');
ylabel('Probability density (1/data units)');
title('Colored noise: amplitude distribution');
legend('Normalized histogram','Gaussian PDF using sample mean and std');
grid on;

%% 6. Numerical checks
% Integrating the one-sided target PSD predicts the variance.
targetVariance = trapz(freqVec,targetOneSided);

% Re-estimate the baseline PSD for a variance comparison
[baselinePSD,baselineFreq] = pwelch( ...
    outNoise,welchWindow,nOverlap,nFFT,sampFreq);

estimatedVariance = trapz(baselineFreq,baselinePSD);

fprintf('\nBaseline record:\n');
fprintf('Samples                      : %d\n',nSamples);
fprintf('Duration                     : %.2f s\n',nSamples/sampFreq);
fprintf('Sample mean                  : %.4f\n',sampleMean);
fprintf('Sample variance              : %.4f\n',var(outNoise,1));
fprintf('Integrated Welch PSD         : %.4f\n',estimatedVariance);
fprintf('Integrated target PSD        : %.4f\n',targetVariance);
fprintf('Welch frequency spacing      : %.2f Hz\n',sampFreq/nFFT);