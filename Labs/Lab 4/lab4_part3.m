%% Lab Topic 2_2 - Part 3: Simulated initial-LIGO noise
clear;
close all;
clc;
rng(5394);

%% 1. Load the initial-LIGO sensitivity curve
% Column 1: frequency (Hz)
% Column 2: ONE-SIDED amplitude spectral density (strain/sqrt(Hz))
sensitivityData = load('iLIGOSensitivity.txt','-ascii');

frequencyInput = sensitivityData(:,1);
asdInput = sensitivityData(:,2);

% Logarithmic interpolation requires positive, finite values.
validRows = isfinite(frequencyInput) & isfinite(asdInput) ...
          & frequencyInput > 0 & asdInput > 0;

frequencyInput = frequencyInput(validRows);
asdInput = asdInput(validRows);

[frequencyInput,sortIndex] = sort(frequencyInput);
asdInput = asdInput(sortIndex);

[frequencyInput,uniqueIndex] = unique(frequencyInput);
asdInput = asdInput(uniqueIndex);

%% 2. Simulation parameters
sampFreq = 2048;             % Hz; Nyquist frequency = 1024 Hz
nSamples = 65536;           % 32 s of retained noise
fltrOrdr = 2000;            % FIR filter order

lowCutoff = 50;             % Hz
highCutoff = 700;           % Hz

assert(frequencyInput(1) <= lowCutoff ...
    && frequencyInput(end) >= highCutoff, ...
    'The sensitivity file must cover 50 to 700 Hz.');

assert(sampFreq/2 > highCutoff, ...
    'Nyquist frequency must exceed 700 Hz.');

timeVec = (0:nSamples-1)/sampFreq;

%% 3. Construct the modified target spectrum
% Uniform design grid including DC and Nyquist.
freqVec = linspace(0,sampFreq/2,8193).';

% Clamp interpolation frequencies:
% f < 50 Hz uses the ASD at 50 Hz;
% f > 700 Hz uses the ASD at 700 Hz.
evaluationFreq = min(max(freqVec,lowCutoff),highCutoff);

% Interpolate in log-frequency and log-ASD.
targetASD = 10.^interp1( ...
    log10(frequencyInput),log10(asdInput), ...
    log10(evaluationFreq),'linear');

% Convert ONE-SIDED ASD to ONE-SIDED PSD.
targetPSDOneSided = targetASD.^2;

% Mohanty's generator expects a TWO-SIDED PSD.
targetPSDTwoSided = targetPSDOneSided/2;

%% 4. Generate the simulated strain noise
% Generate extra samples and remove the FIR startup transient.
noiseWithTransient = statgaussnoisegen( ...
    nSamples + fltrOrdr, ...
    [freqVec,targetPSDTwoSided], ...
    fltrOrdr,sampFreq);

strainNoise = noiseWithTransient(fltrOrdr+1:end);
strainNoise = strainNoise(:);

%% 5. Estimate the one-sided PSD
windowLength = 2048;         % 1-second Welch segments
welchWindow = hamming(windowLength);
nOverlap = windowLength/2;
nFFT = windowLength;

[estimatedPSD,estimatedFreq] = pwelch( ...
    strainNoise,welchWindow,nOverlap,nFFT,sampFreq);

estimatedASD = sqrt(estimatedPSD);

% Exclude DC from logarithmic frequency plots.
positiveFreq = estimatedFreq > 0;
positiveTarget = freqVec > 0;

%% 6. Plot the original and modified sensitivity curves
figure('Name','Initial-LIGO target sensitivity','Color','w');

loglog(frequencyInput,asdInput,'Color',[0.6 0.6 0.6], ...
       'LineWidth',1.2);
hold on;

loglog(freqVec(positiveTarget),targetASD(positiveTarget), ...
       'k','LineWidth',1.7);

xline(lowCutoff,'b--','50 Hz','HandleVisibility','off');
xline(highCutoff,'r--','700 Hz','HandleVisibility','off');

xlabel('Frequency (Hz)');
ylabel('ASD (strain / sqrt(Hz))');
title('Initial-LIGO sensitivity: original and modified targets');
legend('Original design curve','Modified simulation target', ...
       'Location','best');
xlim([10 sampFreq/2]);
grid on;

%% 7. Compare estimated and target PSDs
figure('Name','Simulated LIGO noise: PSD comparison','Color','w');

loglog(estimatedFreq(positiveFreq),estimatedPSD(positiveFreq), ...
       'b','LineWidth',1);
hold on;

loglog(freqVec(positiveTarget), ...
       targetPSDOneSided(positiveTarget), ...
       'k--','LineWidth',1.7);

xlabel('Frequency (Hz)');
ylabel('One-sided PSD (strain^2/Hz)');
title('Simulated initial-LIGO noise: Welch PSD versus target');
legend('Welch estimate','Modified target PSD','Location','best');
xlim([10 sampFreq/2]);
grid on;

%% 8. Compare amplitude spectral densities
figure('Name','Simulated LIGO noise: ASD comparison','Color','w');

loglog(estimatedFreq(positiveFreq),estimatedASD(positiveFreq), ...
       'b','LineWidth',1);
hold on;

loglog(freqVec(positiveTarget),targetASD(positiveTarget), ...
       'k--','LineWidth',1.7);

xlabel('Frequency (Hz)');
ylabel('ASD (strain / sqrt(Hz))');
title('Simulated initial-LIGO noise: ASD comparison');
legend('Estimated ASD','Modified target ASD','Location','best');
xlim([10 sampFreq/2]);
grid on;

%% 9. Plot the simulated strain time series
figure('Name','Simulated LIGO strain noise','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(timeVec,strainNoise,'b');
xlabel('Time (s)');
ylabel('Strain');
title('Simulated initial-LIGO noise: full record');
xlim([0 timeVec(end)]);
grid on;

nexttile;
plot(timeVec,strainNoise,'b');
xlabel('Time (s)');
ylabel('Strain');
title('Simulated initial-LIGO noise: first 0.1 s');
xlim([0 0.1]);
grid on;

%% 10. Numerical checks
targetVariance = trapz(freqVec,targetPSDOneSided);
sampleVariance = var(strainNoise,1);
estimatedVariance = trapz(estimatedFreq,estimatedPSD);

fprintf('\nInitial-LIGO noise simulation:\n');
fprintf('Sampling frequency       : %.0f Hz\n',sampFreq);
fprintf('Number of samples        : %d\n',nSamples);
fprintf('Retained duration        : %.2f s\n',nSamples/sampFreq);
fprintf('FIR filter order         : %d\n',fltrOrdr);
fprintf('Welch frequency spacing  : %.2f Hz\n',sampFreq/nFFT);
fprintf('Sample mean              : %.6e strain\n',mean(strainNoise));
fprintf('Sample standard deviation: %.6e strain\n',std(strainNoise,1));
fprintf('Target variance          : %.6e strain^2\n',targetVariance);
fprintf('Sample variance          : %.6e strain^2\n',sampleVariance);
fprintf('Integrated Welch PSD     : %.6e strain^2\n',estimatedVariance);
fprintf('Variance difference      : %.2f %%\n', ...
        100*(sampleVariance-targetVariance)/targetVariance);