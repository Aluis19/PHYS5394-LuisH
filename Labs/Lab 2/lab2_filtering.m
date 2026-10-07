%% Lab 2 - Filtering: stage 1, construct and verify the input
% Place crcbgensinsig.m on the MATLAB path.
% This stage does not design or apply filters yet.
% No Signal Processing Toolbox required for this stage.

%% Parameters specified in the assignment
nSamples = 2048;
samplFreq = 1024; % Hz
amplitudes = [10, 5, 2.5]; % Physical sinusoidal amplitudes, NOT norm/SNR
frequencies = [100, 200, 300]; % Hz
phases = [0, pi/6, pi/4]; % radians
signalDuration = nSamples/samplFreq; % 2 seconds
nyquistFrequency = samplFreq/2; % 512 Hz
timeVec = (0:nSamples-1)/samplFreq;
freqRes = samplFreq/nSamples; % 0.5 Hz

%% Generate each component using the Lab 1 function
components = zeros(3,nSamples); % rows: components; columns: time samples
for componentIdx = 1:3
    f = frequencies(componentIdx);
    phi = phases(componentIdx);
    rawSignal = sin(2*pi*f*timeVec+phi);
    % The generator fixes norm(sigVec), rather than physical amplitude.
    % Request the norm corresponding to the required physical amplitude.
    requestedNorm = amplitudes(componentIdx)*norm(rawSignal);
    components(componentIdx,:) = crcbgensinsig(timeVec,requestedNorm,[f,phi]);
end
inputSignal = sum(components,1);

%% One-sided rectangular-window power periodogram of the input
fftInput = fft(inputSignal);
kNyq = floor(nSamples/2)+1;
posFreq = (0:kNyq-1)*freqRes;
inputPSD = abs(fftInput(1:kNyq)).^2/(samplFreq*nSamples);
inputPSD(2:end-1) = 2*inputPSD(2:end-1); % N is even; do not double DC/Nyquist

%% Numerical checks against the analytical powers
expectedPowers = amplitudes.^2/2;
measuredPowers = mean(components.^2,2).';
inputMeanSquare = mean(inputSignal.^2);
integratedPSD = sum(inputPSD)*freqRes;
% All tones complete integer cycles, so they are orthogonal in this record.
assert(max(abs(measuredPowers-expectedPowers)) < 1e-10, ...
    'Component powers differ from A^2/2.');
assert(abs(inputMeanSquare-sum(expectedPowers)) < 1e-10, ...
    'Total power differs from the sum of component powers.');
assert(abs(integratedPSD-inputMeanSquare) < 1e-10, ...
    'Periodogram failed the Parseval check.');

% Tone frequencies coincide exactly with DFT bins in this experiment.
toneBins = round(frequencies/freqRes)+1;
peakPSD = inputPSD(toneBins);
toneTable = table(frequencies.',amplitudes.',phases.', ...
    measuredPowers.',peakPSD.', ...
    'VariableNames',{'Frequency_Hz','Amplitude','Phase_rad', ...
    'MeanSquare','PeakPSD'});
disp(toneTable);
fprintf('N = %d, Fs = %g Hz, T = %g s, df = %g Hz\n', ...
    nSamples,samplFreq,signalDuration,freqRes);
fprintf('Nyquist frequency = %g Hz\n',nyquistFrequency);
fprintf('Input mean square = %.6f; integrated PSD = %.6f\n', ...
    inputMeanSquare,integratedPSD);

%% Display components, sum and input periodogram
figure('Name','Lab 2 - Filtering input');
layout = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for componentIdx = 1:3
    nexttile(2*componentIdx-1);
    plot(timeVec,components(componentIdx,:),'.-');
    xlim([0,0.05]); ylim(1.15*[-1,1]*amplitudes(componentIdx));
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Component %d: %g Hz, A = %g',componentIdx, ...
        frequencies(componentIdx),amplitudes(componentIdx)));
    grid on;
end
nexttile(2);
plot(timeVec,inputSignal,'.-');
xlim([0,0.05]); xlabel('Time (s)'); ylabel('Amplitude');
title('Input: sum of the three sinusoids'); grid on;
nexttile(4,[2,1]);
stem(posFreq,inputPSD,'.','MarkerSize',8);
xlim([0,nyquistFrequency]);
xlabel('Frequency (Hz)'); ylabel('PSD (amplitude^2/Hz)');
title('Input power periodogram'); grid on;
sgtitle(layout,'Filtering stage 1: verify the input before designing filters');

%% Next stage: fir1 designs and filter application will be added here

%% Stage 2 - Design three FIR filters (Signal Processing Toolbox required)
filterOrder = 100; % 101 taps, even order for the highpass design
cutoffsHz = {150, [150,250], 250};
filterTypes = {'low','bandpass','high'};
filterNames = {'Lowpass: retain 100 Hz'; ...
    'Bandpass: retain 200 Hz'; 'Highpass: retain 300 Hz'};
filterCoeffs = cell(3,1);
filteredSignals = zeros(3,nSamples);
toneResponse = zeros(3,3); % rows: filters, columns: input tones
responseFrequency = linspace(0,nyquistFrequency,4097);
filterResponses = zeros(3,numel(responseFrequency));

for filterIdx = 1:3
    normalizedCutoff = cutoffsHz{filterIdx}/nyquistFrequency;
    % fir1 uses a Hamming window by default, specified explicitly here.
    filterCoeffs{filterIdx} = fir1(filterOrder,normalizedCutoff, ...
        filterTypes{filterIdx},hamming(filterOrder+1));
    filteredSignals(filterIdx,:) = filter(filterCoeffs{filterIdx},1,inputSignal);
    H = freqz(filterCoeffs{filterIdx},1,responseFrequency,samplFreq);
    filterResponses(filterIdx,:) = H(:).';
    Htones = freqz(filterCoeffs{filterIdx},1,frequencies,samplFreq);
    toneResponse(filterIdx,:) = Htones(:).';
end

groupDelaySamples = filterOrder/2;
fprintf('FIR order = %d, taps = %d, group delay = %d samples (%.6f s)\n', ...
    filterOrder,filterOrder+1,groupDelaySamples,groupDelaySamples/samplFreq);
toneGainDB = 20*log10(max(abs(toneResponse),1e-12));
disp(table(filterNames,toneGainDB(:,1),toneGainDB(:,2),toneGainDB(:,3), ...
    'VariableNames',{'Filter','Gain100Hz_dB','Gain200Hz_dB','Gain300Hz_dB'}));

%% Verify selectivity at the actual tones
% These are our verification criteria, not specifications in the handout.
for filterIdx = 1:3
    unwanted = setdiff(1:3,filterIdx);
    assert(abs(toneGainDB(filterIdx,filterIdx)) < 0.1, ...
        'Desired-tone gain differs from unity by more than 0.1 dB.');
    assert(all(toneGainDB(filterIdx,unwanted) < -40), ...
        'An unwanted tone is attenuated by less than 40 dB.');
end

%% Display filter magnitude responses
figure('Name','Lab 2 - FIR responses');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
for filterIdx = 1:3
    nexttile;
    plot(responseFrequency,20*log10(max(abs(filterResponses(filterIdx,:)),1e-6)));
    hold on;
    plot(frequencies,toneGainDB(filterIdx,:),'o');
    hold off;
    xlim([0,nyquistFrequency]); ylim([-100,5]);
    xlabel('Frequency (Hz)'); ylabel('Gain (dB)');
    title(filterNames{filterIdx}); grid on;
end

%% Periodograms of the full input/output records, including startup
% Zero initial conditions in filter cause a transient during the first
% filterOrder samples. Full-record PSDs therefore include transient leakage.
allRecords = [inputSignal; filteredSignals];
allFFT = fft(allRecords,[],2);
allPSD = abs(allFFT(:,1:kNyq)).^2/(samplFreq*nSamples);
allPSD(:,2:end-1) = 2*allPSD(:,2:end-1);
for recordIdx = 1:4
    assert(abs(sum(allPSD(recordIdx,:))*freqRes ...
        -mean(allRecords(recordIdx,:).^2)) < 1e-10, ...
        'Input/output PSD failed the Parseval check.');
end
figure('Name','Lab 2 - Input and filtered output periodograms');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
recordNames = [{'Input: all three tones'}; filterNames];
for recordIdx = 1:4
    nexttile;
    plot(posFreq,10*log10(max(allPSD(recordIdx,:),1e-12)));
    xlim([0,nyquistFrequency]); ylim([-100,25]);
    xlabel('Frequency (Hz)'); ylabel('PSD (dB re 1 amplitude^2/Hz)');
    title(recordNames{recordIdx}); grid on;
end
sgtitle('Full-record periodograms: startup transient included');

%% Check steady-state output against the analytical filter response
steadyIdx = (filterOrder+1):nSamples;
expectedSteadyPower = sum(abs(toneResponse).^2 .* expectedPowers,2);
steadyFitError = zeros(3,1);
measuredSteadyAmplitudes = zeros(3,3);
% Fit all three known tones; this avoids DFT-bin leakage after trimming.
fitMatrix = zeros(numel(steadyIdx),6);
for componentIdx = 1:3
    fitMatrix(:,2*componentIdx-1) = ...
        sin(2*pi*frequencies(componentIdx)*timeVec(steadyIdx)).';
    fitMatrix(:,2*componentIdx) = ...
        cos(2*pi*frequencies(componentIdx)*timeVec(steadyIdx)).';
end
for filterIdx = 1:3
    predicted = zeros(size(timeVec));
    for componentIdx = 1:3
        Htone = toneResponse(filterIdx,componentIdx);
        predicted = predicted + amplitudes(componentIdx)*abs(Htone)* ...
            sin(2*pi*frequencies(componentIdx)*timeVec ...
            +phases(componentIdx)+angle(Htone));
    end
    steadyFitError(filterIdx) = max(abs( ...
        filteredSignals(filterIdx,steadyIdx)-predicted(steadyIdx)));
    assert(steadyFitError(filterIdx) < 1e-9, ...
        'Steady-state output differs from the analytical response.');
    fitCoeffs = fitMatrix\filteredSignals(filterIdx,steadyIdx).';
    measuredSteadyAmplitudes(filterIdx,:) = ...
        hypot(fitCoeffs(1:2:end),fitCoeffs(2:2:end)).';
end
disp(table(filterNames,measuredSteadyAmplitudes(:,1), ...
    measuredSteadyAmplitudes(:,2),measuredSteadyAmplitudes(:,3), ...
    expectedSteadyPower,steadyFitError, ...
    'VariableNames',{'Filter','Amplitude100Hz','Amplitude200Hz', ...
    'Amplitude300Hz','PredictedSteadyPower','MaxSteadyError'}));
fprintf('Filter selectivity, PSD normalization and steady-state checks passed.\n');

%% Show causal startup and delay for each output
figure('Name','Lab 2 - Filtered signals in time');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
for filterIdx = 1:3
    nexttile;
    plot(timeVec,filteredSignals(filterIdx,:));
    xline(filterOrder/samplFreq,'--r','End of startup');
    xlim([0,0.2]); xlabel('Time (s)'); ylabel('Amplitude');
    title(filterNames{filterIdx}); grid on;
end
