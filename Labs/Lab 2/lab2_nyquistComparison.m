%% Lab 2 - Nyquist comparison: all eight signals
% Keep this script, the parameter script and Lab 1 generators on the path.
% No Signal Processing Toolbox required.
clear;
clc;

lab2_signalParameters;

%% Select the signal and use a common physical amplitude
signalIdx = 8; % 1=sin, 2=LC, 3=QC, 4=transient LC, 5=SG, 6=AM, 7=FM, 8=AM-FM
assert(ismember(signalIdx,1:8), 'Select an integer signalIdx from 1 to 8.');
p = signals(signalIdx).params;
f0 = p(1);
T = signals(signalIdx).duration;
A = 1; % Common waveform amplitude coefficient, NOT discrete norm snr.
rates = signals(signalIdx).samplingRates;
% rawWave matches the unnormalized waveform in each Lab 1 generator.
switch signalIdx
    case 1
        phi0 = p(2);
        rawWave = @(t) sin(2*pi*p(1)*t + phi0);
    case 2
        rawWave = @(t) sin(2*pi*(p(1)*t + 0.5*p(2)*t.^2) + p(3));
    case 3
        % Phase in cycles: a1*t + a2*t^2 + a3*t^3.
        rawWave = @(t) sin(2*pi*(p(1)*t + p(2)*t.^2 + p(3)*t.^3));
    case 4
        % Local time tau=t-ta. The quadratic coefficient has NO factor 1/2.
        % Match the generator's inclusive transient interval [ta, ta+L].
        rawWave = @(t) double(t >= p(1) & t <= p(1)+p(5)) .* ...
            sin(2*pi*(p(2)*(t-p(1)) + p(3)*(t-p(1)).^2) + p(4));
    case 5
        % Gaussian envelope times sinusoidal carrier: [t0,sigma,f0,phi0].
        rawWave = @(t) exp(-(t-p(1)).^2/(2*p(2)^2)) .* ...
            sin(2*pi*p(3)*t + p(4));
    case 6
        % Suppressed-carrier AM: [f0,f1,phi0].
        rawWave = @(t) cos(2*pi*p(2)*t) .* sin(2*pi*p(1)*t+p(3));
    case 7
        % FM: [b,f0,f1], modulation index in radians.
        rawWave = @(t) sin(2*pi*p(2)*t+p(1)*cos(2*pi*p(3)*t));
    case 8
        % Common modulation frequency for amplitude and phase.
        rawWave = @(t) cos(2*pi*p(3)*t) .* ...
            sin(2*pi*p(2)*t+p(1)*cos(2*pi*p(3)*t));
end
referenceTime = linspace(0,T,4001);
referenceSignal = A*rawWave(referenceTime);
results = struct([]);
figure('Name',['Lab 2 - ', signals(signalIdx).name]);
layout = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
psdAxes = gobjects(1,3);
maxPSD = 0;

for rateIdx = 1:numel(rates)
    samplFreq = rates(rateIdx);
    nSamples = round(T*samplFreq);
    timeVec = (0:nSamples-1)/samplFreq;
    rawSignal = rawWave(timeVec);

    % Exactly Nyquist with phi0=0 yields sin(pi*n)=0 analytically.
    % Do not amplify floating-point residuals with norm normalization.
    degenerate = false;
    if signalIdx == 1
        degenerate = abs(samplFreq-2*p(1)) < 100*eps(samplFreq) ...
            && abs(sin(p(2))) < 100*eps;
    elseif signalIdx == 6
        degenerate = abs(samplFreq-2*p(1)) < 100*eps(samplFreq) ...
            && abs(sin(p(3))) < 100*eps;
    elseif signalIdx == 5
        degenerate = abs(samplFreq-2*p(3)) < 100*eps(samplFreq) ...
            && abs(sin(p(4))) < 100*eps;
    end
    if degenerate
        sigVec = zeros(size(timeVec));
        fprintf('Fs = %g Hz: zero-phase carrier samples are zero at Nyquist.\n', ...
            samplFreq);
    else
        assert(norm(rawSignal) > 1e-12, 'Cannot normalize a zero waveform.');
        generator = str2func(signals(signalIdx).generator);
        normalizedSignal = generator(timeVec,signals(signalIdx).snr, ...
            signals(signalIdx).params);
        % Undo the generator's norm scaling to compare fixed peak amplitude.
        sigVec = A*normalizedSignal*norm(rawSignal)/signals(signalIdx).snr;
    end

    %% One-sided power periodogram, rectangular window, no zero padding
    fftSig = fft(sigVec);
    kNyq = floor(nSamples/2)+1;
    posFreq = (0:kNyq-1)*samplFreq/nSamples;
    powerPSD = abs(fftSig(1:kNyq)).^2/(samplFreq*nSamples);
    if rem(nSamples,2) == 0
        powerPSD(2:end-1) = 2*powerPSD(2:end-1);
    else
        powerPSD(2:end) = 2*powerPSD(2:end);
    end

    %% Check that integrated PSD equals mean squared sample amplitude
    freqRes = samplFreq/nSamples;
    meanSquare = mean(sigVec.^2);
    integratedPSD = sum(powerPSD)*freqRes;
    assert(abs(integratedPSD-meanSquare) < 1e-12*max(1,meanSquare), ...
        'Periodogram normalization failed the Parseval check.');
    fprintf('Fs = %g Hz, N = %d, mean square = %.6f, PSD integral = %.6f\n', ...
        samplFreq,nSamples,meanSquare,integratedPSD);

    nexttile;
    plot(referenceTime,referenceSignal,'Color',[0.75 0.75 0.75]);
    hold on;
    stem(timeVec,sigVec,'.','MarkerSize',10);
    hold off;
    xlim([0,T]); ylim([-1.2 1.2]*A);
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Fs = %g Hz, N = %d',samplFreq,nSamples));
    legend('Continuous reference','Samples','Location','southoutside');
    grid on;
    if signalIdx == 4
        xline(p(1),'--r','Start');
        xline(p(1)+p(5),'--r','End');
    end

    psdAxes(rateIdx) = nexttile;
    stem(posFreq,powerPSD,'.','MarkerSize',10);
    xlim([0,1.2*signals(signalIdx).maxInstFreq]);
    maxPSD = max(maxPSD,max(powerPSD));
    xline(samplFreq/2,'--r','Fs/2');
    xlabel('Frequency (Hz)'); ylabel('PSD (amplitude^2/Hz)');
    title(sprintf('%g x Nyquist: power periodogram', ...
        nyquistMultipliers(rateIdx)));
    grid on;

    results(rateIdx).samplFreq = samplFreq;
    results(rateIdx).timeVec = timeVec;
    results(rateIdx).sigVec = sigVec;
    results(rateIdx).posFreq = posFreq;
    results(rateIdx).powerPSD = powerPSD;
    results(rateIdx).degenerate = degenerate;
end
% Common vertical scale for a fair spectral comparison.
for rateIdx = 1:numel(rates)
    ylim(psdAxes(rateIdx),[0,max(1e-6,1.15*maxPSD)]);
end
sgtitle(layout,sprintf('%s: common amplitude coefficient A = %g', ...
    signals(signalIdx).name,A));

%% Sampling interpretation
% Rates follow maximum instantaneous carrier frequency, as requested.
% SG: the Gaussian envelope broadens the spectrum around 30 Hz.
% AM: exact tones at 35 and 45 Hz; Fs=80 aliases 45 Hz onto 35 Hz.
% With phi0=0 these equal-amplitude components cancel after aliasing.
% FM and AM-FM: sidebands extend beyond the 60 Hz instantaneous maximum.
% Oversampling reduces aliasing of these tails; it does not change df=1/T.
