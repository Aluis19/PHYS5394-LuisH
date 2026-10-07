%% Lab 2 - Signal parameters and instantaneous-frequency sampling rates
% Run this script before lab2_nyquistComparison.m.
% Parameters and argument order match the uploaded Lab 1 generators.
% This stage creates configuration and a table; it does not generate signals.
% All frequencies are in Hz, times in seconds, and phases in radians.

%% Common settings
signalDuration = 1.0;
snr = 10; % Requested discrete Euclidean norm, NOT peak amplitude.
nyquistMultipliers = [1, 2, 10];

%% Parameters in the order required by each generator
sinParams  = [10, 0];              % [f0, phi0]
lcParams   = [10, 20, 0];          % [f0, f1, phi0]
qcCoefs    = [10, 3, 3];           % [a1, a2, a3], phase in cycles
ltcParams  = [0.2, 10, 20, 0, 0.5]; % [ta, f0, f1, phi0, L]
sgParams   = [0.5, 0.1, 30, 0];    % [t0, sigma, f0, phi0]
amParams   = [40, 5, 0];           % [f0, f1, phi0]
fmParams   = [4, 40, 5];           % [b, f0, f1]
amfmParams = [4, 40, 5];           % [b, f0, f1]

assert(ltcParams(1) >= 0 && ltcParams(5) > 0 && ...
    ltcParams(1)+ltcParams(5) <= signalDuration, ...
    'The transient must lie inside the observation interval.');
assert(sgParams(2) > 0, 'Gaussian sigma must be positive.');

%% Maximum absolute instantaneous frequency on the observation interval
% Sinusoid: f_inst = f0.
sinMax = abs(sinParams(1));

% Linear chirp: f_inst = f0 + f1*t.
lcMax = max(abs([lcParams(1), ...
    lcParams(1)+lcParams(2)*signalDuration]));

% Quadratic chirp: f_inst = a1 + 2*a2*t + 3*a3*t^2.
% Check endpoints and any interior stationary point of f_inst.
qcTimes = [0, signalDuration];
if qcCoefs(3) ~= 0
    stationaryTime = -qcCoefs(2)/(3*qcCoefs(3));
    if stationaryTime > 0 && stationaryTime < signalDuration
        qcTimes = [qcTimes, stationaryTime];
    end
end
qcMax = max(abs(qcCoefs(1) + 2*qcCoefs(2)*qcTimes ...
    + 3*qcCoefs(3)*qcTimes.^2));

% Transient chirp: tau = t-ta, f_inst = f0 + 2*f1*tau.
ltcMax = max(abs([ltcParams(2), ...
    ltcParams(2)+2*ltcParams(3)*ltcParams(5)]));

% Sine-Gaussian and AM: derivative of the specified carrier phase.
% The envelope introduces spectral width not captured by this derivative.
sgMax = abs(sgParams(3));
amMax = abs(amParams(1));

% FM and AM-FM: f_inst = f0-b*f1*sin(2*pi*f1*t).
% Both present parameter sets cover full modulation cycles in one second.
% Hence the full-cycle maximum absolute frequency is attained.
assert(fmParams(3)*signalDuration >= 1 && ...
    amfmParams(3)*signalDuration >= 1, ...
    'The full-cycle FM bound used here requires at least one cycle.');
fmMax = abs(fmParams(2)) + abs(fmParams(1)*fmParams(3));
amfmMax = abs(amfmParams(2)) + abs(amfmParams(1)*amfmParams(3));

%% Configuration for subsequent scripts
signalNames = {'Sinusoid'; 'Linear chirp'; 'Quadratic chirp'; ...
    'Transient linear chirp'; 'Sine-Gaussian'; 'AM'; 'FM'; 'AM-FM'};
generatorNames = {'crcbgensinsig'; 'crcbgenlcsig'; 'crcbgenqcsig'; ...
    'crcbgenltcsig'; 'crcbgensgsig'; 'crcbgenamsig'; ...
    'crcbgenfmsig'; 'crcbgenamfmsig'};
parameterVectors = {sinParams; lcParams; qcCoefs; ltcParams; ...
    sgParams; amParams; fmParams; amfmParams};
maxInstFreq = [sinMax; lcMax; qcMax; ltcMax; sgMax; amMax; fmMax; amfmMax];
nyquistRate = 2*maxInstFreq;
samplingRates = nyquistRate*nyquistMultipliers;

signals = struct('name', signalNames, 'generator', generatorNames, ...
    'params', parameterVectors);
for signalIdx = 1:numel(signals)
    signals(signalIdx).duration = signalDuration;
    signals(signalIdx).snr = snr;
    signals(signalIdx).maxInstFreq = maxInstFreq(signalIdx);
    signals(signalIdx).samplingRates = samplingRates(signalIdx,:);
end

%% Display the result for this stage
parameterTable = table(signalNames, maxInstFreq, samplingRates(:,1), ...
    samplingRates(:,2), samplingRates(:,3), ...
    'VariableNames', {'Signal','MaxInstFreq_Hz','Fs_Nyquist_Hz', ...
    'Fs_2Nyquist_Hz','Fs_10Nyquist_Hz'});
disp(parameterTable);

%% Interpretation notes for the next stage
% These rates follow the lab's instantaneous-frequency approximation.
% AM has spectral lines at 35 and 45 Hz despite carrier f_inst = 40 Hz.
% FM/AM-FM have infinitely many sidebands; f_inst max is not a band limit.
% Gaussian envelopes and transient gates also broaden spectra.
% A zero-phase 10 Hz sine sampled at exactly 20 Hz is mathematically zero.
% Do not normalize that degenerate sequence: handle it in the next stage.
% Lab 1 plots abs(FFT), a magnitude spectrum, not a power periodogram.
