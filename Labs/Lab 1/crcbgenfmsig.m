function sigVec = crcbgenfmsig(dataX,snr,fmParams)
% Generate a frequency-modulated sinusoidal signal
% S = CRCBGENFMSIG(X,SNR,P)
% Generates an FM sinusoidal signal S. X is the vector of time stamps at
% which the samples of the signal are to be computed. SNR is the matched
% filtering signal-to-noise ratio of S and P is the vector [b,f0,f1].

% b  is the modulation index in radians.
% f0 is the carrier frequency in Hz.
% f1 is the modulation frequency in Hz.

% Luis A. Hernandez, September 2026

b = fmParams(1);
f0 = fmParams(2);
f1 = fmParams(3);

phaseVec = 2*pi*f0*dataX + b*cos(2*pi*f1*dataX);
sigVec = sin(phaseVec);
sigVec = snr*sigVec/norm(sigVec);