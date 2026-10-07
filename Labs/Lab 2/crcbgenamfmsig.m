function sigVec = crcbgenamfmsig(dataX,snr,amfmParams)
% Generate an amplitude- and frequency-modulated sinusoidal signal
% S = CRCBGENAMFMSIG(X,SNR,P)
% Generates an AM-FM sinusoidal signal S. X is the vector of time stamps at
% which the samples of the signal are to be computed. SNR is the matched
% filtering signal-to-noise ratio of S and P is the vector [b,f0,f1].

% b  is the FM modulation index in radians.
% f0 is the carrier frequency in Hz.
% f1 is the common AM and FM modulation frequency in Hz.

% Luis A. Hernandez, September 2026

b = amfmParams(1);
f0 = amfmParams(2);
f1 = amfmParams(3);

amplitudeVec = cos(2*pi*f1*dataX);
phaseVec = 2*pi*f0*dataX + b*cos(2*pi*f1*dataX);
sigVec = amplitudeVec.*sin(phaseVec);
sigVec = snr*sigVec/norm(sigVec);