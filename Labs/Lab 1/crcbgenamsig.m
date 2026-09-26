function sigVec = crcbgenamsig(dataX,snr,amParams)
% Generate an amplitude-modulated sinusoidal signal
% S = CRCBGENAMSIG(X,SNR,P)
% Generates an AM sinusoidal signal S. X is the vector of time stamps at
% which the samples of the signal are to be computed. SNR is the matched
% filtering signal-to-noise ratio of S and P is [f0,f1,phi0].

% f0   is the carrier frequency in Hz.
% f1   is the amplitude-modulation frequency in Hz.
% phi0 is the initial carrier phase in radians.

% Luis A. Hernandez, September 2026

f0 = amParams(1);
f1 = amParams(2);
phi0 = amParams(3);

amplitudeVec = cos(2*pi*f1*dataX);
phaseVec = 2*pi*f0*dataX + phi0;
sigVec = amplitudeVec.*sin(phaseVec);
sigVec = snr*sigVec/norm(sigVec);