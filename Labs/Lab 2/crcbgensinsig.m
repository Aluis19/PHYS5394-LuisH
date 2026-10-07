function sigVec = crcbgensinsig(dataX,snr,sinParams)
% Generate a sinusoidal signal
% S = CRCBGENSINSIG(X,SNR,P)
% Generates a sinusoidal signal S. X is the vector of time stamps at which
% the samples of the signal are to be computed. SNR is the matched filtering
% signal-to-noise ratio of S and P is the vector [f0,phi0].

% Luis A. Hernandez, September 2026

f0 = sinParams(1);
phi0 = sinParams(2);
phaseVec = 2*pi*f0*dataX + phi0;
sigVec = sin(phaseVec);
sigVec = snr*sigVec/norm(sigVec);