function sigVec = crcbgenlcsig(dataX,snr,lcParams)
% Generate a linear chirp signal
% S = CRCBGENLCSIG(X,SNR,P)
% Generates a linear chirp signal S. X is the vector of time stamps at which
% the samples of the signal are to be computed. SNR is the matched filtering
% signal-to-noise ratio of S and P is the vector [f0,f1,phi0].

% f0   is the initial frequency in Hz.
% f1   is the linear frequency-evolution coefficient in Hz/sec.
% phi0 is the initial phase in radians.

% Luis A. Hernandez, September 2026

f0 = lcParams(1);
f1 = lcParams(2);
phi0 = lcParams(3);

phaseVec = 2*pi*(f0*dataX + 0.5*f1*dataX.^2) + phi0;
sigVec = sin(phaseVec);
sigVec = snr*sigVec/norm(sigVec);