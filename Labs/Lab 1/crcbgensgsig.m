function sigVec = crcbgensgsig(dataX,snr,sgParams)
% Generate a sine-Gaussian signal
% S = CRCBGENSGSIG(X,SNR,P)
% Generates a sine-Gaussian signal S. X is the vector of time stamps at
% which the samples of the signal are to be computed. SNR is the matched
% filtering signal-to-noise ratio of S and P is [t0,sigma,f0,phi0].

% t0    is the center time of the Gaussian envelope in seconds.
% sigma is the standard deviation of the envelope in seconds.
% f0    is the sinusoidal frequency in Hz.
% phi0  is the initial phase in radians.

% Luis A. Hernandez, September 2026

t0 = sgParams(1);
sigma = sgParams(2);
f0 = sgParams(3);
phi0 = sgParams(4);

envelopeVec = exp(-(dataX-t0).^2/(2*sigma^2));
phaseVec = 2*pi*f0*dataX + phi0;
sigVec = envelopeVec.*sin(phaseVec);
sigVec = snr*sigVec/norm(sigVec);