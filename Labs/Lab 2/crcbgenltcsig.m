function sigVec = crcbgenltcsig(dataX,snr,ltcParams)
% Generate a linear transient chirp signal
% S = CRCBGENLTCSIG(X,SNR,P)
% Generates a linear transient chirp S. X is the vector of time stamps at
% which the samples of the signal are to be computed. SNR is the matched
% filtering signal-to-noise ratio of S and P is [ta,f0,f1,phi0,L].

% ta   is the starting time of the transient in seconds.
% f0   is the initial frequency in Hz.
% f1   is the quadratic phase coefficient in Hz/sec.
% phi0 is the initial phase in radians.
% L    is the duration of the transient in seconds.

% Luis A. Hernandez, September 2026

ta = ltcParams(1);
f0 = ltcParams(2);
f1 = ltcParams(3);
phi0 = ltcParams(4);
signalDuration = ltcParams(5);


sigVec = zeros(size(dataX));
transientIdx = (dataX >= ta) & (dataX <= ta+signalDuration);
transientTime = dataX(transientIdx)-ta;
phaseVec = 2*pi*(f0*transientTime + f1*transientTime.^2) + phi0;
sigVec(transientIdx) = sin(phaseVec);
sigVec = snr*sigVec/norm(sigVec);

end