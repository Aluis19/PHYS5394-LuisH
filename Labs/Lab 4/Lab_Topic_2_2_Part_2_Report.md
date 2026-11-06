# Part 2: Whitening Colored Noise and Revealing a Mystery Signal

**Author:** Luis Alonso Hernández Galván  
**Course:** PHYS 5394 — Advanced Statistical Methods for Modern Astronomy  
**Instructor:** Prof. Soumya D. Mohanty  
**Laboratory:** Lab Topic 2_2 — Part 2

## 1. Objective

I analyzed `testData.txt`, which contains colored Gaussian noise and a mystery signal added after 5 s. Following the laboratory instructions, I estimated the noise PSD from the signal-free interval, designed an FIR whitening filter, and compared the time series and spectrograms before and after filtering.

## 2. Loading and Preparing the Data

```matlab
data = load('testData.txt');
timeVec = data(:,1);
rawData = data(:,2);

timeSteps = diff(timeVec);
dt = (timeVec(end)-timeVec(1))/(numel(timeVec)-1);
sampFreq = 1/dt;
relativeVariation = max(abs(timeSteps-dt))/dt;

signalStart = 5;
noiseMask = timeVec < signalStart;
nNoise = nnz(noiseMask);

noiseMean = mean(rawData(noiseMask));
inputData = rawData - noiseMean;
noiseData = inputData(noiseMask);
```

I loaded the time stamps and data values from the two columns. I checked that the values were finite and that time increased strictly. I calculated the sampling interval from the full time span to reduce the effect of rounding in individual time stamps.

The inferred sampling frequency was approximately 1024 Hz. The individual intervals ranged from 0.000976 to 0.000977 s, with a maximum relative variation of 0.057597%. These values were consistent with rounded time stamps and passed the script's practical 0.1% screening threshold. This check did not resample the data.

I used only samples with $t<5$ s to estimate the noise mean and PSD. Subtracting this mean removed the estimated offset without using the signal-containing interval. The record contained 16,384 samples, including 5,120 noise-only samples.

## 3. Estimating the Noise PSD

```matlab
windowLength = min(round(0.5*sampFreq),floor(nNoise/4));
windowLength = max(32,windowLength);
welchWindow = hamming(windowLength);
nOverlap = floor(windowLength/2);
nFFT = max(256,2^nextpow2(windowLength));

[noisePSD,freqVec] = pwelch( ...
    noiseData,welchWindow,nOverlap,nFFT,sampFreq);
```

I used Welch's method with 512-sample Hamming windows, 50% overlap, and a 512-point FFT. Each segment covered approximately 0.5 s, and the frequency-grid spacing was 2 Hz. Averaging segment periodograms reduced spectral-estimation fluctuations.

For this real-valued input, `noisePSD` represents the one-sided estimate $\widehat S_n^{(1)}(f)$. I assumed that the noise statistics estimated before 5 s remained representative of the noise later in the record.

## 4. Designing the Whitening Filter

```matlab
psdFloor = 1e-6*max(noisePSD);
regularizedPSD = max(noisePSD,psdFloor);
whitePSDLevel = 2/sampFreq;

whiteningMagnitude = sqrt(whitePSDLevel./regularizedPSD);
whiteningMagnitude([1,end]) = sqrt( ...
    (1/sampFreq)./regularizedPSD([1,end]));

fltrOrdr = min(500,2*floor(nNoise/8));
filterDelay = fltrOrdr/2;
b = fir2(fltrOrdr,freqVec/(sampFreq/2),whiteningMagnitude);
```

I chose the filter magnitude so that

$$
|W(f)|=\sqrt{\frac{2/f_s}{\widehat S_n^{(1)}(f)}}.
$$

Since filtering multiplies the PSD by the squared transfer-function magnitude,

$$
S_{n,\mathrm{white}}^{(1)}(f)
=|W(f)|^2S_n^{(1)}(f)
\approx\frac{2}{f_s}.
$$

This flat target corresponds approximately to unit variance because its integral over $0<f<f_s/2$ is one. I treated DC and Nyquist separately because their one-sided PSD bins are not doubled.

The PSD floor limited excessive gain at frequencies with very small estimated noise power. `fir2` approximated the desired magnitude using an order-500 FIR filter; its frequency input was normalized by the Nyquist frequency.

Statistically, whitening reduces the dominance of noisy frequency bands. It also filters the signal, so the output signal waveform differs from the input waveform.

## 5. Filtering, Delay Compensation, and Boundary Exclusion

```matlab
paddedData = [inputData; zeros(fltrOrdr,1)];
filteredData = fftfilt(b,paddedData);
whiteData = filteredData(filterDelay+(1:numel(inputData)));

validMask = false(numel(inputData),1);
validMask(filterDelay+1:end-filterDelay) = true;
cleanNoiseMask = validMask ...
    & (timeVec < signalStart-filterDelay/sampFreq);
```

I applied the FIR filter with `fftfilt` and compensated its 250-sample group delay, approximately 0.2441 s. Zero padding retained the delayed output through the end of the record, but did not supply actual missing data.

I therefore excluded the affected record edges. Delay compensation makes the aligned filter depend on neighboring samples on both sides of a time stamp. I also excluded the interval immediately before 5 s from the noise-only validation to avoid signal leakage through that filter support.

## 6. PSD, Time-Series, and Spectrogram Comparisons

I estimated the PSD before and after whitening using the same interior noise-only interval and Welch settings. I plotted these estimates with `semilogy` and compared the whitened spectrum with the target $2/f_s$. I also plotted the aligned time series, marked the 5 s signal onset with `xline`, and linked their time axes.

For the spectrograms, I used identical 512-sample Hamming windows, approximately 90% overlap, and 512-point FFTs. The overlap provided closely spaced time views, while the window duration controlled the time-frequency tradeoff. I added the first retained time stamp to the spectrogram time vector to align it with the original record.

```matlab
[~,specFreq,specTime,pAfter] = spectrogram( ...
    whiteData(validMask),specWindow,specOverlap,specNFFT,sampFreq);

validTime = timeVec(validMask);
specTime = specTime + validTime(1);

halfWindowTime = specLength/(2*sampFreq);
noiseColumns = specTime + halfWindowTime ...
             < signalStart-filterDelay/sampFreq;

referencePower = mean(pAfter(:,noiseColumns),2);
relativePowerDB = 10*log10( ...
    max(pAfter,realmin)./max(referencePower,realmin));
```

The initial spectrogram color range obscured weak structure. I therefore divided each frequency bin by its mean power in uncontaminated noise-only windows and expressed the ratio in decibels. The `max(...,realmin)` operations avoided division by zero and logarithms of zero.

```matlab
imagesc(specTime,specFreq,relativePowerDB);
axis xy;
clim([-10 10]);
xlim([4 11]);
ylim([0 100]);
```

I displayed this ratio with frequency increasing upward and focused on 4–11 s and 0–100 Hz. Here 0 dB denotes the mean noise-only power at that frequency, and +10 dB denotes ten times that power. This display normalization did not modify `whiteData`. The original before-and-after spectrograms used different physical units and independently scaled colors, so their colors alone were not a direct power comparison.

## 7. Results and Interpretation

| Diagnostic | Result |
|---|---:|
| Sampling frequency | Approximately 1024 Hz |
| Welch segment length | 512 samples |
| Welch frequency-grid spacing | 2 Hz |
| FIR filter order | 500 |
| Compensated group delay | 0.2441 s |
| Noise variance before whitening | 1.0038 |
| Noise variance after whitening | 1.0334 |
| Approximate target variance | 1 |
| Target one-sided white-noise PSD | $1.9531\times10^{-3}\ \mathrm{Hz}^{-1}$ |

The whitened PSD fluctuated around the expected flat level, and the output variance was approximately 3.34% above the target. The original variance was also near one, demonstrating that variance alone does not establish whiteness: the frequency distribution of power must also be examined.

In the noise-relative spectrogram, I identified a curved upward-frequency track, most clearly visible around 7–9.5 s and reaching approximately 50–60 Hz. I interpreted this as evidence of an upward chirp. These approximate visual estimates do not establish its exact frequency law or true onset; the signal was specified to begin after 5 s, even though its earlier portion was less visible.

## 8. Conclusion

I approximately flattened the colored-noise PSD using a filter designed from the signal-free data. The noise-relative spectrogram made a coherent chirp-like feature easier to identify, while the signal remained difficult to distinguish in the time series. Whitening and improved display contrast both contributed to its visibility.

My PSD check used data from the same interval employed to design the filter, so I treat it as an internal consistency check rather than independent validation. The result also depends on the assumption that the noise spectrum remains representative after 5 s.

## Source Material

I followed the whitening task in Prof. Mohanty's *Lab Topic 2_2* and the FIR-design approach illustrated by `statgaussnoisegen.m`. The PSD floor, delay handling, boundary masks, and noise-relative display were implementation choices in my analysis script.
