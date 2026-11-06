# Part 3: Simulating Initial-LIGO Detector Noise

**Author:** Luis Alonso Hernández Galván  
**Course:** PHYS 5394 — Advanced Statistical Methods for Modern Astronomy  
**Instructor:** Prof. Soumya D. Mohanty  
**Laboratory:** Lab Topic 2_2 — Part 3

## 1. Objective

In this part, I implemented a simulation of stationary Gaussian strain noise using the initial-LIGO design sensitivity curve as the spectral target. I adapted the colored-noise generation procedure from Part 1, converted the supplied amplitude spectral density into the convention required by `statgaussnoisegen.m`, and included comparisons between the estimated and target spectra.

This report explains the implementation. Numerical output and plots from the Part 3 run have not yet been included, so I distinguish the expected behavior from verified simulation results.

## 2. Loading the Sensitivity Curve

```matlab
clear;
close all;
clc;
rng(5394);

sensitivityData = load('iLIGOSensitivity.txt','-ascii');
frequencyInput = sensitivityData(:,1);
asdInput = sensitivityData(:,2);

validRows = isfinite(frequencyInput) & isfinite(asdInput) ...
          & frequencyInput > 0 & asdInput > 0;
frequencyInput = frequencyInput(validRows);
asdInput = asdInput(validRows);

[frequencyInput,sortIndex] = sort(frequencyInput);
asdInput = asdInput(sortIndex);
[frequencyInput,uniqueIndex] = unique(frequencyInput);
asdInput = asdInput(uniqueIndex);
```

I initialized the session and fixed the random-number seed for reproducibility. I loaded the frequency values from the first column and the one-sided amplitude spectral density (ASD) from the second column.

I retained positive, finite values because the interpolation uses logarithms. I also sorted the frequencies and removed duplicates while preserving the corresponding ASD values. These operations prepared an ordered interpolation grid without changing the units.

For strain noise, the supplied curve represents

$$
\mathrm{ASD}(f)=\sqrt{S_n^{(1)}(f)}.
$$

| Quantity | Units |
|---|---|
| Strain | Dimensionless |
| One-sided PSD | $\mathrm{strain}^2/\mathrm{Hz}$ |
| ASD | $\mathrm{strain}/\sqrt{\mathrm{Hz}}$ |

## 3. Simulation Parameters and Spectral Modification

```matlab
sampFreq = 2048;
nSamples = 65536;
fltrOrdr = 2000;
lowCutoff = 50;
highCutoff = 700;

timeVec = (0:nSamples-1)/sampFreq;
freqVec = linspace(0,sampFreq/2,8193).';
evaluationFreq = min(max(freqVec,lowCutoff),highCutoff);

targetASD = 10.^interp1( ...
    log10(frequencyInput),log10(asdInput), ...
    log10(evaluationFreq),'linear');
```

I selected a sampling frequency of 2048 Hz and retained 65,536 samples, corresponding to 32 s. The Nyquist frequency was 1024 Hz, above the upper modification frequency of 700 Hz. The script checks that the input curve covers 50–700 Hz and that Nyquist exceeds 700 Hz before constructing the target.

Following the laboratory instructions, I held the ASD constant below 50 Hz and above 700 Hz:

$$
A_{\mathrm{target}}(f)=
\begin{cases}
A(50\ \mathrm{Hz}),&0\leq f<50\ \mathrm{Hz},\\
A(f),&50\leq f\leq700\ \mathrm{Hz},\\
A(700\ \mathrm{Hz}),&700\ \mathrm{Hz}<f\leq f_s/2.
\end{cases}
$$

Here $A(f)$ denotes the interpolated input ASD. Clamping `evaluationFreq` implements these constant extensions. This simplifies FIR design, especially by avoiding the steep low-frequency rise of the original curve. The modification retains noise outside 50–700 Hz at a constant spectral level; it does not remove those frequencies.

I interpolated linearly in log-frequency and log-ASD to represent the curve across its large dynamic range. The resulting design grid includes DC and Nyquist, which map to the normalized endpoints required by `fir2`. The 700 Hz boundary is a modeling choice for this exercise, not a universal inspiral termination frequency.

## 4. Converting ASD to the Generator's PSD Convention

```matlab
targetPSDOneSided = targetASD.^2;
targetPSDTwoSided = targetPSDOneSided/2;
```

I first squared the ASD to obtain the one-sided PSD. I then divided by two because Mohanty's generator expects the positive-frequency branch of a two-sided PSD:

$$
S_{\mathrm{target}}^{(2)}(f)
=\frac{S_{\mathrm{target}}^{(1)}(f)}{2}
=\frac{A_{\mathrm{target}}(f)^2}{2},
\qquad 0<f<f_s/2.
$$

For real noise, the full two-sided spectrum is even. DC and Nyquist are not doubled when constructing a discrete one-sided spectral estimate, so those endpoint bins require separate interpretation when comparing `pwelch` with the continuous target curve.

## 5. Generating the Strain Noise

```matlab
noiseWithTransient = statgaussnoisegen( ...
    nSamples + fltrOrdr, ...
    [freqVec,targetPSDTwoSided], ...
    fltrOrdr,sampFreq);

strainNoise = noiseWithTransient(fltrOrdr+1:end);
strainNoise = strainNoise(:);
```

I used a 2000th-order FIR filter, containing 2001 coefficients. The generator designs its magnitude to approximate the square root of the two-sided target PSD and applies it to independent Gaussian samples. Its additional $\sqrt{f_s}$ scaling ensures the intended PSD normalization per hertz.

I requested 2000 extra samples and discarded the first 2000 outputs to remove the FIR startup transient. I then converted the retained record into a column vector for subsequent analysis.

Linear filtering preserves Gaussianity while introducing the correlations associated with the target spectrum. I kept the physical strain scale rather than normalizing the output to unit variance, since that scale is determined by the sensitivity curve.

## 6. Estimating and Plotting the Spectrum

```matlab
windowLength = 2048;
welchWindow = hamming(windowLength);
nOverlap = windowLength/2;
nFFT = windowLength;

[estimatedPSD,estimatedFreq] = pwelch( ...
    strainNoise,welchWindow,nOverlap,nFFT,sampFreq);
estimatedASD = sqrt(estimatedPSD);

positiveFreq = estimatedFreq > 0;
positiveTarget = freqVec > 0;
```

I estimated the one-sided PSD with one-second Hamming windows, 50% overlap, and a 2048-point FFT. The returned frequency-grid spacing is

$$
\Delta f=\frac{2048}{2048}=1\ \mathrm{Hz}.
$$

This is the bin spacing; the effective spectral resolution also depends on the window. Taking the square root of the PSD estimate provides an ASD estimate for comparison with the input sensitivity curve.

I included three logarithmic spectral plots: the original versus modified ASD, the estimated versus target one-sided PSD, and the estimated versus target ASD. I excluded DC from these plots because logarithmic frequency axes cannot display zero. Axis labels identify the PSD and ASD units, while reference lines mark 50 and 700 Hz.

I also included full-record and first-0.1-second time-series views. These plots show the simulated strain amplitudes and local temporal structure, but the PSD overlay remains the main check of the prescribed spectral shape.

## 7. Numerical Validation and Expected Behavior

```matlab
targetVariance = trapz(freqVec,targetPSDOneSided);
sampleVariance = var(strainNoise,1);
estimatedVariance = trapz(estimatedFreq,estimatedPSD);

varianceDifference = ...
    100*(sampleVariance-targetVariance)/targetVariance;
```

I included three estimates of the noise power: the integrated target PSD, the sample variance, and the integrated Welch PSD. For the modeled zero-mean process,

$$
\sigma_n^2=\int_0^{f_s/2}S_n^{(1)}(f)\,df.
$$

`trapz` approximates this integral numerically, while `var(strainNoise,1)` calculates the mean-subtracted sample variance with denominator $N$. The printed diagnostics also report the sample mean, standard deviation, duration, filter order, and frequency spacing. Variances have units of strain squared, and the standard deviation has units of strain.

I expect the Welch estimate to follow the modified target, including its sensitivity minimum and constant low- and high-frequency regions. Exact agreement is not expected because of finite-realization fluctuations, Welch smoothing, and the FIR approximation. Persistent spectral distortions would motivate checking the actual filter response and its design parameters rather than attributing every discrepancy to randomness.

## 8. Conclusion

I implemented the conversion from an initial-LIGO ASD curve to a two-sided PSD suitable for colored Gaussian noise generation. The analysis preserves the physical strain scale and provides spectral overlays and variance checks for evaluating the generated realization.

This model represents stationary Gaussian noise with an approximate design spectrum. It does not reproduce the glitches, spectral lines, or changing statistics that can occur in measured detector data. A numerical conclusion about agreement with the target requires the output of the Part 3 simulation.

## Source Material

I followed the “Simulating LIGO noise” section of Prof. Mohanty's *Lab Topic 2_2* and the supplied `colGaussNoiseDemo.m` / `statgaussnoisegen.m` workflow. Logarithmic interpolation, transient removal, and the additional ASD and variance comparisons were implementation choices in my script.
