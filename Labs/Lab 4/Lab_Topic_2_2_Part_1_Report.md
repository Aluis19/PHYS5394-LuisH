# Colored Gaussian Noise: Generation and Spectral Analysis

**Author:** Luis Alonso Hernández Galván  
**Course:** PHYS 5394 — Advanced Statistical Methods for Modern Astronomy  
**Instructor:** Prof. Soumya D. Mohanty  
**Laboratory:** Lab Topic 2_2 — Part 1

## 1. Objective

In this part of the laboratory, I generated colored Gaussian noise with a prescribed power spectral density (PSD). I compared its estimated PSD with the target spectrum, examined the effect of increasing the record length, and inspected both its time series and amplitude histogram. My central objective was to distinguish Gaussianity, which describes the joint probability distribution, from spectral color, which describes how power and correlations are distributed in frequency and time.

I used MATLAB and the supplied function `statgaussnoisegen.m`. The accompanying script extends Mohanty's demonstration by adding labeled plots, comparisons among three record lengths, removal of the filter startup transient, and numerical checks of the PSD normalization.

## 2. Definitions and conventions

I denote the sampling frequency by $f_s$, the number of retained samples by $N$, and the noise samples by $n[k]$. The sample times are

$$
t_k=\frac{k}{f_s},\qquad k=0,\ldots,N-1.
$$

For a real, zero-mean stationary process, I use $S_n^{(2)}(f)$ for the two-sided PSD and $S_n^{(1)}(f)$ for the one-sided PSD. Away from DC and Nyquist,

$$
S_n^{(1)}(f)=2S_n^{(2)}(f),\qquad 0<f<f_s/2.
$$

Both conventions account for the same variance:

$$
\sigma_n^2
=\int_{-f_s/2}^{f_s/2}S_n^{(2)}(f)\,df
=\int_0^{f_s/2}S_n^{(1)}(f)\,df.
$$

The PSD therefore has units of data units squared per hertz. Integrating over frequency gives data units squared, as required for a variance.

## 3. Initialization and simulation parameters

```matlab
clear;
close all;
clc;
rng(5394);

sampFreq = 1024;
sampleCounts = 16384*[1, 2, 4];
fltrOrdr = 500;
```

I used `clear` to remove workspace variables, `close all` to close existing figures, and `clc` to clear the Command Window. These commands prepare the MATLAB session; they do not define the noise statistics.

I set the random-number seed with `rng(5394)` so that I can reproduce the simulation under the same MATLAB random-number configuration. I set the seed once, before the loop. Consequently, the three records use successive random draws rather than being copies or prefixes of a single realization.

My sampling frequency was $1024$ Hz, giving a Nyquist frequency of $512$ Hz. I selected three record lengths:

| Record | Retained samples | Observation duration $N/f_s$ |
|---|---:|---:|
| Baseline | 16,384 | 16 s |
| Twice the baseline | 32,768 | 32 s |
| Four times the baseline | 65,536 | 64 s |

Although the final sample occurs at $(N-1)/f_s$, I use $N/f_s$ as the record duration in the spectral-analysis sense.

The FIR filter order was 500, corresponding to 501 coefficients. This controls the approximation to the desired spectral shape; it is separate from both the record length and the Welch window length.

## 4. Construction of the target PSD

```matlab
targetPSD = @(f) (f >= 100 & f <= 300) ...
                    .* (f-100) .* (300-f)/10000;

freqVec = 0:0.1:sampFreq/2;
psdVec = targetPSD(freqVec);

targetOneSided = 2*psdVec;
targetOneSided([1,end]) = psdVec([1,end]);
```

I specified the positive-frequency branch of the two-sided target as

$$
S_{\mathrm{target}}^{(2)}(f)=
\begin{cases}
\dfrac{(f-100)(300-f)}{10000},&100\leq f\leq300\ \mathrm{Hz},\\
0,&\text{otherwise},
\end{cases}
\qquad f\geq0.
$$

The full two-sided spectrum has an even extension to negative frequencies. The `/10000` factor follows the demo normalization used in my script; the displayed polynomial in the lab slides omits this factor. This changes the power scale but preserves the spectral shape.

The expression `@(f)` defines an anonymous function. The logical condition selects the 100–300 Hz band. Element-wise multiplication, `.*`, evaluates the polynomial separately at every frequency in the input array.

At 200 Hz, the two-sided target reaches its maximum:

$$
S_{\mathrm{target}}^{(2)}(200)
=\frac{100\times100}{10000}=1.
$$

The corresponding one-sided peak is 2. I made this conversion explicitly because `pwelch` returns a one-sided PSD for the real-valued records used here. I preserved the endpoint convention at DC and Nyquist, where bins are not doubled. Both endpoints have zero target power in this example.

The 0.1 Hz grid describes the filter-design target and the plotted reference. It does not imply that my Welch estimate resolves spectral features separated by 0.1 Hz.

## 5. Fixed Welch settings

```matlab
windowLength = 256;
welchWindow = hamming(windowLength);
nOverlap = windowLength/2;
nFFT = windowLength;
noiseRecords = cell(size(sampleCounts));
```

I used 256-sample Hamming windows with 128-sample overlap and a 256-point FFT. Each segment spans $256/1024=0.25$ s. The returned frequency-grid spacing is

$$
\Delta f=\frac{f_s}{N_{\mathrm{FFT}}}
=\frac{1024}{256}=4\ \mathrm{Hz}.
$$

This is the spacing between reported spectral bins. The ability to distinguish nearby spectral features also depends on the window's spectral width; it is not determined by the bin spacing alone.

Welch's method averages windowed segment periodograms. I kept these settings fixed while increasing the record length so that the principal change was the number of averages. The Hamming window reduces spectral leakage, with the usual tradeoff of broadening the spectral response.

I stored each noise realization in a cell array because the three records have different lengths.

## 6. Noise generation and removal of the startup transient

```matlab
figure('Name','PSD comparison','Color','w');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

for k = 1:numel(sampleCounts)
    nSamples = sampleCounts(k);

    noiseWithTransient = statgaussnoisegen( ...
        nSamples + fltrOrdr, ...
        [freqVec(:), psdVec(:)], ...
        fltrOrdr, sampFreq);

    outNoise = noiseWithTransient(fltrOrdr+1:end);
    noiseRecords{k} = outNoise;
```

I generated 500 extra samples for each record and discarded the first 500 filtered outputs. The FIR filter initially lacks the preceding input samples needed for its complete convolution. Removing one filter order of output samples eliminates that startup region, leaving exactly `nSamples` retained samples.

The expressions `freqVec(:)` and `psdVec(:)` convert the two arrays to columns. Combining them produces a two-column matrix: frequency in the first column and target two-sided PSD in the second column.

### 6.1 Operations inside the supplied generator

The essential operations in `statgaussnoisegen.m` are

```matlab
freqVec = psdVals(:,1);
sqrtPSD = sqrt(psdVals(:,2));
b = fir2(fltrOrdr,freqVec/(sampFreq/2),sqrtPSD);
inNoise = randn(1,nSamples);
outNoise = sqrt(sampFreq)*fftfilt(b,inNoise);
```

I interpret these operations as spectral shaping of independent standard Gaussian samples. `fir2` designs a filter whose magnitude approximates the square root of the target PSD. Dividing frequency by $f_s/2$ maps DC and Nyquist to the normalized endpoints 0 and 1 required by `fir2`.

The input samples produced by `randn` have variance one. Their two-sided PSD per hertz is $1/f_s$ across the sampled band. If $B(f)$ denotes the designed FIR frequency response, then the additional amplitude factor $\sqrt{f_s}$ gives

$$
S_{\mathrm{out}}^{(2)}(f)
=f_s|B(f)|^2\frac{1}{f_s}
=|B(f)|^2
\approx S_{\mathrm{target}}^{(2)}(f).
$$

The approximation reflects the finite filter order. `fftfilt` performs the FIR convolution using an FFT-based algorithm.

The filter preserves Gaussianity because its output samples are linear combinations of jointly Gaussian input samples. It introduces correlations because different output samples depend on overlapping sets of inputs.

## 7. PSD estimation and comparison among record lengths

```matlab
    [pxx,f] = pwelch(outNoise,welchWindow, ...
                    nOverlap,nFFT,sampFreq);

    nexttile;
    plot(f,pxx,'b','LineWidth',1);
    hold on;
    plot(freqVec,targetOneSided,'k--','LineWidth',1.5);

    xlabel('Frequency (Hz)');
    ylabel('PSD (data units^2/Hz)');
    title(sprintf('N = %d; duration = %.0f s', ...
                  nSamples,nSamples/sampFreq));
    legend('Welch estimate','One-sided target','Location','northeast');
    xlim([0 sampFreq/2]);
    ylim([0 3.5]);
    grid on;
end
sgtitle('Colored Gaussian noise: effect of record length');
```

I compared the blue Welch estimate with the dashed one-sided target. `nexttile` places each record in a separate panel, while the shared axis limits support visual comparison. `sprintf` inserts the sample count and duration into each title. The vertical limit is a display choice and does not change the estimated PSD.

For segment length $L=256$ and overlap $D=128$, the number of complete segments is

$$
K=1+\left\lfloor\frac{N-L}{L-D}\right\rfloor.
$$

This gives 127, 255, and 511 segments for the three record lengths. Overlapping segments are correlated, so these counts are not counts of independent averages. Nevertheless, longer records provide additional averaging and generally reduce fluctuations in the estimate.

I observed that the estimated spectra followed the target band and became smoother overall for longer records. This improvement concerns statistical variability rather than a change in the 4 Hz spectral grid. Since the records are different realizations, individual frequency bins need not become closer to the target monotonically.

## 8. Examination of the time series

```matlab
outNoise = noiseRecords{1};
nSamples = sampleCounts(1);
timeVec = (0:nSamples-1)/sampFreq;

figure('Name','Noise time series','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(timeVec,outNoise,'b');
xlabel('Time (s)');
ylabel('Amplitude (data units)');
title('Colored Gaussian noise: full baseline record');
xlim([0 timeVec(end)]);
grid on;

nexttile;
plot(timeVec,outNoise,'b.-');
xlabel('Time (s)');
ylabel('Amplitude (data units)');
title('Colored Gaussian noise: first 0.1 s');
xlim([0 0.1]);
grid on;
```

I retrieved the baseline record with braces, `noiseRecords{1}`, which extract the numeric array from the cell. I constructed its time vector by dividing sample indices by the sampling frequency.

The full-record view shows the overall amplitude behavior, while the short view reveals the local oscillatory structure. The `'.-'` style marks individual samples and connects them for readability. Restricting `xlim` changes only the displayed interval.

The oscillatory appearance is consistent with power concentrated between 100 and 300 Hz. I did not inject a deterministic sinusoid: the oscillations belong to the random noise realization.

## 9. Comparison with white Gaussian noise

```matlab
whiteNoise = randn(size(outNoise));
whiteNoise = whiteNoise - mean(whiteNoise);
whiteNoise = whiteNoise/std(whiteNoise)*std(outNoise);

commonLimit = 1.1*max(abs([outNoise(:); whiteNoise(:)]));

figure('Name','White versus colored noise','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(timeVec,whiteNoise,'r.-');
xlabel('Time (s)');
ylabel('Amplitude (data units)');
title('White Gaussian noise: first 0.1 s');
xlim([0 0.1]);
ylim([-commonLimit commonLimit]);
grid on;

nexttile;
plot(timeVec,outNoise,'b.-');
xlabel('Time (s)');
ylabel('Amplitude (data units)');
title('Colored Gaussian noise: first 0.1 s');
xlim([0 0.1]);
ylim([-commonLimit commonLimit]);
grid on;
```

I generated a separate white Gaussian reference, removed its sample mean, and scaled it to have the same sample standard deviation as the colored record. I used this adjustment only to compare temporal structure at comparable amplitude scales.

The reference was originally an independent Gaussian sequence. Centering and scaling with statistics estimated from the same finite record impose constraints on it, so the displayed normalized vector is not strictly an unconstrained iid Gaussian draw. This distinction does not affect its purpose as a visual reference here.

I concatenated both records as columns when computing `commonLimit`, ensuring a scalar maximum. The factor 1.1 provides plotting margin. Identical time and amplitude limits prevent axis scaling from exaggerating the visual difference.

## 10. Amplitude histogram and Gaussian PDF

```matlab
sampleMean = mean(outNoise);
sampleStd = std(outNoise);

amplitudeGrid = linspace(sampleMean-4*sampleStd, ...
                         sampleMean+4*sampleStd,500);

gaussianPDF = exp(-0.5*((amplitudeGrid-sampleMean)/sampleStd).^2) ...
              /(sqrt(2*pi)*sampleStd);

figure('Name','Amplitude distribution','Color','w');
histogram(outNoise,60,'Normalization','pdf', ...
          'FaceColor',[0.3 0.6 0.9],'EdgeColor','none');
hold on;
plot(amplitudeGrid,gaussianPDF,'k','LineWidth',1.8);
xlabel('Noise amplitude (data units)');
ylabel('Probability density (1/data units)');
title('Colored noise: amplitude distribution');
legend('Normalized histogram','Gaussian PDF using sample mean and std');
grid on;
```

I estimated the mean and standard deviation from the baseline realization. MATLAB's default `std` uses the $N-1$ denominator. I then evaluated

$$
p(a)=\frac{1}{\sqrt{2\pi}\,\widehat\sigma}
\exp\left[-\frac{(a-\widehat\mu)^2}{2\widehat\sigma^2}\right]
$$

on 500 amplitude points spanning four measured standard deviations on each side of the measured mean. Here $a$ is an amplitude variable, not a time index.

I selected `'Normalization','pdf'` so that histogram area represents probability and can be compared directly with the Gaussian density. Its vertical units are inverse data units.

The histogram was consistent with a Gaussian marginal distribution. This visual comparison is not a proof of joint Gaussianity or independence. Gaussianity follows from the linear filtering construction; the histogram provides a useful check of its amplitude distribution. A Gaussian-looking histogram alone cannot distinguish white from colored noise.

## 11. Variance checks and numerical output

```matlab
targetVariance = trapz(freqVec,targetOneSided);

[baselinePSD,baselineFreq] = pwelch( ...
    outNoise,welchWindow,nOverlap,nFFT,sampFreq);
estimatedVariance = trapz(baselineFreq,baselinePSD);

fprintf('\nBaseline record:\n');
fprintf('Samples                      : %d\n',nSamples);
fprintf('Duration                     : %.2f s\n',nSamples/sampFreq);
fprintf('Sample mean                  : %.4f\n',sampleMean);
fprintf('Sample variance              : %.4f\n',var(outNoise,1));
fprintf('Integrated Welch PSD         : %.4f\n',estimatedVariance);
fprintf('Integrated target PSD        : %.4f\n',targetVariance);
fprintf('Welch frequency spacing      : %.2f Hz\n',sampFreq/nFFT);
```

I used `trapz` to integrate the target and estimated one-sided PSDs numerically. I recomputed the baseline PSD because the loop variable `pxx` otherwise contains the estimate for the final, longest record.

The expression `var(outNoise,1)` uses denominator $N$ and removes the sample mean. This differs slightly from the default $N-1$ normalization used by `std`. Since the process is modeled as zero mean, its variance equals its mean-square power; the very small measured mean makes the practical distinction negligible here.

The `fprintf` commands report scalar diagnostics. The format `%d` prints an integer, `%.4f` prints four decimal places, and `\n` starts a new line.

### 11.1 Analytical target variance

I verified the numerical target integral analytically:

$$
\sigma_{\mathrm{target}}^2
=\frac{2}{10000}\int_{100}^{300}(f-100)(300-f)\,df.
$$

With $u=f-100$, the integration limits become 0 and 200, and the polynomial becomes $u(200-u)$. Therefore,

$$
\begin{aligned}
\sigma_{\mathrm{target}}^2
&=\frac{2}{10000}\int_0^{200}(200u-u^2)\,du\\
&=\frac{2}{10000}\left[100u^2-\frac{u^3}{3}\right]_0^{200}\\
&=\frac{2}{10000}\left(4{,}000{,}000-\frac{8{,}000{,}000}{3}\right)\\
&=266.6667\ \text{data units}^2.
\end{aligned}
$$

### 11.2 Results from my baseline run

| Quantity | Measured or calculated value |
|---|---:|
| Number of samples | 16,384 |
| Observation duration | 16.00 s |
| Sample mean | −0.0016 |
| Sample variance | 271.6181 |
| Integrated Welch PSD | 272.2903 |
| Integrated target PSD | 266.6666 |
| Welch frequency-grid spacing | 4.00 Hz |

The sample variance exceeded the target by approximately 1.86%, and the integrated Welch PSD exceeded it by approximately 2.11%. The integrated estimate and sample variance differed by approximately 0.25%, relative to the sample variance. The measured standard deviation was approximately 16.48 data units, so the sample mean was small relative to the noise amplitude scale.

I expect differences because a finite random realization does not reproduce ensemble statistics exactly, the FIR filter approximates the target, and the Welch estimate includes windowing and finite frequency sampling. Numerical integration with `trapz` is also an approximate PSD integral. Agreement here supports the scaling without requiring exact equality.

## 12. Discussion and conclusion

I generated colored Gaussian noise whose estimated PSD followed the prescribed 100–300 Hz target. Increasing the record length with fixed Welch settings provided more averages and reduced the overall spectral fluctuations, while the frequency-grid spacing remained unchanged.

The short time-series plots showed temporal structure associated with the selected band, and the amplitude histogram remained consistent with a Gaussian density. I interpret these observations as complementary: the PSD describes second-order spectral behavior, while the histogram describes the marginal amplitude distribution.

Finally, my sample variance agreed closely with the integrated Welch PSD and was within approximately 1.86% of the target variance. This verified the practical importance of converting the generator's two-sided target to the one-sided convention used by `pwelch`.

This exercise establishes the spectral-shaping framework I use in the next part to construct a whitening filter. It also illustrates why Gaussian detector noise can possess a strongly frequency-dependent spectrum.

## Source material

I followed Prof. Mohanty's *Lab Topic 2_2*, particularly its sections on one-sided and two-sided PSDs, colored Gaussian noise generation, and Welch estimation. The implementation follows the supplied `colGaussNoiseDemo.m` and `statgaussnoisegen.m` approach as reproduced in the working Part 1 script. The transient removal, matched-amplitude reference, and numerical variance checks are additions in my analysis script.
