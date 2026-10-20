# Feedback: Lab Topic 1_1

- Lab specification: `LABS/Lab Topic 1_1.pdf`
- Submission folder: `PHYS5394-LuisH`
- Review date: 2026-10-20
- Review method: Static inspection only; student and instructor code were not executed.

## Overview

The Lab 1 folder contains functions and test scripts for all eight assigned signal types. The scripts generate time plots and positive-frequency FFT-magnitude plots, and `Report.pdf` includes descriptive observations for the spectra. The most useful additions are a sampling-rate convergence comparison and the requested sweep of sinusoid frequencies in steps smaller than $1/T$.

## Submission identification and files reviewed

The matching source is in `Labs/Lab 1/`. Reviewed the eight `crcbgen*.m` functions and corresponding `testcrcbgen*.m` scripts for quadratic chirp, sinusoid, linear chirp, sine-Gaussian, FM, AM, AM-FM, and linear transient chirp, plus `Labs/Lab 1/Report.pdf`. The report is a ten-page exported script/report containing signal generation, plots, and explanatory text; its displayed outputs were treated as submitted artifacts and not independently reproduced.

## Instructor reference material

- `DATASCIENCE_COURSE/SIGNALS/crcbgenqcsig.m`: reference quadratic-chirp equation and normalization.
- `DATASCIENCE_COURSE/SIGNALS/testcrcbgenqcsig.m`: reference sampling and periodogram workflow.

## Lab objectives and requirements

| Requirement | Status | Code evidence and feedback |
|---|---|---|
| Generate the eight specified signal forms (Part 1: signal definitions) | Addressed in code | The function files implement each assigned form, including the transient interval and parameterized phases. On inspection, the formulas and SNR normalization appear internally consistent with the descriptions. |
| Choose a sufficiently fine sampling interval by reducing it until plots stop changing (Part 1: “Plots”) | Partially addressed | The test scripts choose a fixed sampling frequency from the maximum-frequency estimate (generally five times the Nyquist frequency); the report uses a fixed multiple of the maximum frequency. No coarse-to-fine comparison is present in the reviewed Lab 1 files. Add a small comparison of sampled plots at increasing sample rates and state why the selected interval is adequate. |
| Plot the magnitude of the FFT against positive-frequency values, using spacing $1/(N\Delta)$ (Part 2: “FFT”) | Addressed in code | The test scripts retain `floor(nSamples/2)+1` FFT values and calculate spacing as `samplFreq/nSamples`, equivalent to `1/(N*Delta)` for the constructed grid. |
| Interpret periodograms, including finite-duration sinusoid behavior and frequency steps below $1/T$ (Part 2: “Frequencies in a DFT”) | Partially addressed | `Report.pdf` provides a short signal-specific interpretation, and the sinusoid test shows one fixed frequency. No sub-$1/T$ frequency sweep was found, and the report describes a narrow sinusoid peak rather than discussing its finite-window sinc-like shape. Add the sweep and a brief comparison of its spectral shape. |

## Code correctness and suggested improvements

- **Location:** `Labs/Lab 1/testcrcbgenqcsig.m`, time-vector setup. **Observation:** This script includes both 0 and the one-second endpoint, unlike the half-open grids used in the other tests. Its FFT axis is computed consistently for its actual sample count, but the extra endpoint means the observation contains $N$ samples over a span of $N-1$ intervals. **Likely consequence:** The quadratic-chirp test has a slightly different DFT spacing and sampling convention than the other one-second tests. **Suggested improvement:** Use a consistent half-open grid, such as `0:samplIntrvl:(signalDuration-samplIntrvl)`, across the test scripts. This is a consistency concern clear from inspection; its effect on plots was not tested.

No other obvious syntax or formula issue was identified in the reviewed implementations. This does not establish runtime or numerical correctness.

## Documentation feedback

### What works well

Function headers document the parameter vectors and units. The test scripts are divided into readable sections, label axes, and explain key maximum-frequency calculations. The report includes an interpretation for each signal type.

### Suggested improvements

Document the sample-count convention once and apply it consistently. Expand the sinusoid interpretation to explain finite-duration spectral leakage and the requested sub-$1/T$ frequency changes.

## Recommended next steps

1. Add a coarse-to-fine sampling comparison and explain the chosen interval.
2. Add a sinusoid sweep with frequency increments below $1/T$, then interpret how its FFT magnitude changes.
3. Make the quadratic-chirp sample grid consistent with the other one-second scripts.

## Review limitations

Student and instructor code were not executed. Runtime behavior and numerical outputs, including the figures embedded in `Report.pdf`, were not independently verified. Conclusions are based on static inspection of the listed source and report.