# Hybrid RIFT-cVEP Flicker Script (2026)

This folder contains a standalone PsychoPy script that ports the temporal flicker logic from `dualFlickerStimulus.m` into Python with explicit mode support.

## Files

- `stimulation_script_hybrid.py`: New experiment script with per-channel flicker modes.
- `metadata/mgold_61_6521.npz`: Local copy of the m-sequence codebook used by Python.
- `metadata/mgold_61_6521.mat`: Local copy of the MATLAB codebook source for traceability.

## Flicker Modes

The script accepts these modes per channel:

- `periodic` (MATLAB alias: `freq`): sinusoidal carrier.
- `alternating` (MATLAB alias: `alt`): frame-wise 0/1 alternation.
- `code`: pure m-sequence modulation.
- `hybrid`: bipolar code multiplied by bipolar periodic carrier.

## Notes

- Codebook loading defaults to `metadata/mgold_61_6521.npz` in this folder.
- If local metadata is missing, the script falls back to `project/stimulus/assets/codes/`.
- Shared temporal parameters (`frames_per_bit`, `ramp_len_frames`, `trial_taper_frames`) are configured under the `EXPERIMENT["flicker"]` section.
