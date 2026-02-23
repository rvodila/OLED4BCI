# Stimulus Folder Layout

This folder is organized by role so files are easier to find and maintain.

## Structure

- `assets/`: non-code resources used by experiments
  - `audio/`: generated or curated sound files
  - `codes/`: code sequences (`.mat`, `.npz`)
  - `images/`: still-image stimuli
  - `videos/`: video stimuli
- `scripts/`: runnable code
  - `matlab/`: Psychtoolbox and MATLAB scripts
  - `python/`: PsychoPy and Python scripts
- `docs/`: protocol notes and integration docs
- `archive/`: legacy artifacts and cached/generated files kept for reference

## Conventions

- Place new media files in `assets/*`, not next to scripts.
- Place new analysis or experiment code in `scripts/*`.
- Keep external/reference examples in `archive/` or `scripts/*/external/`.
- Avoid absolute paths in scripts; build paths relative to the script location.
