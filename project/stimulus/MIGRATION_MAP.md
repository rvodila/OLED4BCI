# Migration Map

This file tracks the cleanup move from the old layout to the new one.

## Directory Moves

- `project/stimulus/audio_files/` -> `project/stimulus/assets/audio/`
- `project/stimulus/codes/` -> `project/stimulus/assets/codes/`
- `project/stimulus/images/` -> `project/stimulus/assets/images/`
- `project/stimulus/videos/` -> `project/stimulus/assets/videos/`
- `project/stimulus/video_tagging/` -> `project/stimulus/scripts/matlab/video_tagging/`
- `project/stimulus/extern/` -> split into:
  - `project/stimulus/scripts/matlab/external/`
  - `project/stimulus/scripts/python/external/`
  - `project/stimulus/archive/pycache/extern__pycache/`
- `project/stimulus/python/` -> `project/stimulus/scripts/python/`

## File Moves

- `project/stimulus/LSL Integration.md` -> `project/stimulus/docs/lsl_integration.md`
- `project/stimulus/dualFlickerStimulus.m` -> `project/stimulus/scripts/matlab/core/dualFlickerStimulus.m`
- `project/stimulus/minimal_viable_flicker.m` -> `project/stimulus/scripts/matlab/core/minimal_viable_flicker.m`
- `project/stimulus/minimal_viable_flicker.asv` -> `project/stimulus/scripts/matlab/legacy/minimal_viable_flicker.asv`
- `project/stimulus/python/2026/stimulation_script.py` -> `project/stimulus/scripts/python/experiments/2026/stimulation_script.py`
- `project/stimulus/python/2025/flicker_grow.py` -> `project/stimulus/scripts/python/experiments/2025/flicker_grow.py`
- `project/stimulus/python/2025/flicker_psp.ipynb` -> `project/stimulus/scripts/python/experiments/2025/flicker_psp.ipynb`
- `waveforms.py` -> `project/stimulus/scripts/python/utils/waveforms.py`

## Notes

- No source files or assets were intentionally deleted.
- Script paths were updated to resolve assets relative to script location.
