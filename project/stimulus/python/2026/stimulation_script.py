#!/usr/bin/env python
# -*- coding: utf-8 -*-

from psychopy import visual, core, event
import math
import random

# =========================
# Hyperparameters
# =========================
SUBJECT_ID = 42
RANDOM_SEED = SUBJECT_ID
SHOW_DIAGNOSTICS = True

EXPERIMENT = {
    "blocks": {"n_rift": 4, "n_localizer": 2},
    "timing": {
        "rift_trial_duration_s": 5.0,
        "localizer_segment_s": 1.0,
        "localizer_n_cycles": 1,
        "isi_s": 0.3,
    },
    "fixation": {"min_s": 0.5, "max_s": 1.5},
    "frequencies": {
        "rift_hz":      [3.0, 4.0, 5.0, 6.0],
        "localizer_hz": [2.0, 2.0, 2.0, 2.0],
    },
    "localizer": {"order": (0, 1, 2, 3), "others": "off"},
}

# Deterministic randomness per subject (fixation jitter + anything else using global random)
random.seed(RANDOM_SEED)

# Pull values from config into simple names (minimal, avoids KeyErrors)
RIFT_TRIAL_DURATION_S = EXPERIMENT["timing"]["rift_trial_duration_s"]
LOCALIZER_SEGMENT_S   = EXPERIMENT["timing"]["localizer_segment_s"]
LOCALIZER_N_CYCLES    = EXPERIMENT["timing"]["localizer_n_cycles"]
ISI_S                 = EXPERIMENT["timing"]["isi_s"]

FIX_MIN_S = EXPERIMENT["fixation"]["min_s"]
FIX_MAX_S = EXPERIMENT["fixation"]["max_s"]

RIFT_FREQS_HZ      = EXPERIMENT["frequencies"]["rift_hz"]
LOCALIZER_FREQS_HZ = EXPERIMENT["frequencies"]["localizer_hz"]

LOCALIZER_ORDER  = EXPERIMENT["localizer"]["order"]
LOCALIZER_OTHERS = EXPERIMENT["localizer"]["others"]

# Window
FULLSCREEN = True
UNITS = "height"
BG_COLOR = [0, 0, 0]

# Layout
OFFSET = 0.25
STIM_SIZE = 0.25
QUAD_POSITIONS = [(-OFFSET, +OFFSET), (+OFFSET, +OFFSET), (-OFFSET, -OFFSET), (+OFFSET, -OFFSET)]

# Rect colors
RECT_BG_COLOR = [1, 1, 1]
RECT_FG_COLOR = [0.2, 0.2, 0.2]
FG_OPACITY_ON = 0.5
FG_OPACITY_OFF = 0.0
OPACITY_MIN_255 = 20
OPACITY_MAX_255 = 230
OPACITY_MIN = OPACITY_MIN_255 / 255.0
OPACITY_MAX = OPACITY_MAX_255 / 255.0

QUIT_KEYS = ["escape"]

print("==== RUN INFO ====")
print(f"Subject ID: {SUBJECT_ID}")
print(f"Random seed: {RANDOM_SEED}")
print("==================")


# =========================
# Helpers
# =========================
def build_block_sequence(exp_cfg, seed):
    rng = random.Random(seed)  # deterministic order per subject
    seq = (["localizer"] * exp_cfg["blocks"]["n_localizer"] +
           ["rift"] * exp_cfg["blocks"]["n_rift"])
    rng.shuffle(seq)
    return seq

def check_quit():
    if event.getKeys(keyList=QUIT_KEYS):
        safe_quit()


def seconds_to_frames(seconds, fps):
    return int(round(seconds * fps))

def sample_fix_seconds():
    return random.uniform(FIX_MIN_S, FIX_MAX_S)

class SineFlicker:
    __slots__ = ("phase", "inc")
    def __init__(self, freq_hz, fps, phase0=0.0):
        self.phase = phase0
        self.inc = (2.0 * math.pi * float(freq_hz)) / float(fps) if freq_hz > 0 else 0.0
    def step_value(self):
        p = self.phase + self.inc
        if p >= (2.0 * math.pi):
            p -= (2.0 * math.pi)
        self.phase = p
        return math.sin(p)
def set_diag(block_idx=None, block_type=None, phase=None, extra=None):
    if not SHOW_DIAGNOSTICS:
        return

    lines = []
    if block_idx is not None:
        lines.append(f"Block: {block_idx}")
    if block_type is not None:
        lines.append(f"Type: {block_type}")
    if phase is not None:
        lines.append(f"Phase: {phase}")
    if extra is not None:
        lines.append(str(extra))

    diag_text.text = "\n".join(lines)


# =========================
# Build Window + Stimuli
# =========================
win = visual.Window(
    fullscr=FULLSCREEN,
    units=UNITS,
    color=BG_COLOR,
    colorSpace="rgb",
    waitBlanking=True,
)

fps = win.getActualFrameRate(nIdentical=20, nMaxFrames=200, nWarmUpFrames=10)
if fps is None:
    fps = 60.0
fps = float(fps)
# assert int(round(fps)) == 480, f"Expected 480 Hz refresh, got {fps}"
print(f"Measured FPS: {fps}")

diag_text = visual.TextStim(
    win,
    text="",
    pos=(0.48, 0.48),   # upper-right in height units
    height=0.03,
    anchorHoriz="right",
    anchorVert="top",
    color=[1, 1, 0],    # yellow for visibility
    colorSpace="rgb",
)
diag_text.autoDraw = True   # <-- IMPORTANT
diag_text.opacity = 1.0
fixation = visual.TextStim(win, text="+", height=0.05, color=[1, 1, 1], colorSpace="rgb")

bg_rects, fg_rects = [], []
for pos in QUAD_POSITIONS:
    bg = visual.Rect(win, width=STIM_SIZE, height=STIM_SIZE, pos=pos, units=UNITS,
                     fillColor=RECT_BG_COLOR, lineColor=RECT_BG_COLOR, opacity=1.0)
    fg = visual.Rect(win, width=STIM_SIZE, height=STIM_SIZE, pos=pos, units=UNITS,
                     fillColor=RECT_FG_COLOR, lineColor=RECT_FG_COLOR, opacity=FG_OPACITY_OFF)

    bg.autoDraw = True
    fg.autoDraw = True
    bg_rects.append(bg)
    fg_rects.append(fg)



# =========================
# Routines
# =========================
def run_fixation_seconds(duration_s):
    # Keep backgrounds visible; turn flicker off during fixation
    for bg in bg_rects: bg.autoDraw = True
    for fg in fg_rects:
        fg.autoDraw = True
        fg.opacity = FG_OPACITY_OFF

    n_frames = seconds_to_frames(duration_s, fps)
    for _ in range(n_frames):
        fixation.draw()
        if SHOW_DIAGNOSTICS:
            diag_text.draw()    
        win.flip()
        check_quit()

def run_isi_seconds(duration_s):
    # Same as fixation but without drawing the "+"
    for bg in bg_rects: bg.autoDraw = True
    for fg in fg_rects:
        fg.autoDraw = True
        fg.opacity = FG_OPACITY_OFF

    n_frames = seconds_to_frames(duration_s, fps)
    for _ in range(n_frames):
        win.flip()
        check_quit()

def run_trial_rift(duration_s, freqs_hz):
    n_frames = seconds_to_frames(duration_s, fps)
    phases = [SineFlicker(freqs_hz[i], fps, phase0=0.0) for i in range(4)]

    for _ in range(n_frames):
        fg_rects[0].opacity = OPACITY_MIN + ((phases[0].step_value() + 1.0) * 0.5) * (OPACITY_MAX - OPACITY_MIN)
        fg_rects[1].opacity = OPACITY_MIN + ((phases[1].step_value() + 1.0) * 0.5) * (OPACITY_MAX - OPACITY_MIN)
        fg_rects[2].opacity = OPACITY_MIN + ((phases[2].step_value() + 1.0) * 0.5) * (OPACITY_MAX - OPACITY_MIN)
        fg_rects[3].opacity = OPACITY_MIN + ((phases[3].step_value() + 1.0) * 0.5) * (OPACITY_MAX - OPACITY_MIN)

        fixation.draw()
        win.flip()
        check_quit()

def run_trial_localizer(segment_s, n_cycles, freqs_hz, order=(0,1,2,3), others="off"):
    seg_frames = seconds_to_frames(segment_s, fps)
    phases = [SineFlicker(freqs_hz[i], fps, phase0=0.0) for i in range(4)]
    non_target_opacity = FG_OPACITY_OFF if others == "off" else FG_OPACITY_ON

    for _c in range(n_cycles):
        for target in order:
            for _ in range(seg_frames):
                fg_rects[0].opacity = non_target_opacity
                fg_rects[1].opacity = non_target_opacity
                fg_rects[2].opacity = non_target_opacity
                fg_rects[3].opacity = non_target_opacity

                fg_rects[target].opacity = OPACITY_MIN + ((phases[target].step_value() + 1.0) * 0.5) * (OPACITY_MAX - OPACITY_MIN)

                fixation.draw()
                win.flip()
                check_quit()

import sys

def safe_quit():
    """Close PsychoPy cleanly so the window doesn't hang and the kernel doesn't crash."""
    try:
        event.clearEvents()
    except Exception:
        pass

    try:
        if win is not None:
            win.close()
    except Exception:
        pass

    try:
        core.quit()   # raises SystemExit
    except Exception:
        raise SystemExit

# =========================
# Main
# =========================
def main(block_sequence):
    # Start screen
    for bg in bg_rects: bg.autoDraw = False
    for fg in fg_rects: fg.autoDraw = False

    start = visual.TextStim(win, text="Press SPACE to start.\nPress ESC to quit.",
                            height=0.045, color=[1,1,1], colorSpace="rgb")
    start.draw()
    win.flip()
    keys = event.waitKeys(keyList=["space"] + QUIT_KEYS)
    if keys and keys[0] in QUIT_KEYS:
        safe_quit()


    # Restore rects
    for bg in bg_rects: bg.autoDraw = True
    for fg in fg_rects:
        fg.autoDraw = True
        fg.opacity = FG_OPACITY_OFF
    win.flip()

    print(f"Block sequence: {block_sequence}")

    for mode in block_sequence:
        run_fixation_seconds(sample_fix_seconds())

        # Reset foregrounds
        for fg in fg_rects:
            fg.opacity = FG_OPACITY_OFF

        if mode == "rift":
            run_trial_rift(RIFT_TRIAL_DURATION_S, RIFT_FREQS_HZ)
        elif mode == "localizer":
            run_trial_localizer(
                segment_s=LOCALIZER_SEGMENT_S,
                n_cycles=LOCALIZER_N_CYCLES,
                freqs_hz=LOCALIZER_FREQS_HZ,
                order=LOCALIZER_ORDER,
                others=LOCALIZER_OTHERS
            )
        else:
            raise ValueError(f"Unknown mode: {mode}")

        if ISI_S > 0:
            run_isi_seconds(ISI_S)

    # End screen
    for bg in bg_rects: bg.autoDraw = False
    for fg in fg_rects: fg.autoDraw = False

    end = visual.TextStim(win, text="Done. Thank you!", height=0.06, color=[1,1,1], colorSpace="rgb")
    end.draw()
    win.flip()
    core.wait(1.0)

    safe_quit()



if __name__ == "__main__":
    BLOCK_SEQUENCE = build_block_sequence(EXPERIMENT, RANDOM_SEED)
    main(BLOCK_SEQUENCE)
