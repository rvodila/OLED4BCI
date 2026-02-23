#!/usr/bin/env python
# -*- coding: utf-8 -*-

from psychopy import visual, core
from psychopy.hardware import keyboard
import math
import random

# =========================
# Hyperparameters

# =========================
DEFAULT_SUBJECT_ID = 42
DEFAULT_SESSION_ID = 1

def prompt_positive_int(label, default_value):
    prompt = f"{label} [{default_value}]: "
    while True:
        try:
            raw = input(prompt).strip()
        except EOFError:
            return int(default_value)
        if raw == "":
            return int(default_value)
        try:
            value = int(raw)
        except ValueError:
            print("Please enter an integer.")
            continue
        if value <= 0:
            print("Please enter a positive integer.")
            continue
        return value

print("==== PARTICIPANT INFO ====")
SUBJECT_ID = prompt_positive_int("Enter subject index", DEFAULT_SUBJECT_ID)
SESSION_ID = prompt_positive_int("Enter session number", DEFAULT_SESSION_ID)
RANDOM_SEED = SUBJECT_ID
print("==========================")

EXPERIMENT = {
    "design": {
        "n_runs": 4,
        "run_duration_min": 2.0,
        "trials_per_run": 8,
    },
    "timing": {
        "rift_trial_duration_s": 15.0,
        "localizer_trial_duration_s": 15.0,
        "localizer_segment_s": 2.0,
        "trial_onset_jitter_min_s": 0.2,
        "trial_onset_jitter_max_s": 0.5,
    },
    "frequencies": {
        "rift_hz":      [1.0, 2.0, 3.0, 4.0],
        "localizer_hz": [2.0, 2.0, 2.0, 2.0],
    },
    "localizer": {"order": (0, 1, 2, 3), "others": "off", "cycles_per_trial": 1},
    "rift": {"targets_per_trial": 3},
    "opto": {
        "enabled": True,
        "size_px": 40,
        "margin_px": 0,
        "flicker_hz": 15.0,
    },
}

# Deterministic randomness per subject (fixation jitter + anything else using global random)
random.seed(RANDOM_SEED)

# Pull values from config into simple names
N_RUNS = int(EXPERIMENT["design"]["n_runs"])
RUN_DURATION_MIN = float(EXPERIMENT["design"]["run_duration_min"])
TRIALS_PER_RUN = int(EXPERIMENT["design"]["trials_per_run"])
RIFT_TRIAL_DURATION_S = EXPERIMENT["timing"]["rift_trial_duration_s"]
LOCALIZER_TRIAL_DURATION_S = EXPERIMENT["timing"]["localizer_trial_duration_s"]
LOCALIZER_SEGMENT_S   = EXPERIMENT["timing"]["localizer_segment_s"]
TRIAL_ONSET_JITTER_MIN_S = float(EXPERIMENT["timing"]["trial_onset_jitter_min_s"])
TRIAL_ONSET_JITTER_MAX_S = float(EXPERIMENT["timing"]["trial_onset_jitter_max_s"])

RIFT_FREQS_HZ      = EXPERIMENT["frequencies"]["rift_hz"]
LOCALIZER_FREQS_HZ = EXPERIMENT["frequencies"]["localizer_hz"]

LOCALIZER_ORDER  = EXPERIMENT["localizer"]["order"]
LOCALIZER_OTHERS = EXPERIMENT["localizer"]["others"]
LOCALIZER_CYCLES_PER_TRIAL = int(EXPERIMENT["localizer"]["cycles_per_trial"])
TARGETS_PER_TRIAL = int(EXPERIMENT["rift"]["targets_per_trial"])
OPTO_ENABLED = EXPERIMENT["opto"]["enabled"]
OPTO_SIZE_PX = int(EXPERIMENT["opto"]["size_px"])
OPTO_MARGIN_PX = int(EXPERIMENT["opto"]["margin_px"])
OPTO_FLICKER_HZ = float(EXPERIMENT["opto"]["flicker_hz"])

if not math.isclose(RIFT_TRIAL_DURATION_S, LOCALIZER_TRIAL_DURATION_S, rel_tol=0.0, abs_tol=1e-9):
    raise ValueError("RIFT and localizer trial durations must be equal for balanced run timing.")
TRIAL_DURATION_S = float(RIFT_TRIAL_DURATION_S)
RUN_DURATION_S = RUN_DURATION_MIN * 60.0
if not math.isclose(TRIALS_PER_RUN * TRIAL_DURATION_S, RUN_DURATION_S, rel_tol=0.0, abs_tol=1e-6):
    raise ValueError("trials_per_run * trial_duration_s must equal run_duration_min * 60.")
if TRIALS_PER_RUN % 2 != 0:
    raise ValueError("trials_per_run must be even to keep equal numbers of trial types.")
if TRIAL_ONSET_JITTER_MIN_S < 0.0 or TRIAL_ONSET_JITTER_MAX_S < TRIAL_ONSET_JITTER_MIN_S:
    raise ValueError("Invalid trial onset jitter range.")
if LOCALIZER_CYCLES_PER_TRIAL <= 0:
    raise ValueError("localizer.cycles_per_trial must be >= 1.")
if TARGETS_PER_TRIAL < 0:
    raise ValueError("rift.targets_per_trial must be >= 0.")

TRIALS_PER_TYPE_PER_RUN = TRIALS_PER_RUN // 2
TOTAL_TRIALS = N_RUNS * TRIALS_PER_RUN
RIFT_TRIALS_TOTAL = N_RUNS * TRIALS_PER_TYPE_PER_RUN
LOCALIZER_TRIALS_TOTAL = N_RUNS * TRIALS_PER_TYPE_PER_RUN
TOTAL_ACTIVE_DURATION_S = N_RUNS * RUN_DURATION_S
TOTAL_DURATION_MIN = TOTAL_ACTIVE_DURATION_S / 60.0

# Window
FULLSCREEN = True
UNITS = "height"
BG_COLOR = [0, 0, 0]

# Layout
OFFSET = 0.25
STIM_SIZE = 0.15
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
RESPONSE_KEYS = ["space"]
FIX_COLOR_BASE = [1, 1, 1]
FIX_COLOR_TARGET = [1, 0, 0]
TARGET_RESPONSE_WINDOW_S = 3.0
TARGET_ONSET_JITTER_S = 2.0
OPTO_COLOR_ON = [1, 1, 1]
OPTO_COLOR_OFF = [-1, -1, -1]
RUN_TRIAL_SEQUENCES = []
SEEN_KEY_EVENTS = set()
INSTRUCTION_TEXT = (
    "Welcome.\n\n"
    "INSTRUCTION TEXT\n"
    "Press ESC at any time to quit.\n\n"
    "Press SPACE to start."
)
POST_EXPERIMENT_TEXT = (
    "This Experiment is complete.\n\n"
    "Thank you for participating.\n\n"
    "Please wait for the experimenter.\n"
    "Press SPACE to close."
)

print("==== RUN INFO ====")
print(f"Subject ID: {SUBJECT_ID}")
print(f"Session ID: {SESSION_ID}")
print(f"Random seed: {RANDOM_SEED}")
print(f"Design: {N_RUNS} runs x {RUN_DURATION_MIN:.1f} min, {TRIALS_PER_RUN} trials/run")
print(
    f"RIFT: {RIFT_TRIALS_TOTAL} trials total "
    f"({TRIALS_PER_TYPE_PER_RUN} per run x {N_RUNS} runs)"
)
print(
    f"LOCALIZER: {LOCALIZER_TRIALS_TOTAL} trials total "
    f"({TRIALS_PER_TYPE_PER_RUN} per run x {N_RUNS} runs)"
)
print(
    f"Localizer cycles per trial: {LOCALIZER_CYCLES_PER_TRIAL}; "
    f"Rift targets per trial: {TARGETS_PER_TRIAL}"
)
print(
    f"Opto square: {OPTO_SIZE_PX}x{OPTO_SIZE_PX}px, "
    f"upper-left margin {OPTO_MARGIN_PX}px, binary flicker {OPTO_FLICKER_HZ:.1f} Hz"
)
print(
    f"Trial onset jitter: {TRIAL_ONSET_JITTER_MIN_S * 1000:.0f}-{TRIAL_ONSET_JITTER_MAX_S * 1000:.0f} ms"
)
print(f"Total active duration: {TOTAL_DURATION_MIN:.1f} min ({TOTAL_TRIALS} trials x {TRIAL_DURATION_S:.1f}s)")
print("==================")


# =========================
# Helpers
# =========================
def build_run_trial_sequences(n_runs, trials_per_run, seed):
    rng = random.Random(seed)
    trials_per_type = trials_per_run // 2
    run_sequences = []
    for _ in range(n_runs):
        run_seq = (["rift"] * trials_per_type) + (["localizer"] * trials_per_type)
        rng.shuffle(run_seq)
        run_sequences.append(run_seq)
    return run_sequences

def seconds_to_frames(seconds, fps):
    return int(round(seconds * fps))

def sample_trial_onset_jitter_seconds():
    return random.uniform(TRIAL_ONSET_JITTER_MIN_S, TRIAL_ONSET_JITTER_MAX_S)

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

def clear_keyboard_events():
    kb.clearEvents()
    SEEN_KEY_EVENTS.clear()

def key_time_s(key_obj):
    key_tdown = getattr(key_obj, "tDown", None)
    if key_tdown is not None:
        key_t_s = float(key_tdown)
        now_s = EXP_CLOCK.getTime()
        # Some PsychoPy backends expose absolute tDown; convert to experiment-relative when detected.
        if key_t_s > (now_s + 5.0):
            try:
                key_t_s = key_t_s - float(EXP_CLOCK.getLastResetTime())
            except Exception:
                pass
        return key_t_s
    return EXP_CLOCK.getTime()

def key_uid(key_obj):
    key_name = getattr(key_obj, "name", str(key_obj))
    key_tdown = getattr(key_obj, "tDown", None)
    return (
        key_name,
        None if key_tdown is None else round(float(key_tdown), 6),
    )

def poll_keyboard(context, run_idx=None, trial_type=None, trial_idx=None, key_list=None):
    if key_list is None:
        key_list = RESPONSE_KEYS + QUIT_KEYS
    keys = kb.getKeys(keyList=key_list, waitRelease=False, clear=False)
    if not keys:
        return []
    new_keys = []
    for key_obj in keys:
        uid = key_uid(key_obj)
        if uid in SEEN_KEY_EVENTS:
            continue
        SEEN_KEY_EVENTS.add(uid)
        key_name = getattr(key_obj, "name", str(key_obj))
        key_t_s = key_time_s(key_obj)
        # TIMESTAMP HOOK: key press detected at `key_t_s`.
        if key_name in QUIT_KEYS:
            # TIMESTAMP HOOK: quit key pressed.
            safe_quit()
        new_keys.append({"name": key_name, "t_s": key_t_s, "obj": key_obj})
    return new_keys

def build_target_onsets(duration_s, n_targets):
    usable_s = max(0.0, duration_s - TARGET_RESPONSE_WINDOW_S)
    if n_targets <= 0 or usable_s <= 0.0:
        return []
    bin_size_s = usable_s / float(n_targets)
    onset_times_s = []
    for idx in range(n_targets):
        left = idx * bin_size_s
        right = (idx + 1) * bin_size_s
        center = 0.5 * (left + right)
        jitter = min(TARGET_ONSET_JITTER_S, 0.45 * bin_size_s)
        onset = random.uniform(center - jitter, center + jitter)
        onset_times_s.append(max(left, min(onset, right)))
    return onset_times_s

def update_opto_square():
    if not OPTO_ENABLED or opto_square is None:
        return
    if OPTO_FLICKER_HZ <= 0.0:
        color = OPTO_COLOR_ON
    else:
        color = OPTO_COLOR_ON if opto_flicker.step_value() >= 0.0 else OPTO_COLOR_OFF
    opto_square.fillColor = color
    opto_square.lineColor = color

def wait_for_space(context, run_idx=None, trial_type=None):
    clear_keyboard_events()
    while True:
        keys = poll_keyboard(context=context, run_idx=run_idx, trial_type=trial_type)
        if any(k["name"] == "space" for k in keys):
            return
        update_opto_square()
        win.flip()

def show_info_screen(text, context, run_idx=None):
    bg_states = [bg.autoDraw for bg in bg_circles]
    fg_states = [fg.autoDraw for fg in fg_circles]
    fixation_state = fixation.autoDraw
    for bg in bg_circles:
        bg.autoDraw = False
    for fg in fg_circles:
        fg.autoDraw = False
    fixation.autoDraw = False

    info = visual.TextStim(
        win,
        text=text,
        height=0.045,
        color=[1, 1, 1],
        colorSpace="rgb",
    )
    info.autoDraw = True
    wait_for_space(context=context, run_idx=run_idx)
    info.autoDraw = False

    for bg, state in zip(bg_circles, bg_states):
        bg.autoDraw = state
    for fg, state in zip(fg_circles, fg_states):
        fg.autoDraw = state
    fixation.autoDraw = fixation_state
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
kb = keyboard.Keyboard()
EXP_CLOCK = kb.clock
# Single experiment timebase for event logs, onsets, and RTs.
EXP_CLOCK.reset()

fixation = visual.TextStim(win, text="+", height=0.05, color=FIX_COLOR_BASE, colorSpace="rgb")
fixation.autoDraw = False

bg_circles, fg_circles = [], []
for pos in QUAD_POSITIONS:
    bg = visual.Circle(win, radius=STIM_SIZE / 2.0, pos=pos, units=UNITS,
                       fillColor=RECT_BG_COLOR,
                       lineColor=RECT_BG_COLOR,
                       colorSpace="rgb",
                       opacity=1.0)
    fg = visual.Circle(win, radius=STIM_SIZE / 2.0, pos=pos, units=UNITS,
                       fillColor=RECT_FG_COLOR,
                       lineColor=RECT_FG_COLOR,
                       colorSpace="rgb",
                       opacity=FG_OPACITY_OFF)

    bg.autoDraw = True
    fg.autoDraw = True
    bg_circles.append(bg)
    fg_circles.append(fg)

opto_square = None
opto_flicker = None
if OPTO_ENABLED:
    win_w_px, win_h_px = win.size
    opto_pos = (
        -0.5 * win_w_px + OPTO_MARGIN_PX + 0.5 * OPTO_SIZE_PX,
        0.5 * win_h_px - OPTO_MARGIN_PX - 0.5 * OPTO_SIZE_PX,
    )
    opto_square = visual.Rect(
        win,
        width=OPTO_SIZE_PX,
        height=OPTO_SIZE_PX,
        pos=opto_pos,
        units="pix",
        fillColor=OPTO_COLOR_OFF,
        lineColor=OPTO_COLOR_OFF,
        colorSpace="rgb",
    )
    opto_square.autoDraw = True
    opto_flicker = SineFlicker(OPTO_FLICKER_HZ, fps, phase0=0.0)



# =========================
# Routines
# =========================
def run_trial_rift(freqs_hz, duration_s, run_idx=None, trial_idx=None):
    phases = [SineFlicker(freqs_hz[i], fps, phase0=0.0) for i in range(4)]
    response_window_frames = seconds_to_frames(TARGET_RESPONSE_WINDOW_S, fps)
    target_onsets_s = build_target_onsets(duration_s, TARGETS_PER_TRIAL)
    next_target_idx = 0
    target_id_counter = 0
    active_target = None
    false_alarm_counter = 0
    n_frames = seconds_to_frames(duration_s, fps)

    for frame_idx in range(n_frames):
        if frame_idx == 0:
            # TIMESTAMP HOOK: trial_start (flip-locked) on first frame flip.
            pass

        fg_circles[0].opacity = OPACITY_MIN + ((phases[0].step_value() + 1.0) * 0.5) * (OPACITY_MAX - OPACITY_MIN)
        fg_circles[1].opacity = OPACITY_MIN + ((phases[1].step_value() + 1.0) * 0.5) * (OPACITY_MAX - OPACITY_MIN)
        fg_circles[2].opacity = OPACITY_MIN + ((phases[2].step_value() + 1.0) * 0.5) * (OPACITY_MAX - OPACITY_MIN)
        fg_circles[3].opacity = OPACITY_MIN + ((phases[3].step_value() + 1.0) * 0.5) * (OPACITY_MAX - OPACITY_MIN)

        can_start_target = frame_idx < max(0, n_frames - response_window_frames)
        now_rel_s = frame_idx / float(fps)
        if (
            (active_target is None)
            and (next_target_idx < len(target_onsets_s))
            and (now_rel_s >= target_onsets_s[next_target_idx])
            and can_start_target
        ):
            target_id_counter += 1
            active_target = {
                "run_id": run_idx,
                "trial_id": trial_idx,
                "target_id": target_id_counter,
                "onset_s": None,
                "deadline_s_abs": None,
            }
            fixation.color = FIX_COLOR_TARGET

            def set_target_onset(target_ref=active_target):
                onset_s = EXP_CLOCK.getTime()
                target_ref["onset_s"] = onset_s
                target_ref["deadline_s_abs"] = onset_s + TARGET_RESPONSE_WINDOW_S

            win.callOnFlip(set_target_onset)
            # TIMESTAMP HOOK: target_onset (flip-locked) on this flip.
            next_target_idx += 1

        if frame_idx == (n_frames - 1):
            # TIMESTAMP HOOK: trial_end (flip-locked) on final frame flip.
            pass

        update_opto_square()
        win.flip()

        keys = poll_keyboard(
            context="rift",
            run_idx=run_idx,
            trial_type="rift",
            trial_idx=trial_idx,
        )
        if keys:
            for key_ev in keys:
                key_name = key_ev["name"]
                key_t_s = key_ev["t_s"]
                if key_name in RESPONSE_KEYS:
                    now_s = EXP_CLOCK.getTime()
                    if (
                        active_target is not None
                        and active_target["onset_s"] is not None
                        and now_s <= active_target["deadline_s_abs"]
                    ):
                        response_t_s = key_t_s
                        if (
                            response_t_s < active_target["onset_s"]
                            or response_t_s > active_target["deadline_s_abs"]
                        ):
                            response_t_s = now_s
                        # TIMESTAMP HOOK: response hit at `response_t_s`; RT = response_t_s - target_onset_s.
                        active_target = None
                        fixation.color = FIX_COLOR_BASE
                    else:
                        false_alarm_counter += 1
                        # TIMESTAMP HOOK: response outside target window / false alarm.

        if (
            active_target is not None
            and active_target["onset_s"] is not None
            and EXP_CLOCK.getTime() >= active_target["deadline_s_abs"]
        ):
            # TIMESTAMP HOOK: target timeout / deadline passed.
            active_target = None
            fixation.color = FIX_COLOR_BASE

    if active_target is not None and active_target["onset_s"] is not None:
        # TIMESTAMP HOOK: trial ended with unresolved active target (timeout at trial end).
        fixation.color = FIX_COLOR_BASE

def run_trial_localizer(duration_s, segment_s, freqs_hz, order=(0,1,2,3), others="off", n_cycles_per_trial=1, run_idx=None, trial_idx=None):
    n_frames = seconds_to_frames(duration_s, fps)
    phases = [SineFlicker(freqs_hz[i], fps, phase0=0.0) for i in range(4)]
    n_segments = max(1, int(n_cycles_per_trial) * len(order))
    base_seg_frames = max(1, n_frames // n_segments)
    extra_frames = max(0, n_frames - (base_seg_frames * n_segments))
    frame_idx_global = 0

    for seg_idx in range(n_segments):
        target = order[seg_idx % len(order)]
        seg_len = base_seg_frames + (1 if seg_idx < extra_frames else 0)
        for seg_frame_idx in range(seg_len):
            for idx in range(4):
                bg_circles[idx].opacity = 1.0 if idx == target else 0.0
                fg_circles[idx].opacity = FG_OPACITY_OFF

            fg_circles[target].opacity = OPACITY_MIN + ((phases[target].step_value() + 1.0) * 0.5) * (OPACITY_MAX - OPACITY_MIN)

            if frame_idx_global == 0:
                # TIMESTAMP HOOK: localizer trial_start (flip-locked) on first frame flip.
                pass

            if seg_frame_idx == 0:
                # TIMESTAMP HOOK: localizer target_onset (flip-locked) for this segment/quadrant.
                pass

            if frame_idx_global == (n_frames - 1):
                # TIMESTAMP HOOK: localizer trial_end (flip-locked) on final frame flip.
                pass

            poll_keyboard(context="localizer", run_idx=run_idx, trial_type="localizer", trial_idx=trial_idx)
            update_opto_square()
            win.flip()
            frame_idx_global += 1

    # Restore defaults so subsequent trials are unaffected by localizer visibility state.
    for bg in bg_circles:
        bg.opacity = 1.0
    for fg in fg_circles:
        fg.opacity = FG_OPACITY_OFF

def run_trial_onset_jitter_seconds(duration_s, run_idx=None, trial_idx=None, trial_type=None):
    for bg in bg_circles:
        bg.opacity = 0.0
    for fg in fg_circles:
        fg.opacity = FG_OPACITY_OFF
    fixation.color = FIX_COLOR_BASE
    n_frames = seconds_to_frames(duration_s, fps)
    for _ in range(n_frames):
        poll_keyboard(context="trial_onset_jitter", run_idx=run_idx, trial_type=trial_type, trial_idx=trial_idx)
        update_opto_square()
        win.flip()

def safe_quit():
    """Close PsychoPy cleanly so the window doesn't hang and the kernel doesn't crash."""
    try:
        kb.clearEvents()
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
def main(run_trial_sequences):
    global RUN_TRIAL_SEQUENCES
    RUN_TRIAL_SEQUENCES = run_trial_sequences

    # Keep stimuli hidden until a run is announced.
    for bg in bg_circles: bg.autoDraw = False
    for fg in fg_circles: fg.autoDraw = False
    fixation.autoDraw = False

    # Instruction screen
    show_info_screen(
        text=INSTRUCTION_TEXT,
        context="instruction_screen",
        run_idx=None,
    )

    for run_idx, run_trials in enumerate(run_trial_sequences, start=1):
        if run_idx >= 2:
            show_info_screen(
                text="Take a break if needed.\nPress SPACE to continue.",
                context="break_screen",
                run_idx=run_idx,
            )
        show_info_screen(
            text=f"Run {run_idx}/{N_RUNS}\nPress SPACE to continue.",
            context="run_intro_screen",
            run_idx=run_idx,
        )

        # Only enable task stimuli after run intro is acknowledged.
        for bg in bg_circles: bg.autoDraw = True
        for fg in fg_circles:
            fg.autoDraw = True
            fg.opacity = FG_OPACITY_OFF
        fixation.autoDraw = True

        n_rift_in_run = run_trials.count("rift")
        n_localizer_in_run = run_trials.count("localizer")
        print(
            f"Run {run_idx}/{N_RUNS}: "
            f"{n_rift_in_run} rift + {n_localizer_in_run} localizer trials"
        )
        # TIMESTAMP HOOK: run_start (after run intro acknowledged).
        run_start_s = EXP_CLOCK.getTime()

        for trial_idx, trial_type in enumerate(run_trials, start=1):
            jitter_s = sample_trial_onset_jitter_seconds()
            # TIMESTAMP HOOK: trial_onset_jitter start (jitter_s begins here).
            run_trial_onset_jitter_seconds(
                jitter_s,
                run_idx=run_idx,
                trial_idx=trial_idx,
                trial_type=trial_type,
            )

            # Reset foregrounds and fixation state at trial start
            if trial_type == "rift":
                for bg in bg_circles:
                    bg.opacity = 1.0
            else:
                for bg in bg_circles:
                    bg.opacity = 0.0
            for fg in fg_circles:
                fg.opacity = FG_OPACITY_OFF
            fixation.color = FIX_COLOR_BASE

            if trial_type == "rift":
                run_trial_rift(
                    RIFT_FREQS_HZ,
                    duration_s=RIFT_TRIAL_DURATION_S,
                    run_idx=run_idx,
                    trial_idx=trial_idx,
                )
            elif trial_type == "localizer":
                run_trial_localizer(
                    duration_s=LOCALIZER_TRIAL_DURATION_S,
                    segment_s=LOCALIZER_SEGMENT_S,
                    freqs_hz=LOCALIZER_FREQS_HZ,
                    order=LOCALIZER_ORDER,
                    others=LOCALIZER_OTHERS,
                    n_cycles_per_trial=LOCALIZER_CYCLES_PER_TRIAL,
                    run_idx=run_idx,
                    trial_idx=trial_idx,
                )
            else:
                raise ValueError(f"Unknown trial type: {trial_type}")

            fixation.color = FIX_COLOR_BASE

        run_duration_s = EXP_CLOCK.getTime() - run_start_s
        # TIMESTAMP HOOK: run_end (duration_s = run_duration_s).

    # End screen
    for bg in bg_circles: bg.autoDraw = False
    for fg in fg_circles: fg.autoDraw = False
    fixation.autoDraw = False

    # TIMESTAMP HOOK: post_experiment_screen_shown.
    show_info_screen(
        text=POST_EXPERIMENT_TEXT,
        context="post_experiment_screen",
        run_idx=None,
    )
    # TIMESTAMP HOOK: post_experiment_screen_continue (space pressed).

    safe_quit()



if __name__ == "__main__":
    RUN_TRIAL_SEQUENCES = build_run_trial_sequences(N_RUNS, TRIALS_PER_RUN, RANDOM_SEED)
    main(RUN_TRIAL_SEQUENCES)
