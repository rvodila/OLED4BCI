#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import math
import random

import numpy as np
from psychopy import core, visual
from psychopy.hardware import keyboard


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
    "flicker": {
        # Aligns with MATLAB dualFlickerStimulus.m defaults.
        "frames_per_bit": 6,
        "ramp_len_frames": 2,
        "trial_taper_frames": 60,
        "opacity_min_255": 20,
        "opacity_max_255": 230,
        "codebook_npz": "metadata/mgold_61_6521.npz",
    },
    "rift": {
        "targets_per_trial": 3,
        # Modes: periodic | alternating | code | hybrid
        # Aliases: freq -> periodic, alt -> alternating
        "channels": [
            {"mode": "hybrid", "carrier_hz": 10.0, "code_row": 0},
            {"mode": "hybrid", "carrier_hz": 15.0, "code_row": 1},
            {"mode": "periodic", "carrier_hz": 3.0},
            {"mode": "alternating", "carrier_hz": 0.0},
        ],
    },
    "localizer": {
        "order": (0, 1, 2, 3),
        "others": "off",
        "cycles_per_trial": 1,
        "channels": [
            {"mode": "periodic", "carrier_hz": 2.0},
            {"mode": "periodic", "carrier_hz": 2.0},
            {"mode": "periodic", "carrier_hz": 2.0},
            {"mode": "periodic", "carrier_hz": 2.0},
        ],
    },
    "opto": {
        "enabled": True,
        "size_px": 40,
        "margin_px": 0,
        "mode": "alternating",
        "carrier_hz": 15.0,
        "code_row": 0,
    },
}


# Deterministic randomness per subject (fixation jitter + anything else using global random)
random.seed(RANDOM_SEED)


MODE_ALIASES = {
    "freq": "periodic",
    "alt": "alternating",
}
ALLOWED_MODES = {"periodic", "alternating", "code", "hybrid"}


def canonical_mode(mode):
    m = str(mode).strip().lower()
    m = MODE_ALIASES.get(m, m)
    if m not in ALLOWED_MODES:
        raise ValueError(f"Unknown flicker mode: {mode}")
    return m


def resolve_codebook_path(relative_or_absolute_path):
    script_dir = Path(__file__).resolve().parent
    direct = Path(relative_or_absolute_path)
    if direct.is_absolute() and direct.is_file():
        return direct

    local = script_dir / relative_or_absolute_path
    if local.is_file():
        return local

    # Fallback to shared assets (project/stimulus/assets/codes).
    fallback = script_dir.parents[4] / "assets" / "codes" / Path(relative_or_absolute_path).name
    if fallback.is_file():
        return fallback

    raise FileNotFoundError(
        "Could not resolve codebook path. Tried: "
        f"{local} and {fallback}"
    )


def load_codebook(path_like):
    resolved = resolve_codebook_path(path_like)
    with np.load(str(resolved)) as data:
        if "codes" not in data:
            raise KeyError(f"Expected key 'codes' in {resolved}")
        codes = np.asarray(data["codes"], dtype=np.float64)
    if codes.ndim != 2 or codes.shape[0] < 1 or codes.shape[1] < 1:
        raise ValueError(f"Codebook must be 2D with non-zero shape, got {codes.shape}")

    # Keep binary semantics but tolerate tiny numeric noise.
    codes = np.where(codes >= 0.5, 1.0, 0.0)
    return codes, resolved


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


def build_trial_envelope(n_frames, taper_frames):
    env = np.ones(n_frames, dtype=np.float64)
    tf = int(max(0, min(int(taper_frames), n_frames // 2)))
    if tf > 0:
        w = np.hanning(2 * tf)
        env[:tf] = w[:tf]
        env[-tf:] = w[tf:]
    return env


def raised_cosine_smooth(x, ramp_len):
    y = np.asarray(x, dtype=np.float64).copy()
    n = y.size
    if ramp_len <= 0 or n < 3:
        return y

    transitions = np.flatnonzero(np.diff(y) != 0.0) + 1
    for idx in transitions:
        old_level = y[max(0, idx - 1)]
        new_level = y[idx]

        left = max(0, idx - ramp_len)
        right = min(n - 1, idx + ramp_len - 1)
        length = right - left + 1
        if length < 2:
            continue

        u = np.linspace(0.0, 1.0, num=length)
        ramp = 0.5 * (1.0 - np.cos(np.pi * u))
        y[left : right + 1] = ((1.0 - ramp) * old_level) + (ramp * new_level)
    return y


def expand_code(code_bits, n_frames, frames_per_bit, ramp_len):
    bits = np.asarray(code_bits, dtype=np.float64).ravel()
    if bits.size == 0:
        return np.array([], dtype=np.float64)

    bits = np.where(bits >= 0.5, 1.0, 0.0)
    fpb = max(1, int(frames_per_bit))
    expanded = np.repeat(bits, fpb)

    reps = int(math.ceil(float(n_frames) / float(expanded.size)))
    code_long = np.tile(expanded, reps)[:n_frames]
    return raised_cosine_smooth(code_long, int(max(0, ramp_len)))


def mode_map01(
    mode,
    n_frames,
    fps,
    carrier_hz,
    code_bits,
    frames_per_bit,
    ramp_len,
    taper_frames,
):
    mode = canonical_mode(mode)
    t = np.arange(n_frames, dtype=np.float64) / float(fps)
    env = build_trial_envelope(n_frames, taper_frames)
    carrier01 = 0.5 + 0.5 * np.sin(2.0 * np.pi * float(carrier_hz) * t)

    code_long = None
    if mode in {"code", "hybrid"}:
        if code_bits is None:
            raise ValueError(f"Mode '{mode}' needs code_bits.")
        code_long = expand_code(code_bits, n_frames, frames_per_bit, ramp_len)

    if mode == "periodic":
        contrast = (carrier01 - 0.5) * env
        map01 = 0.5 + contrast
    elif mode == "alternating":
        alt = (np.arange(n_frames) % 2).astype(np.float64)
        contrast = (alt - 0.5) * env
        map01 = 0.5 + contrast
    elif mode == "code":
        contrast = (code_long - 0.5) * env
        map01 = 0.5 + contrast
    elif mode == "hybrid":
        mod_bipolar = (2.0 * code_long - 1.0) * (2.0 * carrier01 - 1.0)
        mod_bipolar *= env
        map01 = 0.5 * (mod_bipolar + 1.0)
    else:
        raise ValueError(f"Unsupported mode '{mode}'")

    return np.clip(map01, 0.0, 1.0), code_long


def select_code_for_channel(codebook, channel_spec, default_row):
    row = int(channel_spec.get("code_row", default_row))
    row_wrapped = row % int(codebook.shape[0])
    return codebook[row_wrapped, :]


def build_trial_opacity_traces(channel_specs, duration_s, fps, taper_frames):
    n_frames = seconds_to_frames(duration_s, fps)
    traces = np.zeros((len(channel_specs), n_frames), dtype=np.float64)

    for idx, spec in enumerate(channel_specs):
        mode = canonical_mode(spec.get("mode", "periodic"))
        carrier_hz = float(spec.get("carrier_hz", 0.0))

        code_bits = None
        if mode in {"code", "hybrid"}:
            code_bits = select_code_for_channel(CODEBOOK, spec, default_row=idx)

        map01, _ = mode_map01(
            mode=mode,
            n_frames=n_frames,
            fps=fps,
            carrier_hz=carrier_hz,
            code_bits=code_bits,
            frames_per_bit=FRAMES_PER_BIT,
            ramp_len=RAMP_LEN_FRAMES,
            taper_frames=taper_frames,
        )
        traces[idx, :] = OPACITY_MIN + map01 * (OPACITY_MAX - OPACITY_MIN)

    return traces


def sample_mode01_realtime(mode, frame_idx, fps, carrier_hz, code_bits=None, frames_per_bit=1):
    mode = canonical_mode(mode)
    if mode == "periodic":
        val = 0.5 + 0.5 * math.sin((2.0 * math.pi * carrier_hz * frame_idx) / float(fps))
        return max(0.0, min(1.0, val))

    if mode == "alternating":
        return float(frame_idx % 2)

    if code_bits is None:
        raise ValueError(f"Mode '{mode}' requires code_bits for realtime sampling.")

    bits = np.asarray(code_bits, dtype=np.float64).ravel()
    if bits.size == 0:
        return 0.0

    bit_idx = (frame_idx // max(1, int(frames_per_bit))) % bits.size
    bit = 1.0 if bits[bit_idx] >= 0.5 else 0.0

    if mode == "code":
        return bit

    carrier_bipolar = math.sin((2.0 * math.pi * carrier_hz * frame_idx) / float(fps))
    code_bipolar = 1.0 if bit >= 0.5 else -1.0
    return 0.5 * ((carrier_bipolar * code_bipolar) + 1.0)


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


def update_opto_square():
    global OPTO_FRAME_COUNTER
    if not OPTO_ENABLED or opto_square is None:
        return

    map01 = sample_mode01_realtime(
        mode=OPTO_MODE,
        frame_idx=OPTO_FRAME_COUNTER,
        fps=fps,
        carrier_hz=OPTO_CARRIER_HZ,
        code_bits=OPTO_CODE_BITS,
        frames_per_bit=FRAMES_PER_BIT,
    )
    color = OPTO_COLOR_ON if map01 >= 0.5 else OPTO_COLOR_OFF
    opto_square.fillColor = color
    opto_square.lineColor = color
    OPTO_FRAME_COUNTER += 1


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


def describe_channel(idx, spec):
    mode = canonical_mode(spec.get("mode", "periodic"))
    hz = float(spec.get("carrier_hz", 0.0))
    row = spec.get("code_row", None)
    if mode in {"code", "hybrid"}:
        return f"Q{idx+1}: mode={mode}, carrier={hz:.3f} Hz, code_row={int(row) if row is not None else idx}"
    return f"Q{idx+1}: mode={mode}, carrier={hz:.3f} Hz"


# Pull values from config into simple names
N_RUNS = int(EXPERIMENT["design"]["n_runs"])
RUN_DURATION_MIN = float(EXPERIMENT["design"]["run_duration_min"])
TRIALS_PER_RUN = int(EXPERIMENT["design"]["trials_per_run"])
RIFT_TRIAL_DURATION_S = float(EXPERIMENT["timing"]["rift_trial_duration_s"])
LOCALIZER_TRIAL_DURATION_S = float(EXPERIMENT["timing"]["localizer_trial_duration_s"])
LOCALIZER_SEGMENT_S = float(EXPERIMENT["timing"]["localizer_segment_s"])
TRIAL_ONSET_JITTER_MIN_S = float(EXPERIMENT["timing"]["trial_onset_jitter_min_s"])
TRIAL_ONSET_JITTER_MAX_S = float(EXPERIMENT["timing"]["trial_onset_jitter_max_s"])

FRAMES_PER_BIT = int(EXPERIMENT["flicker"]["frames_per_bit"])
RAMP_LEN_FRAMES = int(EXPERIMENT["flicker"]["ramp_len_frames"])
TRIAL_TAPER_FRAMES = int(EXPERIMENT["flicker"]["trial_taper_frames"])
OPACITY_MIN_255 = int(EXPERIMENT["flicker"]["opacity_min_255"])
OPACITY_MAX_255 = int(EXPERIMENT["flicker"]["opacity_max_255"])
CODEBOOK_PATH_SETTING = EXPERIMENT["flicker"]["codebook_npz"]

RIFT_CHANNELS = list(EXPERIMENT["rift"]["channels"])
TARGETS_PER_TRIAL = int(EXPERIMENT["rift"]["targets_per_trial"])

LOCALIZER_ORDER = tuple(EXPERIMENT["localizer"]["order"])
LOCALIZER_OTHERS = str(EXPERIMENT["localizer"]["others"]).lower()
LOCALIZER_CYCLES_PER_TRIAL = int(EXPERIMENT["localizer"]["cycles_per_trial"])
LOCALIZER_CHANNELS = list(EXPERIMENT["localizer"]["channels"])

OPTO_ENABLED = bool(EXPERIMENT["opto"]["enabled"])
OPTO_SIZE_PX = int(EXPERIMENT["opto"]["size_px"])
OPTO_MARGIN_PX = int(EXPERIMENT["opto"]["margin_px"])
OPTO_MODE = canonical_mode(EXPERIMENT["opto"]["mode"])
OPTO_CARRIER_HZ = float(EXPERIMENT["opto"].get("carrier_hz", 0.0))

if len(RIFT_CHANNELS) != 4:
    raise ValueError("rift.channels must contain exactly 4 channel specs.")
if len(LOCALIZER_CHANNELS) != 4:
    raise ValueError("localizer.channels must contain exactly 4 channel specs.")
if any(canonical_mode(spec.get("mode", "periodic")) not in ALLOWED_MODES for spec in RIFT_CHANNELS):
    raise ValueError("One or more RIFT channel modes are invalid.")
if any(canonical_mode(spec.get("mode", "periodic")) not in ALLOWED_MODES for spec in LOCALIZER_CHANNELS):
    raise ValueError("One or more localizer channel modes are invalid.")
if OPACITY_MIN_255 < 0 or OPACITY_MAX_255 > 255 or OPACITY_MIN_255 >= OPACITY_MAX_255:
    raise ValueError("flicker opacity_min_255/opacity_max_255 must satisfy 0 <= min < max <= 255.")
if FRAMES_PER_BIT <= 0:
    raise ValueError("flicker.frames_per_bit must be >= 1.")
if RAMP_LEN_FRAMES < 0:
    raise ValueError("flicker.ramp_len_frames must be >= 0.")
if TARGETS_PER_TRIAL < 0:
    raise ValueError("rift.targets_per_trial must be >= 0.")
if LOCALIZER_CYCLES_PER_TRIAL <= 0:
    raise ValueError("localizer.cycles_per_trial must be >= 1.")
if LOCALIZER_OTHERS not in {"off", "on"}:
    raise ValueError("localizer.others must be 'off' or 'on'.")
if TRIAL_ONSET_JITTER_MIN_S < 0.0 or TRIAL_ONSET_JITTER_MAX_S < TRIAL_ONSET_JITTER_MIN_S:
    raise ValueError("Invalid trial onset jitter range.")

if not math.isclose(RIFT_TRIAL_DURATION_S, LOCALIZER_TRIAL_DURATION_S, rel_tol=0.0, abs_tol=1e-9):
    raise ValueError("RIFT and localizer trial durations must be equal for balanced run timing.")

TRIAL_DURATION_S = float(RIFT_TRIAL_DURATION_S)
RUN_DURATION_S = RUN_DURATION_MIN * 60.0
if not math.isclose(TRIALS_PER_RUN * TRIAL_DURATION_S, RUN_DURATION_S, rel_tol=0.0, abs_tol=1e-6):
    raise ValueError("trials_per_run * trial_duration_s must equal run_duration_min * 60.")
if TRIALS_PER_RUN % 2 != 0:
    raise ValueError("trials_per_run must be even to keep equal numbers of trial types.")

TRIALS_PER_TYPE_PER_RUN = TRIALS_PER_RUN // 2
TOTAL_TRIALS = N_RUNS * TRIALS_PER_RUN
RIFT_TRIALS_TOTAL = N_RUNS * TRIALS_PER_TYPE_PER_RUN
LOCALIZER_TRIALS_TOTAL = N_RUNS * TRIALS_PER_TYPE_PER_RUN
TOTAL_ACTIVE_DURATION_S = N_RUNS * RUN_DURATION_S
TOTAL_DURATION_MIN = TOTAL_ACTIVE_DURATION_S / 60.0

OPACITY_MIN = OPACITY_MIN_255 / 255.0
OPACITY_MAX = OPACITY_MAX_255 / 255.0

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
FG_OPACITY_OFF = 0.0

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
OPTO_FRAME_COUNTER = 0

INSTRUCTION_TEXT = (
    "Welcome.\\n\\n"
    "Attend the target and respond with SPACE when fixation turns red.\\n"
    "Press ESC at any time to quit.\\n\\n"
    "Press SPACE to start."
)
POST_EXPERIMENT_TEXT = (
    "This experiment is complete.\\n\\n"
    "Thank you for participating.\\n\\n"
    "Please wait for the experimenter.\\n"
    "Press SPACE to close."
)

CODEBOOK, CODEBOOK_PATH = load_codebook(CODEBOOK_PATH_SETTING)
OPTO_CODE_BITS = None
if OPTO_MODE in {"code", "hybrid"}:
    OPTO_CODE_BITS = select_code_for_channel(
        CODEBOOK,
        EXPERIMENT["opto"],
        default_row=int(EXPERIMENT["opto"].get("code_row", 0)),
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
    f"RIFT targets per trial: {TARGETS_PER_TRIAL}"
)
print(
    f"Trial onset jitter: {TRIAL_ONSET_JITTER_MIN_S * 1000:.0f}-{TRIAL_ONSET_JITTER_MAX_S * 1000:.0f} ms"
)
print(f"Codebook: {CODEBOOK_PATH} shape={CODEBOOK.shape} (rows, bits)")
print(
    f"Flicker shared params: frames_per_bit={FRAMES_PER_BIT}, "
    f"ramp_len={RAMP_LEN_FRAMES}, taper={TRIAL_TAPER_FRAMES} frames"
)
print("RIFT channel modes:")
for i, spec in enumerate(RIFT_CHANNELS):
    print("  " + describe_channel(i, spec))
print("Localizer channel modes:")
for i, spec in enumerate(LOCALIZER_CHANNELS):
    print("  " + describe_channel(i, spec))
if OPTO_ENABLED:
    opto_row = int(EXPERIMENT["opto"].get("code_row", 0))
    opto_row_text = f", code_row={opto_row}" if OPTO_MODE in {"code", "hybrid"} else ""
    print(
        f"Opto square: {OPTO_SIZE_PX}x{OPTO_SIZE_PX}px, "
        f"upper-left margin {OPTO_MARGIN_PX}px, mode={OPTO_MODE}, carrier={OPTO_CARRIER_HZ:.1f} Hz{opto_row_text}"
    )
else:
    print("Opto square: disabled")
print(f"Total active duration: {TOTAL_DURATION_MIN:.1f} min ({TOTAL_TRIALS} trials x {TRIAL_DURATION_S:.1f}s)")
print("==================")


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
print(f"Measured FPS: {fps:.3f}")

kb = keyboard.Keyboard()
EXP_CLOCK = kb.clock
# Single experiment timebase for event logs, onsets, and RTs.
EXP_CLOCK.reset()

fixation = visual.TextStim(win, text="+", height=0.05, color=FIX_COLOR_BASE, colorSpace="rgb")
fixation.autoDraw = False

bg_circles, fg_circles = [], []
for pos in QUAD_POSITIONS:
    bg = visual.Circle(
        win,
        radius=STIM_SIZE / 2.0,
        pos=pos,
        units=UNITS,
        fillColor=RECT_BG_COLOR,
        lineColor=RECT_BG_COLOR,
        colorSpace="rgb",
        opacity=1.0,
    )
    fg = visual.Circle(
        win,
        radius=STIM_SIZE / 2.0,
        pos=pos,
        units=UNITS,
        fillColor=RECT_FG_COLOR,
        lineColor=RECT_FG_COLOR,
        colorSpace="rgb",
        opacity=FG_OPACITY_OFF,
    )

    bg.autoDraw = True
    fg.autoDraw = True
    bg_circles.append(bg)
    fg_circles.append(fg)

opto_square = None
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


# =========================
# Routines
# =========================
def run_trial_rift(duration_s, run_idx=None, trial_idx=None):
    opacity_traces = build_trial_opacity_traces(
        channel_specs=RIFT_CHANNELS,
        duration_s=duration_s,
        fps=fps,
        taper_frames=TRIAL_TAPER_FRAMES,
    )
    n_frames = opacity_traces.shape[1]

    response_window_frames = seconds_to_frames(TARGET_RESPONSE_WINDOW_S, fps)
    target_onsets_s = build_target_onsets(duration_s, TARGETS_PER_TRIAL)
    next_target_idx = 0
    target_id_counter = 0
    active_target = None
    false_alarm_counter = 0

    for frame_idx in range(n_frames):
        if frame_idx == 0:
            # TIMESTAMP HOOK: trial_start (flip-locked) on first frame flip.
            pass

        for stim_idx in range(4):
            fg_circles[stim_idx].opacity = float(opacity_traces[stim_idx, frame_idx])

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



def run_trial_localizer(
    duration_s,
    segment_s,
    order=(0, 1, 2, 3),
    others="off",
    n_cycles_per_trial=1,
    run_idx=None,
    trial_idx=None,
):
    opacity_traces = build_trial_opacity_traces(
        channel_specs=LOCALIZER_CHANNELS,
        duration_s=duration_s,
        fps=fps,
        taper_frames=TRIAL_TAPER_FRAMES,
    )
    n_frames = opacity_traces.shape[1]

    n_segments = max(1, int(n_cycles_per_trial) * len(order))
    if segment_s > 0:
        suggested_segments = max(1, int(round(duration_s / segment_s)))
        if suggested_segments != n_segments:
            # Keep behavior deterministic while surfacing mismatch in console.
            print(
                "[localizer] segment_s implies "
                f"{suggested_segments} segments but cycles/order implies {n_segments}; using cycles/order."
            )

    base_seg_frames = max(1, n_frames // n_segments)
    extra_frames = max(0, n_frames - (base_seg_frames * n_segments))
    frame_idx_global = 0

    for seg_idx in range(n_segments):
        target = order[seg_idx % len(order)]
        seg_len = base_seg_frames + (1 if seg_idx < extra_frames else 0)
        for seg_frame_idx in range(seg_len):
            for idx in range(4):
                if others == "on":
                    bg_circles[idx].opacity = 1.0
                else:
                    bg_circles[idx].opacity = 1.0 if idx == target else 0.0
                fg_circles[idx].opacity = FG_OPACITY_OFF

            fg_circles[target].opacity = float(opacity_traces[target, frame_idx_global])

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
    """Close PsychoPy cleanly so the window does not hang and the kernel does not crash."""
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
        core.quit()  # raises SystemExit
    except Exception:
        raise SystemExit


# =========================
# Main
# =========================
def main(run_trial_sequences):
    global RUN_TRIAL_SEQUENCES
    RUN_TRIAL_SEQUENCES = run_trial_sequences

    # Keep stimuli hidden until a run is announced.
    for bg in bg_circles:
        bg.autoDraw = False
    for fg in fg_circles:
        fg.autoDraw = False
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
                text="Take a break if needed.\\nPress SPACE to continue.",
                context="break_screen",
                run_idx=run_idx,
            )
        show_info_screen(
            text=f"Run {run_idx}/{N_RUNS}\\nPress SPACE to continue.",
            context="run_intro_screen",
            run_idx=run_idx,
        )

        # Only enable task stimuli after run intro is acknowledged.
        for bg in bg_circles:
            bg.autoDraw = True
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

            # Reset foregrounds and fixation state at trial start.
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
                    duration_s=RIFT_TRIAL_DURATION_S,
                    run_idx=run_idx,
                    trial_idx=trial_idx,
                )
            elif trial_type == "localizer":
                run_trial_localizer(
                    duration_s=LOCALIZER_TRIAL_DURATION_S,
                    segment_s=LOCALIZER_SEGMENT_S,
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
    for bg in bg_circles:
        bg.autoDraw = False
    for fg in fg_circles:
        fg.autoDraw = False
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
