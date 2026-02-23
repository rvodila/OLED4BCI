

# opacity_flicker_incremental_pairs.py
# Photodiode box toggling every frame (upper-right) + growing set of (image + flicker overlay) pairs in the center box.
# Per level: 10s static (no flicker) -> 100s flicker -> 10s baseline (no flicker), then add one more pair, up to 10.
# %%
import os
import sys
import platform
import psutil
import numpy as np
from pathlib import Path
from psychopy import visual, core, event

# %%
# %%
# -------------------- helpers --------------------
def _safe_set_priority():
    """Best-effort process priority tweaks to reduce dropped frames."""
    try:
        p = psutil.Process(os.getpid())
        if hasattr(p, "cpu_affinity"):
            p.cpu_affinity(list(range(psutil.cpu_count())))
        sysname = platform.system()
        if sysname == "Windows":
            p.nice(psutil.HIGH_PRIORITY_CLASS)
        elif sysname in ("Linux", "Darwin"):
            try:
                p.nice(-10)

            except Exception:
                pass
    except Exception:
        pass


def _grid_in_box(n, box_w=600, box_h=800):
    """
    Return exactly n (x,y) positions centered within a box_w × box_h region (pixels).
    Positions are laid out in a roughly square grid (cols × rows) and evenly spaced over the box area.
    """
    cols = int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))
    xs = [0.0] if cols == 1 else list(np.linspace(-box_w / 2.0, box_w / 2.0, cols))
    ys = [0.0] if rows == 1 else list(np.linspace(+box_h / 2.0, -box_h / 2.0, rows))
    pos = []
    k = 0
    for r in range(rows):
        for c in range(cols):
            if k >= n:
                break
            pos.append((xs[c], ys[r]))
            k += 1
    return pos

# PsychoPy Rect can raise a numpy 2.0 compatibility error in some versions.
def _make_rect(win, width, height, pos, fillColor, lineColor, colorSpace, opacity=1.0):
    try:
        return visual.Rect(
            win, width=width, height=height,
            fillColor=fillColor, lineColor=lineColor, colorSpace=colorSpace,
            pos=pos, opacity=opacity
        )
    except ValueError:
        verts = [(-0.5, 0.5), (0.5, 0.5), (0.5, -0.5), (-0.5, -0.5)]
        line_width = 0 if lineColor is None else 1
        return visual.ShapeStim(
            win,
            vertices=verts,
            size=(width, height),
            pos=pos,
            fillColor=fillColor,
            lineColor=lineColor,
            fillColorSpace=colorSpace,
            lineColorSpace=colorSpace,
            opacity=opacity,
            lineWidth=line_width,
            closeShape=True
        )


def _find_stimulus_root(start_dir: Path) -> Path:
    for candidate in [start_dir, *start_dir.parents]:
        if (candidate / "assets" / "images").is_dir():
            return candidate
    raise FileNotFoundError(
        f"Could not locate stimulus root from {start_dir}. Expected assets/images folder."
    )

# %% 
# -------------------- main --------------------
def main():
    # -------------------- user params --------------------
    # Window
    win_size       = [1920, 1080]
    fullscr        = True

    # Path to the image that will be replicated under each overlay
    stimulus_root  = _find_stimulus_root(Path(__file__).resolve().parent)
    img_path       = stimulus_root / "assets" / "images" / "capybara.png"

    # Photodiode box (upper-right, alternates every frame)
    diode_size     = 80     # px
    diode_margin   = 0      # px from window edge

    # Center stimulus region (fixed, centered)
    BOX_W          = 500    # width  in pixels
    BOX_H          = 500    # height in pixels

    # Tile content sizing (image + overlay share the same size)
    tile_w         = 120    # px
    tile_h         = 120    # px

    # Number of tiles (image+overlay pairs) and overlay flicker params
    max_stimuli    = 10
    opacity_min    = 0.00   # 0..1
    opacity_max    = 0.8   # 0..1
    freqs_hz       = np.arange(14,34,2)  # count 10
    randomize_phase = True


    # Block design per level
    static_sec     = 10.0   # no flicker
    flicker_sec    = 200.0
    baseline_sec   = 0.0   # no flicker

    # Quit keys
    quit_keys      = ["escape", "q"]

    # -----------------------------------------------------

    # Validate / normalize frequencies
    if isinstance(freqs_hz, (int, float)):
        freqs_all = np.full(max_stimuli, float(freqs_hz), dtype=float)
    else:
        freqs_all = np.array(freqs_hz, dtype=float)
        if freqs_all.size < max_stimuli:
            freqs_all = np.concatenate([freqs_all, np.full(max_stimuli - freqs_all.size, freqs_all[-1])])

    # Priority (best-effort)
    _safe_set_priority()

    # Window
    win = visual.Window(
        size=win_size, fullscr=fullscr, screen=0,
        color=(200, 200, 200), colorSpace="rgb255",
        units="pix", allowGUI=False
    )
    win.mouseVisible = False

    # Frame tracking
    win.recordFrameIntervals = True
    est_hz = win.getActualFrameRate(nIdentical=20, nMaxFrames=150, nWarmUpFrames=30, threshold=1) or 60.0
    win.refreshThreshold = (1.0 / est_hz) * 1.5

    # Photodiode box position (upper-right)
    half_w, half_h = win.size[0] / 2.0, win.size[1] / 2.0
    diode_x = -half_w + diode_margin + diode_size / 2.0
    diode_y = half_h - diode_margin - diode_size / 2.0
    diode = _make_rect(
        win, width=diode_size, height=diode_size,
        fillColor=(0, 0, 0), lineColor=None, colorSpace="rgb255",
        pos=(diode_x, diode_y), opacity=1.0
    )
    diode_state_white = False  # toggled each frame

    # Positions for up to max_stimuli tiles within the 600×800 box
    positions_all = _grid_in_box(max_stimuli, box_w=BOX_W, box_h=BOX_H)

    # Build tiles: each tile = (ImageStim, Rect overlay)
    tiles = []
    img_exists = img_path.is_file()
    if not img_exists:
        print(f"WARNING: image not found at {img_path}. Using gray fill instead of images.")

    for i in range(max_stimuli):
        pos = positions_all[i]

        if img_exists:
            img = visual.ImageStim(
                win,
                image=str(img_path),
                size=(tile_w, tile_h),
                pos=pos,
                units="pix",
                interpolate=True
            )
        else:
            # Fallback: plain rectangle to mimic an image tile
            img = _make_rect(
                win, width=tile_w, height=tile_h,
                fillColor=(200, 200, 200), lineColor=(80, 80, 80),
                colorSpace="rgb255", pos=pos, opacity=1.0
            )

        overlay = _make_rect(
            win, width=tile_w, height=tile_h,
            fillColor=(0, 0, 0), lineColor=None, colorSpace="rgb255",
            pos=pos, opacity=opacity_min  # start at min opacity
        )
        tiles.append((img, overlay))

    # Optional: a faint guide box to visualize the stimulus region
    # guide_box = _make_rect(win, width=BOX_W, height=BOX_H, fillColor=None, lineColor=(255, 0, 0), colorSpace="rgb255", pos=(0, 0), opacity=1.0)

    # Phases for sine opacity
    phases_all = np.random.uniform(0, 2*np.pi, size=max_stimuli) if randomize_phase else np.zeros(max_stimuli)

    # Precompute constants
    two_pi = 2.0 * np.pi
    op_span = (opacity_max - opacity_min)

    # Timing
    trial_clock = core.Clock()
    frameN = 0
    win.flip()
    trial_clock.reset()

    def draw_tiles(n, mode, t0):
        """
        Draw the first n tiles. For each tile: draw image first, then overlay.
        mode: 'static' | 'flicker' | 'baseline'
        """
        if n <= 0:
            return
        if mode == "flicker":
            t = trial_clock.getTime() - t0
            freqs = freqs_all[:n]
            phases = phases_all[:n]
            op = opacity_min + op_span * (0.5 + 0.5 * np.sin(two_pi * freqs * t + phases))
            for i in range(n):
                img, overlay = tiles[i]
                img.draw()
                overlay.opacity = float(op[i])
                overlay.draw()
        else:
            # static / baseline = overlay at min opacity (i.e., images visible, no flicker)
            for i in range(n):
                img, overlay = tiles[i]
                img.draw()
                overlay.opacity = opacity_min
                overlay.draw()

    def draw_frame(block_mode, active_pairs, t0):
        """
        One frame:
          - draw the first 'active_pairs' tiles (image + overlay)
          - flip the photodiode square
          - win.flip()
        """
        nonlocal frameN, diode_state_white

        # guide_box.draw()  # uncomment to visualize the region
        draw_tiles(active_pairs, block_mode, t0)

        # Photodiode box (toggle every frame)
        diode_state_white = not diode_state_white
        diode.fillColor = (255, 255, 255) if diode_state_white else (0, 0, 0)
        diode.draw()

        win.flip()
        frameN += 1
        if event.getKeys(keyList=quit_keys):
            return False
        return True

    # ------MAIN--------- block loop: grow number of image+overlay pairs ---------------
    levels_completed = 0
    for active_n in range(1, max_stimuli + 1):
        # --- Static (no flicker), show previous count ---
        t0 = trial_clock.getTime()
        while trial_clock.getTime() < t0 + static_sec:
            if not draw_frame("static", active_n - 1, t0):
                elapsed = trial_clock.getTime()
                print(f"Early exit during static. Elapsed: {elapsed:.3f}s")
                win.close(); core.quit(); return

        # --- Flicker with active_n pairs ---
        t0 = trial_clock.getTime()
        while trial_clock.getTime() < t0 + flicker_sec:
            if not draw_frame("flicker", active_n, t0):
                elapsed = trial_clock.getTime()
                print(f"Early exit during flicker. Elapsed: {elapsed:.3f}s")
                win.close(); core.quit(); return

        # --- Baseline (no flicker) with active_n pairs ---
        t0 = trial_clock.getTime()
        while trial_clock.getTime() < t0 + baseline_sec:
            if not draw_frame("baseline", active_n, t0):
                elapsed = trial_clock.getTime()
                print(f"Early exit during baseline. Elapsed: {elapsed:.3f}s")
                win.close(); core.quit(); return

        levels_completed += 1

    # -------------------- reporting --------------------
    elapsed = trial_clock.getTime()
    if elapsed > 0:
        mean_fps = frameN / elapsed
        print(f"Levels: {levels_completed} | Frames: {frameN} | Elapsed: {elapsed:.3f}s | Mean FPS: {mean_fps:.2f}")

    dropped_native = getattr(win, "nDroppedFrames", None)
    if dropped_native is not None:
        print(f"Dropped frames (PsychoPy native): {dropped_native}")

    if getattr(win, "frameIntervals", None):
        intervals = np.array(win.frameIntervals)
        dropped_custom = int(np.sum(intervals > win.refreshThreshold))
        p95 = float(np.percentile(intervals, 95)) if intervals.size else float("nan")
        worst = float(np.max(intervals)) if intervals.size else float("nan")
        print(
            f"Recorded intervals: {len(intervals)} | Custom drops: {dropped_custom} "
            f"| 95th pct: {p95*1000:.2f} ms | Worst: {worst*1000:.2f} ms "
            f"| Threshold: {win.refreshThreshold*1000:.2f} ms"
        )

    # Cleanup
    win.mouseVisible = True
    win.close()
    core.quit()

# %%
if __name__ == "__main__":
    main()

