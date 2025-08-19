import numpy as np
import matplotlib.pyplot as plt

'''
Step 0 — Trial Hann envelope: shows the fade-in and fade-out shape that we’ll apply to contrast only.

Step 1 — Raw sampled sine: the true 64 Hz sine evaluated at each 480 Hz frame (range [-1,1]).

Step 2 — Map to [0,1]: shift/scale into unit range for unified processing.

Step 3 — Contrast around 0.5: subtract 0.5 so the signal is a deviation around mid-gray (range [-0.5,0.5]).

Step 4 — Apply Hann to contrast: multiply the deviation by the envelope. Mean stays at 0.5.

Step 5 — Back to [0,1]: add 0.5 again, preserving the mean while tapering contrast.

Step 6 — Luminance mapping: convert to real display units via L = lb + (hb-lb)*map01.
'''
# --- Parameters (match your PTB pipeline) ---
fs = 480.0                # refresh rate (Hz)
f  = 64.0                 # carrier frequency (Hz) -- change as you like
duration_s = 1.0          # seconds
T = int(round(fs * duration_s))  # total frames
t = np.arange(T) / fs     # frame times

# Luminance bounds (0..255 8-bit)
lb, hb = 50, 195

# Trial-level Hann taper length (frames) on each edge
tf = 96                   # e.g., ~200 ms at 480 Hz
tf = max(0, min(tf, T//2))

# --- Step 0: Build Hann envelope (0->1 fade-in and 1->0 fade-out) ---
env = np.ones(T)
if tf > 0:
    w = np.hanning(2*tf)            # 0..1..0 across 2*tf points
    env[:tf] = w[:tf]               # fade-in
    env[-tf:] = w[tf:]              # fade-out

# --- Step 1: Raw sampled sine in [-1, 1] ---
sine_pm1 = np.sin(2*np.pi*f*t)      # [-1, +1]

# --- Step 2: Map sine to [0, 1] (unit control signal) ---
sine_01 = 0.5 + 0.5*sine_pm1        # [0, 1]

# --- Step 3: Convert to zero-mean contrast around mid-gray ---
contrast = sine_01 - 0.5            # [-0.5, +0.5]

# --- Step 4: Apply trial envelope to contrast only (keep mean) ---
contrast_tapered = contrast * env   # taper the deviation

# --- Step 5: Return to [0, 1] after taper ---
map01 = 0.5 + contrast_tapered      # [0, 1]

# --- Step 6: Map to physical luminance (0..255) ---
L = lb + (hb - lb) * map01          # luminance per frame


# -------- VISUALIZE (one chart per step) --------

# Envelope alone
plt.figure(figsize=(8,3))
plt.plot(env)
plt.title('Step 0 — Trial Hann envelope (fade-in/out)')
plt.xlabel('Frame')
plt.ylabel('Envelope (0..1)')
plt.tight_layout()
plt.show()

# Raw sine
plt.figure(figsize=(8,3))
plt.plot(sine_pm1)
plt.title('Step 1 — Raw sampled sine in [-1, 1]')
plt.xlabel('Frame')
plt.ylabel('Amplitude')
plt.tight_layout()
plt.show()

# Sine mapped to [0,1]
plt.figure(figsize=(8,3))
plt.plot(sine_01)
plt.title('Step 2 — Sine mapped to [0, 1]')
plt.xlabel('Frame')
plt.ylabel('map01')
plt.tight_layout()
plt.show()

# Zero-mean contrast
plt.figure(figsize=(8,3))
plt.plot(contrast)
plt.title('Step 3 — Contrast around 0.5 ([-0.5, 0.5])')
plt.xlabel('Frame')
plt.ylabel('Contrast')
plt.tight_layout()
plt.show()

# Tapered contrast
plt.figure(figsize=(8,3))
plt.plot(contrast_tapered)
plt.title('Step 4 — Contrast with Hann taper applied')
plt.xlabel('Frame')
plt.ylabel('Contrast (tapered)')
plt.tight_layout()
plt.show()

# Back to [0,1] after taper
plt.figure(figsize=(8,3))
plt.plot(map01)
plt.title('Step 5 — map01 after taper (still centered at 0.5)')
plt.xlabel('Frame')
plt.ylabel('map01 (0..1)')
plt.tight_layout()
plt.show()

# Final luminance
plt.figure(figsize=(8,3))
plt.plot(L)
plt.title('Step 6 — Final per-frame luminance')
plt.xlabel('Frame')
plt.ylabel('Luminance (0..255 scale)')
plt.tight_layout()
plt.show()
