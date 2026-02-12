from psychopy import core
from psychopy.hardware import keyboard
from psychopy.visual import Window, Rect
from pylsl import StreamInfo, StreamOutlet
import serial


def send_triggers(value):
    """Send value as hardware trigger and LSL marker."""
    port.write(value)
    outlet.push_sample(value)


BLACK, WHITE = (-1, -1, -1), (1, 1, 1)
defaultKeyboard = keyboard.Keyboard()

# set correct COM port and use baudrate=2000000 for TriggerBox Plus:
port = serial.Serial("COM4",baudrate=2000000)

n_trials = 300

# adjust according to screen resolution and refresh rate:
win_size = 1920, 1080  # window size
stim_size = 100, 100  # stimulus size
n_frames = 30  # number of frames; 60 frames for 60Hz refresh rate is 1s

LOWER_RIGHT = (win_size[0]/2 - stim_size[0]/2, stim_size[1]/2 - win_size[1]/2)
UPPER_LEFT = (stim_size[0]/2 - win_size[0]/2, win_size[1]/2 - stim_size[1]/2)

win = Window(size=win_size, fullscr=True, color=BLACK, colorSpace="rgb", units="pix")

info = StreamInfo(name="LSL_Markers", type="Markers", channel_count=1,
                  channel_format="int32", source_id="LSL_Markers_001")
outlet = StreamOutlet(info)

rect1 = Rect(
    win=win,
    width=stim_size[0],
    height=stim_size[1],
    fillColor=WHITE,
    lineColor=WHITE,
    colorSpace="rgb",
    pos=UPPER_LEFT,
)
rect2 = Rect(
    win=win,
    width=stim_size[0],
    height=stim_size[1],
    fillColor=BLACK,
    lineColor=BLACK,
    colorSpace="rgb",
    pos=UPPER_LEFT,
)

core.wait(10)  # wait before paradigm starts

# send triggers every time an image changes:
for _ in range(n_trials):
    if defaultKeyboard.getKeys(keyList=["escape"]):
        break
    win.callOnFlip(send_triggers, [1])
    for n in range(n_frames):
        rect1.draw()
        win.flip()
    win.callOnFlip(send_triggers, [0])
    for n in range(n_frames):
        rect2.draw()
        win.flip()

port.close()
win.close()
core.quit()
