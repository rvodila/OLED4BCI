

### Create an LSL _Marker_ outlet once (at experiment start)

In PsychoPy you typically create **one outlet** and then push markers during the task. A widely used structure is:

from pylsl import StreamInfo, StreamOutlet
from psychopy import core

%% Create a marker stream (irregular rate => nominal_srate=0)%%
info = StreamInfo(
    name='PsychoPyMarkers', 
    type='Markers',
    channel_count=1,
    nominal_srate=0,
    channel_format='string',
    source_id='my_experiment_unique_id'
)
outlet = StreamOutlet(info)

def lsl_marker(text: str):
    # Use PsychoPy's clock for a stable timestamp in your script
    outlet.push_sample([text], timestamp=core.getTime())``

### Why `callOnFlip` matters

Instead of calling `lsl_marker()` “whenever”, we will do:

`win.callOnFlip(lsl_marker, "something") win.flip()`
This is to schdeule marker, which egts triggered at window flip

## Where to place markers (a practical checklist)

### 1) Experiment start / end

- After the participant presses “start” and right before the first real display flip
    
- At the final “done” flip
    

### 2) Block begin / end

- **Block begin:** on the flip that shows the first fixation/stimulus of the block
    
- **Block end:** on the flip that ends the block (often start of ISI or end screen)
    

### 3) Trial begin / end

- **Trial begin:** on the flip that shows the first trial frame (often fixation or stimulus)
    
- **Trial end:** on the final flip of the trial (or immediately after last frame, aligned to that flip)
    

### 4) Stimulus onset (the “critical” marker)

- Put the marker on the **exact flip** that first displays the stimulus state you care about (e.g., flicker begins, target appears, mask begins).
    

### 5) Dynamic onsets (your “phase-shift / staggered start” case)

If each of 4 stimuli begins flickering at different frames:

- Emit **one marker per stimulus** on _its_ onset flip: `flicker_onset stim=0`, `stim=1`, etc.
    
- Do it when your code first switches that stimulus from “off/static” to “flickering”
    


>There are multiple clocks:

- PsychoPy’s clocks / Python timing (your code execution)  
- LSL’s clock domain and time-correction system
- The physical display flip timing (refresh-locked)

**MVP Integration**
send a brief **marker “ping”** early to confirm the receiving setup is actually capturing markers before you start.  
For true timing verification, you need an independent reference (e.g., photodiode channel or hardware trigger comparison) to measure latency/jitter end-to-end.


Brainproducts Tutorial(https://pressrelease.brainproducts.com/timing-verification.com)

![[Pasted image 20260209204830.png]]