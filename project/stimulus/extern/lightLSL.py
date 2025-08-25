import serial
import struct
from pylsl import StreamInfo, StreamOutlet

# === Configuratie ===
SERIAL_PORT = 'COM3'         # Pas dit aan aan jouw poort
BAUD_RATE = 1000000
SAMPLES_PER_PACKET = 10      # 10 samples per 1 pakket
BYTES_PER_SAMPLE = 2         # int16 = 2 bytes
PACKET_SIZE = SAMPLES_PER_PACKET * BYTES_PER_SAMPLE

# === LSL Stream-instellingen ===
info = StreamInfo(
    name='LightSensor',
    type='Sensor',
    channel_count=1,                # Eén kanaal (lichtsensor)
    nominal_srate=10000,           # 10 kHz sample rate
    channel_format='int16',        # Dataformaat
    source_id='arduino_light_1'
)
outlet = StreamOutlet(info)

# === Seriële poort openen ===
ser = serial.Serial(SERIAL_PORT, BAUD_RATE, timeout=1)

print("Seriële verbinding geopend, streamen naar LSL...")

try:
    while True:
        # Zoek naar startbyte 0xAA
        start = ser.read(1)
        if start != b'\xAA':
            continue

        # Lees 20 bytes = 10 int16 samples
        raw = ser.read(PACKET_SIZE)
        if len(raw) != PACKET_SIZE:
            print("Onvolledig pakket ontvangen.")
            continue

        try:
            # Ontpak naar 10 int16 samples
            samples = struct.unpack('<' + 'h' * SAMPLES_PER_PACKET, raw)

            # Verstuur elk sample apart naar LSL (10 kHz)
            for s in samples:
                outlet.push_sample([s])
        except struct.error as e:
            print("Fout bij ontpakken:", e)

except KeyboardInterrupt:
    print("\nOnderbroken door gebruiker.")
finally:
    ser.close()
    print("Seriële poort gesloten.")
