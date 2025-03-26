import os
import wave

input_folder = "audio/converted0325"  # Change to your folder path
total_duration = 0.0  # Total duration in seconds

# Walk through all subfolders
for root, dirs, files in os.walk(input_folder):
    for file in files:
        if file.lower().endswith(".wav"):
            file_path = os.path.join(root, file)
            try:
                with wave.open(file_path, "rb") as wf:
                    frames = wf.getnframes()
                    rate = wf.getframerate()
                    duration = frames / float(rate)
                    total_duration += duration
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

print("Total duration of audio files (in seconds):", total_duration)
