import os
import wave
import json
from vosk import Model, KaldiRecognizer
import pandas as pd
import time

start_time = time.time()

# Load model
# model = Model("models/vosk-model-en-us-0.22") # US English 2
# model = Model("models/vosk-model-en-in-0.5") # Indian English 1
# model = Model("models/vosk-model-small-en-in-0.4") # Indian English, small 3
model = Model("models/vosk-model-small-en-us-0.15") # US English, small 4

# Root folder containing WAV files (and subfolders)
input_root = "audio/converted0325"

# Lists to accumulate results from all files
recognized_all = []
failed_all = []

# Walk through all subfolders
for root, dirs, files in os.walk(input_root):
    for file in files:
        # Process only .wav files
        if file.lower().endswith(".wav"):
            file_path = os.path.join(root, file)
            print("Processing:", file_path)
            try:
                wf = wave.open(file_path, "rb")
            except Exception as e:
                print("Error opening file:", file_path, e)
                continue

            rec = KaldiRecognizer(model, wf.getframerate())
            rec.SetWords(True)

            recognized_entries = []

            # Process the audio file in chunks
            while True:
                data = wf.readframes(4000)
                if len(data) == 0:
                    break
                if rec.AcceptWaveform(data):
                    res = json.loads(rec.Result())
                    if "result" in res:
                        for word in res["result"]:
                            recognized_entries.append({
                                "file": file_path,
                                "word": word["word"],
                                "start": word["start"],
                                "end": word["end"]
                            })
            # Process any remaining audio
            final_res = json.loads(rec.FinalResult())
            if "result" in final_res:
                for word in final_res["result"]:
                    recognized_entries.append({
                        "file": file_path,
                        "word": word["word"],
                        "start": word["start"],
                        "end": word["end"]
                    })

            # Get total duration of the audio file
            total_duration = wf.getnframes() / wf.getframerate()
            wf.close()

            # Sort recognized entries by start time
            recognized_entries.sort(key=lambda x: x["start"])
            failed_entries = []
            gap_threshold = 0.1  # in seconds

            # Compute gap before the first word
            if recognized_entries and recognized_entries[0]["start"] > gap_threshold:
                failed_entries.append({
                    "file": file_path,
                    "start": 0.0,
                    "end": recognized_entries[0]["start"],
                    "text": ""
                })
            # Compute gaps between recognized words
            for i in range(len(recognized_entries) - 1):
                cur_end = recognized_entries[i]["end"]
                next_start = recognized_entries[i + 1]["start"]
                if next_start - cur_end > gap_threshold:
                    failed_entries.append({
                        "file": file_path,
                        "start": cur_end,
                        "end": next_start,
                        "text": ""
                    })
            # Compute gap after the last word
            if recognized_entries and (total_duration - recognized_entries[-1]["end"]) > gap_threshold:
                failed_entries.append({
                    "file": file_path,
                    "start": recognized_entries[-1]["end"],
                    "end": total_duration,
                    "text": ""
                })

            # If no words recognized at all, mark whole file as failed
            if not recognized_entries:
                failed_entries.append({
                    "file": file_path,
                    "start": 0.0,
                    "end": total_duration,
                    "text": ""
                })

            # Accumulate results for all files
            recognized_all.extend(recognized_entries)
            failed_all.extend(failed_entries)

# Save results to Excel files
df_recognized = pd.DataFrame(recognized_all)
df_failed = pd.DataFrame(failed_all)

df_recognized.to_excel("recognized_transcription_all.xlsx", index=False)
df_failed.to_excel("failed_transcription_all.xlsx", index=False)

print("Processing complete.")
print("Recognized words saved to recognized_transcription_all.xlsx")
print("Failed segments saved to failed_transcription_all.xlsx")

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Elapsed time: {elapsed_time:.4f} seconds")