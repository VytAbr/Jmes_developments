import sys
sys.path.append(r"D:\VGTUstuff\Projektai\Jahan mes\venv\Lib\site-packages")
import os
from moviepy import VideoFileClip

# Defined paths
video_root_folder = "Interviu"     
output_root_folder = "audio/extracted0325"

# Walk all subfolders
for root, dirs, files in os.walk(video_root_folder):
    for file in files:
        if file.lower().endswith(".webm"):
            video_path = os.path.join(root, file)

            # Create mirrored folder structure in output
            relative_path = os.path.relpath(root, video_root_folder)
            output_folder = os.path.join(output_root_folder, relative_path)
            os.makedirs(output_folder, exist_ok=True)

            # Define output audio path
            audio_filename = os.path.splitext(file)[0] + ".wav"  # or ".wav"
            audio_path = os.path.join(output_folder, audio_filename)

            print(f"Extracting audio from: {video_path}")
            try:
                video_clip = VideoFileClip(video_path)
                video_clip.audio.write_audiofile(audio_path)
                video_clip.close()
            except Exception as e:
                print(f"Failed to extract audio from {video_path}: {e}")

print("✅ All audio extracted successfully.")
