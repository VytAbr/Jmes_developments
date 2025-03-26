import os
import librosa
import soundfile as sf

# Input and output root folders
input_root = "audio/extracted0325"
output_root = "audio/converted0325"

# Supported input formats
supported_ext = (".wav", ".mp3", ".flac", ".ogg", ".m4a")

for root, dirs, files in os.walk(input_root):
    for file in files:
        if file.lower().endswith(supported_ext):
            input_path = os.path.join(root, file)

            # Create mirrored output path
            relative_path = os.path.relpath(root, input_root)
            output_dir = os.path.join(output_root, relative_path)
            os.makedirs(output_dir, exist_ok=True)

            output_filename = os.path.splitext(file)[0] + "_converted.wav"
            output_path = os.path.join(output_dir, output_filename)

            try:
                # Load and convert audio
                audio, sr = librosa.load(input_path, sr=16000, mono=True)
                sf.write(output_path, audio, sr)
                print(f"✅ Converted: {input_path} -> {output_path}")
            except Exception as e:
                print(f"❌ Failed to convert {input_path}: {e}")
