import pandas as pd

# Path to the failed transcription Excel file
failed_file = "failed_transcription_all.xlsx"

# Read the Excel file into a DataFrame
df_failed = pd.read_excel(failed_file)

# Check that required columns exist
if "start" in df_failed.columns and "end" in df_failed.columns:
    # Calculate duration for each row (in seconds)
    df_failed["duration"] = df_failed["end"] - df_failed["start"]
    
    # Compute total failed duration
    total_failed_duration = df_failed["duration"].sum()
    print("Total failed duration (in seconds):", total_failed_duration)
else:
    print("Error: 'start' and/or 'end' columns not found in the Excel file.")