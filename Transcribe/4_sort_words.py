import pandas as pd

# Path to the Excel file produced earlier (recognized_transcription_all.xlsx)
input_excel = "recognized_transcription_all.xlsx"

# Read the Excel file into a DataFrame
df = pd.read_excel(input_excel)

# Ensure the "word" column is treated as string and remove any leading/trailing whitespace
df["word"] = df["word"].astype(str).str.strip()

# Aggregate the words and count occurrences
word_counts = df.groupby("word").size().reset_index(name="count")

# Sort the results by count in descending order
word_counts = word_counts.sort_values(by="count", ascending=False)

# Save the aggregated counts to a new Excel file
output_excel = "word_counts.xlsx"
word_counts.to_excel(output_excel, index=False)

print(f"Aggregated word counts saved to {output_excel}")
