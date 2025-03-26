import pandas as pd

# Paths to the Excel files
word_counts_file = "word_counts.xlsx"    # Aggregated word counts from previous step
keywords_file = "keywords.xlsx"          # Excel file containing a list of keywords (with a column 'keyword')

# Read the aggregated word counts
df_words = pd.read_excel(word_counts_file)

# Read the keywords file; assumes the keywords are in a column named 'keyword'
df_keywords = pd.read_excel(keywords_file)

# Convert both columns to lowercase and strip any whitespace for reliable matching
df_words["word"] = df_words["word"].astype(str).str.lower().str.strip()
df_keywords["keyword"] = df_keywords["keyword"].astype(str).str.lower().str.strip()

# Filter: only keep rows where the word is in the keyword list
filtered_df = df_words[df_words["word"].isin(df_keywords["keyword"])]

# Save the filtered word counts to a new Excel file
output_file = "filtered_word_counts.xlsx"
filtered_df.to_excel(output_file, index=False)

print(f"Filtered word counts saved to {output_file}")
