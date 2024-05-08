import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file into a DataFrame
df = pd.read_csv("RGCAULevels.csv")

# Extract chromosome numbers
df['Chromosome'] = df['Gene'].str[1:3]

# Group by chromosome and calculate mean Polar% and Non-Polar%
grouped = df.groupby('Chromosome').agg({'Polar%': 'mean', 'Non-Polar%': 'mean'}).reset_index()

# Create bar graphs and save to PDF
for index, row in grouped.iterrows():
    plt.figure(figsize=(8, 6))
    plt.bar(['Polar%', 'Non-Polar%'], [row['Polar%'], row['Non-Polar%']], color=['blue', 'orange'])
    plt.title(f'Polar% and Non-Polar% for Chromosome {row["Chromosome"]}')
    plt.xlabel('Amino Acid Type')
    plt.ylabel('Percentage')
    plt.ylim(0, 100)
    plt.savefig(f"Chromosome_{row['Chromosome']}_Stats.pdf", format="pdf", bbox_inches="tight")
    plt.show()