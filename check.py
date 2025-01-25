import seaborn as sns
import matplotlib.pyplot as plt

# Load the 'penguins' dataset
df = sns.load_dataset('penguins')

# Create a scatter plot
sns.scatterplot(data=df, x='bill_length_mm', y='bill_depth_mm', hue='species')
plt.title('Scatter Plot of Bill Length vs Bill Depth')
plt.xlabel('Bill Length (mm)')
plt.ylabel('Bill Depth (mm)')
plt.show()