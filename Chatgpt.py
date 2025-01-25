import matplotlib.pyplot as plt
import numpy as np

# Data
city = ['New York', 'Tokyo', 'London', 'Sydney', 'Cape Town', 'Rio de Janeiro', 'Paris', 'Moscow']
country = ['USA', 'Japan', 'UK', 'Australia', 'South Africa', 'Brazil', 'France', 'Russia']
average_temperature = [12.5, 16.3, 11.1, 21.3, 17.2, 23.4, 13.9, 5.2]
annual_rainfall = [1200, 1500, 715, 850, 500, 1200, 640, 300]

# Create a scatter plot with heatmap effect using rainfall as the color
plt.figure(figsize=(10, 6))

# Normalize the rainfall for colormap
norm = plt.Normalize(min(annual_rainfall), max(annual_rainfall))
cmap = plt.cm.Blues  # Use a colormap like Blues or any other you prefer

# Create the scatter plot
scatter = plt.scatter(annual_rainfall, average_temperature, c=annual_rainfall, cmap=cmap, s=100, edgecolors='black', norm=norm)

# Annotate each point with both city and country names
for i, (city_name, country_name) in enumerate(zip(city, country)):
    label = f'{city_name}, {country_name}'
    plt.text(annual_rainfall[i] + 20, average_temperature[i] + 0.2, label, fontsize=9)

# Add labels and title
plt.xlabel('Annual Rainfall (mm)')
plt.ylabel('Average Temperature (°C)')
plt.title('Average Temperature vs Annual Rainfall (Color-Encoded Rainfall)')

# Add colorbar to represent the rainfall scale
plt.colorbar(scatter, label='Annual Rainfall (mm)')

# Show the plot
plt.tight_layout()
plt.savefig('/Users/subasrees/Desktop/fig5.png')

plt.show()
