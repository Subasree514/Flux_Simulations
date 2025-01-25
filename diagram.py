import numpy as np
import matplotlib.pyplot as plt

# Data
barWidth = 0.25
bars1 = [12, 30, 1, 8, 22]
bars2 = [28, 6, 16, 5, 10]
bars3 = [29, 3, 24, 25, 17]

# Bar positions
r = np.arange(len(bars1))
r2 = r + barWidth
r3 = r2 + barWidth

# Plotting
fig, ax = plt.subplots(dpi=300)
ax.bar(r, bars1, color='#7f6d5f', width=barWidth, edgecolor='white', label='var1')
ax.bar(r2, bars2, color='#557f2d', width=barWidth, edgecolor='white', label='var2')
ax.bar(r3, bars3, color='#2d7f5e', width=barWidth, edgecolor='white', label='var3')

# Xticks
ax.set_xlabel('group', fontweight='bold')
ax.set_xticks(r + barWidth)
ax.set_xticklabels(['A', 'B', 'C', 'D', 'E'])

# Legend and show
ax.legend()
plt.show()