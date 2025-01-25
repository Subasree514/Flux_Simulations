import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import Rectangle

# Reactions(source, target, flux, reaction_name)
reactions = [
    ("A", "B",2, "R1"), 
    ("B", "F", 2, "R2"), 
    ("C", "B", 2,"R3"), 
    ("D", "C",12, "R4"),
    ("F", "G",2, "R5"), 
    ("G", "E",20, "R6"),
    ("E", "H",18, "R7"), 
    ("H", "I",10, "R8"),
    ("E", "D",2, "R9"),
]

# Define Custom Positions for Metabolites
positions = {
    "A": (1, 1),  # Place 'A' at grid position (x,y)
    "B": (2, 1),  # Place 'B' at grid position (x,y)
    "C": (3, 0),  # Place 'C' at grid position (x,y)
    "D": (4, 0),   # Place 'D' at grid position (x,y)
    "E": (5, 1),  # Place 'A' at grid position (x,y)
    "F": (6, 0),  # Place 'B' at grid position (x,y)
    "G": (4, 2),  # Place 'C' at grid position (x,y)
    "H": (5, 2),   # Place 'D' at grid position (x,y)
    "I": (6, 2)   # Place 'D' at grid position (x,y)
}

# Create Graph
G = nx.DiGraph()
for source, target, flux, reaction_name in reactions:
    G.add_edge(source, target, weight=flux, reaction_name=reaction_name)

# Normalize Edge Weights for Arrow Sizes
max_flux = max(flux for _, _, flux, _ in reactions)
min_flux = min(flux for _, _, flux, _ in reactions)
for u, v, d in G.edges(data=True):
    d['normalized_weight'] = (d['weight'] - min_flux) / (max_flux - min_flux) + 0.1  # Avoid zero size

# Draw Graph
fig, ax = plt.subplots(figsize=(10, 6))

# Draw nodes 
nx.draw_networkx_nodes(G, positions, node_size=900, node_color="lightblue", ax=ax)

# Draw edges with scaled widths
for idx, (u, v, d) in enumerate(G.edges(data=True)):
    nx.draw_networkx_edges(
        G, positions,
        edgelist=[(u, v)],
        width=d['normalized_weight'] * 5,  # Scale factor for visualization
        edge_color="gray",
        arrowstyle="->",
        arrowsize=20,
        ax=ax
    )
    
    # Get the custom reaction name
    reaction_name = d['reaction_name']

    # Midpoint calculation for placing the label
    mid_point = ((positions[u][0] + positions[v][0]) / 2, (positions[u][1] + positions[v][1]) / 2)
    
    # Add white rectangle background for the label
    rect_width = 0.2  # Width of the rectangle
    rect_height = 0.08  # Height of the rectangle
    ax.add_patch(Rectangle((mid_point[0] - rect_width / 2, mid_point[1] - rect_height / 2),
                           rect_width, rect_height, facecolor='white', edgecolor='black', lw=1))
    
    # Add text label 
    ax.text(mid_point[0], mid_point[1], reaction_name, fontsize=10, color='blue', ha='center', va='center')

# Draw labels at specified positions for metabolites
nx.draw_networkx_labels(G, positions, font_size=12, font_color="black", ax=ax)

plt.title("Metabolic Pathway")
plt.axis("on")  

# Save the figure with a transparent background
plt.savefig('metabolic_pathway.png', transparent=True)

# Show the plot
plt.show()


