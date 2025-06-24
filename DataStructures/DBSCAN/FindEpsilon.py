import pandas as pd

# Load CSV file with no headers
data = pd.read_csv("generated_datasets/varied.csv", header=None)

# Optional: Rename columns
data.columns = ['x', 'y']
from sklearn.neighbors import NearestNeighbors
import numpy as np
import matplotlib.pyplot as plt

k = 5  # Common starting point; also try 4, 6, etc.
neighbors = NearestNeighbors(n_neighbors=k)
neighbors_fit = neighbors.fit(data)
distances, indices = neighbors_fit.kneighbors(data)

# Take the k-th nearest distance for each point
k_distances = distances[:, k-1]
# Sort distances to identify the "elbow"
k_distances = np.sort(k_distances)

plt.figure(figsize=(8, 4))
plt.plot(k_distances)
plt.title(f'K-Distance Graph (k = {k})')
plt.xlabel('Points sorted by distance')
plt.ylabel(f'{k}th Nearest Neighbor Distance')
plt.grid(True)
plt.show()
