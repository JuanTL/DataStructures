# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
from sklearn import datasets
import os

# Output directory
output_dir = "generated_datasets"
os.makedirs(output_dir, exist_ok=True)

n_samples = 1500

# Dataset generators (only features)
datasets_list = {
    "noisy_circles": datasets.make_circles(n_samples=n_samples, factor=0.5, noise=0.05, random_state=170)[0],
    "noisy_moons": datasets.make_moons(n_samples=n_samples, noise=0.05, random_state=170)[0],
    "blobs": datasets.make_blobs(n_samples=n_samples, random_state=170)[0],
    "no_structure": np.random.RandomState(170).rand(n_samples, 2),
    "aniso": np.dot(
        datasets.make_blobs(n_samples=n_samples, random_state=170)[0],
        [[0.6, -0.6], [-0.4, 0.8]]
    ),
    "varied": datasets.make_blobs(
        n_samples=n_samples, cluster_std=[1.0, 2.5, 0.5], random_state=170
    )[0]
}

# Save to CSV files (no labels, no headers)
for name, X in datasets_list.items():
    filepath = os.path.join(output_dir, f"{name}.csv")
    np.savetxt(filepath, X, delimiter=",", fmt="%.10f")
    print(f"Saved {name} to {filepath}")
