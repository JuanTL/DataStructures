import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler # Good practice for distance-based metrics

def visualize_dbscan_results(filename="dbscan_results_brute_force.csv"):
    """
    Reads DBSCAN results from a CSV file, visualizes clusters,
    and calculates the Silhouette Score.
    """
    try:
        # Read the data, skipping the header line (starts with '#')
        df = pd.read_csv(filename, comment='#')
        
        # Assign column names based on the expected format
        df.columns = ['x', 'y', 'cluster_id']

        # --- Silhouette Score Calculation ---
        # Exclude noise points (cluster_id = 0) for silhouette score calculation
        non_noise_df = df[df['cluster_id'] != 0]

        # Check if there are enough points and clusters to calculate Silhouette Score
        if len(non_noise_df) < 2:
            silhouette_avg = "N/A (Not enough non-noise points)"
        else:
            # Extract features (X) and labels (y) for non-noise points
            # It's good practice to scale data for distance-based metrics,
            # especially if dimensions have very different scales.
            # For your current data (e.g., if all dimensions are in similar ranges),
            # StandardScaler might not be strictly necessary, but it's good to be aware.
            X = non_noise_df[['x', 'y']].values 
            labels = non_noise_df['cluster_id'].values

            # The Silhouette Score requires at least 2 clusters to be meaningful
            # after excluding noise. Check unique labels count.
            if len(np.unique(labels)) < 2:
                silhouette_avg = "N/A (Less than 2 non-noise clusters)"
            else:
                # Calculate the Silhouette Score
                # The metric uses Euclidean distance by default.
                silhouette_avg = silhouette_score(X, labels)
        
        print(f"\n--- Results for {filename} ---")
        print(f"Silhouette Score (excluding noise): {silhouette_avg:.4f}" if isinstance(silhouette_avg, float) else silhouette_avg)


        # --- Visualization ---
        unique_cluster_ids = df['cluster_id'].unique()
        # Sort cluster IDs to ensure consistent plotting order, noise (0) usually first.
        unique_cluster_ids.sort() 

        # Count actual clusters (excluding noise)
        num_actual_clusters = len([c for c in unique_cluster_ids if c != 0])

        plt.figure(figsize=(10, 8))

        # Define a colormap. If num_actual_clusters is 0, cmap will be empty, handle gracefully.
        cmap = plt.cm.get_cmap('viridis', num_actual_clusters) if num_actual_clusters > 0 else None

        cluster_color_map = {} # Map cluster ID to a color
        current_color_idx = 0

        # Plot noise points first (cluster_id = 0)
        noise_points = df[df['cluster_id'] == 0]
        if not noise_points.empty:
            plt.scatter(noise_points['x'], noise_points['y'], 
                        color='gray', marker='x', s=50, label='Noise (Cluster 0)', alpha=0.6)

        # Plot each true cluster
        for cluster_id in unique_cluster_ids:
            if cluster_id == 0: # Skip noise, already plotted
                continue
            
            cluster_points = df[df['cluster_id'] == cluster_id]
            if not cluster_points.empty:
                # Assign a color from the colormap
                color = cmap(current_color_idx) if cmap else 'black' # Fallback for 0 clusters
                cluster_color_map[cluster_id] = color
                plt.scatter(cluster_points['x'], cluster_points['y'], 
                            color=color, marker='o', s=70, 
                            label=f'Cluster {cluster_id}', alpha=0.8)
                current_color_idx += 1

        plt.title(f'DBSCAN Clustering Results (File: {filename.split("/")[-1]})')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.show()

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found. Make sure your C++ program generated it.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Visualize and evaluate the KD-Tree accelerated results
    visualize_dbscan_results("dbscan_results_kd_tree-noisy_moons.csv")
    
    # Visualize and evaluate the Brute-Force results
    visualize_dbscan_results("dbscan_results_brute_force-noisy_moons.csv")