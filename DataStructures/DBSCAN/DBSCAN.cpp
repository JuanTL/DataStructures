// DBSCAN.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "KDTree.h"
// --- DBSCAN Helper ---
enum PointType {
    UNCLASSIFIED,
    NOISE,
    CORE,
    BORDER
};
// --- DBSCAN ---
class DBSCAN {
private:
    double eps;
    size_t minPts;
    const KDTree& kdTree; // Reference to the KDTree instance

    std::vector<int> point_clusters; // Cluster ID for each point (-1 for noise, -2 for unclassified)
    std::vector<bool> visited_points; // True if point has been visited

    // Neighbors using KDTree
    std::vector<int> getNeighbors(const Point& p, const std::vector<Point>& all_points_ref) const {
        std::vector<Point> neighbors_in_range = kdTree.findInRange(p, eps);
        std::vector<int> neighbor_indices;

        // Convert Point objects back to their original indices
        // This requires either:
        // 1. Storing original indices in the Point struct
        // 2. Passing the original points vector by reference to DBSCAN and iterating
        //    through it to find matching points. The latter is less efficient.
        // For simplicity, let's assume `Point` could potentially store an `id` or
        // we'll search by value (less efficient but demonstrates the concept).
        // A better approach is to pass indices around, not copies of points.

        // --- BETTER APPROACH: Use indices for DBSCAN points ---
        // Let's modify the Point class slightly or assume points are given with indices
        // For this example, let's assume all_points_ref is the original data
        // and we need to find the index of each neighbor_in_range point.
        // This is inefficient if many points are identical.
        // A more robust solution would be to pass `const std::vector<Point>& data` to DBSCAN
        // and have `kdTree.findInRange` return `std::vector<int>` (indices) instead of `Point` copies.

        // For now, let's just find the indices by value comparison.
        // In a real scenario, you'd want findInRange to return indices.
        for (const auto& neighbor_p : neighbors_in_range) {
            for (size_t i = 0; i < all_points_ref.size(); ++i) {
                if (all_points_ref[i].coords == neighbor_p.coords) { // Assuming Point equality is coord equality
                    neighbor_indices.push_back(i);
                    break;
                }
            }
        }
        return neighbor_indices;
    }

    // Neighbors without KDTree
    std::vector<int> getNeighborsBruteForce(const Point& p_query, const std::vector<Point>& data) const {
        std::vector<int> neighbors_indices;
        double squaredEps = eps * eps;

        for (size_t i = 0; i < data.size(); ++i) {
            // Calculate squared distance from query point to current point
            // This requires access to a squaredDistance function. Let's make it a static helper
            // or put it directly here, assuming Point has coords.
            double dist_sq = 0.0;
            for (size_t d = 0; d < p_query.dimensions(); ++d) {
                double diff = p_query.coords[d] - data[i].coords[d];
                dist_sq += diff * diff;
            }

            if (dist_sq <= squaredEps) {
                neighbors_indices.push_back(static_cast<int>(i));
            }
        }
        return neighbors_indices;
    }

public:
    DBSCAN(double epsilon, size_t min_points, const KDTree& tree)
        : eps(epsilon), minPts(min_points), kdTree(tree) {
        if (minPts < 1) {
            throw std::invalid_argument("minPts must be at least 1.");
        }
        if (eps <= 0) {
            throw std::invalid_argument("eps must be greater than 0.");
        }
    }


    // --- KD-Tree DBSCAN ---
    std::vector<int> run(const std::vector<Point>& data) {
        size_t n = data.size();
        point_clusters.assign(n, -1);
        visited_points.assign(n, false);

        int clusterId = 0;

        for (size_t i = 0; i < n; ++i) {
            if (visited_points[i]) {
                continue; // Already processed
            }
            visited_points[i] = true; // Mark as visited
            // This is where the KDTree is used! Use to get neighbors!
            std::vector<int> neighbors = kdTree.findInRangeIndices(data[i], eps); 
            // Take out the point itself if it's included in neighbors
            neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(), [&](int idx) {
                return data[idx].coords == data[i].coords;
                }), neighbors.end());
            if (neighbors.size() + 1 < minPts) { // +1 to include itself
                point_clusters[i] = 0; // Assign to noise cluster (cluster 0)
            }
            else {
                // Core Point? start a new cluster
                point_clusters[i] = ++clusterId; // Assign new cluster ID
                std::queue<int> q;
                for (int neighbor_idx : neighbors) {
                    q.push(neighbor_idx);
                }

                while (!q.empty()) {
                    int current_point_idx = q.front();
                    q.pop();

                    if (!visited_points[current_point_idx]) {
                        visited_points[current_point_idx] = true;
                        std::vector<int> current_neighbors = kdTree.findInRangeIndices(data[current_point_idx], eps);
                        current_neighbors.erase(std::remove_if(current_neighbors.begin(), current_neighbors.end(), [&](int idx) {
                            return data[idx].coords == data[current_point_idx].coords;
                            }), current_neighbors.end());


                        if (current_neighbors.size() + 1 >= minPts) {
                            // Current point is Core Point? add neighbors to queue
                            for (int nn_idx : current_neighbors) {
                                q.push(nn_idx);
                            }
                        }
                    }

                    // Not classified? Assign to current cluster
                    if (point_clusters[current_point_idx] <= 0) { // If unclassified (-1) or noise (0)
                        point_clusters[current_point_idx] = clusterId;
                    }
                }
            }
        }
        // In this vector each element is the cluster ID for the point it belongs.
        return point_clusters; 
    }

    // --- Brute-Force DBSCAN (without KD-Tree) ---
    std::vector<int> runWithoutKdTree(const std::vector<Point>& data) {
        size_t n = data.size();
        point_clusters.assign(n, -1); // -1 for unclassified, 0 for noise
        visited_points.assign(n, false);

        int clusterId = 0; // Cluster IDs will start from 1

        for (size_t i = 0; i < n; ++i) {
            if (visited_points[i]) {
                continue; // Already processed
            }
            visited_points[i] = true; // Mark as visited
            // Get neighbors
            std::vector<int> neighbors_indices = getNeighborsBruteForce(data[i], data);
			// To not query the point itself, we remove it from neighbors_indices
            neighbors_indices.erase(std::remove_if(neighbors_indices.begin(), neighbors_indices.end(), [&](int idx) {
                return idx == static_cast<int>(i);
                }), neighbors_indices.end());

            if (neighbors_indices.size() + 1 < minPts) { // +1 to include itself
                point_clusters[i] = 0; // Assign to noise cluster (cluster 0)
            }
            else {
                // Core Point? Start a new cluster
                point_clusters[i] = ++clusterId; // Assign new cluster ID
                std::queue<int> q;
                for (int neighbor_idx : neighbors_indices) {
                    q.push(neighbor_idx);
                }

                while (!q.empty()) {
                    int current_point_idx = q.front();
                    q.pop();

                    if (!visited_points[current_point_idx]) {
                        visited_points[current_point_idx] = true;

                        std::vector<int> current_neighbors_indices = getNeighborsBruteForce(data[current_point_idx], data);
                        current_neighbors_indices.erase(std::remove_if(current_neighbors_indices.begin(), current_neighbors_indices.end(), [&](int idx) {
                            return idx == current_point_idx;
                            }), current_neighbors_indices.end());

                        if (current_neighbors_indices.size() + 1 >= minPts) {
                            // Current point is also a Core Point? Add its neighbors to queue
                            for (int nn_idx : current_neighbors_indices) {
                                q.push(nn_idx);
                            }
                        }
                    }
                    // If not yet classified -> assign to current cluster
                    if (point_clusters[current_point_idx] <= 0) { // If unclassified (-1) or noise (0)
                        point_clusters[current_point_idx] = clusterId;
                    }
                }
            }
        }
        return point_clusters;
    }

};

// ------------------------------------------------------------- //

// Function to read points from a CSV file
std::vector<Point> readPointsFromCSV(const std::string& filename, size_t& dimensions, 
    size_t max_points = std::numeric_limits<size_t>::max()) {
    std::vector<Point> points;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    std::string line;
    bool first_line = true;
    size_t count = 0;
    while (std::getline(file, line) && count < max_points) {
        std::vector<double> coords;
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            try {
                coords.push_back(std::stod(value));
            }
            catch (const std::exception& e) {
                throw std::runtime_error("Invalid number in CSV: " + value);
            }
        }

        if (first_line) {
            dimensions = coords.size();
            if (dimensions == 0) {
                throw std::runtime_error("Empty point in CSV");
            }
            first_line = false;
        }
        else if (coords.size() != dimensions) {
            throw std::runtime_error("Inconsistent number of dimensions in CSV");
        }

        points.emplace_back(coords);
        ++count;
    }

    file.close();
    if (points.empty()) {
        throw std::runtime_error("No points found in CSV");
    }

    return points;
}
// Function to save DBSCAN results to a file for Python visualization
// DBSCAN doesn't use centroids, so we remove that parameter.
void saveResultsToFile(const std::string& filename, const std::vector<Point>& points,
    const std::vector<int>& cluster_assignments) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open output file: " + filename);
    }

    // Save points with cluster assignments
    // Format: x,y,cluster_id
    // Assuming 2D points for easy visualization. If >2D, you'd save all dims or pick two.
    file << "# Points: x,y,cluster_id\n";
    for (size_t i = 0; i < points.size(); ++i) {
        // Ensure points are at least 2-dimensional for x,y
        if (points[i].dimensions() < 2) {
            throw std::runtime_error("Points must have at least 2 dimensions for x,y visualization.");
        }
        file << points[i].coords[0] << "," << points[i].coords[1] << ","
            << cluster_assignments[i] << "\n";
    }

    // No centroids to save for DBSCAN

    file.close();
}


#include <chrono>

int main() {
    std::string input_filename = "generated_datasets/noisy_moons.csv";  
    std::string output_filename_kd = "dbscan_results_kd_tree-noisy_moons.csv";
    std::string output_filename_brute = "dbscan_results_brute_force-noisy_moons.csv";

    size_t dimensions = 0;
    std::vector<Point> points;

    // DBSCAN Parameters
    double eps = 0.05;
    size_t minPts = 5;
    
    try {
        std::cout << "Loading points from " << input_filename << "...\n";
        auto start_load = std::chrono::high_resolution_clock::now();
        points = readPointsFromCSV(input_filename, dimensions);
        auto end_load = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> load_duration = end_load - start_load;
        std::cout << "Loaded " << points.size() << " points with " << dimensions << " dimensions in "
            << load_duration.count() << " seconds.\n";

        if (dimensions < 2) {
            std::cerr << "Warning: Points have less than 2 dimensions. Visualization might be problematic." << std::endl;
        }

        // --- DBSCAN with KD-Tree ---
        std::cout << "\n--- Running DBSCAN with KD-Tree Acceleration ---\n";
        std::cout << "Building KD-Tree...\n";
        auto start_build_tree = std::chrono::high_resolution_clock::now();
        KDTree kdTree(points, dimensions);
        auto end_build_tree = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> build_tree_duration = end_build_tree - start_build_tree;
        std::cout << "KD-Tree built in " << build_tree_duration.count() << " seconds.\n";

        std::cout << "Running DBSCAN (K-d Tree) with eps = " << eps << ", minPts = " << minPts << "...\n";
        auto start_dbscan_kd = std::chrono::high_resolution_clock::now();
        DBSCAN dbscan_kd(eps, minPts, kdTree); // Pass the built KDTree
        std::vector<int> cluster_assignments_kd = dbscan_kd.run(points);
        auto end_dbscan_kd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dbscan_kd_duration = end_dbscan_kd - start_dbscan_kd;
        std::cout << "DBSCAN (K-d Tree) completed in " << dbscan_kd_duration.count() << " seconds.\n";

        int num_clusters_kd = 0;
        size_t noise_points_kd = 0;
        for (int cluster_id : cluster_assignments_kd) {
            if (cluster_id > num_clusters_kd) {
                num_clusters_kd = cluster_id;
            }
            if (cluster_id == 0) {
                noise_points_kd++;
            }
        }
        std::cout << "DBSCAN (K-d Tree) found " << num_clusters_kd << " clusters and " << noise_points_kd << " noise points.\n";
        std::cout << "Saving KD-Tree results to " << output_filename_kd << "...\n";
        saveResultsToFile(output_filename_kd, points, cluster_assignments_kd);
        std::cout << "Results saved successfully.\n";

        // --- DBSCAN without KD-Tree (Brute Force) ---
        std::cout << "\n--- Running DBSCAN without KD-Tree (Brute Force) ---\n";
        std::cout << "Running DBSCAN (Brute Force) with eps = " << eps << ", minPts = " << minPts << "...\n";
        auto start_dbscan_brute = std::chrono::high_resolution_clock::now();
        DBSCAN dbscan_brute(eps, minPts, kdTree);
        std::vector<int> cluster_assignments_brute = dbscan_brute.runWithoutKdTree(points);
        auto end_dbscan_brute = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dbscan_brute_duration = end_dbscan_brute - start_dbscan_brute;
        std::cout << "DBSCAN (Brute Force) completed in " << dbscan_brute_duration.count() << " seconds.\n";
        int num_clusters_brute = 0;
        size_t noise_points_brute = 0;
        for (int cluster_id : cluster_assignments_brute) {
            if (cluster_id > num_clusters_brute) {
                num_clusters_brute = cluster_id;
            }
            if (cluster_id == 0) {
                noise_points_brute++;
            }
        }
        std::cout << "DBSCAN (Brute Force) found " << num_clusters_brute << " clusters and " << noise_points_brute << " noise points.\n";
        std::cout << "Saving Brute Force results to " << output_filename_brute << "...\n";
        saveResultsToFile(output_filename_brute, points, cluster_assignments_brute);
        std::cout << "Results saved successfully.\n";

        std::cout << "\n--- Performance Comparison ---\n";
        std::cout << "KD-Tree DBSCAN Time: " << dbscan_kd_duration.count() << " seconds\n";
        std::cout << "Brute Force DBSCAN Time: " << dbscan_brute_duration.count() << " seconds\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
