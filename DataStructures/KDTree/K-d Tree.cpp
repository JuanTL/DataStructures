#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <limits>
#include <queue>

struct Point {
    std::vector<double> coords; // k-dimensional point
    Point(const std::vector<double>& c = {}) : coords(c) {}
    size_t dimensions() const { return coords.size(); }
    /*
    bool operator<(const Point& other) const {
        return coords < other.coords; // relies on vector's lexicographical comparison
    }
    */
};

struct Node {
    Point point;
    std::unique_ptr<Node> left, right;
    Node(const Point& p) : point(p), left(nullptr), right(nullptr) {}
};

class KDTree {
private:
    std::unique_ptr<Node> root;
    size_t k; // Number of dimensions

    // Build k-d tree recursively
    /*
     /--------------------------------------------------------------------/
     * Recursively builds a balanced k-d tree from a list of points.
     *
     * std::vector<Point>& points -  points Vector of all points to be used in tree construction
     * int depth - depth Current depth in the tree (used to determine axis of comparison)
     * int start - Start index of the subarray of points to consider
     * int end - End index (exclusive) of the subarray of points to consider
     * return - pointer to the root of the subtree
     /--------------------------------------------------------------------/
     * The algorithm:
     * 1. Select an axis based on the current tree depth: axis = depth % k
     * 2. Sort the points in the current subarray [start, end) by the selected axis
     * 3. Choose the median point as the root of this subtree
     * 4. Recursively build left and right subtrees from the points left and right of the median
     /--------------------------------------------------------------------/
     */
    std::unique_ptr<Node> buildTree(std::vector<Point>& points, int depth, int start, int end) {
        if (start >= end) return nullptr; //empty subtree

        size_t axis = depth % k;
        std::sort(points.begin() + start, points.begin() + end, // Sort the current range of points based on the axis coordinate
            [axis](const Point& a, const Point& b) {
                return a.coords[axis] < b.coords[axis];
            });
        // Find the middle index to use as root for this subtree
        int mid = start + (end - start) / 2;
        auto node = std::make_unique<Node>(points[mid]);
		// Recursively build left and right subtrees
        node->left = buildTree(points, depth + 1, start, mid);
        node->right = buildTree(points, depth + 1, mid + 1, end);
        return node;
    }

    // Squared Euclidean distance
    double distance(const Point& a, const Point& b) const {
        double dist = 0.0;
        for (size_t i = 0; i < k; ++i) {
            double diff = a.coords[i] - b.coords[i];
            dist += diff * diff;
        }
        return dist;
    }

    // Nearest neighbor search
    /*
     /--------------------------------------------------------------------/
     * Recursively searches for the nearest neighbor to a given target point in a k-d tree.
     *
     * node       Current node in the k-d tree being visited
     * target     The target point we are searching nearest to
     * depth      Current depth in the tree (used to determine which axis to compare)
     * best       Reference to the current best (closest) point found so far
     * bestDist   Reference to the squared distance from target to the current best point
     /--------------------------------------------------------------------/
     * The function updates `best` and `bestDist` as it finds closer points.
     /--------------------------------------------------------------------/
     */
    void nearestNeighbor(const std::unique_ptr<Node>& node, const Point& target,
        int depth, Point& best, double& bestDist) const {
        if (!node) return; //Null?, thereÅfs nothing to check

        double dist = distance(node->point, target);
        // If this point is closer than the best found so far, update the best
        if (dist < bestDist) {
            best = node->point;
            bestDist = dist;
        }

        size_t axis = depth % k; // Determine axis for comparison based on current depth (Here should cycle through dimensions)
        double diff = target.coords[axis] - node->point.coords[axis];
        // Choose which subtree to explore first:(the one in the direction of the target along this axis)
        auto& nearSubtree = (diff < 0) ? node->left : node->right;
        auto& farSubtree = (diff < 0) ? node->right : node->left;

        nearestNeighbor(nearSubtree, target, depth + 1, best, bestDist);

		// If the distance to the splitting plane is less than the best distance found so far,
		// we need to check the other subtree as well A pruning step.	
        if (diff * diff < bestDist) {
            nearestNeighbor(farSubtree, target, depth + 1, best, bestDist);
        }
    }

    // Range search: collect points within radius of target
    /*
     /--------------------------------------------------------------------/
     * Recursively finds all points within a specified radius of a target point.
     *
     * node     Current node in the k-d tree being visited
     * target   The target point (center of the search region)
     * radius   The search radius (points within this distance from the target will be collected)
     * depth    Current depth in the tree (used to determine axis for splitting)
     * results  Output vector where matching points will be collected
     /--------------------------------------------------------------------/
     * Depth-first traversal of the k-d tree.
     /--------------------------------------------------------------------/
     */
    void rangeSearch(const std::unique_ptr<Node>& node, const Point& target,
        double radius, int depth, std::vector<Point>& results) const {
        if (!node) return;

        double dist = std::sqrt(distance(node->point, target)); // Eiclidean distance between target and current.
        // If the point is within the radius, include it in the results
        if (dist <= radius) {
            results.push_back(node->point);
        }

        
        size_t axis = depth % k;
        double diff = target.coords[axis] - node->point.coords[axis];

        // Explore both subtrees if the splitting plane is within radius
        /**
         * If the distance between the target and the current node along the splitting axis
         * is less than or equal to the radius, it means the search sphere crosses the splitting plane.
         * THEN: Explore BOTH subtrees, either side might contain valid points.
         *
         * Otherwise, explore the subtree on the same side as the target.
         */

        if (std::abs(diff) <= radius) {
            rangeSearch(node->left, target, radius, depth + 1, results);
            rangeSearch(node->right, target, radius, depth + 1, results);
        }
        else {
            auto& subtree = (diff < 0) ? node->left : node->right;
            rangeSearch(subtree, target, radius, depth + 1, results);
        }
    }

    // K-nearest neighbors using a max-heap
    using KNNPair = std::pair<double, Point>; // pair: (distance to target, point itself)
    
    /**
     /--------------------------------------------------------------------/
     * Recursively finds the k-nearest neighbors to a target point in a k-d tree.
     *
     * template Heap  A max-heap type (e.g., std::priority_queue<KNNPair>)
     * node   Current node in the k-d tree
     * target The point we want to find nearest neighbors to
     * k      The number of nearest neighbors to return
     * depth  Current depth in the tree (used to select splitting axis)
     * heap   Max-heap that stores the k best points found so far
     /--------------------------------------------------------------------/
     * The heap holds up to k elements, ordered by farthest-first (largest distance on top).
     * This to discard bad candidates during traversal.
     /--------------------------------------------------------------------/
     */
    
    template <typename Heap>
    void kNearestNeighbors(const std::unique_ptr<Node>& node,
        const Point& target,
        int k,
        int depth,
        Heap& heap) const {
        if (!node) return;

        double dist = distance(node->point, target);
        /**
         * Update the max-heap:
         * - If we haven't found k neighbors, add this one.
         * - If we have k neighbors, but this is closer than the current farthest, replace it.
         *
         * The top of the max-heap is the WORST (farthest) candidate.
         */
        if (heap.size() < static_cast<size_t>(k) || dist < heap.top().first) {
            heap.emplace(dist, node->point);
            if (heap.size() > static_cast<size_t>(k)) heap.pop();
        }

        size_t axis = depth % this->k;
        double diff = target.coords[axis] - node->point.coords[axis];
        // Decide which subtree is closer to the target
        auto& nearSubtree = (diff < 0) ? node->left : node->right;
        auto& farSubtree = (diff < 0) ? node->right : node->left;

        // Recursively search the near subtree
        kNearestNeighbors(nearSubtree, target, k, depth + 1, heap);
        /*
         * Only search the far subtree if there's a possibility it contains closer points.
         * TRUE if:
         *   - We havenÅft yet found k points, OR
         *   - The hypersphere intersects the splitting plane: (diff*diff < current-farthest-distance*current-farthest-distance)
         */
        if (heap.size() < static_cast<size_t>(k) || diff * diff < heap.top().first) {
            kNearestNeighbors(farSubtree, target, k, depth + 1, heap);
        }
    }

public:

    /*
     * Constructor: Builds a k-d tree from a given set of points.
     *
     * points      A vector of points to insert into the tree
     * dimensions  The number of dimensions for each point (k)
     */
    KDTree(const std::vector<Point>& points, size_t dimensions) : k(dimensions) {
        if (points.empty() || points[0].dimensions() != k) {
            throw std::invalid_argument("Invalid points or dimensions");
        }
        std::vector<Point> points_copy = points;
        root = buildTree(points_copy, 0, 0, points_copy.size());
    }
    /**
     * Finds the single nearest neighbor to a given target point.
     *
     * target  The point to find the closest neighbor to
     * return        The nearest point in the k-d tree
     *
     * Calls the recursive nearestNeighbor() function.
     */
    Point findNearest(const Point& target) const {
        if (target.dimensions() != k) {
            throw std::invalid_argument("Target point dimension mismatch");
        }
        Point best;
        double bestDist = std::numeric_limits<double>::max();
        nearestNeighbor(root, target, 0, best, bestDist);
        return best;
    }
    /**
     * Finds all points within a specified radius of a target point.
     *
     * target  The center point of the search region
     * radius  The maximum distance from the target to include points
     * return        A vector of all points within the given radius
     *
     * Calls the recursive rangeSearch() function.
     */
    std::vector<Point> findInRange(const Point& target, double radius) const {
        if (target.dimensions() != k) {
            throw std::invalid_argument("Target point dimension mismatch");
        }
        std::vector<Point> results;
        rangeSearch(root, target, radius, 0, results);
        return results;
    }

    /*
     /--------------------------------------------------------------------/
     * Finds the k-nearest neighbors to a target point.
     *
     * target The query point to which we want to find the nearest neighbors
     * k      The number of neighbors to find
     * return       A vector of the k closest points, sorted from closest to farthest
     /--------------------------------------------------------------------/
     * It prepares the heap, verifies dimensions, and formats the result.
     */

    std::vector<Point> findKNearest(const Point& target, int k) const {
		if (target.dimensions() != this->k) { // Check if target point has the same dimension as the tree
            throw std::invalid_argument("Target point dimension mismatch");
        }

        using KNNPair = std::pair<double, Point>;

        // Lambda comparator for max-heap: compare only by distance (first element of pair)
        //Points with larger distances have higher priority, to discard the farthest neighbor when a better one is found.
        auto cmp = [](const KNNPair& a, const KNNPair& b) {
            return a.first < b.first; // max-heap: larger distance has higher priority
            };

        // Declare priority queue with custom comparator
        std::priority_queue<KNNPair, std::vector<KNNPair>, decltype(cmp)> heap(cmp);

        // Perform recursive k-nearest neighbor search
        // Should accept any heap-like structure(hence the templated function)
        kNearestNeighbors(root, target, k, 0, heap);

        // Collect results from heap into a vector (in reverse order)
        std::vector<Point> results;
        results.reserve(k);
        while (!heap.empty()) {
            results.push_back(heap.top().second);
            heap.pop();
        }
        std::reverse(results.begin(), results.end()); // Closest first
        return results;
    }

};

// K-means clustering using k-d tree
class KMeans {
private:
    size_t k_clusters; // Number of clusters
    size_t dimensions;
    std::vector<Point> centroids; // Current centroid positions
    const KDTree& tree; // Reference to k-d tree

    // Initialize centroids randomly
    void initializeCentroids(const std::vector<Point>& points) {
        std::vector<Point> points_copy = points;
        std::random_shuffle(points_copy.begin(), points_copy.end());
        centroids.assign(points_copy.begin(), points_copy.begin() + k_clusters);
    }

    // Compute new centroid as mean of assigned points
    Point computeCentroid(const std::vector<Point>& points) const {
        if (points.empty()) return Point(std::vector<double>(dimensions, 0.0));
        std::vector<double> sum(dimensions, 0.0);
        for (const auto& p : points) {
            for (size_t i = 0; i < dimensions; ++i) {
                sum[i] += p.coords[i];
            }
        }
        for (size_t i = 0; i < dimensions; ++i) {
            sum[i] /= points.size();
        }
        return Point(sum);
    }

public:
    // Constructor: set cluster count, dimension, and initialize centroids
    KMeans(const std::vector<Point>& points, size_t k, const KDTree& kd_tree)
        : k_clusters(k), dimensions(points[0].dimensions()), tree(kd_tree) {
        initializeCentroids(points);
    }

    // Main k-means loop: assign points and update centroids

    /**
     * Runs the K-Means clustering algorithm on the given set of points.
     *
     * points          The dataset to cluster
     * max_iterations  Maximum number of iterations before stopping
	 /--------------------------------------------------------------------/
     * The algorithm proceeds in two main steps, repeated iteratively:
     *
     * 1. Assignment Step:
     *    Each point is assigned to the cluster whose centroid is nearest (using Euclidean distance). 
          Temporary assigns the data in sets.
     * 2. Update Step:
     *    Each centroid is updated to be the mean (average) of all points assigned to its cluster.
	 /--------------------------------------------------------------------/
     * The process repeats until the centroids stop changing (convergence) or until the
        maximum number of iterations is reached.
     */

    void run(const std::vector<Point>& points, int max_iterations = 100) {
        for (int iter = 0; iter < max_iterations; ++iter) {
            // Prepare empty clusters
            std::vector<std::vector<Point>> clusters(k_clusters);
            // --- Assignment Step ---
            // For each point in the dataset:
            for (const auto& point : points) {
                Point nearest = tree.findNearest(point);
                // Find which centroid is closest
                int cluster_id = 0;
                double min_dist = std::numeric_limits<double>::max();
                // Compare the point to all centroids
                for (size_t i = 0; i < k_clusters; ++i) {
                    double dist = 0.0;
                    for (size_t d = 0; d < dimensions; ++d) {
                        double diff = point.coords[d] - centroids[i].coords[d];
                        dist += diff * diff;
                    }
                    // Keep track of the closest centroid
                    if (dist < min_dist) {
                        min_dist = dist;
                        cluster_id = i;
                    }
                }
                clusters[cluster_id].push_back(point);
            }

            // --- Update Step ---
            // Recalculate centroids as the mean of their assigned cluster
            std::vector<Point> new_centroids;
            bool converged = true;
            for (size_t i = 0; i < k_clusters; ++i) {
                Point new_centroid = computeCentroid(clusters[i]);
                // Check if this centroid has changed
                if (new_centroid.coords != centroids[i].coords) {
                    converged = false;
                }
                new_centroids.push_back(new_centroid);
            }
            centroids = new_centroids; // Replace old centroids with the new ones

            if (converged) break; // Stop early if centroids did not change (convergence reached)
        }
    }

    const std::vector<Point>& getCentroids() const { return centroids; }
};


int main() {
    // 3D points
    std::vector<Point> points = {
        Point({2, 3, 1}),
        Point({5, 4, 2}),
        Point({9, 6, 3}),
        Point({4, 7, 4}),
        Point({8, 1, 5}),
        Point({7, 2, 6}),
        Point({1, 1, 1}),
        Point({10, 10, 10})
    };
    size_t dimensions = 3;
    // Build a k-d tree from the points
    KDTree tree(points, dimensions);

    // Range search
    Point target({ 5, 5, 5 });
    double radius = 5.0;
	// Find all points within a given radius of the target point
    std::vector<Point> in_range = tree.findInRange(target, radius);
    
    std::cout << "Points within radius " << radius << " of ";
    for (double c : target.coords) std::cout << c << " ";
    std::cout << ":\n";
    for (const auto& p : in_range) {
        for (double c : p.coords) std::cout << c << " ";
        std::cout << "\n";
    }

    // K-nearest neighbors
    int k_nn = 3;
    std::vector<Point> knn = tree.findKNearest(target, k_nn);
    std::cout << "\n" << k_nn << "-nearest neighbors to ";
    for (double c : target.coords) std::cout << c << " ";
    std::cout << ":\n";
    for (const auto& p : knn) {
        for (double c : p.coords) std::cout << c << " ";
        std::cout << "\n";
    }

    // K-means clustering
    size_t k_clusters = 2;
    KMeans kmeans(points, k_clusters, tree);
    kmeans.run(points);
    std::cout << "\nK-means centroids:\n";
    for (const auto& c : kmeans.getCentroids()) {
        for (double coord : c.coords) std::cout << coord << " ";
        std::cout << "\n";
    }

    return 0;
}