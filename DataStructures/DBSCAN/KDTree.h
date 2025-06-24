#pragma once
#include <vector>
#include <memory>
#include <algorithm>
#include <limits>    
#include <cmath>     
#include <stdexcept> 
#include <queue>     
#include <functional> 
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

// --- Point struct ---
struct Point {
    std::vector<double> coords;

    // Constructor: Initializes a point with given coordinates
    Point(const std::vector<double>& c = {}) : coords(c) {}

    // Get Dimensions of the point
    size_t dimensions() const { return coords.size(); }
};

// --- Node(K-d Tree) ---
struct Node {
    Point point;                          
    std::unique_ptr<Node> left, right;   
	int original_index; // Original index of the point - added for range search with indices(needed for DBSCAN)
    Node(const Point& p) : point(p), left(nullptr), right(nullptr) {}
    Node(const Point& p, int idx) : point(p), original_index(idx), left(nullptr), right(nullptr) {}

};

// --- KDTree ---
class KDTree {
private:
    std::unique_ptr<Node> root; 
    size_t k;                   // dimensions for the points

    std::unique_ptr<Node> buildTree(std::vector<Point>& points, int depth, int start, int end) {
        if (start >= end) {
            return nullptr; // Empty subtree
        }

        size_t axis = depth % k; // Determine the splitting axis for the current depth
        int mid_idx = start + (end - start) / 2;
        std::nth_element(points.begin() + start, points.begin() + mid_idx, points.begin() + end,
            [axis](const Point& a, const Point& b) {
                return a.coords[axis] < b.coords[axis];
            });
        auto node = std::make_unique<Node>(points[mid_idx]);

        // Recursively build the left and right subtrees
        node->left = buildTree(points, depth + 1, start, mid_idx);
        node->right = buildTree(points, depth + 1, mid_idx + 1, end);

        return node;
    }

    std::unique_ptr<Node> buildTree(std::vector<std::pair<Point, int>>& indexed_points, int depth, int start, int end) {
        if (start >= end) {
            return nullptr;
        }
        size_t axis = depth % k;
        int mid_idx = start + (end - start) / 2;
        std::nth_element(indexed_points.begin() + start, indexed_points.begin() + mid_idx, indexed_points.begin() + end,
            [axis](const std::pair<Point, int>& a, const std::pair<Point, int>& b) {
                return a.first.coords[axis] < b.first.coords[axis];
            });

        // Create node with the point and its original index
        auto node = std::make_unique<Node>(indexed_points[mid_idx].first, indexed_points[mid_idx].second);

        node->left = buildTree(indexed_points, depth + 1, start, mid_idx);
        node->right = buildTree(indexed_points, depth + 1, mid_idx + 1, end);

        return node;
    }
    double squaredDistance(const Point& a, const Point& b) const {
        double dist = 0.0;
        for (size_t i = 0; i < k; ++i) {
            double diff = a.coords[i] - b.coords[i];
            dist += diff * diff;
        }
        return dist;
    }

    // Nearest neighbor search
    void nearestNeighbor(const std::unique_ptr<Node>& node, const Point& target,
        int depth, Point& best, double& bestSquaredDist) const {
        if (!node) {
            return;
        }
        double currentSquaredDist = squaredDistance(node->point, target);
        if (currentSquaredDist < bestSquaredDist) {
            best = node->point;
            bestSquaredDist = currentSquaredDist;
        }

        size_t axis = depth % k; // Determine the splitting axis
        double diff = target.coords[axis] - node->point.coords[axis]; // Difference along the current axis

        // Determine which subtree to explore first
        auto& nearSubtree = (diff < 0) ? node->left : node->right;
        auto& farSubtree = (diff < 0) ? node->right : node->left;

        // Recursively search the near subtree
        nearestNeighbor(nearSubtree, target, depth + 1, best, bestSquaredDist);

        if (diff * diff < bestSquaredDist) {
            nearestNeighbor(farSubtree, target, depth + 1, best, bestSquaredDist);
        }
    }

    // Range search - finds/collects points within a radius
    void rangeSearch(const std::unique_ptr<Node>& node, const Point& target,
        double squaredRadius, int depth, std::vector<Point>& results) const {
        if (!node) {
            return;
        }
        double currentSquaredDist = squaredDistance(node->point, target);

        // Point within the radius? add it in the results
        if (currentSquaredDist <= squaredRadius) {
            results.push_back(node->point);
        }

        size_t axis = depth % k;
        double diff = target.coords[axis] - node->point.coords[axis]; // Difference along the current axis

        if (diff * diff <= squaredRadius) {
            rangeSearch(node->left, target, squaredRadius, depth + 1, results);
            rangeSearch(node->right, target, squaredRadius, depth + 1, results);
        }
        else {
            // Only explore the subtree that contains the target along this axis
            auto& subtree = (diff < 0) ? node->left : node->right;
            rangeSearch(subtree, target, squaredRadius, depth + 1, results);
        }
    }
    // Modified rangeSearch to return indices
    void rangeSearchIndices(const std::unique_ptr<Node>& node, const Point& target,
        double squaredRadius, int depth, std::vector<int>& results) const {
        if (!node) return;

        double currentSquaredDist = squaredDistance(node->point, target);
        if (currentSquaredDist <= squaredRadius) {
            results.push_back(node->original_index); // Store index instead of point
        }

        size_t axis = depth % k;
        double diff = target.coords[axis] - node->point.coords[axis];

        if (diff * diff <= squaredRadius) { // Use squared diff
            rangeSearchIndices(node->left, target, squaredRadius, depth + 1, results);
            rangeSearchIndices(node->right, target, squaredRadius, depth + 1, results);
        }
        else {
            auto& subtree = (diff < 0) ? node->left : node->right;
            rangeSearchIndices(subtree, target, squaredRadius, depth + 1, results);
        }
    }
    template <typename Heap>
    void kNearestNeighbors(const std::unique_ptr<Node>& node,
        const Point& target,
        int k,
        int depth,
        Heap& heap) const {
        if (!node) return;

        double dist = squaredDistance(node->point, target);
        /*
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
         *   - We haven�ft yet found k points, OR
         *   - The hypersphere intersects the splitting plane: (diff*diff < current-farthest-distance*current-farthest-distance)
         */
        if (heap.size() < static_cast<size_t>(k) || diff * diff < heap.top().first) {
            kNearestNeighbors(farSubtree, target, k, depth + 1, heap);
        }
    }

public:
    KDTree(const std::vector<Point>& points, size_t dimensions) : k(dimensions) {
        if (points.empty()) {
            root = nullptr;
            return;
        }
        if (points[0].dimensions() != k) {
            throw std::invalid_argument("Dimension mismatch: points must have 'k' dimensions.");
        }

        std::vector<std::pair<Point, int>> indexed_points;
        indexed_points.reserve(points.size());
        for (int i = 0; i < points.size(); ++i) {
            indexed_points.emplace_back(points[i], i);
        }

        root = buildTree(indexed_points, 0, 0, indexed_points.size());
    }
    // Finds the single nearest neighbor to the target point
    Point findNearest(const Point& target) const {
        if (!root) {
            throw std::runtime_error("KDTree is empty. Cannot find nearest neighbor.");
        }
        if (target.dimensions() != k) {
            throw std::invalid_argument("Target point dimension mismatch.");
        }
        Point best_point = root->point; // Initialize with a point from the tree
        double bestSquaredDist = std::numeric_limits<double>::max();
        nearestNeighbor(root, target, 0, best_point, bestSquaredDist);
        return best_point;
    }

    // Finds all points within a radius of the target point
    std::vector<Point> findInRange(const Point& target, double radius) const {
        if (!root) {
            return {};
        }
        if (target.dimensions() != k) {
            throw std::invalid_argument("Target point dimension mismatch.");
        }
        std::vector<Point> results;
        double squaredRadius = radius * radius;
        rangeSearch(root, target, squaredRadius, 0, results);
        return results;
    }
    std::vector<int> findInRangeIndices(const Point& target, double radius) const {
        if (target.dimensions() != k) {
            throw std::invalid_argument("Target point dimension mismatch");
        }
        std::vector<int> results;
        double squaredRadius = radius * radius;
        rangeSearchIndices(root, target, squaredRadius, 0, results);
        return results;
    }
    // Finds the k-nearest neighbors to the target point
    std::vector<Point> findKNearest(const Point& target, int k) const {
        if (target.dimensions() != this->k) {
            throw std::invalid_argument("Target point dimension mismatch");
        }

        using KNNPair = std::pair<double, Point>;
        auto cmp = [](const KNNPair& a, const KNNPair& b) {
            return a.first < b.first;
            };
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