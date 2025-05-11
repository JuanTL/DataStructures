#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

// Data point struct
struct Point {
    double lat, lon;
	std::vector<double> attributes; //All attributes are numerical. Good.
};

// Tree node struct
struct Node {
    std::vector<Point> points;           // Points in this node (leaf-level only for now)
    std::vector<Node*> children;         // Sub-clusters
    double min_lat, max_lat, min_lon, max_lon;  // Spatial bounding box
    std::vector<double> min_attrs, max_attrs;   // Attribute bounding ranges
    int point_count;                     // Number of points in subtree

    Node() : min_lat(std::numeric_limits<double>::max()), max_lat(std::numeric_limits<double>::min()),
        min_lon(std::numeric_limits<double>::max()), max_lon(std::numeric_limits<double>::min()),
        point_count(0) {
    }
};

// Distance functions

// Radius of the Earth in kilometers
const double R = 6371.0;
const double PI = 3.14159265358979323846;
// Convert degrees to radians
double degToRad(double deg) {
    return deg * (PI / 180.0);
}
// Haversine distance between two points
double haversineDistance(double lat1, double lon1, double lat2, double lon2) {
    lat1 = degToRad(lat1);
    lon1 = degToRad(lon1);
    lat2 = degToRad(lat2);
    lon2 = degToRad(lon2);
    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;
    // Haversine formula
    double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
        std::cos(lat1) * std::cos(lat2) *
        std::sin(dlon / 2) * std::sin(dlon / 2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
    // Distance in kilometers
    return R * c;
}

// Spatial distance
double spatial_distance(const Point& p1, const Point& p2) {
	return std::hypot(p1.lat - p2.lat, p1.lon - p2.lon); //Hypotenuse function, still Euclidean distance
}

// Attribute distance
double attribute_distance(const Point& p1, const Point& p2) {
    double sum = 0.0;
    for (size_t i = 0; i < p1.attributes.size(); ++i) {
        sum += (p1.attributes[i] - p2.attributes[i]) * (p1.attributes[i] - p2.attributes[i]);
    }
	return std::sqrt(sum); //Euclidean distance
}

// Combined distance
double combined_distance(const Point& p1, const Point& p2, double alpha, double beta) {
    return alpha * spatial_distance(p1, p2) + beta * attribute_distance(p1, p2);
}

// Build the HGA-Tree
Node* build_tree(std::vector<Point> points, double alpha, double beta) {
    Node* node = new Node();
    node->points = points;
    node->point_count = points.size();

    // Bounding box and attribute ranges
    for (const auto& p : points) {
        node->min_lat = std::min(node->min_lat, p.lat);
        node->max_lat = std::max(node->max_lat, p.lat);
        node->min_lon = std::min(node->min_lon, p.lon);
        node->max_lon = std::max(node->max_lon, p.lon);
        if (node->min_attrs.empty()) {
            node->min_attrs = p.attributes;
            node->max_attrs = p.attributes;
        }
        else {
            for (size_t i = 0; i < p.attributes.size(); ++i) {
                node->min_attrs[i] = std::min(node->min_attrs[i], p.attributes[i]);
                node->max_attrs[i] = std::max(node->max_attrs[i], p.attributes[i]);
            }
        }
    }
    // Split into two children if more than one point
    if (points.size() > 1) {
        std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
            return a.lat < b.lat;
            });
        auto mid = points.begin() + points.size() / 2;
        std::vector<Point> left(points.begin(), mid);
        std::vector<Point> right(mid, points.end());
        node->children.push_back(build_tree(left, alpha, beta));
        node->children.push_back(build_tree(right, alpha, beta));
    }
    return node;
}

// Range query
void range_query(const Node* node, double lat_min, double lat_max, double lon_min, double lon_max,
    std::vector<Point>& result) {
    if (!node) return;
    if (node->max_lat < lat_min || node->min_lat > lat_max || node->max_lon < lon_min || node->min_lon > lon_max) {
        return; // No overlap
    }
    if (node->min_lat >= lat_min && node->max_lat <= lat_max && node->min_lon >= lon_min && node->max_lon <= lon_max) {
        result.insert(result.end(), node->points.begin(), node->points.end()); // Fully inside
    }
    else {
        for (const auto& p : node->points) {
            if (p.lat >= lat_min && p.lat <= lat_max && p.lon >= lon_min && p.lon <= lon_max) {
                result.push_back(p);
            }
        }
        for (auto child : node->children) {
            range_query(child, lat_min, lat_max, lon_min, lon_max, result);
        }
    }
}

// Data comparison: Find points similar to A and B
void find_similar_points(const Node* node, const Point& a, const Point& b, double attr_threshold,
    std::vector<Point>& result) {
    if (!node) return;
    bool possible = true;
    for (size_t i = 0; i < a.attributes.size(); ++i) {
        double min_range = std::min(a.attributes[i], b.attributes[i]) - attr_threshold;
        double max_range = std::max(a.attributes[i], b.attributes[i]) + attr_threshold;
        if (node->max_attrs[i] < min_range || node->min_attrs[i] > max_range) {
            possible = false;
            break;
        }
    }
    if (!possible) return;
    for (const auto& p : node->points) {
        if (attribute_distance(p, a) <= attr_threshold && attribute_distance(p, b) <= attr_threshold) {
            result.push_back(p);
        }
    }
    for (auto child : node->children) {
        find_similar_points(child, a, b, attr_threshold, result);
    }
}

// Region retrieval: Count points and collect clusters
void region_retrieval(const Node* node, double lat_min, double lat_max, double lon_min, double lon_max,
    int& count, std::vector<Node*>& clusters) {
    if (!node) return;
    if (node->max_lat < lat_min || node->min_lat > lat_max || node->max_lon < lon_min || node->min_lon > lon_max) {
        return; // No overlap
    }
    if (node->min_lat >= lat_min && node->max_lat <= lat_max && node->min_lon >= lon_min && node->max_lon <= lon_max) {
        count += node->point_count;
        if (!node->children.empty()) clusters.push_back(const_cast<Node*>(node));
    }
    else {
        for (const auto& p : node->points) {
            if (p.lat >= lat_min && p.lat <= lat_max && p.lon >= lon_min && p.lon <= lon_max) {
                count++;
            }
        }
        for (auto child : node->children) {
            region_retrieval(child, lat_min, lat_max, lon_min, lon_max, count, clusters);
        }
    }
}

int main() {
    std::vector<Point> points = {
		{40.7128, -74.0060, {1.0, 2.0}},  // New York coordinates
		{37.7749, -122.4194, {2.0, 3.0}}, // San Francisco coordinates
		{41.8781, -87.6298, {5.0, 6.0}},   // Chicago coordinates
    };
    /*
    Some coordinates to work with later:
        {34.0522, -118.2437, {3.0, 4.0}}, // Los Angeles coordinates
		{51.5074, -0.1278, {4.0, 5.0}},   // London coordinates
		{48.8566, 2.3522, {5.0, 6.0}},    // Paris coordinates
		{35.6895, 139.6917, {6.0, 7.0}},  // Tokyo coordinates
		{55.7558, 37.6173, {7.0, 8.0}},   // Moscow coordinates
    */

    double alpha = 1.0, beta = 1.0;
    Node* root = build_tree(points, alpha, beta);

    // Range query
    std::vector<Point> range_result;
    range_query(root, 30.0, 45.0, -120.0, -70.0, range_result);
    std::sort(range_result.begin(), range_result.end(), [](const Point& a, const Point& b) {
        return a.lat < b.lat;
        });
    std::cout << "Range Query Results:\n";
    for (const auto& p : range_result) {
        std::cout << "Lat: " << p.lat << ", Lon: " << p.lon << "\n";
    }

    // Data comparison
    Point a = { 40.7128, -74.0060, {1.0, 2.0} };
    Point b = { 41.8781, -87.6298, {5.0, 6.0} };
    std::vector<Point> similar_points;
    find_similar_points(root, a, b, 5.0, similar_points);
    std::cout << "\nSimilar Points to A and B:\n";
    for (const auto& p : similar_points) {
        std::cout << "Lat: " << p.lat << ", Lon: " << p.lon << "\n";
    }

    // Region retrieval
    int count = 0;
    std::vector<Node*> clusters;
    region_retrieval(root, 30.0, 45.0, -120.0, -70.0, count, clusters);
    std::cout << "\nRegion Retrieval:\n";
    std::cout << "Number of points: " << count << "\n";
    std::cout << "Number of clusters: " << clusters.size() << "\n";

    return 0;
}
