#include <string>
#include <vector>
#include <iostream>
#include "RTree.h"

using namespace std;

#include <fstream>
#include<random>

void print(const int& command, vector<vector<vector<pair<int, int>>>>& objects_n, string& output, int x= 0) {
	output.resize(0);
	output = to_string(command);
	
	// Open file for rectangles
	ofstream rectFile("rectangles"+ to_string(x) + ".txt", ios::out | ios::trunc);
	if (!rectFile.is_open()) {
		cerr << "Error opening rectangles.txt" << endl;
		return;
	}

	switch (command) {
	case 0: // OBJECTS
		output += "|" + to_string(objects_n[0].size());
		for (auto& polygon : objects_n[0]) {
			output += "|" + to_string(polygon.size());

			// Calculate bounding rectangle for the polygon
			if (!polygon.empty()) {
				int minX = polygon[0].first, minY = polygon[0].second;
				int maxX = polygon[0].first, maxY = polygon[0].second;
				for (auto& point : polygon) {
					minX = min(minX, point.first);
					minY = min(minY, point.second);
					maxX = max(maxX, point.first);
					maxY = max(maxY, point.second);
					output += "|" + to_string(point.first) + "|" + to_string(point.second);
				}
				// Write rectangle to file (minX, minY, maxX, maxY)
				rectFile << minX << " " << minY << " " << maxX << " " << maxY << "\n";
			}
		}
		break;

	case 1: // MBRS
		output += "|" + to_string(objects_n.size());
		for (auto& objects : objects_n) {
			output += "|" + to_string(objects.size());
			for (auto& polygon : objects) {
				output += "|" + to_string(polygon.size());
				if (!polygon.empty()) {
					int minX = polygon[0].first, minY = polygon[0].second;
					int maxX = polygon[0].first, maxY = polygon[0].second;
					for (auto& point : polygon) {
						minX = min(minX, point.first);
						minY = min(minY, point.second);
						maxX = max(maxX, point.first);
						maxY = max(maxY, point.second);
						output += "|" + to_string(point.first) + "|" + to_string(point.second);
					}
					// Write rectangle to file (minX, minY, maxX, maxY)
					rectFile << minX << " " << minY << " " << maxX << " " << maxY << "\n";
				}
			}
		}
		break;

	default:
		output += "|0";
		break;
	}

	output += "|END";
	rectFile.close();
	cout << output << endl;
}
void print_pair(vector<pair<int, int>> output, int x = 0) {
	// Open file for points
	ofstream pointFile("points" + to_string(x) + ".txt", ios::out | ios::trunc);
	if (!pointFile.is_open()) {
		cerr << "Error opening points.txt" << endl;
		return;
	}

	for (auto& x : output) {
		pointFile << x.first << " " << x.second << "\n";
		cout << x.first << " ; " << x.second << " - ";
	}
	pointFile.close();
}
 

#include<chrono>

vector<pair<int, int>> generateRandomPoints(mt19937& gen, int numPoints, int minCoord, int maxCoord) {
	uniform_int_distribution<int> dist(minCoord, maxCoord);
	vector<pair<int, int>> points;

	for (int i = 0; i < numPoints; ++i) {
		int x = dist(gen);
		int y = dist(gen);
		points.emplace_back(x, y);
	}
	return points;
}

// Same as the already existing MBR function...
Rect computeBoundingRect(const vector<pair<int, int>>& points) {
	assert(!points.empty());
	int minX = points[0].first, minY = points[0].second;
	int maxX = points[0].first, maxY = points[0].second;

	for (const auto& point : points) {
		minX = min(minX, point.first);
		minY = min(minY, point.second);
		maxX = max(maxX, point.first);
		maxY = max(maxY, point.second);
	}

	return Rect(minX, minY, maxX, maxY);
}

// Function to compare two point sets (order-independent)
bool comparePointSets(const vector<pair<int, int>>& set1, const vector<pair<int, int>>& set2) {
	if (set1.size() != set2.size()) return false;
	vector<pair<int, int>> sorted1 = set1, sorted2 = set2;
	sort(sorted1.begin(), sorted1.end());
	sort(sorted2.begin(), sorted2.end());
	return sorted1 == sorted2;
}

void runTest(int numEntries, int maxNodes, mt19937& gen, int testNum = 0) {
	// Set up random number generator
	uniform_int_distribution<int> sizeDist(2, 3); // Randomly choose 2 or 3 points per set

	// Create R-tree (MAXNODES is set via preprocessor or global define)
	#define MAXNODES maxNodes
	#define MINNODES ((maxNodes + 1) / 2) // Ensure MINNODES is roughly half of MAXNODES
	RTree rtree;

	// Store inserted entries for search verification
	vector<pair<Rect, vector<pair<int, int>>>> insertedEntries;

	// Generate and insert entries
	vector<vector<pair<int, int>>> allPointSets;
	for (int i = 0; i < numEntries; ++i) {
		// Generate random point set (2 or 3 points)
		int numPoints = sizeDist(gen);
		vector<pair<int, int>> points = generateRandomPoints(gen, numPoints, 0, 100); // Coordinates in [0, 100]
		// Compute bounding rectangle
		Rect rect = computeBoundingRect(points);

		// Insert into R-tree
		int min[2] = { rect.m_min[0], rect.m_min[1] };
		int max[2] = { rect.m_max[0], rect.m_max[1] };
		rtree.Insert(min, max, points);

		// Store for search verification
		insertedEntries.emplace_back(rect, points);

		// Store point set for printing
		allPointSets.push_back(points);
	}
	// Search
	bool allSearchesSuccessful = true;
	for (size_t i = 0; i < insertedEntries.size(); ++i) {
		const Rect& rect = insertedEntries[i].first;
		const vector<pair<int, int>>& expectedPoints = insertedEntries[i].second;

		// Prepare and search
		int min[2] = { rect.m_min[0], rect.m_min[1] };
		int max[2] = { rect.m_max[0], rect.m_max[1] };
		vector<vector<pair<int, int>>> results = rtree.Search(min, max);

		// Check if expected point set is in results
		bool found = false;
		for (const auto& result : results) {
			if (comparePointSets(result, expectedPoints)) {
				found = true;
				break;
			}
		}
		if (!found) {
			allSearchesSuccessful = false;
			cout << "Search failed for entry " << i << ": Expected point set not found" << endl;
			cout << "Query rectangle: (" << min[0] << ", " << min[1] << ", " << max[0] << ", " << max[1] << ")" << endl;
			cout << "Expected points: ";
			for (const auto& p : expectedPoints) {
				cout << "(" << p.first << ", " << p.second << ") ";
			}
			cout << endl;
		}
	}
	vector<vector<vector<pair<int, int>>>> objects_n = { allPointSets };
	string output;

	// Print to files (rectangles.txt and points.txt)
	print(0, objects_n, output, testNum);

	// Print one set of points to points.txt for verification
	if (!allPointSets.empty()) {
		print_pair(allPointSets[0],testNum);
	}
	// Output
	cout << "Test with " << numEntries << " entries, MAXNODES=" << maxNodes << " completed." << endl;
	cout << "Search verification: " << (allSearchesSuccessful ? "All searches successful" : "Some searches failed") << endl;
	cout << "Rectangles written to rectangles.txt, first point set to points.txt" << endl;

	// Undefine to allow different MAXNODES in next test
	#undef MAXNODES
	#undef MINNODES
}



int main(int argc, char* argv[])
{
	//AUTOMATIC TESTS
	//*
	vector<int> numEntriesList = { 5, 12 }; // Test with 5 and 12 entries
	vector<int> maxNodesList = { 3, 5 };    // Test with MAXNODES = 3 and 5
	
	int x = 0;
	random_device rd;
	for (int numEntries : numEntriesList) {
		for (int maxNodes : maxNodesList) {
			cout << "\nRunning test: Entries=" << numEntries << ", MAXNODES=" << maxNodes << endl;
			x++;
			mt19937 gen(rd());
			runTest(numEntries, maxNodes, gen, x);
		}
	}

	return 0;
	//*/
	/*
	vector<vector<pair<int, int>>> vpoints;
	vector<float> coord;

	coord.push_back(20);  coord.push_back(59); coord.push_back(20);  coord.push_back(43);
	coord.push_back(50);  coord.push_back(58); coord.push_back(48);  coord.push_back(67);
	coord.push_back(105); coord.push_back(68); coord.push_back(74);  coord.push_back(64);
	coord.push_back(83);  coord.push_back(40); coord.push_back(104); coord.push_back(54);

	coord.push_back(20);  coord.push_back(59); coord.push_back(20);  coord.push_back(43);
	coord.push_back(48);  coord.push_back(67); coord.push_back(105); coord.push_back(68);

	vector<pair<int, int>> points;
	// Store inserted entries for search verification
	vector<pair<Rect, vector<pair<int, int>>>> insertedEntries;

	for (int i = 0; i < coord.size(); i += 2) {
		pair<int, int> A;
		A.first = coord[i];
		A.second = coord[i + 1];
		points.push_back(A);
	}

	for (unsigned int i = 0; i < points.size(); i += 2) {
		vector<pair<int, int>> sub1;
		// Fix out of bounds bug
		if (i + 1 < points.size()) {
			sub1.push_back(points[i]);
			sub1.push_back(points[i + 1]);
			vpoints.push_back(sub1);
			cout << i << " >> " << sub1[0].first << " " << sub1[0].second << endl;
		}
	}

	RTree rtree;
	string output;
	vector<vector<pair<int, int>>> objects;
	vector<vector<vector<pair<int, int>>>> objects_n;

	for (auto& x : vpoints)
	{
		cout << "inserting " << x.size() << ": ";
		print_pair(x);
		Rect rect = rtree.MBR(x);
		rtree.Insert(rect.m_min, rect.m_max, x);

		cout << endl;
	}
	rtree.getMBRs(objects_n);

	// Search with fixed rectangle (0,0) to (25,25)
	int min[2] = { 0, 0 };
	int max[2] = { 25, 25 };

	cout << "\nSearching for rectangle: (0,0) to (25,25)" << endl;
	vector<vector<pair<int, int>>> results = rtree.Search(min, max);
	cout << "\nResults from search:\n";
	for (const auto& r : results) {
		cout << "Result MBR: [("
			<< r[0].first << ", " << r[0].second << ") -> ("
			<< r[1].first << ", " << r[1].second << ")]\n";
	}

	//print(1, objects_n, output);
	return 0;
	*/ 
}