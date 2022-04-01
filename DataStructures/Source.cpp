#include<iostream>
#include<random>
#include<vector>
#include<math.h>
#include<string>
#include<sstream>
#include<fstream>
#include<utility>
using namespace std;
//Dimensions 1 to 10, then 20, 30, 40...100
//Data set size = 20k
//Data range = 1 to 100

//uniform_int_distribution<> distrib(1, 100);

double EuclideanDistance(vector<int> dataPointA, vector<int> dataPointB)
{
	int sum = 0, v = 0;
	for(int i= 0; i < dataPointA.size(); i++)
	{
		v = dataPointA[i] - dataPointB[i];
		sum += v * v;
	}
	return sqrt(sum);
}
template<int dimension,int size>
struct RandomCoincidence
{
	//Parameters
	const string savePath = "NumPair/";
	//fill numbers
	vector<vector<int>> Set;			//Set of Points
	vector<int> Point;					//A single Point of dimension k
	//vector<double> DistSet;				//Distances
	vector<pair<double, int>> DCount;   //Distance and times it was repeated
	RandomCoincidence()
	{
		GenerateRandom();
		FindDistances();
		//CountCoincidences();
		Write();
	}
	void GenerateRandom()
	{
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<> distrib(1, 100);
		for (int i = 0; i < size; i++)
		{
			vector<int> Point;
			for (int j = 0; j < dimension; j++)
			{
				Point.push_back(distrib(gen));
			}
			Set.push_back(Point);
		}
	}
	void FindDistances()
	{
		double distance = 0;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				bool inArray = false;
				int pos = 0;
				distance = EuclideanDistance(Set[i], Set[j]);
				for(int k = 0; k<DCount.size();k++)
				{
					if (DCount[k].first == distance)
					{
						inArray = true;
						pos = k;
						break;
					}
				}
				if(inArray == true)
				{
					DCount[pos].second = DCount[pos].second+1;
				}
				else
				{
					pair<double, int> tmp;
					tmp.first = distance;
					tmp.second = 1;
					DCount.push_back(tmp);
				}
				//DistSet.push_back(distance);
			}
		}
		return;
	}
	/*
	void CountCoincidences()
	{
		for (int i = 0; i < DistSet.size(); i++)
		{
			int coincidences = 0;
			double actual = DistSet[i];
			for (int j = i; j < DistSet.size(); j++)
			{
				if (DistSet[j] == actual)
				{
					coincidences++;
					DistSet.erase(DistSet.begin()+ j);
				}
					
			}
			pair<double, int> result;
			result.first = actual;
			result.second = coincidences;
			DCount.push_back(result);
		}
	}
	*/
	void Write()
	{
		string detail = savePath + to_string(dimension) + ".txt";
		string content = "";
		for (int i = 0; i < DCount.size(); i++)
		{
			content += to_string(DCount[i].first) + " " + to_string(DCount[i].second) + "\n";
		}
		ofstream output(detail);
		output << content;
		output.close();
	}

};

/*
void DimensionTest(int dimension)
{
	//fill numbers
	vector<vector<int>> Set;
	vector<int> Point;
	//Random
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> distrib(1, 100);
	//Fill Random
	for (int i = 0; i < 20000; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			Point.push_back(distrib(gen));
		}
		Set.push_back(Point);
		Point.clear();
	}
	//Calculate Distances
	double distance = 0;
	vector<double> DistSet;
	for (int i = 0; i < 20000/2; i++)
	{
		for (int j = 0; j < 20000; j++)
		{
			distance = EuclideanDistance(Set[i], Set[j]);
			DistSet.push_back(distance);
		}
	}
	//Cleaning
	vector<pair<double, int>> DCount;
	for (int i = 0; i < DistSet.size(); i++)
	{
		int coincidences = 0;
		double actual = DistSet[i];
		for (int j = i; j < DistSet.size(); j++)
		{
			if (DistSet[j] == actual)
				coincidences++;
		}
		pair<double, int> result;
		result.first = actual;
		result.second = coincidences;
		DCount.push_back(result);
	}
	//Write in Disk
	string content;
}
*/
int main() 
{
	//RandomCoincidence<1, 2> A;
	
	RandomCoincidence<1, 20000> A;
	RandomCoincidence<2, 20000> B;
	RandomCoincidence<3, 20000> C;
	RandomCoincidence<4, 20000> D;
	RandomCoincidence<5, 20000> E;
	RandomCoincidence<6, 20000> F;
	RandomCoincidence<7, 20000> G;
	RandomCoincidence<8, 20000> H;
	RandomCoincidence<9, 20000> I;
	RandomCoincidence<10, 20000> J;
	
}