#include<cmath>
#include<vector>
#include<string>
#include<sstream>
#include<fstream>
using namespace std;

int main() {
  ifstream scanFile; scanFile.open("metadata.csv"); string line; getline(scanFile, line);
  double binSize=1.5; int negativeBinNumber=2; double maxEnergy=249.89;
  vector<unsigned int> distribution(log(maxEnergy)/log(binSize)+negativeBinNumber+1);
  while (getline(scanFile, line)) {
    string token; double energy;
    stringstream ss; ss.str(line);
    getline(ss, token, ','); getline(ss, token, ','); getline(ss, token, ','); getline(ss, token, ','); getline(ss, token, ','); getline(ss, token, ','); getline(ss, token, ',');
    getline(ss, token, ','); energy=stod(token);
    ++distribution[log(energy)/log(binSize)+negativeBinNumber];
  } scanFile.close();
  
  ofstream printFileDistribution("Distribution.txt");
  for (int exponent=-negativeBinNumber; exponent<(int)distribution.size()-negativeBinNumber; ++exponent) {
    printFileDistribution << pow(binSize, exponent+0.5) << "\t" << distribution[exponent+negativeBinNumber] << "\n";
  }
  printFileDistribution.close();
}
