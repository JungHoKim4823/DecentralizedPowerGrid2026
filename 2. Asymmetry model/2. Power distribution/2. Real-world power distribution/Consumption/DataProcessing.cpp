#include<random>
#include<vector>
#include<set>
#include<map>
#include<string>
#include<sstream>
#include<fstream>
using namespace std;

int main() {
  ifstream scanFile; scanFile.open("CC_LCL-FullData.csv"); string line; getline(scanFile, line);
  map<unsigned int, vector<double>> dateTimeMap; set<unsigned int> dateTimeSet; vector<unsigned int> dateTimeVector; double maxEnergy=22;
  while (getline(scanFile, line)) {
    string token; unsigned int houseID, date, time; double energy;
    stringstream ss; ss.str(line);
    getline(ss, token, ','); houseID=stoi(token.substr(3, 6));
    getline(ss, token, ',');
    getline(ss, token, ','); date=stoul(token.substr(0, 4)+token.substr(5, 2)+token.substr(8, 2)); time=2*stoul(token.substr(11, 2))+(stoul(token.substr(14,2))==30 ? 1 : 0);
    getline(ss, token); energy=(token=="Null\r" ? -1 : stod(token)*2);
    if (energy>-1e-9) {
      unsigned int dateTime=date*100+time; dateTimeMap[dateTime].emplace_back(energy);
      if (dateTimeSet.emplace(dateTime).second) {dateTimeVector.emplace_back(dateTime);}
    }
  } scanFile.close();
  
  random_device rd; mt19937 gen(rd());
  for (unsigned int sampleIndex=0; sampleIndex<10; ++sampleIndex) {
    unsigned int dateTime=dateTimeVector[gen()%dateTimeVector.size()];
    ofstream printFileDistribution("Distribution"+to_string(sampleIndex)+"_"+to_string(dateTime)+".txt");
    unsigned int binNumber=44; vector<unsigned int> distribution(binNumber+1);
    for (auto &energy : dateTimeMap[dateTime]) {++distribution[binNumber*energy/maxEnergy];}
    for (unsigned int binIndex=0; binIndex<=binNumber; ++binIndex) {printFileDistribution << ((binIndex+0.5)/binNumber)*maxEnergy << "\t" << distribution[binIndex] << "\n";}
    printFileDistribution.close();
  }
}
