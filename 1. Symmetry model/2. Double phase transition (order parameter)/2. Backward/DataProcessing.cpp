#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<random>
#include<iomanip>
#include<sstream>
using namespace std;

int main() {
  unsigned int N=256, bootstrapNumber=1e4;
  ofstream frequencyOrderParameterFile("FrequencyOrderParameter.txt");
  frequencyOrderParameterFile << "lambda" << "\t" << "F" << "\t" << "0.1Quantile" << "\t" << "0.9Quantile" << "\n";
  ofstream susceptibilityFile("Susceptibility.txt");
  susceptibilityFile << "lambda" << "\t" << "chi" << "\t" << "0.1Quantile" << "\t" << "0.9Quantile" << "\n";
  random_device rd; mt19937 gen(rd());
  for (unsigned int dataIndex=0; dataIndex<=100; ++dataIndex) {
    double lambda=0.2*dataIndex;
    string fileName="FrequencyOrderParameter_lambda="+to_string(lambda)+".txt";
    ifstream inputFile(fileName); vector<double> datumVector; datumVector.reserve(1e3);
    double datum; while (inputFile >> datum) {datumVector.push_back(datum);} inputFile.close();
    unsigned int datumNumber=datumVector.size();
    //F
    double average=0, squareAverage=0;
    for (auto &datum : datumVector) {average+=datum; squareAverage+=datum*datum;}
    average/=datumNumber; squareAverage/=datumNumber;
    sort(datumVector.begin(), datumVector.end());
    frequencyOrderParameterFile << lambda << "\t" << average << "\t" << datumVector[0.1*datumNumber+0.5] << "\t" << datumVector[0.9*datumNumber+0.5] << "\n";
    //chi
    double susceptibility=N*(squareAverage-average*average);
    vector<double> bootstrapSusceptibilityVector(bootstrapNumber);
    for (unsigned int bootstrapIndex=0; bootstrapIndex<bootstrapNumber; ++bootstrapIndex) {
      double bootstrapAverage=0, bootstrapSquareAverage=0;
      for (unsigned int i=0; i<datumNumber; ++i) {
        double bootstrapDatum=datumVector[gen()%datumNumber];
        bootstrapAverage+=bootstrapDatum;
        bootstrapSquareAverage+=bootstrapDatum*bootstrapDatum;
      } bootstrapAverage/=datumNumber; bootstrapSquareAverage/=datumNumber;
      bootstrapSusceptibilityVector[bootstrapIndex]=N*(bootstrapSquareAverage-bootstrapAverage*bootstrapAverage);
    }
    sort(bootstrapSusceptibilityVector.begin(), bootstrapSusceptibilityVector.end());
    susceptibilityFile << lambda << "\t" << susceptibility << "\t" << bootstrapSusceptibilityVector[0.1*bootstrapNumber+0.5] << "\t" << bootstrapSusceptibilityVector[0.9*bootstrapNumber+0.5] << "\n";
  }
  frequencyOrderParameterFile.close(); susceptibilityFile.close();
}
