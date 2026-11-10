#include<cmath>
#include<vector>
#include<fstream>
using namespace std;

extern const unsigned int N, ensembleNumber;
extern const double alpha, beta, ip, dp, fp, ilambda, dlambda, flambda, it, dt, rt, ft, synchronizationThreshold;
extern const string suffix;

class Datum {
public:
  vector<double> numberVector;
  Datum() {numberVector.reserve(ensembleNumber);}
  void Add(double number) {numberVector.emplace_back(number);}
  unsigned int SampleNumber() {
    unsigned int sampleNumber=0;
    for (auto &number : numberVector) {if (number>synchronizationThreshold) {++sampleNumber;}}
    return sampleNumber;
  }
  double Average(unsigned int exponent) {
    double sum=0;
    for (auto &number : numberVector) {sum+=pow(number, exponent);}
    return sum/numberVector.size();
  }
};

class Data {
public:
  vector<Datum> datumVector=vector<Datum>((unsigned int)(((fp-ip)/dp)+1.5)*(unsigned int)(((flambda-ilambda)/dlambda)+1.5));
  void Add(unsigned int datumIndex, double number) {datumVector[datumIndex].Add(number);}
  void Print(string filename) {
    ofstream printFile(filename+".txt"); printFile.precision(9); unsigned int lambdaPointNumber=(((flambda-ilambda)/dlambda)+1.5);
    printFile << "N=" << N << "\talpha=" << alpha << "\tbeta=" << beta << "\tit=" << it << "\tdt=" << dt << "\trt=" << rt << "\tft=" << ft << "\tSynchronization threshold=" << synchronizationThreshold << "\n";
    printFile << "Ensembles" << "\tSamples" << "\tp" << "\tlambda" << "\t<X>" << "\t<X^2>" << "\t<X^4>" << "\n";
    for (auto &datum : datumVector) {
      unsigned int datumIndex=&datum-&datumVector.front();
      printFile << datum.numberVector.size() << "\t" << datum.SampleNumber() << "\t" << (ip+dp*(datumIndex/lambdaPointNumber)) << "\t" << (ilambda+dlambda*(datumIndex%lambdaPointNumber)) << "\t";
      printFile << datum.Average(1) << "\t"  << datum.Average(2) << "\t" << datum.Average(4) << "\n";
    } printFile.close();
  }
};
