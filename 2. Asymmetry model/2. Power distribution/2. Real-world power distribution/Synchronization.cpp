#include"Data.hpp"
#include<cmath>
#include<random>
#include<array>
#include<vector>
#include<set>
#include<map>
#include<algorithm>
#include<sstream>
#include<fstream>
using namespace std;

const unsigned int N=256;
const double alpha_c=0.05, beta_c=0.02, alpha_g=0.5, beta_g=0.1;
const double ilambda=0, dlambda=3, flambda=300;
const double it=0, dt=0.001, rt=900, ft=1000;
const double lowerBoundI=0.002, lowerBoundD=0.001, synchronizationThreshold=1e-12;
const unsigned int ensembleNumber=1000;

const double ip=1/32.0, dp=1/16.0, fp=ip; const string suffix="Part1";
//const double ip=3/32.0, dp=1/16.0, fp=ip; const string suffix="Part2";
//const double ip=5/32.0, dp=1/16.0, fp=ip; const string suffix="Part3";
//const double ip=7/32.0, dp=1/16.0, fp=ip; const string suffix="Part4";
//const double ip=9/32.0, dp=1/16.0, fp=ip; const string suffix="Part5";
//const double ip=11/32.0, dp=1/16.0, fp=ip; const string suffix="Part6";
//const double ip=13/32.0, dp=1/16.0, fp=ip; const string suffix="Part7";
//const double ip=15/32.0, dp=1/16.0, fp=ip; const string suffix="Part8";
//const double ip=17/32.0, dp=1/16.0, fp=ip; const string suffix="Part9";
//const double ip=19/32.0, dp=1/16.0, fp=ip; const string suffix="Part10";
//const double ip=21/32.0, dp=1/16.0, fp=ip; const string suffix="Part11";
//const double ip=23/32.0, dp=1/16.0, fp=ip; const string suffix="Part12";
//const double ip=25/32.0, dp=1/16.0, fp=ip; const string suffix="Part13";
//const double ip=27/32.0, dp=1/16.0, fp=ip; const string suffix="Part14";
//const double ip=29/32.0, dp=1/16.0, fp=ip; const string suffix="Part15";
//const double ip=31/32.0, dp=1/16.0, fp=ip; const string suffix="Part16";

class Node {public: double I=0, D=0, P=0, inverseI, DInverseI, PInverseI, omega=0, theta=0, frequency; unsigned int generator=false;};
class Network {
public:
  vector<Node> nodeVector; vector<pair<unsigned int, unsigned int>> linkVector;
  double frequencyOrderParameter;
  
  Network();
  void Derivatives(double lambda, array<double, N>& thetaVector, array<double, N>& omegaVector, array<double, N>& dTheta, array<double, N>& dOmega);
};

Network::Network() {
  random_device rd; mt19937 gen(rd());
  array<double, N> RK4TempThetaArray, k1Array, k2Array, k3Array, k4Array;
  array<double, N> RK4TempOmegaArray, l1Array, l2Array, l3Array, l4Array;
  Data frequencyOrderParameterData=Data();
  //Consumption
  ifstream scanFileConsumption; scanFileConsumption.open("Consumption/CC_LCL-FullData.csv"); string line; getline(scanFileConsumption, line);
  map<unsigned int, vector<double>> dateTimeConsumptionMap; set<unsigned int> dateTimeSet; vector<unsigned int> dateTimeVector;
  while (getline(scanFileConsumption, line)) {
    string token; unsigned int houseID, date, time; double energy;
    stringstream ss; ss.str(line);
    getline(ss, token, ','); houseID=stoi(token.substr(3, 6));
    getline(ss, token, ',');
    getline(ss, token, ','); date=stoul(token.substr(0, 4)+token.substr(5, 2)+token.substr(8, 2)); time=2*stoul(token.substr(11, 2))+(stoul(token.substr(14,2))==30 ? 1 : 0);
    getline(ss, token); energy=(token=="Null\r" ? -1 : stod(token)*2);
    if (energy>-1e-9) {
      unsigned int dateTime=date*100+time; dateTimeConsumptionMap[dateTime].emplace_back(energy);
      if (dateTimeSet.emplace(dateTime).second) {dateTimeVector.emplace_back(dateTime);}
    }
  } scanFileConsumption.close();
  //Generation
  ifstream scanFilePVGeneration; scanFilePVGeneration.open("Generation/metadata.csv"); getline(scanFilePVGeneration, line);
  vector<double> generationVector;
  while (getline(scanFilePVGeneration, line)) {
    string token; double energy;
    stringstream ss; ss.str(line);
    getline(ss, token, ','); getline(ss, token, ','); getline(ss, token, ','); getline(ss, token, ','); getline(ss, token, ','); getline(ss, token, ','); getline(ss, token, ',');
    getline(ss, token, ','); energy=stod(token);
    generationVector.emplace_back(energy);
  } scanFilePVGeneration.close();
  for (unsigned int ensembleIndex=1; ensembleIndex<=ensembleNumber; ++ensembleIndex) {
    unsigned int datumIndex=0;
    for (double p=ip; p<fp+0.5*dp; p+=dp) {
      nodeVector=vector<Node>(N); linkVector.clear(); linkVector.reserve(N*3/2);
      //Consumer
      const unsigned int dateTime=dateTimeVector[gen()%dateTimeVector.size()];
      vector<double> consumptionSampleVector(N); double sumConsumptionSample=0;
      for (auto &consumptionSample : consumptionSampleVector) {
        consumptionSample=dateTimeConsumptionMap[dateTime][gen()%dateTimeConsumptionMap[dateTime].size()];
        sumConsumptionSample+=consumptionSample;
      }
      for (auto &consumptionSample : consumptionSampleVector) {consumptionSample*=N/sumConsumptionSample;}
      for (auto &node : nodeVector) {
        double consumption=-consumptionSampleVector[&node-&nodeVector.front()];
        node.P+=consumption; node.I+=abs(consumption)*alpha_c; node.D+=abs(consumption)*beta_c;
      }
      shuffle(nodeVector.begin(), nodeVector.end(), gen);
      //Generator
      unsigned int N_g=N*p+0.5; vector<double> generationSampleVector(N_g); double sumGenerationSample=0;
      for (auto &generationSample : generationSampleVector) {
        generationSample=generationVector[gen()%generationVector.size()];
        sumGenerationSample+=generationSample;
      }
      for (auto &generationSample : generationSampleVector) {generationSample*=N/sumGenerationSample;}
      for (auto &node : nodeVector) {
        if (&node-&nodeVector.front()==N_g) {break;}
        double generation=generationSampleVector[&node-&nodeVector.front()];
        node.P+=generation; node.I+=abs(generation)*alpha_g; node.D+=abs(generation)*beta_g; node.generator=true;
      }
      shuffle(nodeVector.begin(), nodeVector.end(), gen);
      //Lower bound
      for (auto &node : nodeVector) {if (node.I<lowerBoundI) {node.I=lowerBoundI;} if (node.D<lowerBoundD) {node.D=lowerBoundD;}}
      //Division
      for (auto &node : nodeVector) {node.inverseI=1/node.I; node.PInverseI=node.P/node.I; node.DInverseI=node.D/node.I;}
      //Network construction
      unsigned int L=sqrt(N)+0.5;
      for (unsigned int node1Index=0; node1Index<N; ++node1Index) {
        unsigned int node2Index, node3Index; unsigned int connect=node1Index/L+node1Index%L;
        if (connect%2==0) {node2Index=(node1Index+1)%L+(node1Index/L)*L; linkVector.emplace_back(pair<unsigned int, unsigned int>(node1Index, node2Index));}
        node3Index=(node1Index+L)%N; linkVector.emplace_back(pair<unsigned int, unsigned int>(node1Index, node3Index));
      }
      //Kuramoto model
      for (double lambda=ilambda; flambda>=ilambda ? lambda<flambda+0.5*dlambda : lambda>flambda+0.5*dlambda; lambda+=dlambda, ++datumIndex) {
        frequencyOrderParameter=0;
        for (double t=it; t<ft+0.5*dt; t+=dt) {
          //4th Order Runge-Kutta
          for (unsigned int i=0; i<N; ++i) {RK4TempThetaArray[i]=nodeVector[i].theta; RK4TempOmegaArray[i]=nodeVector[i].omega;}
          Derivatives(lambda, RK4TempThetaArray, RK4TempOmegaArray, k1Array, l1Array);
          for (unsigned int i=0; i<N; ++i) {RK4TempThetaArray[i]=nodeVector[i].theta+dt*k1Array[i]/2; RK4TempOmegaArray[i]=nodeVector[i].omega+dt*l1Array[i]/2;}
          Derivatives(lambda, RK4TempThetaArray, RK4TempOmegaArray, k2Array, l2Array);
          for (unsigned int i=0; i<N; ++i) {RK4TempThetaArray[i]=nodeVector[i].theta+dt*k2Array[i]/2; RK4TempOmegaArray[i]=nodeVector[i].omega+dt*l2Array[i]/2;}
          Derivatives(lambda, RK4TempThetaArray, RK4TempOmegaArray, k3Array, l3Array);
          for (unsigned int i=0; i<N; ++i) {RK4TempThetaArray[i]=nodeVector[i].theta+dt*k3Array[i]; RK4TempOmegaArray[i]=nodeVector[i].omega+dt*l3Array[i];}
          Derivatives(lambda, RK4TempThetaArray, RK4TempOmegaArray, k4Array, l4Array);
          for (unsigned int i=0; i<N; ++i) {
            auto &node=nodeVector[i];
            double deltaTheta=dt*(k1Array[i]+2*k2Array[i]+2*k3Array[i]+k4Array[i])/6, deltaOmega=dt*(l1Array[i]+2*l2Array[i]+2*l3Array[i]+l4Array[i])/6;
            node.frequency=deltaTheta/dt; node.theta+=deltaTheta; node.omega+=deltaOmega;
            while (node.theta>M_PI) {node.theta-=2*M_PI;} while (node.theta<-M_PI) {node.theta+=2*M_PI;}
          }
          //Results
          if (t>rt+0.5*dt) {
            double averageFrequency=0; for (auto &node : nodeVector) {averageFrequency+=node.frequency;} averageFrequency/=N;
            for (auto &node : nodeVector) {frequencyOrderParameter+=pow(node.frequency-averageFrequency, 2);}
          }
        }
        frequencyOrderParameter/=N*(ft-rt)/dt; frequencyOrderParameterData.Add(datumIndex, frequencyOrderParameter);
      }
    }
    if (ensembleIndex%(unsigned int)pow(10, (unsigned int)log10(ensembleIndex))==0) {frequencyOrderParameterData.Print("FrequencyOrderParameter"+suffix+"_"+to_string(ensembleIndex));}
  }
}

void Network::Derivatives(double lambda, array<double, N>& thetaVector, array<double, N>& omegaVector, array<double, N>& dTheta, array<double, N>& dOmega) {
  for (unsigned int i=0; i<N; ++i) {
    dTheta[i]=omegaVector[i];
    dOmega[i]=nodeVector[i].PInverseI-nodeVector[i].DInverseI*omegaVector[i];
  }
  for (auto &link : linkVector) {
    const unsigned int node1Index=link.first, node2Index=link.second;
    const double coupling=lambda*sin(thetaVector[node2Index]-thetaVector[node1Index]);
    dOmega[node1Index]+=coupling*nodeVector[node1Index].inverseI; dOmega[node2Index]-=coupling*nodeVector[node2Index].inverseI;
  }
}

int main() {Network();}
