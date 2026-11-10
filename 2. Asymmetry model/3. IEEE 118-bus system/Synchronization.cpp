#include"Data.hpp"
#include<cmath>
#include<random>
#include<array>
#include<vector>
#include<tuple>
#include<algorithm>
#include<string>
#include<fstream>
#include<iostream>
#include<stdexcept>
using namespace std;

const unsigned int N=118; const double p0=19.0/N;
const double alpha_c=0.05, beta_c=0.02, alpha_g=0.5, beta_g=0.1;
const double ilambda=0, dlambda=0.2, flambda=20;
const double it=0, dt=0.001, rt=900, ft=1000;
const double lowerBoundI=0.002, lowerBoundD=0.001, synchronizationThreshold=1e-12;
const unsigned int ensembleNumber=200;
double ip, dp, fp; string suffix;

class Node {public: double I=0, D=0, P=0, inverseI, DInverseI, PInverseI, omega=0, theta=0, frequency; unsigned int generator=false;};
class Network {
public:
  vector<Node> nodeVector; vector<tuple<unsigned int, unsigned int, double>> linkVector;
  double frequencyOrderParameter;
  
  Network();
  void Derivatives(double lambda, array<double, N>& thetaVector, array<double, N>& omegaVector, array<double, N>& dTheta, array<double, N>& dOmega);
};

Network::Network() {
  random_device rd; mt19937 gen(rd());
  array<double, N> RK4TempThetaArray, k1Array, k2Array, k3Array, k4Array;
  array<double, N> RK4TempOmegaArray, l1Array, l2Array, l3Array, l4Array;
  Data frequencyOrderParameterData=Data();
  
  for (unsigned int ensembleIndex=1; ensembleIndex<=ensembleNumber; ++ensembleIndex) {
    array<double, N> rescaledConsumptionArray{}, rescaledGenerationArray{};
    vector<unsigned int> originalGeneratorIndexVector, candidateGeneratorIndexVector;
    ifstream scanFileNode("ieee118cdf.txt");
    string line; getline(scanFileNode, line); getline(scanFileNode, line);
    double totalLoad=0, totalGeneration=0;
    while (getline(scanFileNode, line)) {
      if (line.find("-999")!=string::npos) {break;}
      const unsigned int busNumber=stoi(line.substr(0, 4));
      const unsigned int nodeIndex=busNumber-1;
      const double consumption=stod(line.substr(40, 9)), generation=stod(line.substr(59, 8));
      rescaledConsumptionArray[nodeIndex]=consumption+max(-generation, 0.0);
      rescaledGenerationArray[nodeIndex]=max(generation, 0.0);
      totalLoad+=rescaledConsumptionArray[nodeIndex]; totalGeneration+=rescaledGenerationArray[nodeIndex];
    }
    scanFileNode.close();
    for (unsigned int nodeIndex=0; nodeIndex<N; ++nodeIndex) {
      rescaledConsumptionArray[nodeIndex]*=N/totalLoad;
      rescaledGenerationArray[nodeIndex]*=N/totalGeneration;
      if (rescaledGenerationArray[nodeIndex]>0) {originalGeneratorIndexVector.emplace_back(nodeIndex);}
      else {candidateGeneratorIndexVector.emplace_back(nodeIndex);}
    }
    shuffle(candidateGeneratorIndexVector.begin(), candidateGeneratorIndexVector.end(), gen);
    unsigned int datumIndex=0;
    for (double p=ip; p<fp+0.5*dp; p+=dp) {
      nodeVector=vector<Node>(N);
      //Consumer
      for (auto &node : nodeVector) {
        const double consumption=-rescaledConsumptionArray[&node-&nodeVector.front()];
        node.P+=consumption; node.I+=abs(consumption)*alpha_c; node.D+=abs(consumption)*beta_c;
      }
      //Generator
      const unsigned int N_g=(unsigned int)(N*p+0.5); const double rho=(p-p0)/(1-p0);
      const unsigned int additionalGeneratorNumber=N_g-originalGeneratorIndexVector.size();;
      for (auto &nodeIndex : originalGeneratorIndexVector) {nodeVector[nodeIndex].generator=true;}
      for (unsigned int nodeIndex=0; nodeIndex<additionalGeneratorNumber; ++nodeIndex) {nodeVector[candidateGeneratorIndexVector[nodeIndex]].generator=true;}
      for (auto &node : nodeVector) {
        const double distributedGeneration=node.generator ? rho*N/N_g : 0;
        const double generation=(1-rho)*rescaledGenerationArray[&node-&nodeVector.front()]+distributedGeneration;
        node.P+=generation; node.I+=generation*alpha_g; node.D+=generation*beta_g;
      }
      //Lower bound
      for (auto &node : nodeVector) {if (node.I<lowerBoundI) {node.I=lowerBoundI;} if (node.D<lowerBoundD) {node.D=lowerBoundD;}}
      //Division
      for (auto &node : nodeVector) {node.inverseI=1/node.I; node.PInverseI=node.P/node.I; node.DInverseI=node.D/node.I;}
      //Network construction
      linkVector.clear(); linkVector.reserve(186);
      ifstream scanFileLink("ieee118cdf.txt");
      while (getline(scanFileLink, line) && line.find("BRANCH DATA FOLLOWS")==string::npos) {}
      while (getline(scanFileLink, line)) {
        if (line.find("-999")!=string::npos) {break;}
        const unsigned int node1Index=stoi(line.substr(0, 4))-1, node2Index=stoi(line.substr(5, 4))-1;
        const unsigned int type=stoi(line.substr(18, 1));
        const double X=stod(line.substr(29, 11)), CDFTap=stod(line.substr(76, 6));
        double tap=1; if (type!=0 && CDFTap>0) {tap=CDFTap;}
        linkVector.emplace_back(node1Index, node2Index, 1/(tap*X));
      }
      scanFileLink.close();
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
    if (ensembleIndex%(unsigned int)pow(10, (unsigned int)log10(ensembleIndex))==0 && ensembleIndex>=100) {frequencyOrderParameterData.Print("FrequencyOrderParameter"+suffix+"_"+to_string(ensembleIndex));}
  }
}

void Network::Derivatives(double lambda, array<double, N>& thetaVector, array<double, N>& omegaVector, array<double, N>& dTheta, array<double, N>& dOmega) {
  for (unsigned int i=0; i<N; ++i) {
    dTheta[i]=omegaVector[i];
    dOmega[i]=nodeVector[i].PInverseI-nodeVector[i].DInverseI*omegaVector[i];
  }
  for (auto &link : linkVector) {
    const unsigned int node1Index=get<0>(link), node2Index=get<1>(link);
    const double coupling=lambda*get<2>(link)*sin(thetaVector[node2Index]-thetaVector[node1Index]);
    dOmega[node1Index]+=coupling*nodeVector[node1Index].inverseI; dOmega[node2Index]-=coupling*nodeVector[node2Index].inverseI;
  }
}

int main(int argc, char* argv[]) {
  int arg=stoi(argv[1]);
  if (arg==0) {ip=7/32.0; dp=1/16.0; fp=ip; suffix="_Part04";}
  else if (arg==1) {ip=9/32.0; dp=1/16.0; fp=ip; suffix="_Part05";}
  else if (arg==2) {ip=11/32.0; dp=1/16.0; fp=ip; suffix="_Part06";}
  else if (arg==3) {ip=13/32.0; dp=1/16.0; fp=ip; suffix="_Part07";}
  else if (arg==4) {ip=15/32.0; dp=1/16.0; fp=ip; suffix="_Part08";}
  else if (arg==5) {ip=17/32.0; dp=1/16.0; fp=ip; suffix="_Part09";}
  else if (arg==6) {ip=19/32.0; dp=1/16.0; fp=ip; suffix="_Part10";}
  else if (arg==7) {ip=21/32.0; dp=1/16.0; fp=ip; suffix="_Part11";}
  else if (arg==8) {ip=23/32.0; dp=1/16.0; fp=ip; suffix="_Part12";}
  else if (arg==9) {ip=25/32.0; dp=1/16.0; fp=ip; suffix="_Part13";}
  else if (arg==10) {ip=27/32.0; dp=1/16.0; fp=ip; suffix="_Part14";}
  else if (arg==11) {ip=29/32.0; dp=1/16.0; fp=ip; suffix="_Part15";}
  else if (arg==12) {ip=31/32.0; dp=1/16.0; fp=ip; suffix="_Part16";}
  Network();
}
