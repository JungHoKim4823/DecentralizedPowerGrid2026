#include<cmath>
#include<random>
#include<array>
#include<vector>
#include<algorithm>
#include<fstream>
using namespace std;

const unsigned int N=256;
const double alpha=0.5, beta=0.1;
const double ip=19/32.0, dp=1/16.0, fp=ip;
const double ilambda=0, dlambda=0.2, flambda=20;
const double it=0, dt=0.01, rt=900, ft=1000;
const double synchronizationThreshold=1e-12;
const unsigned int ensembleNumber=1000;

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
  for (unsigned int ensembleIndex=1; ensembleIndex<=ensembleNumber; ++ensembleIndex) {
    for (double p=ip; p<fp+0.5*dp; p+=dp) {
      nodeVector=vector<Node>(N); linkVector.clear(); linkVector.reserve(N*3/2);
      //Consumer
      for (auto &node : nodeVector) {
        double consumption=-1;
        node.P+=consumption; node.I+=abs(consumption)*alpha; node.D+=abs(consumption)*beta;
      }
      //Generator
      unsigned int N_g=N*p+0.5;
      for (auto &node : nodeVector) {
        if (&node-&nodeVector.front()==N_g) {break;}
        double generation=(double)N/N_g;
        node.P+=generation; node.I+=abs(generation)*alpha; node.D+=abs(generation)*beta; node.generator=true;
      }
      shuffle(nodeVector.begin(), nodeVector.end(), gen);
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
      for (double lambda=ilambda; flambda>=ilambda ? lambda<flambda+0.5*dlambda : lambda>flambda+0.5*dlambda; lambda+=dlambda) {
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
        frequencyOrderParameter/=N*(ft-rt)/dt;
        ofstream printFile("FrequencyOrderParameter_lambda="+to_string(lambda)+".txt", ios::app); printFile.precision(9);
        printFile << frequencyOrderParameter << "\n";
        printFile.close();
      }
    }
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
