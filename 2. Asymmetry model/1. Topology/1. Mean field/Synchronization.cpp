#include"Data.hpp"
#include<cmath>
#include<random>
#include<array>
#include<vector>
#include<queue>
#include<algorithm>
using namespace std;

const unsigned int N=256;
const double alpha_c=0.05, beta_c=0.02, alpha_g=0.5, beta_g=0.1;
const double ilambda=0, dlambda=3, flambda=1200;
const double it=0, dt=0.01, rt=900, ft=1000;
const double synchronizationThreshold=1e-12;
const unsigned int ensembleNumber=1000;

const double ip=1/32.0, dp=1/16.0, fp=7/32.0; const string suffix="Part1";
//const double ip=9/32.0, dp=1/16.0, fp=15/32.0; const string suffix="Part2";
//const double ip=17/32.0, dp=1/16.0, fp=23/32.0; const string suffix="Part3";
//const double ip=25/32.0, dp=1/16.0, fp=31/32.0; const string suffix="Part4";

class Node {public: double I=0, D=0, P=0, inverseI, omega=0, theta, frequency; unsigned int generator=false;};
class Network {
public:
  vector<Node> nodeVector;
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
    unsigned int datumIndex=0;
    for (double p=ip; p<fp+0.5*dp; p+=dp) {
      nodeVector=vector<Node>(N);
      uniform_real_distribution<double> urdTheta(-M_PI, M_PI);
      //Consumer
      for (auto &node : nodeVector) {
        double consumption=-1;
        node.P+=consumption; node.I+=abs(consumption)*alpha_c; node.D+=abs(consumption)*beta_c; node.theta=urdTheta(gen);
      }
      //Generator
      unsigned int N_g=N*p+0.5;
      for (auto &node : nodeVector) {
        if (&node-&nodeVector.front()==N_g) {break;}
        double generation=(double)N/N_g;
        node.P+=generation; node.I+=abs(generation)*alpha_g; node.D+=abs(generation)*beta_g; node.generator=true;
      }
      //Division
      for (auto &node : nodeVector) {node.inverseI=1/node.I;}
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
  double real=0, imaginary=0, radius, phase;
  for (unsigned int i=0; i<N; ++i) {real+=cos(thetaVector[i]); imaginary+=sin(thetaVector[i]);} real/=N; imaginary/=N;
  radius=sqrt(real*real+imaginary*imaginary); phase=atan2(imaginary, real);
  for (unsigned int i=0; i<N; ++i) {
    dTheta[i]=omegaVector[i];
    dOmega[i]=(nodeVector[i].P-nodeVector[i].D*omegaVector[i]+lambda*radius*sin(phase-thetaVector[i]))*nodeVector[i].inverseI;
  }
}

int main() {Network();}
