#include"Data.hpp"
#include<cmath>
#include<random>
#include<array>
#include<vector>
#include<queue>
#include<algorithm>
using namespace std;

const unsigned int N=256; const double f=0.05;
const double alpha_c=0.05, beta_c=0.02, alpha_g=0.5, beta_g=0.1;
const double ilambda=0, dlambda=3, flambda=300;
const double it=0, dt=0.01, rt=900, ft=1000;
const double synchronizationThreshold=1e-12;
const unsigned int ensembleNumber=1000;

const double ip=1/32.0, dp=1/16.0, fp=7/32.0; const string suffix="Part1";
//const double ip=9/32.0, dp=1/16.0, fp=15/32.0; const string suffix="Part2";
//const double ip=17/32.0, dp=1/16.0, fp=23/32.0; const string suffix="Part3";
//const double ip=25/32.0, dp=1/16.0, fp=31/32.0; const string suffix="Part4";

class Node {
public:
  double I=0, D=0, P=0, inverseI, DInverseI, PInverseI, omega=0, theta=0, frequency; unsigned int generator=false;
  unsigned int componentIndex; unsigned int visit, giantComponent;
  vector<Node*> neighborNodeVector;
};
class Network {
public:
  vector<Node> nodeVector; vector<pair<unsigned int, unsigned int>> linkVector;
  unsigned int giantComponentSize;
  double frequencyOrderParameter;
  
  Network();
  void Derivatives(double lambda, array<double, N>& thetaVector, array<double, N>& omegaVector, array<double, N>& dTheta, array<double, N>& dOmega);
  void FindComponents();
};

Network::Network() {
  random_device rd; mt19937 gen(rd()); uniform_real_distribution<double> urd(0, 1);
  array<double, N> RK4TempThetaArray, k1Array, k2Array, k3Array, k4Array;
  array<double, N> RK4TempOmegaArray, l1Array, l2Array, l3Array, l4Array;
  Data frequencyOrderParameterData=Data();
  for (unsigned int ensembleIndex=1; ensembleIndex<=ensembleNumber; ++ensembleIndex) {
    unsigned int datumIndex=0;
    for (double p=ip; p<fp+0.5*dp; p+=dp) {
      nodeVector=vector<Node>(N); linkVector.clear(); linkVector.reserve(N*3);
      //Consumer
      for (auto &node : nodeVector) {
        double consumption=-1;
        node.P+=consumption; node.I+=abs(consumption)*alpha_c; node.D+=abs(consumption)*beta_c;
      }
      //Generator
      unsigned int N_g=N*p+0.5;
      for (auto &node : nodeVector) {
        if (&node-&nodeVector.front()==N_g) {break;}
        double generation=(double)N/N_g;
        node.P+=generation; node.I+=abs(generation)*alpha_g; node.D+=abs(generation)*beta_g; node.generator=true;
      }
      shuffle(nodeVector.begin(), nodeVector.end(), gen);
      //Division
      for (auto &node : nodeVector) {node.inverseI=1/node.I; node.PInverseI=node.P/node.I; node.DInverseI=node.D/node.I;}
      //Network construction
      unsigned int L=sqrt(N)+0.5;
      for (unsigned int node1Index=0; node1Index<N; ++node1Index) {
        Node *node1, *node2, *node3; unsigned int connect=node1Index/L+node1Index%L;
        if (connect%2==0) {
          node1=&nodeVector[node1Index]; node2=&nodeVector[(node1Index+1)%L+(node1Index/L)*L];
          if (f>urd(gen)) {
            if (gen()%2==0) {do {node2=&nodeVector[gen()%N];} while (node1==node2 || find(node1->neighborNodeVector.begin(), node1->neighborNodeVector.end(), node2)!=node1->neighborNodeVector.end());}
            else {do {node1=&nodeVector[gen()%N];} while (node1==node2 || find(node2->neighborNodeVector.begin(), node2->neighborNodeVector.end(), node1)!=node2->neighborNodeVector.end());}
          }
          node1->neighborNodeVector.emplace_back(node2); node2->neighborNodeVector.emplace_back(node1);
          linkVector.emplace_back(pair<unsigned int, unsigned int>(node1-&nodeVector.front(), node2-&nodeVector.front()));
        }
        node1=&nodeVector[node1Index], node3=&nodeVector[(node1Index+L)%N];
        if (f>urd(gen)) {
          if (gen()%2==0) {do {node3=&nodeVector[gen()%N];} while (node1==node3 || find(node1->neighborNodeVector.begin(), node1->neighborNodeVector.end(), node3)!=node1->neighborNodeVector.end());}
          else {do {node1=&nodeVector[gen()%N];} while (node1==node3 || find(node3->neighborNodeVector.begin(), node3->neighborNodeVector.end(), node1)!=node3->neighborNodeVector.end());}
        }
        node1->neighborNodeVector.emplace_back(node3); node3->neighborNodeVector.emplace_back(node1);
        linkVector.emplace_back(pair<unsigned int, unsigned int>(node1-&nodeVector.front(), node3-&nodeVector.front()));
      }
      //Find the giant component
      giantComponentSize=0; FindComponents();
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
            double averageFrequency=0; for (auto &node : nodeVector) {if (node.giantComponent) {averageFrequency+=node.frequency;}} averageFrequency/=giantComponentSize;
            for (auto &node : nodeVector) {if (node.giantComponent) {frequencyOrderParameter+=pow(node.frequency-averageFrequency, 2);}}
          }
        }
        frequencyOrderParameter/=giantComponentSize*(ft-rt)/dt; frequencyOrderParameterData.Add(datumIndex, frequencyOrderParameter);
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

void Network::FindComponents() {
  giantComponentSize=0; unsigned int componentIndex=0, giantComponentIndex;
  for (auto &node : nodeVector) {node.visit=false; node.giantComponent=false;}
  for (auto &node : nodeVector) {
    if (node.visit) {continue;}
    unsigned int componentSize=0; ++componentIndex; queue<Node*> BFSAlgorithm;
    BFSAlgorithm.emplace(&node); node.visit=true; node.componentIndex=componentIndex; ++componentSize;
    while (!BFSAlgorithm.empty()) {
      Node* parentNode=BFSAlgorithm.front();
      for (auto &childNode : parentNode->neighborNodeVector) {
        if (childNode->visit) {continue;}
        BFSAlgorithm.emplace(childNode); childNode->visit=true; childNode->componentIndex=componentIndex; ++componentSize;
      }
      BFSAlgorithm.pop();
    }
    if (giantComponentSize<componentSize) {giantComponentSize=componentSize; giantComponentIndex=componentIndex;}
  }
  for (auto &node : nodeVector) {if (node.componentIndex==giantComponentIndex) {node.giantComponent=true;}}
}

int main() {Network();}
