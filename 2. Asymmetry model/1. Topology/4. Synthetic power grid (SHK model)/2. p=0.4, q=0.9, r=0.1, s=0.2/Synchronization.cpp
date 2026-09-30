#include"Data.hpp"
#include<cmath>
#include<random>
#include<array>
#include<vector>
#include<queue>
#include<algorithm>
#include<string>
#include<sstream>
#include<fstream>
using namespace std;

const unsigned int N=256;
const double alpha_c=0.05, beta_c=0.02, alpha_g=0.5, beta_g=0.1;
const double p_SHK=0.4, q_SHK=0.9, r_SHK=0.1, s_SHK=0.2;
const double ilambda=0, dlambda=3, flambda=300;
const double it=0, dt=0.005, rt=900, ft=1000;
const double synchronizationThreshold=1e-12;
const unsigned int ensembleNumber=500;
double ip, dp, fp; string suffix="";

class Node {public: double x, y, I=0, D=0, P=0, inverseI, DInverseI, PInverseI, omega=0, theta=0, frequency; unsigned int generator=false, visit; vector<Node*> nodeVector;};
class Network {
public:
  vector<Node> nodeVector; vector<pair<unsigned int, unsigned int>> linkVector;
  double frequencyOrderParameter;
  
  Network();
  void Derivatives(double lambda, array<double, N>& thetaVector, array<double, N>& omegaVector, array<double, N>& dTheta, array<double, N>& dOmega);
};

Network::Network() {
  random_device rd; mt19937 gen(rd()); uniform_real_distribution<double> urd(0, 1);
  array<double, N> RK4TempThetaArray, k1Array, k2Array, k3Array, k4Array;
  array<double, N> RK4TempOmegaArray, l1Array, l2Array, l3Array, l4Array;
  Data frequencyOrderParameterData=Data();
  for (unsigned int ensembleIndex=1; ensembleIndex<=ensembleNumber; ++ensembleIndex) {
    unsigned int datumIndex=0;
    for (double p=ip; p<fp+0.5*dp; p+=dp) {
      nodeVector=vector<Node>(N); linkVector.clear(); linkVector.reserve(N*3/2);
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
      unsigned int node1Index=0, node2Index=1; Node *node1=&nodeVector[node1Index], *node2=&nodeVector[node2Index]; unsigned int E=0;
      node1->x=urd(gen); node1->y=urd(gen); node2->x=urd(gen); node2->y=urd(gen);
      node1->nodeVector.emplace_back(node2); node2->nodeVector.emplace_back(node1); linkVector.emplace_back(pair<unsigned int, unsigned int>(node1Index, node2Index)); ++E;
      for (unsigned int node_iIndex=2; node_iIndex<N; ++node_iIndex) {
        Node* node_i=&nodeVector[node_iIndex];
        if (s_SHK>urd(gen)) {
          //G5
          pair<unsigned int, unsigned int>* link=&linkVector[gen()%E]; unsigned int node_aIndex=link->first, node_bIndex=link->second; Node *node_a=&nodeVector[node_aIndex], *node_b=&nodeVector[node_bIndex];
          node_i->x=(node_a->x+node_b->x)/2; node_i->y=(node_a->y+node_b->y)/2;
          node_i->nodeVector.emplace_back(node_a);
          *find(node_a->nodeVector.begin(), node_a->nodeVector.end(), node_b)=node_i;
          link->second=node_iIndex;
          node_i->nodeVector.emplace_back(node_b);
          *find(node_b->nodeVector.begin(), node_b->nodeVector.end(), node_a)=node_i;
          linkVector.emplace_back(pair<unsigned int, unsigned int>(node_iIndex, node_bIndex)); ++E;
        } else {
          //G1
          node_i->x=urd(gen); node_i->y=urd(gen);
          //G2
          double minimumDistance2=2; unsigned int node_jIndex; Node* node_j;
          for (unsigned int nodeIndex=0; nodeIndex<node_iIndex; ++nodeIndex) {
            Node* node=&nodeVector[nodeIndex];
            double spatialDistance2=(node->x-node_i->x)*(node->x-node_i->x)+(node->y-node_i->y)*(node->y-node_i->y);
            if (minimumDistance2>spatialDistance2) {minimumDistance2=spatialDistance2; node_jIndex=nodeIndex; node_j=node;}
          }
          node_i->nodeVector.emplace_back(node_j); node_j->nodeVector.emplace_back(node_i); linkVector.emplace_back(pair<unsigned int, unsigned int>(node_iIndex, node_jIndex)); ++E;
          if (p_SHK>urd(gen)) {
            //G3
            Node* node_l=nullptr; double maximumTargetFunction=0;
            for (auto &node : nodeVector) {node.visit=false;}
            queue<Node*> BFSQueue; unsigned int distance=1; Node indicator=Node();
            BFSQueue.emplace(node_i); node_i->visit=true; BFSQueue.emplace(&indicator);
            while (!BFSQueue.empty()) {
              Node* parentNode=BFSQueue.front();
              for (auto &childNode : parentNode->nodeVector) {
                if (childNode->visit) {continue;}
                BFSQueue.emplace(childNode); childNode->visit=true;
                if (childNode==node_j) {continue;}
                double targetFunction=pow(distance+1, r_SHK)/sqrt((childNode->x-node_i->x)*(childNode->x-node_i->x)+(childNode->y-node_i->y)*(childNode->y-node_i->y));
                if (maximumTargetFunction<targetFunction) {maximumTargetFunction=targetFunction; node_l=childNode;}
              }
              if (BFSQueue.front()==&indicator) {++distance; if (BFSQueue.size()!=1) {BFSQueue.emplace(&indicator);}}
              BFSQueue.pop();
            }
            if (node_l!=nullptr) {
              unsigned int node_lIndex=node_l-&nodeVector.front();
              node_i->nodeVector.emplace_back(node_l); node_l->nodeVector.emplace_back(node_i); linkVector.emplace_back(pair<unsigned int, unsigned int>(node_iIndex, node_lIndex)); ++E;
            }
          }
          if (q_SHK>urd(gen)) {
            //G4
            unsigned int node_ipIndex=gen()%node_iIndex; Node *node_ip=&nodeVector[node_ipIndex], *node_lp=nullptr; double maximumTargetFunction=0;
            for (auto &node : nodeVector) {node.visit=false;}
            queue<Node*> BFSQueue; unsigned int distance=1; Node indicator=Node();
            BFSQueue.emplace(node_ip); node_ip->visit=true; BFSQueue.emplace(&indicator);
            while (!BFSQueue.empty()) {
              Node* parentNode=BFSQueue.front();
              for (auto &childNode : parentNode->nodeVector) {
                if (childNode->visit) {continue;}
                BFSQueue.emplace(childNode); childNode->visit=true;
                if (distance==1) {continue;}
                if (childNode==node_i) {continue;}
                double targetFunction=pow(distance+1, r_SHK)/sqrt((childNode->x-node_ip->x)*(childNode->x-node_ip->x)+(childNode->y-node_ip->y)*(childNode->y-node_ip->y));
                if (maximumTargetFunction<targetFunction) {maximumTargetFunction=targetFunction; node_lp=childNode;}
              }
              if (BFSQueue.front()==&indicator) {++distance; if (BFSQueue.size()!=1) {BFSQueue.emplace(&indicator);}}
              BFSQueue.pop();
            }
            if (node_lp!=nullptr) {
              unsigned int node_lpIndex=node_lp-&nodeVector.front();
              node_ip->nodeVector.emplace_back(node_lp); node_lp->nodeVector.emplace_back(node_ip); linkVector.emplace_back(pair<unsigned int, unsigned int>(node_ipIndex, node_lpIndex)); ++E;
            }
          }
        }
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
    if (ensembleIndex%(unsigned int)pow(10, (unsigned int)log10(ensembleIndex))==0 && ensembleIndex>=100) {frequencyOrderParameterData.Print("FrequencyOrderParameter"+suffix+"_"+to_string(ensembleIndex));}
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

int main(int argc, char* argv[]) {
  int arg=stoi(argv[1]);
  if (arg==0) {ip=1/32.0; dp=1/16.0; fp=ip; suffix="_Part01";}
  else if (arg==1) {ip=3/32.0; dp=1/16.0; fp=ip; suffix="_Part02";}
  else if (arg==2) {ip=5/32.0; dp=1/16.0; fp=ip; suffix="_Part03";}
  else if (arg==3) {ip=7/32.0; dp=1/16.0; fp=ip; suffix="_Part04";}
  else if (arg==4) {ip=9/32.0; dp=1/16.0; fp=ip; suffix="_Part05";}
  else if (arg==5) {ip=11/32.0; dp=1/16.0; fp=ip; suffix="_Part06";}
  else if (arg==6) {ip=13/32.0; dp=1/16.0; fp=ip; suffix="_Part07";}
  else if (arg==7) {ip=15/32.0; dp=1/16.0; fp=ip; suffix="_Part08";}
  else if (arg==8) {ip=17/32.0; dp=1/16.0; fp=ip; suffix="_Part09";}
  else if (arg==9) {ip=19/32.0; dp=1/16.0; fp=ip; suffix="_Part10";}
  else if (arg==10) {ip=21/32.0; dp=1/16.0; fp=ip; suffix="_Part11";}
  else if (arg==11) {ip=23/32.0; dp=1/16.0; fp=ip; suffix="_Part12";}
  else if (arg==12) {ip=25/32.0; dp=1/16.0; fp=ip; suffix="_Part13";}
  else if (arg==13) {ip=27/32.0; dp=1/16.0; fp=ip; suffix="_Part14";}
  else if (arg==14) {ip=29/32.0; dp=1/16.0; fp=ip; suffix="_Part15";}
  else if (arg==15) {ip=31/32.0; dp=1/16.0; fp=ip; suffix="_Part16";}
  Network();
}

