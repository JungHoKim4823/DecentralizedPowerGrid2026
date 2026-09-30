#include"Data.hpp"
#include<cmath>
#include<random>
#include<array>
#include<vector>
#include<algorithm>
using namespace std;

const unsigned int N=256; const double RSD=0.2;
const double alpha_c=0.05, beta_c=0.02, alpha_g=0.5, beta_g=0.1;
const double ilambda=0, dlambda=3, flambda=300;
const double it=0, dt=0.01, rt=900, ft=1000;
const double synchronizationThreshold=1e-12;
const unsigned int ensembleNumber=1000;

const double ip=1/32.0, dp=1/16.0, fp=7/32.0; const string suffix="Part1";
//const double ip=9/32.0, dp=1/16.0, fp=15/32.0; const string suffix="Part2";
//const double ip=17/32.0, dp=1/16.0, fp=23/32.0; const string suffix="Part3";
//const double ip=25/32.0, dp=1/16.0, fp=31/32.0; const string suffix="Part4";

class Node {public: double I=0, D=0, P=0, inverseI, DInverseI, PInverseI, omega=0, theta=0, frequency; unsigned int generator=false;};
class Network {
public:
  vector<Node> nodeVector; vector<pair<unsigned int, unsigned int>> linkVector;
  double frequencyOrderParameter;
  
  Network();
  void Derivatives(double lambda, array<double, N>& thetaVector, array<double, N>& omegaVector, array<double, N>& dTheta, array<double, N>& dOmega);
  double inverseNormalCDF(double p, double mu, double sigma);
  vector<double> generateRegularGaussianSample(unsigned int N, double mu, double sigma);
};

Network::Network() {
  random_device rd; mt19937 gen(rd());
  array<double, N> RK4TempThetaArray, k1Array, k2Array, k3Array, k4Array;
  array<double, N> RK4TempOmegaArray, l1Array, l2Array, l3Array, l4Array;
  Data frequencyOrderParameterData=Data();
  for (unsigned int ensembleIndex=1; ensembleIndex<=ensembleNumber; ++ensembleIndex) {
    unsigned int datumIndex=0;
    for (double p=ip; p<fp+0.5*dp; p+=dp) {
      nodeVector=vector<Node>(N); linkVector.clear(); linkVector.reserve(N*3/2);
      //Consumer
      vector<double> regularGaussianSample=generateRegularGaussianSample(N, -1, RSD);
      for (auto &node : nodeVector) {
        unsigned int nodeIndex=&node-&nodeVector.front();
        double consumption=regularGaussianSample[nodeIndex];
        node.P+=consumption; node.I+=abs(consumption)*alpha_c; node.D+=abs(consumption)*beta_c;
      }
      shuffle(nodeVector.begin(), nodeVector.end(), gen);
      //Generator
      unsigned int N_g=N*p+0.5;
      regularGaussianSample=generateRegularGaussianSample(N_g, (double)N/N_g, ((double)N/N_g)*RSD);
      for (auto &node : nodeVector) {
        unsigned nodeIndex=&node-&nodeVector.front();
        if (nodeIndex==N_g) {break;}
        double generation=regularGaussianSample[nodeIndex];
        node.P+=generation; node.I+=abs(generation)*alpha_g; node.D+=abs(generation)*beta_g; node.generator=true;
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

double Network::inverseNormalCDF(double p, double mu, double sigma) {
  //assert(p>0.0 && p<1.0); // p는 (0,1) 사이 여야 함
  // 사용된 근사: probit 함수 (Acklam approximation)
  const double a[]={-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  const double b[]={-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
  const double c[]={-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00,  2.938163982698783e+00};
  const double d[]={7.784695709041462e-03,  3.224671290700398e-01, 2.445134137142996e+00,  3.754408661907416e+00};
  double q, r;
  if (p<0.02425) {
    // Lower region
    q=sqrt(-2*log(p));
    return mu+sigma*(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p>1-0.02425) {
    // Upper region
    q=sqrt(-2*log(1-p));
    return mu-sigma*(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else {
    // Central region
    q=p-0.5; r=q*q;
    return mu+sigma*(((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  }
}

vector<double> Network::generateRegularGaussianSample(unsigned int N, double mu, double sigma) {
  vector<double> regularGaussianSample(N);
  for (unsigned int i=0; i<N; ++i) {
    double p=(i+0.5)/N;  // 균등하게 분할된 CDF 값 (avoid 0 or 1)
    regularGaussianSample[i]=inverseNormalCDF(p, mu, sigma);
  }
  return regularGaussianSample;
}

int main() {Network();}
