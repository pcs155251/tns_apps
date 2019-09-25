#include "channel.h"

using namespace std;
using namespace uni10;

int main(){
  const int dimD = 2;
  const int dimMPS = 12;
  const int iterMax = 10000;
  const double canonicalErr = 10e-14;
  const double arpErr = 10e-14;
  const string rootFolder = "./tensors/";
  vector<Bond> ansazBds( 4, Bond( BD_OUT, dimD ) );
  ansazBds.insert( ansazBds.begin(), Bond( BD_IN, 2 ) );
  UniTensor<double> ansaz( ansazBds );
  ansaz.Load("./A.ten");

  UniTensor<double> hamiltonian = transverseIsing( 0.5, 3.1, false );
  channel TIM(dimMPS, ansaz, hamiltonian);

  TIM.LoadEnvTensors( rootFolder );
  TIM.measureEnergy( true );

  UniTensor<complex<double>> gradient = TIM.findGradient();
  cout<<gradient<<endl;
  gradient.Save("./gradient.ten");

  return 0;
}
