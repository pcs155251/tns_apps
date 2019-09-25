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
  UniTensor<double> temp;
  temp.Load("./A.ten");
  UniTensor<complex<double>> ansaz( temp );
  cout<<ansaz<<endl;

  UniTensor<double> hamiltonian = transverseIsing( 0.5, 3.1, false );
  channel TIM(dimMPS, ansaz, hamiltonian);

  TIM.findAllEnv( iterMax, canonicalErr, arpErr );
  TIM.measureEnergy( true );
  TIM.SaveEnvTensors( rootFolder );

  return 0;
}
