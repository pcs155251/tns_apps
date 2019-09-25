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

  UniTensor<double> hamiltonian = transverseIsing( 0.5, 3.1, false );
  //channel TIM(dimMPS, ansaz, hamiltonian);
  channel TIM(dimMPS, temp, hamiltonian);

  TIM.findAllEnv( iterMax, canonicalErr, arpErr );
  TIM.measureEnergy( true );

  /*
  UniTensor<double> numerGrad = TIM.numericalGradient( 1.0e-3, iterMax, canonicalErr, arpErr );
  numerGrad.Save("./numerGrad3.ten");
  numerGrad = TIM.numericalGradient( 1.0e-4, iterMax, canonicalErr, arpErr );
  numerGrad.Save("./numerGrad4.ten");
  numerGrad = TIM.numericalGradient( 1.0e-5, iterMax, canonicalErr, arpErr );
  numerGrad.Save("./numerGrad5.ten");
  */

  UniTensor<double> deltaTenOld( ansaz.bond() );
  UniTensor<double> gradientOld( ansaz.bond() );
  deltaTenOld.Zeros();
  gradientOld.Zeros();
  double alpha = 0.1;
  for (int i=0; i!=10000; i++){
    UniTensor<complex<double>> cgradient = TIM.findGradient();
    cgradient *= 1.0/cgradient[0];
    UniTensor<double> gradient = cgradient;
    double beta = (i==0)? 0:pow( Norm( gradient.GetBlock() )/Norm( gradientOld.GetBlock() ) , 2);
    UniTensor<double> deltaTen = -1.0*gradient+beta*deltaTenOld;
    deltaTenOld = deltaTen;
    gradientOld = gradient;

    TIM.addTensor( alpha*deltaTen );
    TIM.findAllEnv( iterMax, canonicalErr, arpErr );
    TIM.measureEnergy( true );
  }

  return 0;
}
