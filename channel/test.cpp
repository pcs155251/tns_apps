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
  vector<double> tempElem( temp.ElemNum(), -1.0 );
  tempElem[0] = 1.0;
  temp.SetRawElem(tempElem);
  UniTensor<complex<double>> ansaz( temp );

  UniTensor<double> hamiltonian = transverseIsing( 0.5, 3.1, false );
  channel TIM(dimMPS, temp, hamiltonian);

//  TIM.findAllEnv( iterMax, canonicalErr, arpErr );
  //double eold = TIM.measureEnergy( true );
  /*
  UniTensor<complex<double>> gradient0 = TIM.findGradient();
  cout<<gradient0<<endl;
  gradient0 *= 1.0/ maxAbsElem( gradient0 ) ;
  cout<<gradient0<<endl;
  */

  /*
  //cg tensor 
  UniTensor<double> deltaTenOld( ansaz.bond() );
  UniTensor<double> gradientOld( ansaz.bond() );
  deltaTenOld.Zeros();
  gradientOld.Zeros();
  double alpha = 0.1;
  for (int i=0; i!=10000; i++){
    //UniTensor<complex<double>> cgradient = TIM.findGradient();
    //cgradient *= 1.0/cgradient[0];
    UniTensor<double> gradient = cgradient;
    double beta = (i==0)? 0:pow( Norm( gradient.GetBlock() )/Norm( gradientOld.GetBlock() ) , 2);
    UniTensor<double> deltaTen = -1.0*gradient+beta*deltaTenOld;
    deltaTenOld = deltaTen;
    gradientOld = gradient;

    TIM.addTensor( alpha*deltaTen );
    TIM.findAllEnv( iterMax, canonicalErr, arpErr );
    TIM.measureEnergy( true );
  }
  */

  //cg numerical
  UniTensor<double> deltaTenOld( ansaz.bond() );
  UniTensor<double> gradientOld( ansaz.bond() );
  deltaTenOld.Zeros();
  gradientOld.Zeros();
  double m1=0.2;
  double alpha0=10.0;
  for (int i=0; i!=10000; i++){
    TIM.findAllEnv( iterMax, canonicalErr, arpErr );
    double eold = TIM.measureEnergy( true ).real();
    UniTensor<double> gradient = TIM.numericalGradientReal( 1.0e-6, iterMax, canonicalErr, arpErr );
    double beta = (i==0)? 0:pow( Norm( gradient.GetBlock() )/Norm( gradientOld.GetBlock() ) , 2);
    UniTensor<double> deltaTen = -1.0*gradient+beta*deltaTenOld;

    //line search
    gradient.SetLabel({0,1,2,3,4});
    deltaTen.SetLabel({0,1,2,3,4});
    double deltaDer = Contract( gradient, deltaTen, false )[0];
    double alphaR = alpha0;
    TIM.addTensor( alphaR*deltaTen );

    while (1){
      FILE *result = fopen( "result.dat", "a");
      fprintf( result, "line search " );
      fflush( result );
      fclose( result );

      TIM.addTensor( -0.5*alphaR*deltaTen );
      TIM.findAllEnv( iterMax, canonicalErr, arpErr );
      double etemp = TIM.measureEnergy( false ).real();
      if (etemp>(eold+m1*0.5*alphaR*deltaDer)){
        alphaR=0.5*alphaR;
      }
      else {
        TIM.normalize();
        FILE *result = fopen( "result.dat", "a");
        fprintf( result, "done \n" );
        fflush( result );
        fclose( result );
        break;
      }
    }
      

    deltaTenOld = deltaTen;
    gradientOld = gradient;
  }

  return 0;
}
