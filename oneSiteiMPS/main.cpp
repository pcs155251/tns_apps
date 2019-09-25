#include "oneSiteiMPS.h"

using namespace std;
using namespace uni10;

int main(){
  const int phyDim = 2;
  const int virDim = 8;
  const int iterMax = 1000;
  const double errTol = 1.0e-14;
  const double arpTol = 1.0e-14;

  UniTensor<double> ansaz;
  ansaz.Load( "./A.ten" );
  cout<<ansaz<<endl;
  UniTensor<complex<double>> cAnsaz = ansaz;
  cout<<cAnsaz<<endl;
  oneSiteiMPS<double> channel( virDim, cAnsaz );
  channel.fixedPoint( iterMax, errTol, arpTol );
  channel.testCanonical();
  channel.printTensors();


  return 0;
}
