#include "hamiltonian.h"

using namespace std;
using namespace uni10;

int main(){

  //UniTensor<double> ham_XXZ = XXZ(1., 1.);
  //cout << ham_XXZ;

  //UniTensor<double> ham_hei = Heisenberg();
  //cout << ham_hei;

  cout << transverseIsing(0.5, 2. * 1.05, true);

  return 0;
}
