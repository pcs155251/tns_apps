#include "uni10.hpp"
#include "paraIpeps.h"
#include "ipeps.h"
#include "aux.h"
#include "../nsy_ham/hamiltonian.h"
#include <cstdio>

using namespace std;
using namespace uni10;

int main( int argc, char* argv[] ){
  paraIpeps input( argc, argv );
  ipeps<double> tim( input.virDim ); 
  UniTensor<double> ham = transverseIsing( 0.5, input.fields[0], false );
  
  for ( int i=0; i!=input.fields.size(); ++i )
  {
    ham = transverseIsing( 0.5, input.fields[i], false );

    if (input.ifsu)
    {
      for ( int ii=0; ii!=input.taus.size(); ++ii)
      {
        tim.simpUpdate( ham, input.taus[ii], input.suStop, input.suIter );
        meaSimpTimQuants( tim, ham, input.fields[i], input.dataFolder  );
      }
      tim.saveTensors( input.tensorFolder, input.fields[i] );
    } else {}

    if (input.iffu)
    {
      tim.loadEn( input.enLoadFolder, input.fieldLoad, input.edgeDimL );
      tim.changeEdgeDim( input.edgeDimU );
      tim.computeEn( input.edgeDimU, input.ctmIter, input.ctmStop, ham);
      for ( int ii=0; ii!=input.taus.size(); ++ii)
      {
        tim.fullUpdateGauge( ham, input.taus[ii], input.fuStop, input.fuIter, input.ctmIter, input.ctmStop, input.edgeDimU );
        meaFullTimQuants( tim, ham, input.fields[i], input.dataFolder );
      }
      tim.saveTensors( input.tensorFolder, input.fields[i] );
      tim.saveEn( input.tensorFolder, input.fields[i] );

    } else {}
  }

  return 0;
}

