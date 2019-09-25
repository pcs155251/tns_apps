#include "pess.h"
#include "paraPess.h"
#include "quant.h"
#include "time.h"

using namespace std;
using namespace uni10;

int main( int argc, char* argv[]){
  paraPess paras( argc, argv );
  quantities result;

  //construct Hamiltonian
  UniTensor<double> ham = Heisenberg( 0.5, 1.0 );
  ham.PutBlock( ham.GetBlock() + paras.field*(DMinteraction().GetBlock() ) );
  ham = periodicHamiltonian( 3, ham );

  for ( int id=0; id!=paras.dimDs.size(); ++id )
  {
    int dimVir = paras.dimDs[id];
    int dimCTM = ceil(paras.start*dimVir*dimVir);
    //pess<double,enPess> kagome( dimVir );
    pess<double,enRedPess> kagome( dimVir );

    if (paras.ifGs)
    {
      for (int it=0; it!=paras.taus.size(); ++it)
      {
        kagome.optimize2ndOrder( ham, paras.taus[it], paras.gsIter, paras.gsErr );
      }
      kagome.saveWv( paras.tensorFolder, paras.field );
    }
    else
    {
      kagome.loadWv( paras.loadFolder, paras.field );
    }

    if (paras.ifEn)
    {
      kagome.initEn( paras.start*dimVir*dimVir );

      double eold = 0;
      result.writeColName( paras.dataFolder, paras.field, dimVir );
      for ( int iedge=0; iedge!=paras.edgeDimN ; ++iedge )
      {
        time_t start = clock();
        int dimCTM = ceil((paras.start+iedge*paras.slope)*dimVir*dimVir);
        kagome.changeEdgeDim( dimCTM );

        kagome.computeEn( paras.ctmStep, 0.1*paras.edgeErr, ham, paras.ifsvd ) ;
        kagome.meaQuantities( ham, result.energyUp, result.energyDn, result.spinComp );
        time_t end = clock();
        result.elapsedT = (end-start)/(double)(CLOCKS_PER_SEC);
        result.printAll( paras.dataFolder, paras.field, dimVir, dimCTM );
        double ediff = fabs( eold - result.ePerSite );

        if (ediff<paras.edgeErr){
          kagome.saveEn( paras.enviFolder, paras.field  );
          break;
        }
        else {
          eold = result.ePerSite;
        }
      }
    }
  }

  return 0;
}

