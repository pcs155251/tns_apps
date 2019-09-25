#include "enPess.h"

#define IF_DEBUG false

using namespace std;
using namespace uni10;

template<typename T>
UniTensor<T> constructBten( const vector<UniTensor<T>> &coreTensors, const vector<UniTensor<T>> &projTensors, const bool ifUp )
{
  uni10::UniTensor<double> out;
  if (ifUp){
    Network_dev out_net( "./pessNets/aUp.net" );
    ContractArgs( out, out_net, coreTensors[0], coreTensors[1], projTensors[0], projTensors[1], projTensors[2] );//8, 7, -5, -4; 1, 2, 3
  } else {
    Network_dev out_net( "./pessNets/aDn.net" );
    ContractArgs( out, out_net, coreTensors[0], coreTensors[1], projTensors[0], projTensors[1], projTensors[2] );//5, 4, 8, 7; 1, 2, 3
  }
  out.SetLabel( vector<int> {1,2,3,4,5,6,7} );

  return out;
}

template<typename T>
UniTensor<T> getaten( UniTensor<T> &bten )
{
  bten.SetLabel( vector<int> {1,2,3,4,5,6,7} );
  UniTensor<double> bDag = Conj( bten );
  bDag.SetLabel( vector<int> {-1,-2,-3,-4,5,6,7} );
  UniTensor<double> aten = Contract( bten, bDag, false );
  combineTwoLayer( aten );
  return aten;
}

template<typename T>
enPess<T>::enPess(  vector<UniTensor<T>>& coreTensors, vector<UniTensor<T>>& projTensors, const int edgeDimIn, bool ifUp ):
ctmBase<T> ( 2, 2, edgeDimIn )
{
  UniTensor<T> bten = constructBten( coreTensors, projTensors, ifUp );
  UniTensor<T> aten = getaten( bten );
  vector<UniTensor<T>> unit( 4, aten );
  ctmBase<T>::setGroups( unit );

  auxTens = vector<UniTensor<T>> { bten };
}

template<typename T>
double enPess<T>::ctmIteration( const int niter, const double errTol, UniTensor<T> &triSiteHam, bool ifsvd ){
  if (ifsvd)
  {
    cout<<endl<<"SVD edge dim = "<<ctmBase<T>::ctmDim<<endl<<flush;
  }
  else 
  {
    cout<<endl<<"QR edge dim = "<<ctmBase<T>::ctmDim<<endl<<flush;
  }
  double energyOld = 0;
  for (int i=0; i!=niter; i++){
    ctmBase<T>::ctmrgAstep( ifsvd ); 
    double energy = measureTriSite( triSiteHam )/measureNorm( );
    cout<<i<<" "<<setprecision(10)<<energy<<" "<<flush<<endl;
    if (fabs(energy-energyOld)<errTol){
      return energy;
    }
    else {
      energyOld = energy;
    }
  }
  return energyOld;
}

template<typename T>
double enPess<T>::measureQuantities( UniTensor<T> &ham, double &energy, vector<vector<double>> &spinComp  )
{
  cout<<"measureQuantities"<<endl<<flush;

  double normVal = measureNorm( );
  energy = measureTriSite( ham )/normVal;
  UniTensor<double> spinX( matSx(0.5) );
  UniTensor<double> spinZ( matSz(0.5) );
  vector<double> staggerM(3);
  for (int isite=0; isite!=3; isite++){
    spinComp[isite][0] = measureOneSite( spinX, isite )/normVal;
    spinComp[isite][1] = measureOneSite( spinZ, isite )/normVal;

    staggerM[isite] = sqrt( spinComp[isite][0]*spinComp[isite][0] + spinComp[isite][1]*spinComp[isite][1] );
  }


  return (2.0/3.0)*energy;
}

template<typename T>
double enPess<T>::measureNorm( ){
  const int posit = 0;
  UniTensor<T> norm;
  Network_dev measureNorm_net( "./enPessNets/measureNorm.net" );
  ContractArgs( norm, measureNorm_net, 
                ctmBase<T>::groups.at(posit).corners.at(0), ctmBase<T>::groups.at(posit).corners.at(1), ctmBase<T>::groups.at(posit).corners.at(2), ctmBase<T>::groups.at(posit).corners.at(3), 
                ctmBase<T>::groups.at(posit).edges.at(0)  , ctmBase<T>::groups.at(posit).edges.at(1)  , ctmBase<T>::groups.at(posit).edges.at(2)  , ctmBase<T>::groups.at(posit).edges.at(3)  , 
                ctmBase<T>::groups.at(posit).core );
  return norm[0];
}

template<typename T>
double enPess<T>::measureOneSite( UniTensor<T> &oneSiteOp, const int isite ){
  assert( (isite>=0)&&(isite<3) );
  const int posit = 0;
  auxTens[0].SetLabel( vector<int> {1,2,3,4,5,6,7} );
  switch(isite){
    case 0: oneSiteOp.SetLabel( vector<int> { 8,5} ); break;
    case 1: oneSiteOp.SetLabel( vector<int> { 9,6} ); break; 
    case 2: oneSiteOp.SetLabel( vector<int> {10,7} ); break;
  }
  UniTensor<T> bop = Contract( auxTens[0], oneSiteOp, false );//1,2,3,4,8,9,10
  UniTensor<T> bdag = Conj( auxTens[0] );
  bdag.SetLabel( vector<int> {-1,-2,-3,-4,5,6,7} );
  switch(isite){
    case 0: bop.SetLabel ( vector<int> {1,2,3,4,6,7,5} ); break;
    case 1: bop.SetLabel ( vector<int> {1,2,3,4,5,7,6} ); break;
    case 2: bop.SetLabel ( vector<int> {1,2,3,4,5,6,7} ); break;
  }
  UniTensor<T> bOpBdag = Contract( bop, bdag, false );
  combineTwoLayer( bOpBdag );

  UniTensor<T> expec;
  Network_dev measureTriSite_net( "./enPessNets/measureNorm.net" );
  ContractArgs( expec, measureTriSite_net, 
                ctmBase<T>::groups.at(posit).corners.at(0), ctmBase<T>::groups.at(posit).corners.at(1), ctmBase<T>::groups.at(posit).corners.at(2), ctmBase<T>::groups.at(posit).corners.at(3), 
                ctmBase<T>::groups.at(posit).edges.at(0)  , ctmBase<T>::groups.at(posit).edges.at(1)  , ctmBase<T>::groups.at(posit).edges.at(2)  , ctmBase<T>::groups.at(posit).edges.at(3)  , 
                bOpBdag );
  return expec[0];
}

template<typename T>
double enPess<T>::measureTriSite( UniTensor<T> &threeSiteOp ){
  const int posit = 0;
  UniTensor<T> bDag = Conj(auxTens[0]);
  auxTens[0].SetLabel( vector<int> {1,2,3,4,5,6,7} );
  threeSiteOp.SetLabel( vector<int> {8,9,10,5,6,7} );
  UniTensor<T> bop = Contract( auxTens[0], threeSiteOp, false );//1,2,3,4,8,9,10
  UniTensor<T> bdag = Conj( auxTens[0] );
  bdag.SetLabel( vector<int> {-1,-2,-3,-4,8,9,10} );
  UniTensor<T> bOpBdag = Contract( bop, bdag, false );
  combineTwoLayer( bOpBdag );
  UniTensor<T> expec;
  Network_dev measureTriSite_net( "./enPessNets/measureNorm.net" );
  ContractArgs( expec, measureTriSite_net, 
                ctmBase<T>::groups.at(posit).corners.at(0), ctmBase<T>::groups.at(posit).corners.at(1), ctmBase<T>::groups.at(posit).corners.at(2), ctmBase<T>::groups.at(posit).corners.at(3), 
                ctmBase<T>::groups.at(posit).edges.at(0)  , ctmBase<T>::groups.at(posit).edges.at(1)  , ctmBase<T>::groups.at(posit).edges.at(2)  , ctmBase<T>::groups.at(posit).edges.at(3)  , 
                bOpBdag );
  return expec[0];
}

template class enPess<double>;
