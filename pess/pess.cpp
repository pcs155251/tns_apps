#include "pess.h"

using namespace std;
using namespace uni10;

/*
template<typename T, template<typename T> typename enType>
UniTensor<double> pess<T,enType>::constructBten( const bool ifUp )
{
  uni10::UniTensor<double> out;
  if (ifUp){
    Network out_net( "./pessNets/aUp.net" );
    ContractArgs( out, out_net, coreTensors[0], coreTensors[1], projTensors[0], projTensors[1], projTensors[2] );//8, 7, -5, -4; 1, 2, 3
  } else {
    Network out_net( "./pessNets/aDn.net" );
    ContractArgs( out, out_net, coreTensors[0], coreTensors[1], projTensors[0], projTensors[1], projTensors[2] );//5, 4, 8, 7; 1, 2, 3
  }
  out.SetLabel( vector<int> {1,2,3,4,5,6,7} );

  return out;
}

template<typename T, template<typename T> class enType>
UniTensor<double> pess<T,enType>::getaten( UniTensor<double> &bten )
{
  bten.SetLabel( vector<int> {1,2,3,4,5,6,7} );
  UniTensor<double> bDag = Conj( bten );
  bDag.SetLabel( vector<int> {-1,-2,-3,-4,5,6,7} );
  UniTensor<double> aten = Contract( bten, bDag );
  combineTwoLayer( aten );
  //aten = Permute( aten, vector<int> {1,-1,2,-2,3,-3,4,-4}, 8 );
  return aten;
}
*/

template<typename T, template<typename T> typename enType>
pess<T,enType>::pess( const int dimVirIn ): dimVir{dimVirIn}, thetaLabForHosvd( vector<vector<int>> ( 2, vector<int> (2*nBdsOfCore)))
{

  vector<Bond> coreBds( nBdsOfCore, Bond( BD_IN, dimVir) );
  UniTensor<T> oneCore( coreBds );
  oneCore.Randomize();
  coreTensors = vector<UniTensor<T>> ( nCore, oneCore );

  vector<Bond> projBds( nCore, Bond( BD_OUT, dimVir ) );
  projBds.insert( projBds.begin(), Bond( BD_IN, 2 ) );
  UniTensor<T> oneProj( projBds );
  oneProj.Randomize();
  projTensors = vector<UniTensor<T>> ( nBdsOfCore, oneProj );

  Matrix<T> oneBondVector( dimVir, dimVir, true );
  uni10_rand( oneBondVector, uni10_mt19937, uni10_uniform_real, 0, 1, uni10_clock);
  bondVectors = vector<vector<Matrix<T>>> ( nCore, vector<Matrix<T>> (nBdsOfCore, oneBondVector) );

  //set core tensor labels
  vector<int> coreLabUp(nBdsOfCore);
  vector<int> coreLabDn(nBdsOfCore);
  for (int j=0; j!=nBdsOfCore; j++){
    coreLabUp.at(j) =   nBdsOfCore+1+j;
    coreLabDn.at(j) = -(nBdsOfCore+1+j);
  }
  coreTensors.at(0).SetLabel( coreLabUp );
  coreTensors.at(1).SetLabel( coreLabDn );
  //set proj tensor labels
  vector<int> projLab(3);
  for (int i=0; i!=nBdsOfCore; i++){
    projLab.at(0) = i+1;
    projLab.at(1) =   nBdsOfCore+1+i;
    projLab.at(2) = -(nBdsOfCore+1+i);
    projTensors.at(i).SetLabel(projLab);
  }

  //set thetaLabForHosvd
  for (int iCore=0; iCore!=2; iCore++){
    for (int i=0; i!=nBdsOfCore; i++){
      thetaLabForHosvd.at(iCore).at(2*i)   = -(i+1);
      thetaLabForHosvd.at(iCore).at(2*i+1) = (iCore)? (i+nBdsOfCore+1):-(i+nBdsOfCore+1);
    }
  }
}

template<typename T, template<typename T> class enType>
void pess<T,enType>::setHamLab( UniTensor<T> &ham )
{
  vector<int> hamLab(2*nBdsOfCore);
  for (int i=0; i!=nBdsOfCore; i++){
    hamLab.at(i)            = -(i+1);
    hamLab.at(nBdsOfCore+i) =   i+1;
  }
  ham.SetLabel( hamLab );
}

template<typename T, template<typename T> class enType>
double pess<T,enType>::meaTriSiteLocal( const UniTensor<T> &oper, const int iCore ){
  UniTensor<T> operCopy = oper;
  UniTensor<T> theta  = ContractCore(iCore);
  UniTensor<T> thetaT = Conj( theta );
  UniTensor<T> value = Contract( theta, thetaT );

  //vector<int> tempLab = {1,2,3,-1,-2,-3};
  vector<int> tempLab = {-1,-2,-3,1,2,3};
  operCopy.SetLabel(tempLab);
  UniTensor<T> thetaH = Contract( theta, operCopy );
  vector<int> newLab = thetaH.label();
  for (int i=0; i!=nBdsOfCore; i++){
    newLab.at( nBdsOfCore+i ) *= -1;
  }
  thetaH.SetLabel( newLab );
  UniTensor<T> expec = Contract( thetaH, thetaT );
  T energy = expec[0]/value[0];
  complex<double> complexEnergy = (complex<double> (0,0))+energy;
  assert( fabs(complexEnergy.imag())<1.0e-14 );
  return complexEnergy.real();
}

template<typename T, template<typename T> class enType>
int pess<T,enType>::optimize1stOrder( UniTensor<T> ham, const double tau, const int maxIter, const double criterion ){
  setHamLab( ham );
  UniTensor<T> evoOperator = getEvoOperator(ham, tau);

  printf( "\nfirst order evolution\n" );
  printf( "%8s%12s%12s\n", "istep", "diff", "energy" );
  vector<vector<Matrix<T>>> bondVectorOlds = bondVectors;
  double energyOld = 0;
  for ( int istep=0; istep!=maxIter; istep++){
    timeEvoAstep( evoOperator, 0 );
    timeEvoAstep( evoOperator, 1 );
  
    vector<vector<double>> diff( 2, vector<double> (nBdsOfCore) );
    double energy = ( meaTriSiteLocal( ham, 0 ) + meaTriSiteLocal( ham, 1 ) )/3.0;
    double energyDiff = abs(energyOld-energy);
   
    if (istep%100==0){
      printf( "%8d%12.4e%12.8f\n", istep, energyDiff, energy );
    } else { }

    if ( energyDiff>criterion ){
      energyOld = energy;
      bondVectorOlds = bondVectors;
    }
    else {
      double energy = ( meaTriSiteLocal( ham, 0 ) + meaTriSiteLocal( ham, 1 ) )/3.0;
      printf( "%8d%12.4e%12.8f\n", istep, energyDiff, energy );
      printf( "evolution converges\n" );
      return istep;
    }
  }
  printf( "\ncannot converge...\n" );
  return -1;
}

template<typename T, template<typename T> class enType>
int pess<T,enType>::optimize2ndOrder( UniTensor<T> ham, const double tau, const int maxIter, const double criterion ){
  setHamLab( ham );
  UniTensor<T> evoOperatorFull = getEvoOperator(ham,     tau);
  UniTensor<T> evoOperatorHalf = getEvoOperator(ham, 0.5*tau);

  printf( "\nsecond order evolution\n" );
  printf( "%8s%12s%12s\n", "istep", "diff", "energy" );
  vector<vector<Matrix<T>>> bondVectorOlds = bondVectors;
  double energyOld = 0;
  for ( int istep=0; istep!=maxIter; istep++){
    timeEvoAstep( evoOperatorHalf, 0 );
    timeEvoAstep( evoOperatorFull, 1 );
    timeEvoAstep( evoOperatorHalf, 0 );
  
    vector<vector<double>> diff( 2, vector<double> (nBdsOfCore) );
    double energy = ( meaTriSiteLocal( ham, 0 ) + meaTriSiteLocal( ham, 1 ) )/3.0;
    double energyDiff = abs(energyOld-energy);
   
    if (istep%100==0){
      printf( "%8d%12.4e%12.8f\n", istep, energyDiff, energy );
    } else { }

    if ( energyDiff>criterion ){
      energyOld = energy;
      bondVectorOlds = bondVectors;
    }
    else {
      double energy = ( meaTriSiteLocal( ham, 0 ) + meaTriSiteLocal( ham, 1 ) )/3.0;
      printf( "%8d%12.4e%12.8f\n", istep, energyDiff, energy );
      printf( "evolution converges\n" );
      return istep;
    }
  }
  printf( "\ncannot converge...\n" );
  return -1;
}

template<typename T, template<typename T> class enType>
UniTensor<T> pess<T,enType>::ContractCore( const int iCore ){
  assert( iCore==0||iCore==1);
  //bondcat all bond vectors
  vector<UniTensor<T>> projCopy = projTensors; 
  for ( int i=0; i!=nBdsOfCore; i++ ){
    bondcat( projCopy.at(i), bondVectors.at(!iCore).at(i), (!iCore)+1 );
  }
  //Contract theta
  UniTensor<T> theta = Contract( coreTensors.at(iCore), projCopy.at(0) );
  for ( int i=1; i!=nBdsOfCore; i++){
    theta = Contract( theta, projCopy.at(i) );
  }
  return theta;
}

template<typename T, template<typename T> class enType>
double pess<T,enType>::timeEvoAstep( UniTensor<T> &evoOperator, const int iCore ){
  //act evoOperator and prepare for Hosvd
  UniTensor<T> theta = ContractCore(iCore);
  theta = Contract( theta, evoOperator );
  theta = Permute( theta, thetaLabForHosvd.at(iCore), 2*nBdsOfCore );

  //Hosvd
  vector<int> groups( nBdsOfCore, 2);
  UniTensor<T> core;
  vector<UniTensor<T>> unitary;
  Hosvd( theta, thetaLabForHosvd.at(iCore), groups, unitary, core, bondVectors.at(iCore), INPLACE );

  //truncate and update tensors;
  core = Permute( core, nBdsOfCore );
  vector<int> uniOldLab = { -1, (!iCore), iCore };
  vector<int> uniNewLab = { -1, 0, 1 };
  double trunErr = 0;
  for (int i=0; i!=nBdsOfCore; i++){
    //truncate
    truncateOneBond( core, i, dimVir );
    truncateOneBond( unitary.at(i), 2, dimVir );
    bondVectors.at(iCore).at(i) *= 1.0/Norm( bondVectors.at(iCore).at(i) );
    Resize( bondVectors.at(iCore).at(i), dimVir, dimVir, INPLACE );
    trunErr += (1.0-pow( Norm(bondVectors.at(iCore).at(i)), 2));
    //Permute
    unitary.at(i).SetLabel( uniOldLab );
    projTensors.at(i).PutBlock( Permute( unitary.at(i), uniNewLab, 1 ).GetBlock() );
    //Normalize and cat inverse lambda
    bondVectors.at(iCore).at(i) *= 1.0/Norm( bondVectors.at(iCore).at(i) );
    bondrm( projTensors.at(i), bondVectors.at(!iCore).at(i), (!iCore)+1 );
  }
  trunErr *= 1.0/nBdsOfCore;
  coreTensors.at(iCore).PutBlock( core.GetBlock() );
  coreTensors.at(iCore) *= 1.0/Norm( coreTensors.at(iCore).GetBlock() );
  return trunErr;
}

template<typename T, template<typename T> class enType>
void pess<T,enType>::Randomize( ){
  for (int i=0; i!=2; i++){
    coreTensors.at(i).Randomize();
    for (int j=0; j!=3; j++){
      Matrix<T> temp = bondVectors.at(i).at(j);
      uni10_rand( temp, uni10_mt19937, uni10_uniform_real, 0, 1, uni10_clock);
      bondVectors.at(i).at(j) = temp;
    }
  }
  for (int j=0; j!=3; j++){
    projTensors.at(j).Randomize();
  }
}

template<typename T, template<typename T> class enType>
void pess<T,enType>::saveWv( const string &folder, double field ){
  cout<<folder<<endl;
  char buffer[16];
  sprintf( buffer, "J%04i/", int(round(1000*field)) );
  string path = folder + string(buffer);
  mkdir( path.c_str(), 0755 );
  
  sprintf( buffer, "D%03i/", dimVir );
  path = path + string( buffer);
  mkdir( path.c_str(), 0755 );

  SaveVecOfTensors( coreTensors, path, "C" );
  SaveVecOfTensors( projTensors, path, "U" );
  SaveVecOfMatrices( bondVectors.at(0), path, "L0" );
  SaveVecOfMatrices( bondVectors.at(1), path, "L1" );
}

template<typename T, template<typename T> class enType>
void pess<T,enType>::loadWv( const string &folder, double field ){
  char buffer[16];
  sprintf( buffer, "J%04i/", int(round(1000*field)) );
  string path = folder + string(buffer);

  sprintf( buffer, "D%03i/", dimVir );
  path = path + string( buffer);

  LoadVecOfTensors( coreTensors, path, "C" );
  LoadVecOfTensors( projTensors, path, "U" );
  LoadVecOfMatrices( bondVectors.at(0), path, "L0" );
  LoadVecOfMatrices( bondVectors.at(1), path, "L1" );
}

template<typename T, template<typename T> class enType>
void pess<T,enType>::computeEn( int niter, double errTol, UniTensor<T>& ham, bool ifsvd ) 
{
  assert( enUp!=0 );
  assert( enDn!=0 );
  enUp->ctmIteration( niter, errTol, ham, ifsvd );
  enDn->ctmIteration( niter, errTol, ham, ifsvd );
}

template<typename T, template<typename T> class enType>
void pess<T,enType>::changeEdgeDim( int edgeDim )
{
  assert( enUp!=0 );
  assert( enDn!=0 );
  enUp->changeEdgeDim( edgeDim );
  enDn->changeEdgeDim( edgeDim );
}

template<typename T, template<typename T> class enType>
void pess<T,enType>::meaQuantities( uni10::UniTensor<T> ham, double &energyUp, double &energyDn, std::vector<std::vector<double>> &spinComp )
{
  assert( enUp!=0 );
  assert( enDn!=0 );
  enUp->measureQuantities( ham, energyUp, spinComp );
  enDn->measureQuantities( ham, energyDn, spinComp );
}

template<typename T, template<typename T> class enType>
void pess<T,enType>::initEn( int edgeDim ) 
{
  enUp = new enType<T> ( coreTensors, projTensors, edgeDim, true );
  enDn = new enType<T> ( coreTensors, projTensors, edgeDim, false );
}

template<typename T, template<typename T> class enType>
void pess<T,enType>::releaseEn() 
{
  if (enUp!=0)
  {
    delete enUp; 
    enUp=0;
  } else {}
  if (enDn!=0)
  {
    delete enDn; 
    enDn=0;
  } else {}
}

template<typename T, template<typename T> class enType>
void pess<T,enType>::saveEn( string folder, double field ) 
{
  char buffer[16];
  sprintf( buffer, "J%04i/", int(round(1000*field)) );
  string path = folder + string(buffer);
  mkdir( path.c_str(), 0755 );

  sprintf( buffer, "D%03i/", dimVir );
  path = path + string(buffer);
  mkdir( path.c_str(), 0755 );
  if (enUp!=0)
  {
    enUp->saveTensors( path + "_up" );
  } else  {}
  if (enDn!=0)
  {
    enDn->saveTensors( path + "_dn" );
  } else  {}
}

template<typename T, template<typename T> class enType>
void pess<T,enType>::loadEn( string folder, double field, int edgeDim ) 
{
  char buffer[16];
  sprintf( buffer, "J%04i/", int(round(1000*field)) );
  string path = folder + string(buffer);

  sprintf( buffer, "D%03i/", dimVir );
  path = path + string(buffer);
  if (enUp==0)
  {
    initEn( edgeDim );
  } else {};
  enUp->loadTensors( path + "_up" );
  enDn->loadTensors( path + "_dn" );
}

template class pess<double,enPess>;
template class pess<double,enRedPess>;
