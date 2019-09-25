#include "enIpeps.h"

using namespace std;
using namespace uni10;

#define IF_DEBUG false

template<typename T>
enIpeps<T>::enIpeps( const vector<UniTensor<T>>& unit , const int edgeDimIn ):
ctmBase<T> ( 2, 2, unit, edgeDimIn )
{
}

template<typename T>
double enIpeps<T>::meaOneSiteNorm( const vector<UniTensor<T>>& gammas, int isite ) const
{
  assert( 0<=isite );
  assert( isite<=3 );
  UniTensor<T> theta;
  if (isite==1 || isite==2 )
  {
    theta = gammas[1];
  }
  else
  {
    theta = gammas[0];
  }
  theta.SetLabel( vector<int> {-1,1,2,3,4} );
  UniTensor<T> thetaDag = Conj( theta );
  thetaDag.SetLabel( vector<int> {-1,11,12,13,14} );
  theta = Contract( theta, thetaDag );
  combineTwoLayer( theta );

  Network meaNorm_net( "./ctmIpepsNets/meaOneSiteNorm.net" );
  UniTensor<T> norm;
  uni10::ContractArgs
  (
    norm, meaNorm_net,
    ctmBase<T>::groups[isite].corners[0], ctmBase<T>::groups[isite].corners[1], 
    ctmBase<T>::groups[isite].corners[2], ctmBase<T>::groups[isite].corners[3], 
    ctmBase<T>::groups[isite].edges[0],   ctmBase<T>::groups[isite].edges[1], 
    ctmBase<T>::groups[isite].edges[2],   ctmBase<T>::groups[isite].edges[3], 
    theta 
  );
  return norm[0];
}

template<typename T>
double enIpeps<T>::meaTwoSiteNorm( const vector<UniTensor<T>>& gammas, int ibond ) const
{
  assert( 0<=ibond );
  assert( ibond<=3 );
  vector<int> alab = {-1,1,2,3,4};
  vector<int> blab = {-2,5,6,7,8}; 
  alab[1+ibond] = 9;
  blab[1+(ibond+2)%4] = 9;
  vector<UniTensor<T>> gammaTmp = gammas;
  gammaTmp[0].SetLabel( alab );
  gammaTmp[1].SetLabel( blab );

  for (int i=0; i!=4; ++i)
  {
    alab[1+i%4] = (ibond+i)%4+1;
    blab[1+(i+2)%4] = (ibond+i)%4+5;
  }
  alab[1] = 9;
  blab[1] = 9;
  gammaTmp[0] = Permute( gammaTmp[0], alab, 1 );
  gammaTmp[1] = Permute( gammaTmp[1], blab, 1 );

  UniTensor<T> theta;
  if (ibond%2)
  {
    theta = Contract( gammaTmp[1], gammaTmp[0] ); 
  }
  else
  {
    theta = Contract( gammaTmp[0], gammaTmp[1] ); 
  }

  UniTensor<T> thetaDag = Conj( theta );
  theta.SetLabel( vector<int> {-1,-11,-12,-13,-2,-14,-15,-16} );
  thetaDag.SetLabel( vector<int> {-1,11,12,13,-2,14,15,16} );
  theta = Contract( theta, thetaDag );
  combineTwoLayer( theta );

  Network meaNorm_net( "./ctmIpepsNets/meaTwoSiteNorm.net" );
  UniTensor<T> norm;

  int il, ir, lu0;
  switch( ibond )
  {
    case 0:
      il=3; ir=2; lu0=2;
      break;
    case 1:
      il=1; ir=3; lu0=1;
      break;
    case 2:
      il=0; ir=1; lu0=0;
      break;
    case 3:
      il=2; ir=0; lu0=3;
      break;
  }
  int lu1 = (lu0+1)%4;
  int lu2 = (lu0+2)%4;
  int lu3 = (lu0+3)%4;

  uni10::ContractArgs
  (
    norm, meaNorm_net,
    ctmBase<T>::groups[il].corners[lu0], ctmBase<T>::groups[il].corners[lu3], 
    ctmBase<T>::groups[il].edges[lu0]  , ctmBase<T>::groups[il].edges[lu2]  , ctmBase<T>::groups[il].edges[lu3],
    ctmBase<T>::groups[ir].corners[lu1], ctmBase<T>::groups[ir].corners[lu2], 
    ctmBase<T>::groups[ir].edges[lu0]  , ctmBase<T>::groups[ir].edges[lu1]  , ctmBase<T>::groups[ir].edges[lu2],
    theta 
  );
  return norm[0];
}

template<typename T>
double enIpeps<T>::meaOneSiteExp ( const vector<UniTensor<T>>& gammas, uni10::UniTensor<T> &oneSiteOp, int isite ) const
{
  assert( 0<=isite );
  assert( isite<=3 );
  UniTensor<T> theta;
  if (isite==1 || isite==2 )
  {
    theta = gammas[1];
  }
  else
  {
    theta = gammas[0];
  }
  theta.SetLabel( vector<int> {-1,1,2,3,4} );
  UniTensor<T> thetaDag = Conj( theta );
  oneSiteOp.SetLabel( vector<int> {-2,-1} );
  theta = Contract( oneSiteOp, theta );
  thetaDag.SetLabel( vector<int> {-2,11,12,13,14} );
  theta = Contract( theta, thetaDag );
  combineTwoLayer( theta );

  Network meaNorm_net( "./ctmIpepsNets/meaOneSiteNorm.net" );
  UniTensor<T> exp;
  uni10::ContractArgs
  (
    exp, meaNorm_net,
    ctmBase<T>::groups[isite].corners[0], ctmBase<T>::groups[isite].corners[1], 
    ctmBase<T>::groups[isite].corners[2], ctmBase<T>::groups[isite].corners[3], 
    ctmBase<T>::groups[isite].edges[0],   ctmBase<T>::groups[isite].edges[1], 
    ctmBase<T>::groups[isite].edges[2],   ctmBase<T>::groups[isite].edges[3], 
    theta 
  );
  return exp[0];
}

template<typename T>
double enIpeps<T>::meaTwoSiteExp ( const vector<UniTensor<T>>& gammas, uni10::UniTensor<T> &twoSiteOp, int ibond ) const
{
  assert( 0<=ibond );
  assert( ibond<=3 );
  vector<int> alab = {-1,1,2,3,4};
  vector<int> blab = {-2,5,6,7,8}; 
  alab[1+ibond] = 9;
  blab[1+(ibond+2)%4] = 9;
  vector<UniTensor<T>> gammaTmp = gammas;
  gammaTmp[0].SetLabel( alab );
  gammaTmp[1].SetLabel( blab );

  for (int i=0; i!=4; ++i)
  {
    alab[1+i%4] = (ibond+i)%4+1;
    blab[1+(i+2)%4] = (ibond+i)%4+5;
  }
  alab[1] = 9;
  blab[1] = 9;
  gammaTmp[0] = Permute( gammaTmp[0], alab, 1 );
  gammaTmp[1] = Permute( gammaTmp[1], blab, 1 );

  UniTensor<T> theta;
  if (ibond%2)
  { 
    theta = Contract( gammaTmp[1], gammaTmp[0] ); 
  }
  else
  {
    theta = Contract( gammaTmp[0], gammaTmp[1] ); 
  }

  UniTensor<T> thetaDag = Conj( theta );
  theta.SetLabel( vector<int> {-1,1,2,3,-2,4,5,6} );
  twoSiteOp.SetLabel( vector<int> {-3,-4,-1,-2} );
  theta = Contract( theta, twoSiteOp );//1,2,3,4,5,6;-3,-4
  thetaDag.SetLabel( vector<int> {-3,11,12,13,-4,14,15,16} );
  thetaDag = Permute( thetaDag, vector<int> {-3,-4,11,12,13,14,15,16}, 2 );
  theta = Contract( theta, thetaDag );
  combineTwoLayer( theta );

  int il, ir, lu0;
  switch( ibond )
  {
    case 0:
      il=3; ir=2; lu0=2;
      break;
    case 1:
      il=1; ir=3; lu0=1;
      break;
    case 2:
      il=0; ir=1; lu0=0;
      break;
    case 3:
      il=2; ir=0; lu0=3;
      break;
  }
  int lu1 = (lu0+1)%4;
  int lu2 = (lu0+2)%4;
  int lu3 = (lu0+3)%4;

  Network meaNorm_net( "./ctmIpepsNets/meaTwoSiteNorm.net" );
  UniTensor<T> exp;
  uni10::ContractArgs
  (
    exp, meaNorm_net,
    ctmBase<T>::groups[il].corners[lu0], ctmBase<T>::groups[il].corners[lu3], 
    ctmBase<T>::groups[il].edges[lu0],   ctmBase<T>::groups[il].edges[lu2],   ctmBase<T>::groups[il].edges[lu3],
    ctmBase<T>::groups[ir].corners[lu1], ctmBase<T>::groups[ir].corners[lu2], 
    ctmBase<T>::groups[ir].edges[lu0],   ctmBase<T>::groups[ir].edges[lu1],   ctmBase<T>::groups[ir].edges[lu2],
    theta 
  );
  return exp[0];
}

template<typename T>
bool enIpeps<T>::ctmIteration( const vector<UniTensor<T>>& gammas, int niter, double errTol, UniTensor<T> &ham, double &energy, bool ifprint, bool ifsvd )
{
  vector<double> ebondOld(4,0);
  for (int i=0; i!=niter; i++){
    ctmBase<T>::ctmrgAstep( ifsvd ); 

    vector<double> ebond(4);
    energy = 0;
    double ediff = 0;
    for (int ibond=0; ibond!=4; ++ibond)
    {
      double normVal = meaTwoSiteNorm( gammas, ibond );
      ebond[ibond] = meaTwoSiteExp( gammas, ham, ibond )/normVal;
      energy += 0.25*ebond[ibond];
      ediff += 0.25*abs( ebond[ibond] - ebondOld[ibond]);
    }

    if (ifprint)
    {
      printf( "%4i", i );
      fflush( stdout);
      for (int ibond=0; ibond!=4; ++ibond)
      {
        printf( "%17.12f", ebond[ibond] );
        fflush( stdout);
      }
      printf( "%17.12f%17.12f\n", energy, ediff );
      fflush( stdout);
    } else {}

    if (ediff<errTol)
    {
      return true;
    }
    else 
    {
      ebondOld = ebond;
    }
  }
  return false;
}

template class enIpeps<double>;
