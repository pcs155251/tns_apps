#include "ipeps.h"

using namespace std;
using namespace uni10;

#define if_verbose true  

///public member functions
template<typename T>
ipeps<T>::ipeps( int virDim ):
_virDim( virDim ),
gammas( vector<UniTensor<T>> (2,UniTensor<T> ( vector<Bond> { Bond(BD_IN,_phyDim), Bond(BD_OUT,virDim), Bond(BD_OUT,virDim), Bond(BD_OUT,virDim), Bond(BD_OUT,virDim) } ) ) ),
lambdas( vector<Matrix<T>> (4,Matrix<T> (virDim,virDim,true) ) )
{
  init();
}

template<typename T>
void ipeps<T>::init( )
{
  gammas[0].Randomize();
  gammas[1].Randomize();
  for (int i=0; i!=lambdas.size(); ++i )
  {
    lambdas[i].Identity();
  }
}

template<typename T>
void ipeps<T>::absrobSimpEn( int idir )
{
  int il=0, ir=1;
  for (int i=0; i!=4; ++i)
  {
    if (i!=idir)
    {
      bondsqrtcat( gammas[il], lambdas[i], i+1);
      bondsqrtcat( gammas[ir], lambdas[i], (i+2)%4+1 );
    } else {}
  }
}

template<typename T>
void ipeps<T>::releaseSimpEn( int idir )
{
  int il=0, ir=1;
  for (int i=0; i!=4; ++i)
  {
    if (i!=idir)
    {
      bondsqrtrm( gammas[il], lambdas[i], i+1);
      bondsqrtrm( gammas[ir], lambdas[i], (i+2)%4+1 );
    } else {}
  }
}


template<typename T>
UniTensor<T> ipeps<T>::contractAB( int idir )
{
  assert( 0<=idir );
  assert( idir<=3 );
  vector<int> alab = {-1,1,2,3,4};
  vector<int> blab = {-2,5,6,7,8}; 
  alab[1+idir] = 9;
  blab[1+(idir+2)%4] = 9;
  gammas[0].SetLabel( alab );
  gammas[1].SetLabel( blab );
  UniTensor<T> out = Contract( gammas[0], gammas[1] ); 
  return out; 
}

template<typename T>
UniTensor<T> ipeps<T>::actTwoSiteOp( UniTensor<T>& twoSiteOp, int idir )
{
  UniTensor<T> out = contractAB( idir );
  twoSiteOp.SetLabel( vector<int> {-3,-4,-1,-2} );
  out = Contract( out, twoSiteOp ); 
  vector<int> newLab = out.label();
  newLab.erase( newLab.begin()+6, newLab.end() );
  newLab.push_back( -1 );
  newLab.push_back( -2 );
  out.SetLabel( newLab );
  return out; 
}

template<typename T>
double ipeps<T>::simpUpdateAstep( UniTensor<T>& evolOp, int idir )
{
  absrobSimpEn( idir );
  UniTensor<T> out = actTwoSiteOp( evolOp, idir );//the bond of out: virtual bond of aten, virtual bond of bten, physical bonds
  vector<int> newLab = out.label();
  newLab.erase( newLab.begin()+6, newLab.end() );
  newLab.insert( newLab.begin()  , -1 );
  newLab.insert( newLab.begin()+4, -2 );
  out = Permute( out, newLab, 4 );

  double trunErr = svdSeperateGeneral( out, gammas[0], gammas[1], lambdas[idir], _virDim );

  vector<int> llab( newLab.begin(), newLab.begin()+4 );
  vector<int> rlab( newLab.begin()+4, newLab.end() );
  llab.push_back( 1+idir );
  rlab.insert( rlab.begin(), 1+(idir+2)%4+4 );
  gammas[0].SetLabel( llab );
  gammas[1].SetLabel( rlab );

  gammas[0] = Permute( gammas[0], vector<int> {-1,1,2,3,4}, 1);
  gammas[1] = Permute( gammas[1], vector<int> {-2,5,6,7,8}, 1);

  releaseSimpEn( idir );
  bondsqrtcat( gammas[0], lambdas[idir],  idir+1 );
  bondsqrtcat( gammas[1], lambdas[idir], (idir+2)%4+1 );
  return trunErr;
}

template<typename T>
UniTensor<T> ipeps<T>::contractSkel( int ibond, const UniTensor<T> &xten, const UniTensor<T>& yten ) const
{
  UniTensor<T> out;
  string netFile;
  if (ibond%2)
  { 
    netFile = "./fullUpdateNets/contractSkel1.net";
  }
  else
  {
    netFile = "./fullUpdateNets/contractSkel0.net";
  }
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

  UniTensor<T> xdagten = Conj( xten );
  UniTensor<T> ydagten = Conj( yten );
  Network contractSkel_net( netFile );
  uni10::ContractArgs
  (
    out, contractSkel_net,
    en->ctmBase<T>::groups[il].corners[lu0], en->ctmBase<T>::groups[il].corners[lu3], 
    en->ctmBase<T>::groups[il].edges[lu0],   en->ctmBase<T>::groups[il].edges[lu2],   en->ctmBase<T>::groups[il].edges[lu3],
    en->ctmBase<T>::groups[ir].corners[lu1], en->ctmBase<T>::groups[ir].corners[lu2], 
    en->ctmBase<T>::groups[ir].edges[lu0],   en->ctmBase<T>::groups[ir].edges[lu1],   en->ctmBase<T>::groups[ir].edges[lu2],
    xten, xdagten, yten, ydagten 
  );
  return out;
}

template<typename T>
UniTensor<T> ipeps<T>::contractT( const UniTensor<T> &skel, const UniTensor<T> &thetaPrime ) const
{
  UniTensor<T> thetaPrimeDag = Conj( thetaPrime );
  string netFile;
  netFile = "./fullUpdateNets/contractT.net";
  Network contractT_net( netFile );
  UniTensor<T> out;
  uni10::ContractArgs
  (
    out, contractT_net,
    skel, thetaPrime, thetaPrimeDag 
  );
  return out;
}

template<typename T>
UniTensor<T> ipeps<T>::contractR( const UniTensor<T> &skel, const UniTensor<T> &subTen, string type ) const
{
  UniTensor<T> subTenDag = Conj( subTen );
  string netFile;
  if (type=="bl")
  {
    netFile = "./fullUpdateNets/contractR_bl.net";
  }
  else if (type=="ar")
  {
    netFile = "./fullUpdateNets/contractR_ar.net";
  }
  else
  {
    cerr<<"undefined type for contracting R: type = "<<type<<endl;
  }
  Network contractR_net( netFile );
  UniTensor<T> out;
  uni10::ContractArgs
  (
    out, contractR_net,
    skel, subTen, subTenDag 
  );
  return out;
}

template<typename T>
UniTensor<T> ipeps<T>::contractS( const UniTensor<T> &skel, const UniTensor<T> &subTen, const UniTensor<T>& theta, string type ) const
{
  UniTensor<T> subTenDag = Conj( subTen );
  string netFile;
  if (type=="bl")
  {
    netFile = "./fullUpdateNets/contractS_bl.net";
  }
  else if (type=="ar")
  {
    netFile = "./fullUpdateNets/contractS_ar.net";
  }
  else
  {
    cerr<<"undefined type for contracting S: type = "<<type<<endl;
  }
  Network contractS_net( netFile );
  UniTensor<T> out;
  uni10::ContractArgs
  (
    out, contractS_net,
    skel, theta, subTenDag 
  );
  return out;
}

template<typename T>
UniTensor<T> ipeps<T>::contractRInv( const UniTensor<T>& rten, const UniTensor<T>& sten, string type )
{
  UniTensor<T> rinv = rten;
  {
    Matrix<T> rinvMat = Inverse( rten.GetBlock() );
    //Matrix<T> rinvMat = pseudoInverse( rten.GetBlock() );
    rinv.PutBlock( rinvMat );
  }

  string netFile;
  if (type=="bl")
  {
    netFile = "./fullUpdateNets/contractRinv_bl.net";
  }
  else if (type=="ar")
  {
    netFile = "./fullUpdateNets/contractRinv_ar.net";
  }
  else
  {
    cerr<<"undefined type for new Subten: type = "<<type<<endl;
  }

  Network contractRinv_net( netFile );
  UniTensor<T> out;
  uni10::ContractArgs
  (
    out, contractRinv_net,
    rinv, sten
  );
  return out;
}

template<typename T> 
double ipeps<T>::costFunction( UniTensor<T> subten, UniTensor<T> rten, UniTensor<T> sten, const UniTensor<T>& tten, string type )
{
  UniTensor<T> subDag = Conj( subten );
  if (type=="ar")
  {
    subten.SetLabel( vector<int> {1,0,2} );
    subDag.SetLabel( vector<int> {-1,0,-2} );
  }
  else if (type=="bl")
  {
    subten.SetLabel( vector<int> {0,1,2} );
    subDag.SetLabel( vector<int> {0,-1,-2} );
  } 
  else 
  {
    cerr<<"undefined type for input subten : type = "<<type<<endl;
  }
  rten.SetLabel( vector<int> {-1,-2,1,2} );
  UniTensor<T> aRa = Contract( subten, rten );
  aRa = Contract( aRa, subDag );
  UniTensor<T> sdag = Conj( sten );
  sten.SetLabel( vector<int> {-1,-2,0 } );
  sdag.SetLabel( vector<int> {1,2,0 } );
  UniTensor<T> aS = Contract( subDag, sten );
  UniTensor<T> Sa = Contract( sdag, subten );
  return aRa[0]-aS[0]-Sa[0]+tten[0];
}

template<typename T>
bool ipeps<T>::fullUpdateAstep( UniTensor<T>& ham, UniTensor<T>& evolOp, double &energy, int idir, int costIter, double costErr, int iterEn, double errEn, bool ifprint )
{
  UniTensor<T> xten, arten, yten, blten;
  getSubTens( idir, xten, arten, yten, blten );

  
  arten.SetLabel( vector<int> {1,-1,2} );
  blten.SetLabel( vector<int> {-2,2,3} );
  evolOp.SetLabel( vector<int> {-3,-4,-1,-2} );
  UniTensor<T> thetaP = Contract( arten, blten );
  thetaP = Contract( thetaP, evolOp );
  UniTensor<T> skel = contractSkel( idir, xten, yten );
  skel *= 1.0/abs(maxAbsElem(skel));
  double cost = 100;
  for (int i=0; i!=costIter; ++i )
  {
    UniTensor<T> tten = contractT( skel, thetaP );
    UniTensor<T> rten = contractR( skel, blten, "bl" );
    UniTensor<T> sten = contractS( skel, blten, thetaP, "bl" );
    arten = contractRInv( rten, sten, "ar" );

    tten = contractT( skel, thetaP );
    rten = contractR( skel, arten, "ar" );
    sten = contractS( skel, arten, thetaP, "ar" );
    blten = contractRInv( rten, sten, "bl" );
    double costB = costFunction( blten, rten, sten, tten, "bl" );
    if ( if_verbose)
    {
      printf( "cost %i%14.7e%6s%14.7e\n", i, costB, "diff", abs((cost-costB)/cost) );
      fflush( stdout);
    } else {}
    if ( abs((cost-costB)/cost)<costErr )
    {
      if ( if_verbose)
      {
        printf( "\n" );
        fflush( stdout);
      } else {}
      break;
    }
    else 
    {
      cost = costB;
    }
  }

  updateSubTens( idir, xten, arten, yten, blten );

  normalizeGamma( );
  en->updateUnit( getUnit() );
  return en->ctmIteration( gammas, iterEn, errEn, ham, energy, ifprint, true );
}

template<typename T>
//void herPosApproxGauge( UniTensor<T> &skel, UniTensor<T> &arten, UniTensor<T> &blten )
void ipeps<T>::herPosApproxGauge( UniTensor<T> &skel, UniTensor<T> &arten, UniTensor<T> &blten, Matrix<T> &rminv, Matrix<T> &lminv )
{
  bool iftest = false;
  UniTensor<T> test, result;
  if (iftest)
  { 
    test = skel;
  } else {}

  ///take hermitian
  Matrix<T> mat = skel.GetBlock();
  Matrix<T> dag = Dagger( mat );
  mat = 0.5*(mat+dag);
  skel.PutBlock( mat );
  ///positive approximate
  /// A = du1Dag * du0 * du1
  vector<Matrix<T>> du = EigH( mat );
  cout<<du[0]<<endl;

  /// testing scheme
  int nelem = du[0].ElemNum();
  if (du[0][0]<0.0)
  {
    if (du[0][nelem-1]<0.0)
    {
      for (int i=0; i!=nelem; ++i )
      {
        du[0][i] = sqrt( abs(du[0][i]) );
      }
    }
    else
    {
      if (abs(du[0][0])>abs(du[0][nelem-1]))
      {
        du[0] = -1.0*du[0];
        for (int i=0; i!=nelem; ++i )
        {
          if (du[0][i]>=0)
          {
            du[0][i] = sqrt( du[0][i] );
          }
          else
          {
            du[0][i] = 0;
          }
        }
      }
      else
      {
        for (int i=0; i!=nelem; ++i )
        {
          if (du[0][i]<0)
          {
            du[0][i] = 0;
          }
          else
          {
            du[0][i] = sqrt( du[0][i] );
          }
        }
      }
    }
  }
  else 
  {
    for (int i=0; i!=nelem; ++i )
    {
      du[0][i] = sqrt( du[0][i] );
    }
  }

  /*
  ///original scheme
  for (int i=0; i!=du[0].ElemNum(); ++i )
  {
    if (du[0][i]>=0)
    {
      du[0][i] = sqrt( du[0][i] );
    } 
    else 
    {
      du[0][i] = sqrt( abs(du[0][i]) );
    }
  }
  */



  UniTensor<T> zten( vector<Bond> { skel.bond(0), skel.bond(1), Bond(BD_OUT,du[0].row() ) } );
  zten.PutBlock( Dagger(du[1]) );
  bondcat( zten, du[0], 2 );

  if (iftest)
  {
    zten.SetLabel( vector<int> {1,2,0} );
    UniTensor<T> zdag = Conj( zten );
    zdag.SetLabel( vector<int> {-1,-2,0} );
    UniTensor<T> result = Contract( zten, zdag );
    cout<<setprecision(10)<<Norm( result.GetBlock() - test.GetBlock() )<<endl;
  } else {}
  
  //seperate
  UniTensor<T> rten, lten;
  UniTensor<T> temp;
  zten.SetLabel( vector<int> {1,2,0} );
  zten = Permute( zten, vector<int> {1,0,2}, 1);
  LQSeperate( zten, lten, temp );
  zten = Permute( zten, vector<int> {0,1,2}, 2);
  QRSeperate( zten, temp, rten );
  Matrix<T> rmat = rten.GetBlock();
  Matrix<T> lmat = lten.GetBlock();
  rminv = Inverse( rmat );
  lminv = Inverse( lmat );
  //rminv = pseudoInverse( rmat );
  //lminv = pseudoInverse( lmat );

  absorbMatrixtoTensor( arten, lmat, 0, true );
  absorbMatrixtoTensor( blten, rmat, 2, false );
  absorbMatrixtoTensor( zten, lminv, 1, false );
  absorbMatrixtoTensor( zten, rminv, 2, true );

  UniTensor<T> zdag = Conj( zten );
  zdag.SetLabel( vector<int> {0,-1,-2} );
  skel = Contract( zdag, zten );
}

template<typename T>
bool ipeps<T>::fullUpdateAstepGauge( UniTensor<T>& ham, UniTensor<T>& evolOp, double &energy, int idir, int costIter, double costErr, int iterEn, double errEn, bool ifprint )
{
  if ( ifprint )
  {
    cout<<"idir = "<<idir<<endl;
  } else {}

  UniTensor<T> xten, arten, yten, blten;
  getSubTens( idir, xten, arten, yten, blten );

  //gauge fixing
  UniTensor<T> skel = contractSkel( idir, xten, yten );
  skel *= 1.0/abs(maxAbsElem(skel));
  Matrix<T> rminv, lminv;
  herPosApproxGauge( skel, arten, blten, rminv, lminv );
  
  arten.SetLabel( vector<int> {1,-1,2} );
  blten.SetLabel( vector<int> {-2,2,3} );
  evolOp.SetLabel( vector<int> {-3,-4,-1,-2} );
  UniTensor<T> thetaP = Contract( arten, blten );
  thetaP = Contract( thetaP, evolOp );

  double cost = 100;
  for (int i=0; i!=costIter; ++i )
  {
    UniTensor<T> tten = contractT( skel, thetaP );
    UniTensor<T> rten = contractR( skel, blten, "bl" );
    UniTensor<T> sten = contractS( skel, blten, thetaP, "bl" );
    arten = contractRInv( rten, sten, "ar" );

    tten = contractT( skel, thetaP );
    rten = contractR( skel, arten, "ar" );
    sten = contractS( skel, arten, thetaP, "ar" );
    blten = contractRInv( rten, sten, "bl" );
    double costB = costFunction( blten, rten, sten, tten, "bl" );
    if ( if_verbose)
    {
      printf( "cost %i%14.7e%6s%14.7e\n", i, costB, "diff", abs((cost-costB)/cost) );
      fflush( stdout);
    } else {}
    if ( abs((cost-costB)/cost)<costErr )
    {
      if ( if_verbose)
      {
        printf( "\n" );
        fflush( stdout);
      } else {}
      break;
    }
    else 
    {
      cost = costB;
    }
  }
  absorbMatrixtoTensor( xten, lminv, 3, false );
  absorbMatrixtoTensor( yten, rminv, 0, true );

  updateSubTens( idir, xten, arten, yten, blten );

  normalizeGamma( );
  en->updateUnit( getUnit() );
  return en->ctmIteration( gammas, iterEn, errEn, ham, energy, ifprint, true );
}

template<typename T>
bool ipeps<T>::fullUpdate( UniTensor<T> &ham, double tau, double evoStop, int maxIter, int iterEn, double errEn, int edgeDim )
{
  computeEn( edgeDim, iterEn, errEn, ham );
  UniTensor<T> evolOp = getEvoOperator( ham, tau);
  double eOld = 0;
  int costIter = 100;
  double costErr = 1.0e-8;
  for (int i=0; i!=maxIter; ++i )
  {
    cout<<i<<endl;
    double eAver;
    bool ifCtmConverge = true;
    ifCtmConverge = ifCtmConverge && fullUpdateAstep( ham, evolOp, eAver, 0, costIter, costErr, iterEn, errEn, if_verbose );
    ifCtmConverge = ifCtmConverge && fullUpdateAstep( ham, evolOp, eAver, 2, costIter, costErr, iterEn, errEn, if_verbose );
    ifCtmConverge = ifCtmConverge && fullUpdateAstep( ham, evolOp, eAver, 1, costIter, costErr, iterEn, errEn, if_verbose );
    ifCtmConverge = ifCtmConverge && fullUpdateAstep( ham, evolOp, eAver, 3, costIter, costErr, iterEn, errEn, true );

    if ( abs( eAver - eOld) <evoStop )
    {
      return true;
    }
    else if (!ifCtmConverge)
    {
      return true;
    }
    else
    { 
      eOld = eAver;
      //cout<<i<<" "<<setprecision(8)<<eAver<<endl;
    }
  }
  return false;
}

template<typename T>
bool ipeps<T>::fullUpdateGauge( UniTensor<T> &ham, double tau, double evoStop, int maxIter, int iterEn, double errEn, int edgeDim )
{
  //measure site
  int isite = 1;
    
  computeEn( edgeDim, iterEn, errEn, ham );
  UniTensor<T> evolOp = getEvoOperator( ham, tau);
  double mzOld = 0, eOld = 0;
  int costIter = 100;
  double costErr = 1.0e-8;
  for (int i=0; i!=maxIter; ++i )
  {
    double eAver;
    bool ifCtmConverge = true;
    ifCtmConverge = ifCtmConverge && fullUpdateAstepGauge( ham, evolOp, eAver, 0, costIter, costErr, iterEn, errEn, if_verbose );
    ifCtmConverge = ifCtmConverge && fullUpdateAstepGauge( ham, evolOp, eAver, 2, costIter, costErr, iterEn, errEn, if_verbose );
    ifCtmConverge = ifCtmConverge && fullUpdateAstepGauge( ham, evolOp, eAver, 1, costIter, costErr, iterEn, errEn, if_verbose );
    ifCtmConverge = ifCtmConverge && fullUpdateAstepGauge( ham, evolOp, eAver, 3, costIter, costErr, iterEn, errEn, if_verbose );

    Matrix<double> sigmaz = 2.0*matSz();
    UniTensor<double> tenZ( sigmaz );
    T normVal = this->meaOneSiteNorm( gammas, isite );
    T mz = this->meaOneSiteExp( gammas, tenZ, isite )/normVal;
    cout<<"check "<<i<<" "<<setprecision(10)<<eAver<<" "<<setprecision(10)<<mz<<endl;

    if ( abs( eAver - eOld) <evoStop )
    {
      return true;
    }
    else if ( !ifCtmConverge )
    {
      return true;
    }
    else
    { 
      eOld = eAver;
    }

    //if ( abs( mz - mzOld) <evoStop )
    //{
    //  return true;
    //}
    //else
    //{ 
    //  mzOld = mzOld;
    //}
  }
  return false;
}

/*
template<typename T>
void ipeps<T>::test( UniTensor<T>& evolOp, int idir )
{
  UniTensor<T> xten, arten, yten, blten;
  getSubTens( idir, xten, arten, yten, blten );
  arten.SetLabel( vector<int> {1,-1,2} );
  blten.SetLabel( vector<int> {-2,2,3} );
  evolOp.SetLabel( vector<int> {-3,-4,-1,-2} );
  UniTensor<T> thetaP = Contract( arten, blten );
  thetaP = Contract( thetaP, evolOp );
  UniTensor<T> skel = contractSkel( idir, xten, yten );
  cout<<"b"<<endl;
  UniTensor<T> tten = contractT( skel, thetaP );
  cout<<"t"<<endl;
  UniTensor<T> rten = contractR( skel, blten, "bl" );
  cout<<"r"<<endl;
  UniTensor<T> sten = contractS( skel, blten, thetaP, "bl" );
  cout<<"s"<<endl;

  tten = contractT( skel, thetaP );
  cout<<"t"<<endl;
  rten = contractR( skel, arten, "ar" );
  cout<<"r"<<endl;
  sten = contractS( skel, arten, thetaP, "ar" );
  cout<<"s"<<endl;
}
*/

template<typename T>
void ipeps<T>::getSubTens( int idir, UniTensor<T> &xten, UniTensor<T>& arten, UniTensor<T>& yten, UniTensor<T>& blten ) const
{
  int posit0 = 1+idir;
  int posit1 = 1+(idir+2)%4;
  vector<int> lab0 = {-1,1,2,3,4};
  vector<int> lab1 = {-2,5,6,7,8};
  vector<UniTensor<T>> gammaTmps = gammas;
  gammaTmps[0].SetLabel( lab0 );
  gammaTmps[1].SetLabel( lab1 );
  lab0 = vector<int> {(idir+1)%4+1, (idir+2)%4+1, (idir+3)%4+1, -1, idir+1};
  lab1 = vector<int> {-2, 5+(idir+2)%4, 5+(idir+3)%4, 5+(idir)%4, 5+(idir+1)%4 };

  gammaTmps[0] = Permute( gammaTmps[0], lab0, 3 );
  gammaTmps[1] = Permute( gammaTmps[1], lab1, 2 );

  QRSeperate( gammaTmps[0], xten, arten );
  LQSeperate( gammaTmps[1], blten, yten );
}

template<typename T>
void ipeps<T>::updateSubTens( int idir, UniTensor<T> &xten, UniTensor<T>& arten, UniTensor<T>& yten, UniTensor<T>& blten )
{
  int posit0 = 1+idir;
  int posit1 = 1+(idir+2)%4;
  vector<int> lab0 = {-1,1,2,3,4};
  vector<int> lab1 = {-2,5,6,7,8};

  xten.SetLabel( vector<int> {0,1,2,3} );
  arten.SetLabel( vector<int> {3,4,5} );
  blten.SetLabel( vector<int> {0,1,2} );
  yten.SetLabel( vector<int> {2,3,4,5} );

  lab0 = vector<int> {(idir+1)%4+1, (idir+2)%4+1, (idir+3)%4+1, -1, idir+1};
  lab1 = vector<int> {-2, 5+(idir+2)%4, 5+(idir+3)%4, 5+(idir)%4, 5+(idir+1)%4 };

  gammas[0] = Contract( xten, arten );
  gammas[1] = Contract( blten, yten );
  gammas[0].SetLabel( lab0 );
  gammas[1].SetLabel( lab1 );
  gammas[0] = Permute( gammas[0], vector<int> {-1,1,2,3,4}, 1 );
  gammas[1] = Permute( gammas[1], vector<int> {-2,5,6,7,8}, 1 );
}

template<typename T>
double ipeps<T>::simpUpdateAstepReduced( UniTensor<T>& evolOp, int idir )
{
  absrobSimpEn( idir );
  UniTensor<T> xten, arten, yten, blten;
  getSubTens( idir, xten, arten, yten, blten );

  arten.SetLabel( vector<int> {1,-1,2} );
  blten.SetLabel( vector<int> {-2,2,3} );
  evolOp.SetLabel( vector<int> {-3,-4,-1,-2} );
  UniTensor<T> theta = Contract( arten, blten );
  theta = Contract( theta, evolOp );
  theta = Permute( theta, vector<int> {-3, 1, -4, 3}, 2 );
  double trunErr = svdSeperateGeneral( theta, arten, blten, lambdas[idir], _virDim );
  arten.SetLabel( vector<int> {0,1,2} );
  arten = Permute( arten, vector<int> {1,0,2}, 2 );
  blten.SetLabel( vector<int> {1,0,2} );
  blten = Permute( blten, vector<int> {0,1,2}, 2 );
  updateSubTens( idir, xten, arten, yten, blten );

  releaseSimpEn( idir );
  bondsqrtcat( gammas[0], lambdas[idir],  idir+1 );
  bondsqrtcat( gammas[1], lambdas[idir], (idir+2)%4+1 );
  return trunErr;
}

template<typename T>
void ipeps<T>::simpUpdate( UniTensor<T> &ham, double tau, double evoStop, int maxIter )
{
  UniTensor<T> evolOp = getEvoOperator( ham, tau);
  T eOld = 0;
  for (int i=0; i!=maxIter; ++i )
  {
    for (int idir=0; idir!=4; ++idir )
    {
      simpUpdateAstepReduced(evolOp,idir);
    }
    T eAver = 0;
    for (int idir=0; idir!=4; ++idir)
    {
      eAver += meaSimpTwoSiteOp( ham, idir );
    }
    eAver *= 0.25;
    if ( abs( eAver - eOld) <evoStop )
    {
      break;
    }
    else
    { 
      eOld = eAver;
      if (if_verbose)
      {
        printf( "%6d%14.7f\n", i, eAver );
        fflush( stdout);
      } else {}
    }
  }
}

template<typename T>
T ipeps<T>::meaSimpNorm( int idir )
{
  absrobSimpEn( idir );
  UniTensor<T> tmp = contractAB( idir );
  UniTensor<T> tmpDag = Conj( tmp );
  UniTensor<T> net = Contract( tmp, tmpDag );
  releaseSimpEn( idir );
  return net[0];
}

template<typename T>
T ipeps<T>::meaSimpTwoSiteOp( UniTensor<T>& twoSiteOp, int idir )
{
  assert( 0<=idir );
  assert( idir<=3 );
  absrobSimpEn( idir );
  UniTensor<T> tmp = contractAB( idir );
  UniTensor<T> tmpDag = Conj( tmp );
  UniTensor<T> norm = Contract( tmp, tmpDag );
  UniTensor<T> theta = actTwoSiteOp( twoSiteOp, idir );
  UniTensor<T> exp = Contract( theta, tmpDag );
  releaseSimpEn( idir );
  return exp[0]/norm[0];
}

template<typename T>
T ipeps<T>::meaSimpOneSiteOp( UniTensor<T>& oneSiteOp, int isite )
{
  assert( isite==0 || isite==1 );
  int idir = 0;
  absrobSimpEn( idir );

  UniTensor<T> tmp = contractAB( idir );
  UniTensor<T> tmpDag = Conj( tmp );
  UniTensor<T> norm = Contract( tmp, tmpDag );
  oneSiteOp.SetLabel( vector<int> {-1-int(isite)-2, -1-int(isite)} );
  tmp = Contract( oneSiteOp, tmp );
  vector<int> newLab = tmp.label();
  newLab[0] = -1-isite;
  tmp.SetLabel( newLab );
  UniTensor<T> exp = Contract( tmp, tmpDag );

  releaseSimpEn( idir );
  return abs( exp[0]/norm[0] );
}

/*
template<typename T>
vector<UniTensor<T>> ipeps<T>::getGammas()
{
  vector<UniTensor<T>> out = gammas;
  return out;
}
*/

template<typename T>
vector<UniTensor<T>> ipeps<T>::getUnit()
{
  vector<UniTensor<T>> gammaTmp = gammas;

  gammaTmp[0].SetLabel( vector<int> {-1,1,2,3,4} );
  UniTensor<T> gammaDag = Conj( gammaTmp[0] );
  gammaDag.SetLabel( vector<int> {-1,5,6,7,8} );
  UniTensor<T> aten = Contract( gammaTmp[0], gammaDag );
  combineTwoLayer( aten );

  gammaTmp[1].SetLabel( vector<int> {-2,5,6,7,8} );
  gammaDag = Conj( gammaTmp[1] );
  gammaDag.SetLabel( vector<int> {-2,1,2,3,4} );
  UniTensor<T> bten = Contract( gammaTmp[1], gammaDag );
  combineTwoLayer( bten );

  return vector<UniTensor<T>> { aten, bten, bten, aten };
}

template<typename T>
void ipeps<T>::initEn( int edgeDim )
{
  vector<UniTensor<double>> unitTmp = getUnit();
  en = new enIpeps<T> ( unitTmp, edgeDim );
}

template<typename T>
void ipeps<T>::changeEdgeDim( int edgeDim )
{
  if (en==0)
  {
    initEn(edgeDim);
  } 
  else if ( en->ctmDim != edgeDim )
  {
    en->changeEdgeDim(edgeDim);
  } else {}

}

template<typename T>
void ipeps<T>::computeEn( int edgeDim, int niter, double errTol, UniTensor<T>& ham )
{
  if (en==0)
  {
    initEn(edgeDim);
  } 
  else if ( en->ctmDim != edgeDim )
  {
    en->changeEdgeDim(edgeDim);
  } else {}
  double energy;
  en->ctmIteration( gammas, niter, errTol, ham, energy, true, true );
}

template<typename T>
void ipeps<T>::saveTensors( string folder, double field )
{
  char buffer[16];
  sprintf( buffer, "h%04i/", int(round(1000*field)) );
  string path = folder + string(buffer);
  mkdir( path.c_str(), 0755 );
  
  sprintf( buffer, "D%02i/", _virDim );
  path = path + string( buffer);
  mkdir( path.c_str(), 0755 );
  
  SaveVecOfTensors( gammas, path, "G" );
  SaveVecOfMatrices( lambdas, path, "L" );
}

template<typename T>
bool ipeps<T>::loadTensors( string folder, double field ) 
{
  char buffer[16];
  sprintf( buffer, "h%04i/", int(round(1000*field)) );
  string path = folder + string( buffer );

  sprintf( buffer, "D%02i/", _virDim );
  path = path + string(buffer);

  bool ifG = LoadVecOfTensors( gammas, path, "G" );
  bool ifM = LoadVecOfMatrices( lambdas, path, "L" );
  return ifG&&ifM;
}

template<typename T>
void ipeps<T>::normalizeGamma( )
{
  gammas[0] *= 1.0/maxAbsElem( gammas[0] );
  gammas[1] *= 1.0/maxAbsElem( gammas[1] );
}

template<typename T>
void ipeps<T>::saveEn( string folder, double field )
{
  if (en!=0)
  {
    char buffer[16];
    sprintf( buffer, "h%04i/", int(round(1000*field)) );
    string path = folder + string(buffer);
    mkdir( path.c_str(), 0755 );

    sprintf( buffer, "D%02i/", _virDim );
    path = path + string(buffer);
    mkdir( path.c_str(), 0755 );

    en->saveTensors( path );
  } else  {}
}

template<typename T>
void ipeps<T>::loadEn( string folder, double field, int edgeDim )
{
  if (en==0)
  {
    initEn( edgeDim );
  } else {}
  char buffer[16];
  sprintf( buffer, "h%04i/", int(round(1000*field)) );
  string path = folder + string(buffer);
  sprintf( buffer, "D%02i/", _virDim );
  path = path + string(buffer);
  en->loadTensors( path );
}

template<typename T>
ipeps<T>::~ipeps() 
{ 
  if (en!=0)
  {
    releaseEn(); 
  } else {}
}

template class ipeps<double>;
//template class ipeps<complex<double>>;
