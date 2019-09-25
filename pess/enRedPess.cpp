#include "enRedPess.h"

#define IF_DEBUG false

void fromPessToCTM( const vector<UniTensor<double>> &coreTensors, const vector<UniTensor<double>> &projTensors, UniTensor<double> &cross, UniTensor<double> &crossDag, UniTensor<double> &ket, UniTensor<double> &bra, const bool ifUp )
{
  if (ifUp){
    Network_dev cross_net( "./enRedPessNets/cross.net" );
    ContractArgs( cross, cross_net, coreTensors[0], coreTensors[1], projTensors[0], projTensors[1], projTensors[2] );//-8,3,-7,-5,-4,1,2
  } else {
    Network_dev cross_net( "./enRedPessNets/crossDn.net" );
    ContractArgs( cross, cross_net, coreTensors[0], coreTensors[1], projTensors[0], projTensors[1], projTensors[2] );//-8,3,-7,-5,-4,1,2
  }
  cross.SetLabel( vector<int> {-8,3,-7,-5,-4,1,2} );
  crossDag = Conj( cross );
  Permute( crossDag, vector<int> {-8,-7,3,-5,1,2,-4}, 7, INPLACE );
  cross *= 1.0/maxAbsElem( cross );
  crossDag *= 1.0/maxAbsElem( crossDag );
  ket = cross;
  bra = crossDag;

  cross.CombineBond( vector<int> {-8,3} );
  cross.CombineBond( vector<int> {-4,1,2} );
  crossDag.CombineBond( vector<int> {-7,3} );
  crossDag.CombineBond( vector<int> {-5,1,2} );
}

void getIdens( int dimVir, UniTensor<double> &threeIden, UniTensor<double> &fourIden, UniTensor<double> &pureIden )
{
  vector<Bond> idenBds( 2, Bond(BD_IN,dimVir) );
  idenBds.push_back( Bond(BD_IN,2));
  idenBds.push_back( Bond(BD_OUT,dimVir));
  idenBds.push_back( Bond(BD_OUT,dimVir));
  idenBds.push_back( Bond(BD_OUT,2));
  threeIden = UniTensor<double> ( idenBds );
  threeIden.Identity();
  threeIden.SetLabel( vector<int> {1,2,3,-1,-2,-3} );
  Permute( threeIden, vector<int> {1,2,-1,3,-2,-3}, 6, INPLACE );
  threeIden.CombineBond( vector<int> {-1,3} );
  threeIden.CombineBond( vector<int> {-2,-3} );

  idenBds.insert( idenBds.begin()+3, Bond(BD_IN,2) );
  idenBds.push_back( Bond(BD_OUT,2) );
  fourIden = UniTensor<double> ( idenBds );
  fourIden.Identity();
  fourIden.SetLabel( vector<int> {1,2,3,4,-1,-2,-3,-4} );
  Permute( fourIden, vector<int> {1,3,4,2,-3,-4,-1,-2}, 8, INPLACE );
  fourIden.CombineBond( vector<int> {1,3,4} );
  fourIden.CombineBond( vector<int> {2,-3,-4} );

  vector<Bond> pureIdenBds = {Bond(BD_IN,dimVir), Bond(BD_IN,dimVir), Bond(BD_OUT,dimVir), Bond(BD_OUT,dimVir)};
  pureIden = UniTensor<double> ( pureIdenBds ); 
  pureIden.Identity();
  Permute( pureIden, 4, INPLACE );
}

template<typename T>
enRedPess<T>::enRedPess(  vector<UniTensor<T>>& coreTensors, vector<uni10::UniTensor<T>>& projTensors, const int edgeDimIn, bool ifUp ):
ctmBase<T> ( 2, 2, edgeDimIn )
{
  UniTensor<T> threeIden, fourIden, pureIden;
  int dimVir = coreTensors[0].bond(0).dim();
  getIdens( dimVir, threeIden, fourIden, pureIden );

  UniTensor<T> cross, crossDag;
  UniTensor<double> bra, ket;
  fromPessToCTM( coreTensors, projTensors, cross, crossDag, ket, bra, ifUp );
  vector<UniTensor<double>> unit = {threeIden, cross, crossDag, fourIden};
  ctmBase<T>::setGroups( unit );
  auxTens = { pureIden, ket, bra };
}

/*
template<typename T>
enRedPess<T>::enRedPess( const vector<UniTensor<T>> &coresIn, const int edgeDimIn, vector<UniTensor<T>> auxTensIn ):
ctmBase<T> ( 2, 2, coresIn, edgeDimIn )
,auxTens( auxTensIn )
{
}
*/

template<typename T>
double enRedPess<T>::measureNorm( const int isite ){
  Network_dev measureNorm_net( "./enRedPessNets/measureNorm.net" );
  UniTensor<T> norm;
  ContractArgs( norm, measureNorm_net, 
                ctmBase<T>::groups.at(0).core, ctmBase<T>::groups.at(0).corners.at(0), 
                ctmBase<T>::groups.at(0).edges.at(0), ctmBase<T>::groups.at(0).edges.at(3),
                ctmBase<T>::groups.at(1).core, ctmBase<T>::groups.at(1).corners.at(1), 
                ctmBase<T>::groups.at(1).edges.at(1), ctmBase<T>::groups.at(1).edges.at(0),
                ctmBase<T>::groups.at(3).core, ctmBase<T>::groups.at(3).corners.at(2), 
                ctmBase<T>::groups.at(3).edges.at(2), ctmBase<T>::groups.at(3).edges.at(1),
                ctmBase<T>::groups.at(2).core, ctmBase<T>::groups.at(2).corners.at(3), 
                ctmBase<T>::groups.at(2).edges.at(3), ctmBase<T>::groups.at(2).edges.at(2));

  return norm[0];
}

template<typename T>
double enRedPess<T>::measureOneSite( const UniTensor<T> &oneSiteOp, const int iSite, const int posit )
{
  assert( (iSite>=0) && (iSite<=2) );
  char buffer[16];
  sprintf( buffer, "%i", iSite );
  string netString = string("./enRedPessNets/measureOneExpec") + string(buffer) + string(".net");
  Network_dev measureExpec_net( netString );
  UniTensor<T> expec;
  ContractArgs( expec, measureExpec_net, 
                auxTens[0], ctmBase<T>::groups.at(0).corners.at(0), ctmBase<T>::groups.at(0).edges.at(0), ctmBase<T>::groups.at(0).edges.at(3),
                auxTens[1] , ctmBase<T>::groups.at(1).corners.at(1), ctmBase<T>::groups.at(1).edges.at(1), ctmBase<T>::groups.at(1).edges.at(0),
                auxTens[0], ctmBase<T>::groups.at(3).corners.at(2), ctmBase<T>::groups.at(3).edges.at(2), ctmBase<T>::groups.at(3).edges.at(1),
                auxTens[2] , ctmBase<T>::groups.at(2).corners.at(3), ctmBase<T>::groups.at(2).edges.at(3), ctmBase<T>::groups.at(2).edges.at(2),
                oneSiteOp);
  return expec[0];
}

template<typename T>
double enRedPess<T>::measureTriSite( const UniTensor<T> &threeSiteOp, const int posit )
{
  Network_dev measureExpec_net( "./enRedPessNets/measureTriExpec.net" );
  UniTensor<T> expec;
  ContractArgs( expec, measureExpec_net, 
                auxTens[0], ctmBase<T>::groups.at(0).corners.at(0), ctmBase<T>::groups.at(0).edges.at(0), ctmBase<T>::groups.at(0).edges.at(3),
                auxTens[1] , ctmBase<T>::groups.at(1).corners.at(1), ctmBase<T>::groups.at(1).edges.at(1), ctmBase<T>::groups.at(1).edges.at(0),
                auxTens[0], ctmBase<T>::groups.at(3).corners.at(2), ctmBase<T>::groups.at(3).edges.at(2), ctmBase<T>::groups.at(3).edges.at(1),
                auxTens[2] , ctmBase<T>::groups.at(2).corners.at(3), ctmBase<T>::groups.at(2).edges.at(3), ctmBase<T>::groups.at(2).edges.at(2),
                threeSiteOp);
  return expec[0];
}

template<typename T>
double enRedPess<T>::measureQuantities( const UniTensor<T> &ham, double &energy, vector<vector<double>> &spinComp  )
{
  //isite dummy parameter
  int isite = 0;
  cout<<"measureQuantities"<<endl<<flush;
  double normVal = measureNorm( isite );
  energy = measureTriSite( ham, isite )/normVal;
  spinComp = vector<vector<double>> ( 3, vector<double> (2,0) );
  UniTensor<double> spinX( matSx(0.5) );
  UniTensor<double> spinZ( matSz(0.5) );
  vector<double> staggerM(3);
  for (int i=0; i!=3; i++){
    spinComp[i][0] = measureOneSite( spinX, i, isite )/normVal;
    spinComp[i][1] = measureOneSite( spinZ, i, isite )/normVal;

    staggerM[i] = sqrt( spinComp[i][0]*spinComp[i][0] + spinComp[i][1]*spinComp[i][1] );
  }
  return (2.0/3.0)*energy;
}

template<typename T>
double enRedPess<T>::ctmIteration( const int niter, const double errTol, const UniTensor<T> &upHam, bool ifsvd )
{
  if (ifsvd)
  {
    cout<<endl<<"(SVD) dimCTM=";
  }
  else
  {
    cout<<endl<<"(QR) dimCTM=";
  }
  cout<<ctmBase<T>::ctmDim<<endl<<flush;
  //dummy parameter isite
  int isite = 0;
  double energyOld = 0;
  for (int i=0; i!=niter; i++){
    ctmBase<T>::ctmrgAstep( ifsvd ); 
    double normVal = measureNorm( isite );
    double energy = (2.0/3.0)*measureTriSite( upHam, isite )/normVal;
    cout<<i<<" "<<setprecision(10)<<energy<<" "<<normVal<<flush<<endl;
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
void enRedPess<T>::changeEdgeDim( const int newDim )
{
  ctmBase<T>::changeEdgeDim( newDim );
}

template class enRedPess<double>;
