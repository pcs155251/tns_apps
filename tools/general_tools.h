#ifndef __GENERAL_TOOLS_H__
#define __GENERAL_TOOLS_H__
#include "uni10.hpp"
#include <stdio.h>

using namespace std;
using namespace uni10;

template<typename T>
vector<Matrix<T>> QrSdd( const Matrix<T> &Amat );

template<typename T1, typename T2>
void absorbMatrixtoTensor( UniTensor<T1> &ten, const Matrix<T2> &mat, const int posit, const bool tensorFirst );

template<typename T>
Matrix<T> pseudoInverse( const Matrix<T> &mat );

template<typename T>
T maxAbsElem(const UniTensor<T>& Tin);

template<typename T>
double maxAbsElem(const Matrix<T>& Tin);

template<typename T>
vector<Matrix<T>> rsvd( const Matrix<T> &amat, const int nRankIn, double &trunRsvd );

template<typename T>
void truncateOneBond( UniTensor<T> &ten, const int posit, const int newDim );

//! \ingroup aux
//! Change one bond dimension of a tensor.
/*! 
 *  Truncate original elements/ filled with zeros if new bond dimension is smaller/ larger than original one.
 *  \param[in] ten Tensor being changed.
 *  \param[in] posit Bond being changed.
 *  \param[in] newDim New bond dimension.
 */
template<typename T>
void changeOneBond( UniTensor<T> &ten, const int posit, const int newDim );

//! \ingroup aux
//! Get imaginary time evolution operator from given two site hamiltonain.
/*! 
 *  \param[in] hamiltonian 
 *  \param[in] tau Imaginary time step.
 */
UniTensor<double> getEvoOperator( const UniTensor<double> &hamiltonian, const double tau );

UniTensor<complex<double>> getEvoOperator( const UniTensor<complex<double>> &hamiltonian, const double tau );

template<typename T>
void seperBond( UniTensor<T> &tensor, const int posit, const vector<int> newDim );

template<typename T>
void SaveVecOfTensors( const vector<UniTensor<T>> &tensors, const string &SaveFolder, const string &prefix );

template<typename T>
void SaveVecOfMatrices( const vector<Matrix<T>> &matrices, const string &SaveFolder, const string &prefix );

template<typename T>
bool LoadVecOfTensors( vector<UniTensor<T>> &tensors, const string &LoadFolder, const string &prefix );

template<typename T>
bool LoadVecOfMatrices( vector<Matrix<T>> &matrices, const string &LoadFolder, const string &prefix );

template<typename T>
double svdSeperate( const UniTensor<T> &tin, vector<UniTensor<T>> &result, Matrix<T> &singVal, int chi );

template<typename T>
double svdSeperateNoNormalize( const UniTensor<T> &tin, vector<UniTensor<T>> &result, Matrix<T> &singVal, int chi, const double eps=0 );
//seperate with svd into two tensor u*sqrt(s), sqrt(s)*vT, with repect to inbond and outbond return singular value, doesn't do rescale
//may set eps, if not will use chi directly
//return truncation error

template<typename T>
double svdSeperateNoLambda( const UniTensor<T> &tin, UniTensor<T> &resultL, UniTensor<T> &resultR, Matrix<T> &singVal, int chi );

template<typename T>
double svdSeperateSymNoLambda( const UniTensor<T> &tin, UniTensor<T> &resultL, UniTensor<T> &resultR, Matrix<T> &singVal, int chi );

template<typename T>
void QRSeperate( const UniTensor<T> &tin, UniTensor<T> &qten, UniTensor<T> &rten );

template<typename T>
void LQSeperate( const UniTensor<T> &tin, UniTensor<T> &lten, UniTensor<T> &qten );

template<typename T>
double svdSeperateGeneral( const UniTensor<T> &tin, UniTensor<T> &resultL, UniTensor<T> &resultR, Matrix<T> &singVal, int chi );

template<typename T>
double tebdHastings( const UniTensor<T> &theta, UniTensor<T> &tenL, UniTensor<T> &tenR, Matrix<T> &matL, Matrix<T> &matC, const int chi );

template<typename T1, typename T2>
void bondcat(uni10::UniTensor<T1>& Tout, const uni10::Matrix<T2>& L, uni10_int32 bidx);

template<typename T1, typename T2>
void bondrm(uni10::UniTensor<T1>& Tout, const uni10::Matrix<T2>& L, uni10_int32 bidx);

template<typename T>
uni10::UniTensor<T> ContractSelf( const uni10::UniTensor<T> &tin, const std::vector<int> &posits);

Matrix<double> matId( const int dim );

Matrix<uni10_complex128> cMatrix( const Matrix<double> &test );

Matrix<uni10_double64> EignLeft( const Matrix<uni10_double64> &test, uni10_double64 &lambda, const double cutoff );

Matrix<uni10_double64> EignRigh( const Matrix<uni10_double64> &test, uni10_double64 &lambda, const double cutoff );

template<typename T>
T maxAbsElem(const UniTensor<T>& Tin){
  int nelem = Tin.ElemNum();
  T maxAbsVal = Tin[0];
  for (int i=0; i!=nelem; i++ ){
    if ( fabs(maxAbsVal) < fabs(Tin[i]) ){
      maxAbsVal = Tin[i];
    } else {};
  }
  return maxAbsVal;
}

template<typename T>
double maxAbsElem(const Matrix<T>& Tin){
  int nelem = Tin.ElemNum();
  double maxAbsVal=0;
  for (int i=0; i!=nelem; i++ ){
    maxAbsVal = max( maxAbsVal, fabs(Tin[i]) );
  }
  return maxAbsVal;
}

template<typename T>
Matrix<T> sqrtDiagMat( const Matrix<T>& L ){
  Matrix<T> sqrtL = L;
  for(uni10_uint64 i=0; i!=L.col(); i++){
    sqrtL[i] = sqrt( sqrtL[i] );
  }
  return sqrtL;
}

template<typename T1, typename T2>
void bondcat(UniTensor<T1>& Tout, const Matrix<T2>& L, uni10_int32 bidx){

  uni10_int32 InBondNum = Tout.InBondNum();
  vector<uni10_int32> labels = Tout.label();
  vector<uni10_int32> per_labels = labels;
  uni10_int32 l = labels[bidx];
  per_labels.erase(per_labels.begin() + bidx);
  per_labels.insert(per_labels.begin(), l);

  UniTensor<T1> T_c = Permute(Tout, per_labels, 1);

  T_c.PutBlock((Dot(L, T_c.GetBlock())));
  //debug
  //Matrix<T1> result = Dot(L,T_c.GetBlock());
  //T_c.PutBlock(result);

  Tout = Permute( T_c, labels, InBondNum);
}

template<typename T1, typename T2>
void bondsqrtcat(UniTensor<T1>& Tout, const Matrix<T2>& L, uni10_int32 bidx){
  Matrix<T2> sqrtL = sqrtDiagMat( L );
  bondcat( Tout, sqrtL, bidx);
}

template<typename T1, typename T2>
void bondrm(UniTensor<T1>& Tout, const Matrix<T2>& L, uni10_int32 bidx){

  Matrix<T2> invL = L;
  for(uni10_uint64 i=0; i!=L.col(); i++){
    invL[i] = invL[i] == 0.0 ? 0.0 : ( 1.0 / invL[i]);
  }
  bondcat(Tout, invL, bidx);
}

template<typename T1, typename T2>
void bondsqrtrm(UniTensor<T1>& Tout, const Matrix<T2>& L, uni10_int32 bidx){

  Matrix<T2> invL = L;
  for(uni10_uint64 i=0; i!=L.col(); i++){
    invL[i] = invL[i] == 0.0 ? 0.0 : ( 1.0 / sqrt(invL[i]) );
  }
  bondcat(Tout, invL, bidx);
}

template<typename T>
void mergeMatBeforeTen( UniTensor<T> &tenIn, const Matrix<T> &matIn, const int posit ){
  int inBdN = tenIn.InBondNum();
  vector<int> oldLab = tenIn.label();
  vector<int> newLab = oldLab;
  int moveLab = newLab.at( posit );
  newLab.erase( newLab.begin()+posit );
  newLab.insert( newLab.begin(), moveLab );
  tenIn = Permute( tenIn, newLab, 1 );
  Matrix<T> matOut = Dot( matIn, tenIn.GetBlock() );
  tenIn.PutBlock( matOut );
  tenIn = Permute( tenIn, oldLab, inBdN );
}

template<typename T>
void mergeMatAfterTen( UniTensor<T> &tenIn, const Matrix<T> &matIn, const int posit ){
  int bdN = tenIn.BondNum();
  int inBdN = tenIn.InBondNum();
  vector<int> oldLab = tenIn.label();
  vector<int> newLab = oldLab;
  int moveLab = newLab.at( posit );
  newLab.erase( newLab.begin()+posit );
  newLab.push_back( moveLab );
  tenIn = Permute( tenIn, newLab, bdN-1 );
  Matrix<T> matOut = Dot( tenIn.GetBlock(), matIn );
  tenIn.PutBlock( matOut );
  tenIn = Permute( tenIn, oldLab, inBdN );
}

template<typename T>
UniTensor<T> ContractSelf( UniTensor<T> &tin, const vector<int> &posits){
  vector<int> tcopyLab = tin.label();
  int max=*max_element( tcopyLab.begin(), tcopyLab.end());
  for ( int i=0; i!=tcopyLab.size(); i++){
    tcopyLab.at(i) = max+i+1;
  }
  for ( int i=0; i!=posits.size(); i++){
    tcopyLab.at( posits.at(i) ) = tin.label().at( posits.at(i) );
  }
  //UniTensor<T> tcopy = Dagger( tin );
  //tcopy = Transpose( tin );
  UniTensor<T> tcopy = tin;
  tcopy.SetLabel( tcopyLab );
  UniTensor<T> tout;
  tout = Contract( tin, tcopy, false );
  return tout;
}

template<typename T>
UniTensor<T> ContractDagger( UniTensor<T> &tin, const vector<int> &posits){
  assert( posits.size() <= tin.label().size() );
  vector<int> tcopyLab = tin.label();
  int max=*max_element( tcopyLab.begin(), tcopyLab.end());
  for ( int i=0; i!=tcopyLab.size(); i++){
    tcopyLab.at(i) = max+i+1;
  }
  for ( int i=0; i!=posits.size(); i++){
    tcopyLab.at( posits.at(i) ) = tin.label().at( posits.at(i) );
  }
  UniTensor<T> tcopy = Conj( tin );
  tcopy.SetLabel( tcopyLab );
  UniTensor<T> tout;
  //note that tin is first input here, in my convention tin should be ket state tensors
  tout = Contract( tin, tcopy, false );
  return tout;
}

template<typename T>
void combineTwoLayer( UniTensor<T> &tensor ){
  const int bdN = tensor.BondNum();
  assert( bdN%2==0 );
  int newbdN = bdN/2;
  vector<int> oldLab = tensor.label();
  vector<int> combineLab(2);
  for (int i=0; i!=newbdN; i++){
    combineLab.at(0) = oldLab.at(i);
    combineLab.at(1) = oldLab.at(newbdN+i);
    tensor.CombineBond( combineLab );
  }
}

template<typename T>
void truncateOneBond( UniTensor<T> &tin, const int posit, const int newDim ){
  int oldDim = tin.bond(posit).dim();
  assert( newDim <= oldDim );
  //Permute tin
  int inBdN = tin.InBondNum();
  vector<int> oldLab = tin.label();
  vector<int> newLab = oldLab;
  int moveLab = newLab.at( posit );
  newLab.erase( newLab.begin()+posit );
  newLab.insert( newLab.begin(), moveLab );
  tin = Permute( tin, newLab, 1 );
  //get matrix and Resize
  Matrix<T> matOut = tin.GetBlock();
  Resize( matOut, newDim, matOut.col(), INPLACE );
  //new tin and put matrix
  vector<Bond> newBds = tin.bond();
  newBds.erase( newBds.begin() );
  newBds.insert( newBds.begin(), Bond( BD_IN, newDim ) );
  tin = UniTensor<T> (newBds);
  tin.SetLabel( newLab );
  tin.PutBlock( matOut );
  tin = Permute( tin, oldLab, inBdN );
}

template<typename T>
void changeOneBond( UniTensor<T> &ten, const int posit, const int newDim )
{
  int oldDim = ten.bond(posit).dim();
  if (newDim==oldDim){
  }
  else {
    //Permute ten
    int inBdN = ten.InBondNum();
    vector<int> oldLab = ten.label();
    vector<int> newLab = oldLab;
    int moveLab = newLab.at( posit );
    newLab.erase( newLab.begin()+posit );
    newLab.insert( newLab.begin(), moveLab );
    ten = Permute( ten, newLab, 1 );
    //get matrix and Resize
    Matrix<T> matOut = ten.GetBlock();
    Resize( matOut, newDim, matOut.col(), INPLACE );
    //new ten and put matrix
    vector<Bond> newBds = ten.bond();
    newBds.erase( newBds.begin() );
    newBds.insert( newBds.begin(), Bond( BD_IN, newDim ) );
    ten = UniTensor<T> (newBds);
    ten.SetLabel( newLab );
    ten.PutBlock( matOut );
    ten = Permute( ten, oldLab, inBdN );
  }
}

template<typename type>
void convertVecToMat( Matrix<type> &EignVector, const int firstDim ){
  //this function reshape a row vector or a column vector into a matrix with row number = firstDim
  assert( (EignVector.col()==1 )||( EignVector.row()==1) );
  if (EignVector.row()==1){
    //left EignVector
    int totalDim = EignVector.col();
    assert( totalDim%firstDim == 0 );
    int secondDim = totalDim/firstDim;
    Matrix<type> copy = EignVector;
    Resize( EignVector, firstDim, secondDim, INPLACE );
    EignVector.SetElem( copy.GetElem() );
  }
  else if (EignVector.col()==1){
    //right EignVector
    int totalDim = EignVector.row();
    assert( totalDim%firstDim == 0 );
    int secondDim = totalDim/firstDim;
    Matrix<type> copy = EignVector;
    Resize( EignVector, firstDim, secondDim, INPLACE );
    EignVector.SetElem( copy.GetElem() );
  }
}

template<typename T>
void SaveVecOfTensors( const vector<UniTensor<T>> &tensors, const string &SaveFolder, const string &prefix ){
  for ( int i=0; i<tensors.size(); i++ ){
    char buffer[32];
    sprintf( buffer, "%s%d.ten", prefix.c_str(), i );
    string SavePath = SaveFolder + string(buffer);
    tensors.at(i).Save(SavePath);
  }
}

template<typename T>
void SaveVecOfMatrices( const vector<Matrix<T>> &matrices, const string &SaveFolder, const string &prefix ){
  for ( int i=0; i<matrices.size(); i++ ){
    char buffer[32];
    sprintf( buffer, "%s%d.mat", prefix.c_str(), i );
    string SavePath = SaveFolder + string(buffer);
    matrices.at(i).Save(SavePath);
  }
}

template<typename T>
void removeVecOfTensors( const vector<UniTensor<T>> &tensors, const string &SaveFolder, const string &prefix ){
  for ( int i=0; i<tensors.size(); i++ ){
    char buffer[32];
    sprintf( buffer, "%s%d.ten", prefix.c_str(), i );
    string SavePath = SaveFolder + string(buffer);
    remove(SavePath.c_str());
  }
}

template<typename T>
bool LoadVecOfTensors( vector<UniTensor<T>> &tensors, const string &LoadFolder, const string &prefix ){
  struct stat check;   
  for ( int i=0; i<tensors.size(); i++ ){
    char buffer[32];
    sprintf( buffer, "%s%d.ten", prefix.c_str(), i );
    string LoadPath = LoadFolder + string(buffer);
    if (stat(LoadPath.c_str(),&check)==0){
      cout<<"Load "<<LoadPath<<endl;
      tensors.at(i).Load(LoadPath);
    }
    else {
      cout<<"didn't Load "<<LoadPath<<endl;
      return false;
    }
  }
  return true;
}

template<typename T>
bool LoadVecOfMatrices( vector<Matrix<T>> &matrices, const string &LoadFolder, const string &prefix ){
  struct stat check;   
  for ( int i=0; i<matrices.size(); i++ ){
    char buffer[32];
    sprintf( buffer, "%s%d.mat", prefix.c_str(), i );
    string LoadPath = LoadFolder + string(buffer);
    if (stat(LoadPath.c_str(),&check)==0){
      cout<<"Load "<<LoadPath<<endl;
      matrices.at(i).Load(LoadPath);
    }
    else {
      cout<<"didn't Load "<<LoadPath<<endl;
      return false;
    }
  }
  return true;
}

/*
template<typename T>
void seperBond( UniTensor<T> &tensor, const int posit, const int newDim1 ){
  int oldDim = tensor.bond(posit).dim();
  assert( oldDim%newDim1 == 0 );
  assert( posit>=0&&posit<tensor.bond().size() );
  int newDim2 = oldDim/newDim1;
  vector<Bond> oldBds = tensor.bond();
  Bond newBd1( oldBds.at(posit).type(), newDim1 );
  Bond newBd2( oldBds.at(posit).type(), newDim2 );
  oldBds.erase( oldBds.begin()+posit );
  oldBds.insert( oldBds.begin()+posit, newBd1 );
  oldBds.insert( oldBds.begin()+posit, newBd2 );
  UniTensor<T> outTen( oldBds );
  outTen.PutBlock( tensor.GetBlock() );
  tensor = outTen;
}
*/

template<typename T>
void seperBond( UniTensor<T> &tensor, const int posit, const vector<int> newDim ){
  //check position and bond dimension
  assert( posit>=0&&posit<tensor.bond().size() );
  int oldDim = tensor.bond(posit).dim();
  int factor = 1;
  for (int i=0; i!=newDim.size(); i++){
    factor *= newDim[i];
  }
  assert( factor == oldDim );
  //declare new tensor
  vector<Bond> oldBds = tensor.bond();
  Bond oldBond = oldBds.at(posit);
  oldBds.erase( oldBds.begin()+posit );
  for (int i=0; i!=newDim.size(); i++){
    oldBds.insert( oldBds.begin()+posit+i, Bond( oldBond.type(), newDim[i] ));
  }
  UniTensor<T> outTen( oldBds );
  outTen.PutBlock( tensor.GetBlock() );
  tensor = outTen;
}


template<typename T>
double svdSeperate( const UniTensor<T> &tin, vector<UniTensor<T>> &result, Matrix<T> &singVal, int chi ){
  //svd, truncate, Normalize singular matrix
  vector<Matrix<T> > usv = Sdd(tin.GetBlock());

  if (chi>usv[1].col()){
    chi=usv[1].col();
  } else {}

  double trunErr = 0;
  for ( int i=chi; i!=usv[1].col(); i++ ){
    trunErr += pow( fabs(usv[1][i]), 2);
  }
  trunErr = trunErr/pow( Norm(usv[1]), 2);

  Resize( usv[0], usv[0].row(),          chi, INPLACE);
  Resize( usv[1],          chi,          chi, INPLACE);
  Resize( usv[2],          chi, usv[2].col(), INPLACE);
  usv[1] *= 1.0/Norm( usv[1] );
  
  //allocate result tensors
  result.clear();
  int bdN = tin.BondNum();
  int inBdN = tin.InBondNum();
  vector<Bond> tBonds = tin.bond();

  vector<Bond> usBonds( tBonds.begin(), tBonds.begin()+inBdN );
  usBonds.push_back( Bond( BD_OUT, chi) );
  result.push_back( UniTensor<T> (usBonds) );

  vector<Bond> svBonds( tBonds.begin()+inBdN, tBonds.end());
  svBonds.insert( svBonds.begin(), Bond( BD_IN, chi) );
  result.push_back( UniTensor<T> (svBonds) );


  //assign results
  singVal = usv[1];

  for (int i=0; i!=chi; i++){
    usv[1][i] = sqrt( usv[1][i] );
  }
  result[0].PutBlock( usv[0] );
  bondcat( result[0], usv[1], inBdN );
  result[1].PutBlock( usv[2] );
  bondcat( result[1], usv[1], 0 );

  return trunErr;
}

template<typename T>
double svdSeperateNoNormalize( const UniTensor<T> &tin, vector<UniTensor<T>> &result, Matrix<T> &singVal, int chi, const double eps ){
  //svd, truncate, return two tensors u*sqrt(s), sqrt(s)*vT, without rescale s
  vector<Matrix<T> > usv = Sdd(tin.GetBlock());

  double trunErr = 0;
  double normSqr = pow( Norm( usv[1] ), 2 );
  if (eps==0){
    if (chi>usv[1].col()){
      chi=usv[1].col();
    } else {}
    for ( int i=chi; i!=usv[1].col(); i++ ){
      trunErr += pow( fabs(usv[1][i]), 2);
    }
    trunErr = trunErr/normSqr;
  } 
  else {
    //decide chi by eps
    for ( int i=usv[1].col()-1; i!=-1; i-- ){
      trunErr += pow( fabs(usv[1][i]), 2)/normSqr;
      if (trunErr>eps){
        trunErr -= pow( fabs(usv[1][i]), 2)/normSqr;
        chi = i+1;
        break;
      } else {}
    }
  }

  Resize( usv[0], usv[0].row(),          chi, INPLACE);
  Resize( usv[1],          chi,          chi, INPLACE);
  Resize( usv[2],          chi, usv[2].col(), INPLACE);
  
  //allocate result tensors
  result.clear();
  int bdN = tin.BondNum();
  int inBdN = tin.InBondNum();
  vector<Bond> tBonds = tin.bond();

  vector<Bond> usBonds( tBonds.begin(), tBonds.begin()+inBdN );
  usBonds.push_back( Bond( BD_OUT, chi) );
  result.push_back( UniTensor<T> (usBonds) );

  vector<Bond> svBonds( tBonds.begin()+inBdN, tBonds.end());
  svBonds.insert( svBonds.begin(), Bond( BD_IN, chi) );
  result.push_back( UniTensor<T> (svBonds) );


  //assign results
  singVal = usv[1];

  for (int i=0; i!=chi; i++){
    usv[1][i] = sqrt( usv[1][i] );
  }
  result[0].PutBlock( usv[0] );
  bondcat( result[0], usv[1], inBdN );
  result[1].PutBlock( usv[2] );
  bondcat( result[1], usv[1], 0 );

  return trunErr;
}

template<typename T>
double svdSeperateNoLambda( const UniTensor<T> &tin, UniTensor<T> &resultL, UniTensor<T> &resultR, Matrix<T> &singVal, int chi ){
  //svd, truncate, Normalize singular matrix
  //vector<Matrix<T> > usv = Svd(tin.GetBlock());
  vector<Matrix<T> > usv = QrSdd(tin.GetBlock());

  if (chi>usv[1].col()){
    chi=usv[1].col();
  } else {}

  double trunErr = 0;
  for ( int i=chi; i!=usv[1].col(); i++ ){
    trunErr += pow( fabs(usv[1][i]), 2);
  }
  double sumSquSing = Norm( usv[1] );
  trunErr = trunErr/(sumSquSing*sumSquSing);


  //dynamical truncation
  /*
  int chiPrime = chi;
  for ( int i=0; i!=chi; i++){
    if ( fabs(usv[1].At(i,i)/usv[1].At(0,0))<1.0e-8){
      chiPrime = i;
      //cout<<chiPrime<<endl;
      break;
    } else {}
  }
  //Resize( usv[0], usv[0].row(),     chiPrime, INPLACE);
  Resize( usv[1],     chiPrime,     chiPrime, INPLACE);
  //Resize( usv[2],     chiPrime, usv[2].col(), INPLACE);
  */

  Resize( usv[0], usv[0].row(),          chi, INPLACE);
  Resize( usv[1],          chi,          chi, INPLACE);
  Resize( usv[2],          chi, usv[2].col(), INPLACE);

  usv[1] *= 1.0/Norm( usv[1] );
  
  //allocate result tensors
  int bdN = tin.BondNum();
  int inBdN = tin.InBondNum();
  vector<Bond> tBonds = tin.bond();

  vector<Bond> usBonds( tBonds.begin(), tBonds.begin()+inBdN );
  usBonds.push_back( Bond( BD_OUT, chi) );
  resultL = UniTensor<T> (usBonds);

  vector<Bond> svBonds( tBonds.begin()+inBdN, tBonds.end());
  svBonds.insert( svBonds.begin(), Bond( BD_IN, chi) );
  resultR = UniTensor<T> (svBonds);


  //assign results
  singVal = usv[1];

  resultL.PutBlock( usv[0] );
  resultR.PutBlock( usv[2] );

  return trunErr;
}


template<typename T>
double randSvdSeperateNoLambda( const UniTensor<T> &tin, UniTensor<T> &resultL, UniTensor<T> &resultR, Matrix<T> &singVal, int chi ){
  //svd, truncate, Normalize singular matrix
  Matrix<T> amat = tin.GetBlock();
  complex<double> factor = Trace( Dot(amat, Dagger(amat)) );
  assert(factor.imag()<1.0e-14);
  amat *= 1.0/sqrt(factor.real());
  double trunErr = 0;
  vector<Matrix<T> > usv = rsvd(amat, chi, trunErr);
  usv[1] *= 1.0/Norm( usv[1] );
  
  //allocate result tensors
  int bdN = tin.BondNum();
  int inBdN = tin.InBondNum();
  vector<Bond> tBonds = tin.bond();

  vector<Bond> usBonds( tBonds.begin(), tBonds.begin()+inBdN );
  usBonds.push_back( Bond( BD_OUT, chi) );
  resultL = UniTensor<T> (usBonds);

  vector<Bond> svBonds( tBonds.begin()+inBdN, tBonds.end());
  svBonds.insert( svBonds.begin(), Bond( BD_IN, chi) );
  resultR = UniTensor<T> (svBonds);


  //assign results
  singVal = usv[1];

  resultL.PutBlock( usv[0] );
  resultR.PutBlock( usv[2] );

  return trunErr;
}

template<typename T>
double svdSeperateSymNoLambda( const UniTensor<T> &tin, UniTensor<T> &resultL, UniTensor<T> &resultR, Matrix<T> &singVal, int chi, double cutDiff ){
  //svd, truncate, Normalize singular matrix
  vector<Matrix<T> > usv = Svd(tin.GetBlock());

  if (chi>usv[1].col()){
    chi=usv[1].col();
  } else {}

  T lastS = usv[1][chi-1];
  for ( int i=chi; i!=usv[1].col(); i++ ){
    if (fabs(lastS-usv[1][i])<cutDiff){
      chi++;
    }
    else {
      break;
    }
  }

  double trunErr = 0;
  for ( int i=chi; i!=usv[1].col(); i++ ){
    trunErr += pow( fabs(usv[1][i]), 2);
  }
  double sumSquSing = Norm( usv[1] );
  trunErr = trunErr/(sumSquSing*sumSquSing);

  Resize( usv[0], usv[0].row(),          chi, INPLACE);
  Resize( usv[1],          chi,          chi, INPLACE);
  Resize( usv[2],          chi, usv[2].col(), INPLACE);
  usv[1] *= 1.0/Norm( usv[1] );
  
  //allocate result tensors
  int bdN = tin.BondNum();
  int inBdN = tin.InBondNum();
  vector<Bond> tBonds = tin.bond();

  vector<Bond> usBonds( tBonds.begin(), tBonds.begin()+inBdN );
  usBonds.push_back( Bond( BD_OUT, chi) );
  resultL = UniTensor<T> (usBonds);

  vector<Bond> svBonds( tBonds.begin()+inBdN, tBonds.end());
  svBonds.insert( svBonds.begin(), Bond( BD_IN, chi) );
  resultR = UniTensor<T> (svBonds);


  //assign results
  singVal = usv[1];

  resultL.PutBlock( usv[0] );
  resultR.PutBlock( usv[2] );

  return trunErr;
}

template<typename T>
void QRSeperate( const UniTensor<T> &tin, UniTensor<T> &qten, UniTensor<T> &rten ){
  vector<Matrix<T>> qrResult = Qr( tin.GetBlock() );
  int chi = qrResult[0].col();

  //allocate result tensors
  int bdN = tin.BondNum();
  int inBdN = tin.InBondNum();
  vector<Bond> tBonds = tin.bond();

  vector<int> labs(tin.bond().size()+1);
  for (unsigned i=0; i!=labs.size(); ++i )
  {
    labs[i] = i;
  }

  vector<Bond> usBonds( tBonds.begin(), tBonds.begin()+inBdN );
  usBonds.push_back( Bond( BD_OUT, chi) );
  qten = UniTensor<T> (usBonds);
  qten.PutBlock( qrResult[0] );
  qten.SetLabel( vector<int> (labs.begin(), labs.begin()+inBdN+1 ) );

  vector<Bond> svBonds( tBonds.begin()+inBdN, tBonds.end());
  svBonds.insert( svBonds.begin(), Bond( BD_IN, chi) );
  rten = UniTensor<T> (svBonds);
  rten.PutBlock( qrResult[1] );
  rten.SetLabel( vector<int> (labs.begin()+inBdN, labs.end() ) );
}

template<typename T>
void LQSeperate( const UniTensor<T> &tin, UniTensor<T> &lten, UniTensor<T> &qten ){
  vector<Matrix<T>> lqResult = Lq( tin.GetBlock() );
  int chi = lqResult[0].col();

  //allocate result tensors
  int bdN = tin.BondNum();
  int inBdN = tin.InBondNum();
  vector<Bond> tBonds = tin.bond();

  vector<int> labs(tin.bond().size()+1);
  for (unsigned i=0; i!=labs.size(); ++i )
  {
    labs[i] = i;
  }

  vector<Bond> usBonds( tBonds.begin(), tBonds.begin()+inBdN );
  usBonds.push_back( Bond( BD_OUT, chi) );
  lten = UniTensor<T> (usBonds);
  lten.PutBlock( lqResult[0] );
  lten.SetLabel( vector<int> (labs.begin(), labs.begin()+inBdN+1 ) );

  vector<Bond> svBonds( tBonds.begin()+inBdN, tBonds.end());
  svBonds.insert( svBonds.begin(), Bond( BD_IN, chi) );
  qten = UniTensor<T> (svBonds);
  qten.PutBlock( lqResult[1] );
  qten.SetLabel( vector<int> (labs.begin()+inBdN, labs.end() ) );
}

//sperate tensors using svd with intter connecting bond truncate into chi(if chi=0 don't truncate)
template<typename T>
double svdSeperateGeneral( const UniTensor<T> &tin, UniTensor<T> &resultL, UniTensor<T> &resultR, Matrix<T> &singVal, int chi )
{
  vector<Matrix<T>> svdResult = Svd( tin.GetBlock() );
  if ( (chi == 0) || (chi>svdResult[0].col()) )
  {
    chi = svdResult[0].col();
  } else {}

  //Resize
  Resize( svdResult[0], svdResult[0].row(), chi, INPLACE );
  Resize( svdResult[1], chi,                chi, INPLACE );
  Resize( svdResult[2], chi, svdResult[2].col(), INPLACE );

  //update
  singVal = svdResult[1];
  singVal *= 1.0/Norm( singVal );

  int bdN = tin.BondNum();
  int inBdN = tin.InBondNum();

  vector<Bond> tBonds = tin.bond();
  vector<Bond> usBonds( tBonds.begin(), tBonds.begin()+inBdN );
  usBonds.push_back( Bond( BD_OUT, chi) );
  resultL = UniTensor<T> (usBonds);
  resultL.PutBlock( svdResult[0] );

  vector<Bond> svBonds( tBonds.begin()+inBdN, tBonds.end());
  svBonds.insert( svBonds.begin(), Bond( BD_IN, chi) );
  resultR = UniTensor<T> (svBonds);
  resultR.PutBlock( svdResult[2] );

  double trunErr = 0;
  return trunErr;
}

template<typename T>
vector<Matrix<T>> rsvd( const Matrix<T> &amat, const int nRankIn, double &trunRsvd ){
  int nRow = amat.row();
  int nCol = amat.col();
  int nRank = min(nRankIn,min(nRow,nCol));

  Matrix<T> omat(nRow,nRank);
  omat.Randomize('N', 0, 1, uni10_clock);

  Matrix<T> ymat = Dot( Dagger(amat),omat );
  vector<Matrix<T>> yqr = Qr( ymat );
  ymat = yqr[0];
  
  Matrix<T> bmat = Dot(amat,ymat);
  Matrix<T> pmat(nRank,nRank);
  pmat.Randomize('N', 0, 1, uni10_clock);
  Matrix<T> zmat = Dot(bmat,pmat);
  vector<Matrix<T>> zqr = Qr( zmat );
  zmat = zqr[0];

  Matrix<T> cmat = Dot(Dagger(zmat),bmat);
  vector<Matrix<T>> usv = Sdd(cmat);

  vector<Matrix<T>> result;
  result.push_back( Dot(zmat,usv[0]) );
  result.push_back( usv[1] );
  result.push_back( Dot(usv[2],Dagger(ymat)) );

  trunRsvd = 1.0-pow( Norm( usv[1] ), 2);
  return result;
}


template<typename T1, typename T2>
void absorbMatrixtoTensor( UniTensor<T1> &ten, const Matrix<T2> &mat, const int posit, const bool tensorFirst ){
  assert( (posit>=0)&&(posit<ten.BondNum()) );
  int connDim;
  if (tensorFirst){
    connDim = mat.row();
  }
  else {
    connDim = mat.col();
  }
  assert( connDim==ten.bond(posit).dim());
  //record
  int allbdN = ten.BondNum();
  int inbdN  = ten.InBondNum();
  vector<int> oldlab = ten.label();
  //assign lab and contract
  vector<int> newlab(allbdN);
  for (int i=0; i!=allbdN; i++){
    newlab[i] = i+1;
  }
  ten.SetLabel( newlab );
  UniTensor<T2> matT(mat);
  if (tensorFirst){
    matT.SetLabel(vector<int> { (posit+1),-(posit+1)});
  }
  else {
    matT.SetLabel(vector<int> {-(posit+1), (posit+1)});
  }
  ten = Contract( ten, matT, false );
  //permute back
  ten.SetLabel( posit+1, allbdN-1 );
  ten = Permute( ten, newlab, inbdN );
  ten.SetLabel(oldlab);
}

template<typename T>
double tebdHastings( UniTensor<T> &theta, UniTensor<T> &tenL, UniTensor<T> &tenR, Matrix<T> &matL, Matrix<T> &matC, const int chi ){
  int inBdN = theta.InBondNum();
  vector<Bond> thetaBds = theta.bond();
  assert( thetaBds.at(inBdN-1).dim()==thetaBds.back().dim() );
  vector<int> oldLab = theta.label();
  UniTensor<T> thetaAll = theta;
  bondcat( thetaAll, matL, inBdN-1 );
  vector<Matrix<T>> usv = QrSdd( thetaAll.GetBlock() );
  usv[1] *= 1.0/Norm(usv[1]);
  Resize( usv[0], usv[0].row(),          chi, INPLACE);
  Resize( usv[1],          chi,          chi, INPLACE);
  Resize( usv[2],          chi, usv[2].col(), INPLACE);
  vector<Bond> vbds( thetaBds.begin()+inBdN, thetaBds.end() );
  vbds.insert( vbds.begin(), Bond(BD_IN,chi) );
  //vector<Bond> vbds = { Bond(BD_IN,chi), Bond(BD_OUT, chi), Bond(BD_OUT, dimD) };
  tenR = UniTensor<T> (vbds);
  tenR.PutBlock( usv[2] );
  vector<int> newThetaLab;
  for (int i=0; i!=thetaBds.size(); i++){
    newThetaLab.push_back(i);
  }
  theta.SetLabel( newThetaLab );
  vector<int> thetaRLab( newThetaLab.begin()+inBdN, newThetaLab.end() );
  thetaRLab.insert( thetaRLab.begin(), thetaBds.size() );
  tenR.SetLabel( thetaRLab );
  UniTensor<T> tenRc=  Conj(tenR);
  tenL = Contract( theta, tenRc, false );
  double trunErr = 1.0-Norm(usv[1]);
  matC = usv[1];
  matC *= 1.0/( Norm(matC) );
  theta.SetLabel(oldLab);
  return trunErr;
}

template<typename T>
vector<Matrix<T>> cholesky_decomposition( Matrix<T> &Eigh_mat ){
  bool if_debug = false;
  ///return U_Dagger*s^(0.5), s^(-0.5)*U;
  vector<Matrix<T>> sU = EigH( Eigh_mat );
  Matrix<T> U = sU.at(1);
  int Dbond = Eigh_mat.col();
  Matrix<T> sqrt_lambdas( Dbond, Dbond, true);
  Matrix<T> sqrt_lambdas_inverse( Dbond, Dbond, true);
  complex<double> test = ( fabs(sU.at(0).At(0,0)) > fabs(sU.at(0).At(Dbond-1,Dbond-1)) )? sU.at(0).At(0,0):sU.at(0).At(Dbond-1,Dbond-1);
  double factor = (test.real()>0)? 1.0:-1.0;///prevent input vector take an overall -1.0 factor
  if (factor>0){
    for (int i=0; i!=Dbond; i++){
      if ( fabs( sU.at(0).At( i, i)) > 1.0e-15){
        sqrt_lambdas[Dbond-i-1] = sqrt(sU.at(0).At( i, i));
        sqrt_lambdas_inverse[Dbond-i-1] = 1.0/sqrt_lambdas.At( Dbond-i-1, Dbond-i-1);
      }
      else {
        sqrt_lambdas[Dbond-i-1] = 0;
        sqrt_lambdas_inverse[Dbond-i-1] = 0;
      }
    }
    for (int i=0; i!=Dbond; i++){
      for (int j=0; j!=Dbond; j++){
        U[i*Dbond+j]=sU.at(1)[(Dbond-i-1)*Dbond+j];
      }
    }
  }
  else {
    for (int i=0; i!=Dbond; i++){
      if ( fabs( sU.at(0).At( i, i)) > 1.0e-15){
        sqrt_lambdas[i] = sqrt(-sU.at(0).At( i, i));
        sqrt_lambdas_inverse[i] = 1.0/sqrt_lambdas.At(i,i); 
      }
      else {
        sqrt_lambdas[i] = 0;
        sqrt_lambdas_inverse[i] = 0;
      }
    }
  }

  Matrix<T> U_Dagger = Transpose( U );

  vector<Matrix<T>> output;
  output.reserve(2);
  output.push_back( Dot( U_Dagger, sqrt_lambdas) ); 
  output.push_back( Dot( sqrt_lambdas_inverse, U ) ); 
  if (if_debug){
    cout<<sU[0]<<endl;
    cout<<output[0]<<endl;
    cout<<output[1]<<endl;
  } else {}
  return output;
}

template<typename T>
vector<Matrix<T>> QrSdd( const Matrix<T> &Amat ){
  vector<Matrix<T>> qdr = QdrColPivot(Amat);
  Matrix<T> rp = Dot(qdr[1],qdr[2]);
  vector<Matrix<T>> qrResult = Svd(rp);
  qrResult[0] = Dot(qdr[0],qrResult[0]);
  //vector<Matrix<T>> qrResult = Svd(Amat);
  return qrResult;
}

template<typename T>
vector<Matrix<T>> Qr_nm( const Matrix<T> &Amat ){
  vector<Matrix<T>> result;
  if (Amat.row()<Amat.col()){
    Matrix<T> Atemp;
    Resize( Atemp, Amat, Amat.row(), Amat.row(), INPLACE );
    result = Qr( Atemp );
    Resize( result[0], Amat.row(), Amat.row(), INPLACE );
    Resize( result[1], Amat.row(), Amat.col(), INPLACE );
  }
  else {
    result = Qr( Amat );
  }
  return result;
}

template<typename T>
Matrix<T> pseudoInverse( const Matrix<T> &mat )
{
  assert( mat.col()==mat.row() );
  vector<Matrix<T>> usv = Sdd( mat );
  Matrix<T> dinv = usv[1];
  for (int i=0; i!=usv[1].col(); ++i)
  {
    if ( usv[1][i]>1.0e-14 )
    {
      dinv[i] = 1.0/usv[1][i];
    }
    else 
    {
      dinv[i] = 0;
    }
  }
  //cout<<setprecision(8)<<usv[1]<<endl;
  Matrix<T> out = Dot( Dagger(usv[2]), Dot(dinv,Dagger(usv[0])));
  return out;
}

#endif
