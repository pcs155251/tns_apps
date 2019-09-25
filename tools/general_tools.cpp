#include "general_tools.h"

using namespace std;
using namespace uni10;

Matrix<double> matId( const int dim ){
  Matrix<double> id( dim, dim, true );
  id.Identity();
  return id;
}

Matrix<uni10_complex128> cMatrix( const Matrix<double> &test ){
  Matrix<uni10_complex128> copy( test.row(), test.col(), false );
  for ( int i=0; i!=test.row(); i++ ){
    for ( int j=0; j!=test.col(); j++ ){
      copy[ i*test.col()+j ].real( test.At( i, j) );
      copy[ i*test.col()+j ].imag( 0 );
    }
  }
  return copy;
}

Matrix<uni10_double64> EignLeft( const Matrix<uni10_double64> &test, uni10_double64 &lambda, const double cutoff ){
  assert( test.row() == test.col() );
  int dim = test.row();
  Matrix<uni10_complex128> copy = cMatrix( test );
  //Eig
  vector<Matrix<uni10_complex128>> EignDecomp = Eig( copy );
  //find max abs Eigenvalue
  int iabsMax = 0;
  double absMax = 0;
  for ( int i=0; i!=dim; i++ ){
    double absNow = abs( EignDecomp.at(0)[i] );
    if ( absNow > absMax ){
      iabsMax = i;
      absMax = absNow;
    } else {}
  }
  //assign lambda and Eigenvector
  assert( (EignDecomp.at(0)[iabsMax].imag())<cutoff );
  lambda = EignDecomp.at(0)[iabsMax].real();
  Matrix<uni10_double64> EignVector( 1, dim, false );
  for ( int i=0; i!=dim; i++){
    assert( abs(EignDecomp.at(1)[ iabsMax*dim+i ].imag())<cutoff );
    EignVector[i] = EignDecomp.at(1)[ iabsMax*dim+i ].real();
  }
  return EignVector;
}

Matrix<uni10_double64> EignRigh( const Matrix<uni10_double64> &test, uni10_double64 &lambda, const double cutoff ){
  Matrix<uni10_double64> trans = Transpose( test );
  Matrix<uni10_double64> EignVector = EignLeft( trans, lambda, cutoff );
  return Transpose( EignVector ); 
}

UniTensor<double> getEvoOperator( const UniTensor<double> &hamiltonian, const double tau ){
  UniTensor<double> evoOperator( hamiltonian.bond() );
  evoOperator.PutBlock( ExpH( -tau, hamiltonian.GetBlock() ) );
  evoOperator.SetLabel(hamiltonian.label());
  return evoOperator;
}

UniTensor<complex<double>> getEvoOperator( const UniTensor<complex<double>> &hamiltonian, const double tau ){
  UniTensor<complex<double>> evoOperator( hamiltonian.bond() );
  evoOperator.PutBlock( ExpH( -tau, hamiltonian.GetBlock() ) );
  evoOperator.SetLabel(hamiltonian.label());
  return evoOperator;
}

void truncateOneBond( UniTensor<double> &ten, const int posit, const int newDim ){
  const vector<int> oldLab = ten.label();
  const int nIn = ten.InBondNum();
  //construct new label
  vector<int> newLab = oldLab;
  const int moveIndex = newLab.at(posit);
  newLab.erase( newLab.begin()+posit );
  newLab.insert( newLab.begin(), moveIndex );
  //Permute and trucnate
  ten = Permute( ten, newLab, 1 );
  Matrix<double> temp = ten.GetBlock();
  Resize( temp, newDim, temp.col(), INPLACE );
  //construct new tensor, put, Permute back
  vector<Bond> newBonds = ten.bond();
  newBonds.erase( newBonds.begin() );
  newBonds.insert( newBonds.begin(), Bond( BD_IN, newDim ) );
  ten = UniTensor<double> ( newBonds );
  ten.PutBlock( temp );
  ten.SetLabel( newLab );
  ten = Permute( ten, oldLab, nIn );
}

