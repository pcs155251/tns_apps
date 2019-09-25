#ifndef __ONESITEIMPS_H_INCLUDED__
#define __ONESITEIMPS_H_INCLUDED__

#include <sys/stat.h> //for mkdir
#include "../arpack_wrapper/arpack_wrapper.h"
#include "../nsy_ham/hamiltonian.h"
#include "../tools/general_tools.h"

template<typename T> 
void reverseCol( uni10::Matrix<T> &matIn, const int iCol ){
  const int nRow = matIn.row();
  const int nCol = matIn.col();
  for ( int i=0; i!=nRow; i++ ){
    matIn[ i*nCol+iCol ]*= -1.0;
  }
}

template<typename T> 
void reverseRow( uni10::Matrix<T> &matIn, const int iRow ){
  const int nCol = matIn.col();
  for ( int i=0; i!=nCol; i++ ){
    matIn[ iRow*nCol+i ]*= -1.0;
  }
}

template<typename T> 
std::vector<uni10::Matrix<T>> qrPositveR( const uni10::Matrix<T> &matIn ){
  vector<Matrix<T>> output = Qr( matIn );
  const int dim = output.at(0).col();
  for (int i=0; i!=dim; i++){
    complex<double> val = output.at(1).At(i,i);
    if ( val.real() < 0 ){
      reverseCol( output.at(0), i );
      reverseRow( output.at(1), i );
    } else {}
  }
  return output;
}

template<typename T> 
std::vector<uni10::Matrix<T>> lqPositveL( const uni10::Matrix<T> &matIn ){
  vector<Matrix<T>> output = Lq( matIn );
  const int dim = output.at(0).row();
  for (int i=0; i!=dim; i++){
    complex<double> val = output.at(0).At(i,i);
    if ( val.real() < 0 ){
      reverseCol( output.at(0), i );
      reverseRow( output.at(1), i );
    } else {}
  }
  return output;
}

//! A class describing infinite MPS with single site unit cell.
/*!
 *  oneSiteiMPS is used for the channel environment.
 *  The method described in PRB 94, 155123 is implemented to find the canonical form of the MPS.
 */
template<typename T> 
class oneSiteiMPS{
  public:
    oneSiteiMPS( );
    //! Constructor
    /*!
     *  \param[in] horizDimIn Boundary bond dimension of MPS.
     *  \param[in] ansatzIn Tensor appears successively and infinitely on the MPS.
     */
    oneSiteiMPS( const int horizDimIn, const uni10::UniTensor<complex<double>> &ansazIn );
    //! Copy constructor.
    oneSiteiMPS( const oneSiteiMPS &input );
    //oneSiteiMPS( const std::string &rootFolder, const std::string &subFolder );
    //! Copy constructor.
    oneSiteiMPS* operator*=( const complex<double> rhs ){
      this->aC *= rhs;
      this->aL *= rhs;
      this->aR *= rhs;
      this->cmat *= rhs;
      this->lambdaL *= pow(fabs(rhs),2);
      this->lambdaR *= pow(fabs(rhs),2);
    }
    
    //get data member
    const uni10::UniTensor<complex<double>>& getAnsaz() const { return ansaz; }
    const uni10::UniTensor<complex<double>>& getEvolOperator() const { return evolComplex; }
    const uni10::UniTensor<complex<double>>& getIsometryL() const { return isometryLeft; }
    const uni10::UniTensor<complex<double>>& getIsometryR() const { return isometryRight; }
    const uni10::UniTensor<complex<double>>& getAC() const { return aC; }
    const uni10::UniTensor<complex<double>>& getAL() const { return aL; }
    const uni10::UniTensor<complex<double>>& getAR() const { return aR; }
    //interface
    void setEvolOperator( const uni10::UniTensor<double> &evolOperatorIn );
    void canonical( const int maxIter, const double errTor );
    void fixedPoint( const int maxIter, const double errTor, const double arpackErr );
    void printTensors( );
    //test
    void testCanonical( ); //test three conditions: left canonical, right canonical, aL*C = C*aR
    void testFixedPoint( );
    //auxiliary functions
    void SaveTensors( const std::string &rootFolder );
    void LoadTensors( const std::string &rootFolder );

  private:
    //tensors
    uni10::UniTensor<complex<double>> aC;
    uni10::UniTensor<complex<double>> aL;
    uni10::UniTensor<complex<double>> aR;
    uni10::UniTensor<complex<double>> isometryLeft;
    uni10::UniTensor<complex<double>> isometryRight;
    uni10::Matrix<complex<double>> cmat;
    uni10::UniTensor<complex<double>> ansaz;
    uni10::UniTensor<complex<double>> evolComplex;
    complex<double> lambdaL;
    complex<double> lambdaR;
    complex<double> lambdaAC;
    complex<double> lambdaC;
    //data
    int horizDim;
    int vertiDim;
    //std::vector<int> evolOperLab;
    std::vector<int> gammaLab;
    //implementation
    void isometries( const double arpackErr );
    void Normalize( );
    void improveAcC( const int arpackIter, const double arpackErr );
};

#endif
