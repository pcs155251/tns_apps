#ifndef __CHANNEL_H_INCLUDED__
#define __CHANNEL_H_INCLUDED__

#include <sys/stat.h> //for mkdir
#include "../oneSiteiMPS/oneSiteiMPS.h"

//! A class defined for channel environment of 2D iPEPS ansatz.
/*!
 *  2D infinite iPEPS ansatz with single tensor unit cell is assumed.
 *  The method in this class is based on PRB 94, 155123.
 *  The computation of channel environment is done while variational optimization is under development.
 */
class channel{
  public:
    //! Constructor
    /*!
     *  \param[in] dimMPS Boundary bond dimension of channel environment.
     *  \param[in] ansazIn iPEPS ansatz with single tensor unit cell. The bond start from left bond and continue counterclockwisely.
     *  \param[in] hamiltonianIn Two site hamiltonian which the ground state is optimized for.
     */
    channel( const int dimMPSIn, const UniTensor<complex<double>> &ansazIn, const UniTensor<double> &hamiltonianIn );
    
    //! Constructor
    /*!
     *  \param[in] iterMax Maximum iteration step for finding channel environment.
     *  \param[in] canonicalErr Error criterion for canonical form.:
     *  \param[in] arpErr Error criterion for finding dominant eigen value and eigen vector by arnoldi method.
     */
    void findAllEnv( const int iterMax, const double canonicalErr, const double arpErr );

    //! Find gradient by the method proposed in the paper.
    UniTensor<complex<double>> findGradient();
    //! Measure energy with current environment
    /*  
     *  \param[in] ifShow Whether showing detailed information or not.
     */
    complex<double> measureEnergy( const bool ifShow) const ;
    //! Save environment tensors.
    /*!
     *  \param[in] rootFolder Folder where tensors are saved.
     */
    void SaveEnvTensors( const string &rootFolder );
    //! Load environment tensors.
    /*!
     *  \param[in] rootFolder Folder where tensors are loaded from.
     */
    void LoadEnvTensors( const string &rootFolder );
    void addTensor( const UniTensor<complex<double>> &addOn );
    //! Measure expectation of twoSite operator.
    /*!
     *  \param[in] iDir Specify which direction of bond is measured. 0-> left, 1-> up, 2->right, 3->bottom.
     */
    complex<double> twoSiteOpExpec( const UniTensor<double> &twoSiteOp, const int iDir ) const;

    //!  Normalize tensors such that its L2 norm is identity.
    void normalize();
    void test();
    UniTensor<double> numericalGradientReal( const double delta, const int iterMax, const double canonicalErr, const double arpErr );
    //! Find gradient numerically.
    UniTensor<complex<double>> numericalGradientComplex( const double delta, const int iterMax, const double canonicalErr, const double arpErr );

    vector<uni10::UniTensor<complex<double>>> sups;
    vector<uni10::UniTensor<complex<double>>> bluesI;
    vector<uni10::UniTensor<complex<double>>> bluesII;
    vector<uni10::UniTensor<complex<double>>> bluecaps;
    vector<uni10::UniTensor<complex<double>>> reds;
  private:
    //tensors
    const int dimPhy;
    int dimMPS;
    int dimD;
    uni10::UniTensor<complex<double>> ansaz;
    uni10::UniTensor<double> hamiltonian;
    vector<oneSiteiMPS<complex<double>>> edges; 
    vector<uni10::UniTensor<complex<double>>> caps; 
    vector<uni10::UniTensor<complex<double>>> corners;
    //implementation
    //environment
    void findOneCorner( const int iDir, unsigned int max_iter, double err_tol );//checked
    void findOneCap( const int iDir, unsigned int max_iter, double err_tol, complex<double> &EignValue );//checked
    //structure factor (to be checked)
    void findOneSup( const int iDir, const double q_y);//checked
    //gradient
    void findOneBlueCap( const int iDir );//checked
    void findOneBlueI( const int iDir );//checked
    void findOneBlueII( const int iDir );//checked
    void findOneRed( const int iDir );//checked
   
    uni10::UniTensor<complex<double>> grad0( const int iDir );//checked
    uni10::UniTensor<complex<double>> grad1( const int iDir );//checked 
    uni10::UniTensor<complex<double>> grad2( const int iDir );//checked
    uni10::UniTensor<complex<double>> grad3( const int iDir );//checked
    uni10::UniTensor<complex<double>> grad4( const int iDir );//checked
    uni10::UniTensor<complex<double>> grad5( const int iDir );//checked
    //double check
    UniTensor<double> twoSiteCluster( UniTensor<double> twoSiteOp, const int iDir ) const;
    complex<double> findNorm( );
};     

#endif
