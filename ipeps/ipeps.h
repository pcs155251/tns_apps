#ifndef __TWODIMITEBD_H_INCLUDED__
#define __TWODIMITEBD_H_INCLUDED__

#include <sys/stat.h> //for mkdir
#include <ctime>
#include <complex.h> 
#include <fstream>
#include "uni10.hpp"
#include "../tools/general_tools.h"
#include "../nsy_opt/operator.h"
#include "enIpeps.h"

//! An ipeps ansatz on square lattice. Optimization method such as full update and simple update is implemented. 2 by 2 unit cell with checkerboard pattern is assumed.
template<typename T>
class ipeps
{
  public:
    //! Constructor
    /*!
     *  \param[in] virDim Virtual bond dimension of ansatz.
     */
    ipeps( int virDim );
    ~ipeps();
    //! Perform simple update.
    /*!
     *  \param[in] ham Two site operator by which evolution operator is constructed.
     *  \param[in] tau Time step of time evolution.
     *  \param[in] evoStop Convergence criterion for time evolution. When energy difference between two steps is smaller than evoStop, tha state
     *  is considered converged the iteration stop.
     *  \param[in] maxIter Maximum iteration step.
     */
    void simpUpdate( UniTensor<T> &ham, double tau, double evoStop, int maxIter );
    //! Perform full udpate.
    /*!
     *  \param[in] ham Two site operator by which evolution operator is constructed.
     *  \param[in] tau Time step of time evolution.
     *  \param[in] evoStop Convergence criterion for time evolution. When energy difference between two steps is smaller than evoStop, tha state
     *  is considered converged the iteration stop.
     *  \param[in] maxIter Maximum iteration step.
     *  \param[in] iterEn Maximum iteration step for environment convergence.
     *  \param[in] errEn Convergence critertion for computing environment by CTM.
     *  \param[in] edgeDim Boundary bond dimension for environment (CTM).
     *  \return Whether the environment is converged by CTM of not after the wavefunction is updated.
     */
    bool fullUpdate( UniTensor<T> &ham, double tau, double evoStop, int maxIter, int iterEn, double errEn, int edgeDim );
    //! Perform full udpate with gauge fixing.
    /*!
     *  \param[in] ham Two site operator by which evolution operator is constructed.
     *  \param[in] tau Time step of time evolution.
     *  \param[in] evoStop Convergence criterion for time evolution. When energy difference between two steps is smaller than evoStop, tha state
     *  is considered converged the iteration stop.
     *  \param[in] maxIter Maximum iteration step.
     *  \param[in] iterEn Maximum iteration step for environment convergence.
     *  \param[in] errEn Convergence critertion for computing environment by CTM.
     *  \param[in] edgeDim Boundary bond dimension for environment (CTM).
     *  \return Whether the environment is converged by CTM of not after the wavefunction is updated.
     */
    bool fullUpdateGauge( UniTensor<T> &ham, double tau, double evoStop, int maxIter, int iterEn, double errEn, int edgeDim );
    //! Change boundary bond dimension of environment.
    /*!
     *  \param[in] edgeDim New boundary bond dimension.
     */
    void changeEdgeDim( int edgeDim );
    //! Measure norm using simple environment (singular matrix).
    /*!
     *  \param[in] idir Specify which direction of bond is being measured. Start from left -> up -> right -> dn. 0 <= idir <= 3.
     *  \return Norm value.
     */
    T meaSimpNorm( int idir );
    //! Measure two site operator using simple environment (singular matrix).
    /*!
     *  \param[in] twoSiteOp Two site operator being measured.
     *  \param[in] idir Specify which direction of bond is being measured. Start from left -> up -> right -> dn. 0 <= idir <= 3.
     *  \return two site expecatation value.
     */
    T meaSimpTwoSiteOp( UniTensor<T>& twoSiteOp, int idir );
    //! Measure one site operator using simple environment (singular matrix).
    /*!
     *  \param[in] oneSiteOp One site operator being measured.
     *  \param[in] idir Specify which direction of bond is being measured. Start from left -> up -> right -> dn. 0 <= idir <= 3.
     *  \return One site expecatation value.
     */
    T meaSimpOneSiteOp( UniTensor<T>& oneSiteOp, int isite );
    //! Measure one site norm using full environment (obtained by CTM).
    /*!
     *  \param[in] gammas iPEPS ansatz being measured.
     *  \param[in] isite Specify which site is being measured. left site -> 0, right site -> 1.
     *  \return One site norm value.
     */
    double meaOneSiteNorm( const std::vector<uni10::UniTensor<T>>& gammas, int isite ) const 
    { 
      assert( en!=0 );
      return en->meaOneSiteNorm( gammas, isite ); 
    }
    //! Measure one site norm using full environment (obtained by CTM).
    /*!
     *  \param[in] gammas iPEPS ansatz being measured.
     *  \param[in] idir Specify which direction of bond is being measured. Start from left -> up -> right -> dn. 0 <= idir <= 3.
     *  \return One site norm value.
     */
    double meaTwoSiteNorm( const std::vector<uni10::UniTensor<T>>& gammas, int ibond ) const 
    { 
      assert( en!=0 );
      return en->meaTwoSiteNorm( gammas, ibond ); 
    }
    //! Measure one site expecatation using full environment (obtained by CTM).
    /*!
     *  \param[in] gammas iPEPS ansatz being measured.
     *  \param[in] isite Specify which site is being measured. left site -> 0, right site -> 1.
     *  \return One site expecatation value.
     */
    double meaOneSiteExp ( const std::vector<uni10::UniTensor<T>>& gammas, uni10::UniTensor<T> &oneSiteOp, int isite ) const 
    { 
      assert( en!=0 );
      return en->meaOneSiteExp( gammas, oneSiteOp, isite );
    }
    //! Measure two site expecatation using full environment (obtained by CTM).
    /*!
     *  \param[in] gammas iPEPS ansatz being measured.
     *  \param[in] idir Specify which direction of bond is being measured. Start from left -> up -> right -> dn. 0 <= idir <= 3.
     *  \return two site expecatation value.
     */
    double meaTwoSiteExp ( const std::vector<uni10::UniTensor<T>>& gammas, uni10::UniTensor<T> &twoSiteOp, int ibond ) const
    { 
      assert( en!=0 );
      return en->meaTwoSiteExp( gammas, twoSiteOp, ibond );
    }
    //! Get virtual bond dimension of current ansatz.
    /*!
     *  \return current virtual bond dimension.
     */
    int getVirDim( ) const { return _virDim; }
    //! Get reduced tensors within unit cell.
    /*!
     *  \return current reduced tensors within unit cell.
     */
    vector<UniTensor<T>> getUnit();
    //! Initialize the environment tensors for CTM.
    /*!
     *  \param[in] edgeDim Boundary bond dimension used for CTM.
     */
    void initEn( int edgeDim );
    //! Compute the environment by CTM.
    /*!
     *  \param[in] edgeDim Boundary bond dimension used for CTM.
     *  \param[in] niter Maximum iteration step of CTM.
     *  \param[in] errTol Convergence criterion for CTM.
     *  \param[in] ham Two site hamiltonain. Used to measure the energy and monitor the convergence.
     */
    void computeEn( int edgeDim, int niter, double errTol, UniTensor<T>& ham );
    //! Normalize gammas matrix (singular matrix obtained from svd of simple update).
    /*!
     *  The matrix are normalized such that it has L2 norm equal to one.
     */
    void normalizeGamma( );
    //! Release the memory used to stored the environment tensors.
    void releaseEn() { delete en; en=0;}
    //void test( UniTensor<T>& evolOp, int idir );
    //! Get iPEPS ansatz tensors (A and B).
    /*!
     *  \return iPEPS ansatz tensors (A and B).
     */
    const vector<UniTensor<T>>& getGammas(){return gammas;}

    //! Save iPEPS ansatz tensors into given folder and field.
    /*!
     *  \param[in] folder The folder where tensors are saved.
     *  \param[in] field A parameter in hamiltonian. 
     *  A corresponding subfolder will be created to save the wave function evolved by the hamiltonian with specified field.
     */
    void saveTensors( string folder, double field );
    
    //! Load iPEPS ansatz tensors from given folder and field.
    /*!
     *  \param[in] folder The folder where tensors are loaded.
     *  \param[in] field A parameter in hamiltonian. 
     */
    bool loadTensors( string folder, double field );
    
    //! Save envionment tensors into given folder and field.
    /*!
     *  \param[in] folder The folder where environment tensors are saved.
     *  \param[in] field A parameter of hamiltonian. 
     *  A corresponding subfolder will be created. 
     *  A subsubfolder indicating currentent boundary bond dimension will be created to save the environment tensors.
     */
    void saveEn( string folder, double field);
    //! Load envionment tensors from given folder and field.
    /*!
     *  \param[in] folder The folder where environment tensors are loaded.
     *  \param[in] field A parameter of hamiltonian. 
     *  \param[in] edgeDim Environment tensors are loaded from folder/fieldSubFolder/edgeDimSubFolder/
     */
    void loadEn( string folder, double field, int edgeDim );

  private: 
    int _phyDim = 2;
    int _virDim;
    vector<UniTensor<T>> gammas;
    vector<Matrix<T>> lambdas;
    void bondrmSqrAll( );
    void absrobSimpEn( int idir );
    void releaseSimpEn( int idir );

    //implementation
    void init();
    UniTensor<T> contractAB( int idir ); //contract AB 
    UniTensor<T> actTwoSiteOp( UniTensor<T>& twoSiteOp, int idir ); //contractAB and act two site operator upon
    double simpUpdateAstep( UniTensor<T>& evolOp, int idir );
    double simpUpdateAstepReduced( UniTensor<T>& evolOp, int idir );
    bool fullUpdateAstep( UniTensor<T>& ham, UniTensor<T>& evolOp, double &energy, int idir, int costIter, double costErr, int iterEn, double errEn, bool ifprint );
    bool fullUpdateAstepGauge( UniTensor<T>& ham, UniTensor<T>& evolOp, double &energy, int idir, int costIter, double costErr, int iterEn, double errEn, bool ifprint );
    void getSubTens( int idir, UniTensor<T> &xten, UniTensor<T>& arten, UniTensor<T>& yten, UniTensor<T>& blten ) const;
    void updateSubTens( int idir, UniTensor<T> &xten, UniTensor<T>& arten, UniTensor<T>& yten, UniTensor<T>& blten );
    UniTensor<T> contractSkel( int idir, const UniTensor<T> &xten, const UniTensor<T>& yten ) const;
    UniTensor<T> contractT( const UniTensor<T> &skel, const UniTensor<T> &thetaPrime ) const;
    UniTensor<T> contractR( const UniTensor<T> &skel, const UniTensor<T> &subTen, string type ) const;
    UniTensor<T> contractS( const UniTensor<T> &skel, const UniTensor<T> &subTen, const UniTensor<T>& thetaPrime, string type ) const;
    UniTensor<T> contractRInv( const UniTensor<T>& rten, const UniTensor<T>& sten, string type );
    double costFunction( UniTensor<T> subten, UniTensor<T> rten, UniTensor<T> sten, const UniTensor<T>& tten, string type );
    void herPosApproxGauge( UniTensor<T> &skel, UniTensor<T> &arten, UniTensor<T> &blten, Matrix<T> &rminv, Matrix<T> &lminv );
    enIpeps<T>* en = 0;
};

#endif
