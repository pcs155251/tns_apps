#ifndef __PESS_HPP_INCLUDED__
#define __PESS_HPP_INCLUDED__

#include <sys/stat.h> //for mkdir
#include <complex.h> 
#include <fstream>
#include <vector>
#include <string>
#include "uni10.hpp"
#include "../tools/general_tools.h"
#include "enPess.h"
#include "enRedPess.h"

//! A 3-PESS ansatz on square kagome lattice. Simple update and measurement considering full environment by CTM is implemented. 
/*!
 *  There are two conventional CTM approach is implemented as enPess and dimension reduction CTM approach is implemneted as enRedPess.
 *  One should declare the type of environment as second argument of the template class pess.
 */
template<typename T, template<typename T> typename enType>
class pess{
  public:
    //! Constructor.
    pess() {}
    //! Constructor.
    /*!
     *  \param[in] dimVirIn Virtual bond dimension of 3-PESS ansatz.
     */
    pess( const int dimVirIn );
    //! Randomize the 3-PESS ansatz.
    void Randomize( );
    //! Perform simple update iterativley with first order Suzuki-Trotter decomposition.
    /*!
     *  \param[in] ham Three site hamiltonian acting on up or down triangle of Husimi lattice.
     *  it is used to create the time evolution operator and measure energy.
     *  \param[in] tau Time step for time evolution.
     *  \param[in] maxIter Maximum iteration step of simple update.
     *  \param[in] criterion Convergence criterion of simple update.
     *  If energy difference between two time steps is smaller than criterion, the ansatz is considered converged.
     *  \return Number of steps to reach convergence. 
     *  If the ansatz is not converge within maximum iteration step, returns -1.
     */
    int optimize1stOrder( uni10::UniTensor<T> ham, const double tau, const int maxIter, const double criterion );
    //! Perform simple update iterativley with second order Suzuki-Trotter decomposition.
    /*!
     *  \param[in] ham Three site hamiltonian acting on up or down triangle of Husimi lattice.
     *  it is used to create the time evolution operator and measure energy.
     *  \param[in] tau Time step for time evolution.
     *  \param[in] maxIter Maximum iteration step of simple update.
     *  \param[in] criterion Convergence criterion of simple update.
     *  If energy difference between two time steps is smaller than criterion, the ansatz is considered converged.
     *  \return Number of steps to reach convergence. 
     *  If the ansatz is not converge within maximum iteration step, returns -1.
     */
    int optimize2ndOrder( uni10::UniTensor<T> ham, const double tau, const int maxIter, const double criterion );
    //! Measure expectation value of three site operator on up triangle or down triangle with simple environment.
    /*! 
     * \param[in] op Three site operator.
     * \param[in] iCore iCore=0 -> up triangle, iCore=! down triangle.
     */
    double meaTriSiteLocal( const uni10::UniTensor<T> &op, const int iCore );
    //! Get virtual bond dimension of 3-PESS ansatz.
    /*! 
     * \return Virtual bond dimension.
     */
    int getVirDim( ) const { return dimVir; }
    
    //! Save iPEPS ansatz tensors into given folder and field.
    /*!
     *  \param[in] folder The folder where tensors are saved.
     *  \param[in] field A parameter in hamiltonian. 
     *  A corresponding subfolder will be created to save the wave function evolved by the hamiltonian with specified field.
     */
    void saveWv( const std::string &tensorFolder, double field ); 
    //! Load iPEPS ansatz tensors from given folder and field.
    /*!
     *  \param[in] folder The folder where tensors are loaded.
     *  \param[in] field A parameter in hamiltonian. 
     */
    void loadWv( const std::string &tensorFolder, double field ); 

    //! Measure important physical quantities.
    /*!
     *  \param[in] ham
     *  \param[out] energyUp Three site energy of up triangle.
     *  \param[out] energyDn Three site energy of down triangle.
     *  \param[out] spinComp Spin component of each site. spinComp = vector<vector<double>> { S0, S1, S2}, where Si = vector<double> { Six, Siy, Siz}.
     */
    void meaQuantities( uni10::UniTensor<T> ham, double &energyUp, double &energyDn, std::vector<std::vector<double>> &spinComp );

    //! Initialized environment.
    /*!
     *  \param[in] edgeDim Boundary bond dimension of envionment.
     */
    void initEn( int edgeDim );

    //! Change boundary bond dimension of envionment.
    /*!
     *  \param[in] edgeDim New boundary bond dimension of envionment.
     */
    void changeEdgeDim( int edgeDim );

    //! Compute the environment by CTM.
    /*!
     *  \param[in] niter Maximum iteration step of CTM.
     *  \param[in] ham Three site hamiltonian. It is used to measure energy to monitor the convergence of environment.
     *  \param[in] ifsvd Whether SVD is employed to find isometries in CTM. It it is set to false, QR decomposition is employed.
     */
    void computeEn( int niter, double errTol, UniTensor<T>& ham, bool ifsvd ) ;
    
    //! Release the memory for storing envionment tensors.
    void releaseEn();

    //! Save envionment tensors into given folder and field.
    /*!
     *  \param[in] folder The folder where environment tensors are saved.
     *  \param[in] field A parameter of hamiltonian. 
     *  A corresponding subfolder will be created. 
     *  A subsubfolder indicating currentent boundary bond dimension and up environment or down environment will be created to save the environment tensors.
     */
    void saveEn( string folder, double field);
    //! Load envionment tensors from given folder and field.
    /*!
     *  \param[in] folder The folder where environment tensors are loaded.
     *  \param[in] field A parameter of hamiltonian. 
     *  \param[in] edgeDim Environment tensors are loaded from folder/fieldSubFolder/edgeDimSubFolder/. 
     *  Environment tensors of up and down triangle will be load automatically.
     */
    void loadEn( string folder, double field, int edgeDim );
  private:
    //default settings
    const int nCore = 2;
    const int nBdsOfCore = 3;
    int dimVir;
    std::vector<std::vector<int>> thetaLabForHosvd;
    //tensors
    std::vector<uni10::UniTensor<T>> coreTensors;
    std::vector<uni10::UniTensor<T>> projTensors;
    std::vector<std::vector<uni10::Matrix<T>>> bondVectors;
    //implementation
    uni10::UniTensor<T> ContractCore( const int iCore ); //contract one of the coreTensors with all projTensors with corresponding bondVectors
    uni10::UniTensor<double> constructBten( const bool ifUp );
    uni10::UniTensor<double> getaten( UniTensor<double> &bten );
    double timeEvoAstep( uni10::UniTensor<T> &evoOperator, const int iCore);
    void setHamLab( uni10::UniTensor<T> &ham);
    enType<T>* enUp = 0;
    enType<T>* enDn = 0;
};

#endif
