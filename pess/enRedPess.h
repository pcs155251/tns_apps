#ifndef __ENREDPESS_H_INCLUDED__
#define __ENREDPESS_H_INCLUDED__

#include <sys/stat.h> //for mkdir
#include <ctime>
#include <complex.h> 
#include <fstream>
#include "../ctmBase/ctmBase.h"
#include "../nsy_ham/hamiltonian.h"
#include "../nsy_opt/operator.h"
#include "../tools/general_tools.h"
#include "../arpack_wrapper/arpack_wrapper.h"

//! A derived class of ctmBase describing the environment with dimension reduction technique of 3-PESS on kagome lattice.
/*!
 *  The 3-PESS ansatz with translational invariance is contracted to a 2D infinite TNS with 2 by 2 unit cell.
 *  The dimension reduction technique is empolyee such that the tensors within unit cell are b, i0, i1, b^{\dagger}.
 *  All the measurement is performed within the unit cell.
 *  enRedPess is uesd in pess to measure the physical quantites.
 *  Methods defined in ctmBase are employee to find the converged environment.
 */
template<typename T>
class enRedPess:
public ctmBase<T>
{
  public:

    //! Constructor
    /*! 
     *  The constructor initialized the environment by the parameters below.
     *  \param[in] coreTensors Core tensors of 3-PESS ansatz. 
     *  \param[in] projTensors Projected tensors of 3-PESS ansatz.
     *  \param[in] edgeDimIn Boundary bond dimension for the CTM environment.
     *  \param[in] ifUp Whether the environment is for the up triangle or down triangle.
     *  Since the unit cell can be defined by up triangle or down triangle, the environment for the two definition are different.
     */
    enRedPess(  std::vector<uni10::UniTensor<T>>& coreTensors, std::vector<uni10::UniTensor<T>>& projTensors, const int edgeDimIn, bool ifUp );
    
    //! Find the approximated environment iteratively by CTM technique.
    /*!
     *  ctmrgAstep() in ctmBase is employee to perform the CTM procedure.
     *  \param[in] niter  Maximum iteration step.
     *  \param[in] errTol Convergence criterion for CTM iteration.
     *  \param[in] triSiteHam Three site opertor acting on up triangle or down triangle. Used to compute the energy.
     *  \param[in] ifsvd Whether SVD is employee when finding isometries for CTM. 
     *  If it is set to false, QR decomposition is employee.
     *  \return Whether the CTM iteration is converged or not.
     */
    double ctmIteration( const int niter, const double errTol, const UniTensor<T> &upHam, bool ifsvd );

    //! Measure physical quantities.
    /*!
     *  \param[in] ham Three site hamiltonian acting on up or down triangle. 
     *  \param[out] energy Three site energy expectation value.
     *  \param[out] spinComp Spin component of each site. spinComp = vector<vector<double>> { S0, S1, S2}, where Si = vector<double> { Six, Siy, Siz}.
     *  \return Energy per-site.
     */
    double measureQuantities( const UniTensor<T> &ham, double &energy, vector<vector<double>> &spinComp  );

    //! Measure norm value.
    /*!
     *  \return norm value.
     */
    double measureNorm( const int isite );
    
    //! Measure one site expectation value.
    /*!
     *  \param[in] isite Site being measured. 0<=isite<=2.
     *  \return expectation value.
     */
    double measureOneSite( const UniTensor<T> &oneSiteOp, const int iSite, const int posit );

    //! Measure three site expectation value.
    /*!
     *  \param[in] isite Site being measured. 0<=isite<=2.
     *  \return expectation value.
     */
    double measureTriSite( const UniTensor<T> &threeSiteOp, const int posit  );

    //! Change the boundary bond dimension of CTM.
    /*!
     *  \param[in] newDim New boundary bond dimension.
     */
    void changeEdgeDim( const int newDim );
  private:
    //iden, ket, bra
    std::vector<uni10::UniTensor<T>> auxTens;
};

#endif
