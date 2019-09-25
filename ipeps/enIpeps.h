#ifndef __ENIPEPS_H_INCLUDED__
#define __ENIPEPS_H_INCLUDED__

#include <sys/stat.h> //for mkdir
#include <ctime>
#include <complex.h> 
#include <fstream>
#include "uni10.hpp"
#include "../ctmBase/ctmBase.h"
#include "../tools/general_tools.h"
#include "../nsy_opt/operator.h"

template<typename T>
class ipeps;

//! A derived class of ctmBase describing the environment of 2D iPEPS on square lattice.
/*!
 *  By default the wave function has 2 by 2 unit cell with checkerboard pattern. 
 *  enIpeps is used in ipeps to measure physical quantities and update the wave function.
 *  Methods defined in ctmBase are employee to find the converged environment.
 */
template<typename T>
class enIpeps: public ctmBase<T>
{
  public:
    //! Constructor
    /*! \param[in] unit Tensors consist of unit cell. By default we assume 2 by 2 unit cell. The tensors should be stroed from left to right, top to bottom.
     *  \param[in] edgeDimIn Boundary bond dimension.
     */  
    enIpeps( const std::vector<uni10::UniTensor<T>>& unit , const int edgeDimIn );
    //! Measure norm of one site environment
    /*!
     *  \param[in] gammas iPEPS wave function ansatz.
     *  \param[in] isite Site being measured (from 0~3). 
     *  \return norm value.
     */
    double meaOneSiteNorm( const std::vector<uni10::UniTensor<T>>& gammas, int isite ) const;
    //! Measure norm of two site environment
    /*!
     *  \param[in] gammas iPEPS wave function ansatz.
     *  \param[in] ibond Bond being measured (from 0~3).
     *  \return norm value.
     */
    double meaTwoSiteNorm( const std::vector<uni10::UniTensor<T>>& gammas, int ibond ) const;
    //! Measure expecatation value of two one site operator.
    /*!
     *  \param[in] gammas iPEPS wave function ansatz.
     *  \param[in] oneSiteOp One site operator
     *  \param[in] isite Site being measured (from 0~3).
     *  \return expectation value.
     */
    double meaOneSiteExp ( const std::vector<uni10::UniTensor<T>>& gammas, uni10::UniTensor<T> &oneSiteOp, int isite ) const;
    //! Measure expecatation value of two site environment
    /*!
     *  \param[in] gammas iPEPS wave function ansatz.
     *  \param[in] twoSiteOp Two site operator
     *  \param[in] ibond Bond being measured (from 0~3).
     *  \return expectation value.
     */
    double meaTwoSiteExp ( const std::vector<uni10::UniTensor<T>>& gammas, uni10::UniTensor<T> &twoSiteOp, int ibond ) const;
    //! Iteratively find the environment which can approximate the 2D infinite TNS. 
    //  ctmrgAstep() in ctmBase is employee to perform each RG step.
    /*!
     *  \param[in] gammas iPEPS wave function ansatz.
     *  \param[in] niter Maximum iteration step.
     *  \param[in] errTol Convergence criterion. If the measured energy between two step is smaller than errTol, the environment is considered converged.
     *  \param[in] ham Two site hamiltonian. It is used to create the evolution operator and measure energy.
     *  \param[out] energy Energy measured in the final step.
     *  \param[in] ifprint Whether detailed data is printed or not.
     *  \param[in] ifsvd Whether svd is employee to find the isometries in CTM procedure or not.
     *  \return Whether the iteration is converged or not.
     */
    bool ctmIteration( const std::vector<uni10::UniTensor<T>>& gammas, int niter, double errTol, uni10::UniTensor<T> &ham, double &energy, bool ifprint, bool ifsvd );
    //void getOneSiteEn( int isite ) const;
    friend class ipeps<T>;
  private:
};

#endif
