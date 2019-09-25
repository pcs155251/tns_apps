#ifndef __CTMBASE_H_INCLUDED__
#define __CTMBASE_H_INCLUDED__

#include <sys/stat.h> //for mkdir
#include <ctime>
#include <complex.h> 
#include <fstream>
#include "uni10.hpp"
#include "../tools/general_tools.h"
#include "../nsy_opt/operator.h"

template<typename T>
class ipeps;

//!A data structure to stored a tensor and its environemnt ( 4 corners and 4 edges ).
/*!
 * This class is used in ctmBase to stored the a tensor obtained from tracing out the physical bond (reduced tensor) of TNS wave function.
 * There are 9 tensors in total, the reduced tensor, and 4 corner tensors (start form upper left), and 4 edge tensors (start upper edge)
 * The accuracy of environment is tuned by the boundary bond dimension edgeDim.
*/
template<typename T>
struct nineTens{
  //! Constructor.
  /*!
   *  \param[in] coreIn Core tensor. The core tensor is obtain by tracing out the physical bond of TNS wave function.
   *  \param[in] edgeDimIn Boundary bond dimension.
  */
  nineTens( const uni10::UniTensor<T> &coreIn, const int edgeDimIn );
  //! A tensor within unit cell. 
  /*! 
   * Obtained by tracing out the physical bond of TNS wave function.
  */
  uni10::UniTensor<T> core;//bond start from left
  //! Corner tensors (C) in ref. 
  /*!
   *  There are 4 corner tensors, starting from upper left corner and oriented in counterclockwise order.
  */
  std::vector<uni10::UniTensor<T>> corners;//tensor start from upper left
  //! Edge tensors (T) in ref.
  /*!
   *  There are 4 edge tensors, starting from upper edge and oriented in counterclockwise order.
  */
  std::vector<uni10::UniTensor<T>> edges;//tensor start from top
  //! Normalize all tensors.
  /*! 
   *  Each tensors is normalized such that the maximum element is unity.
   */
  void normalize();
  //! Print tensor diagrams of all tensors.
  void printAll() const;
  //! Change boundary bond dimension
  /*! 
   *  If newDim is smaller than the original boundary bond dimension, the excessive element will be truncated. 
   *  If newDim is larger than the original one, the empty elements will be filled by zeros.
   */  
  void changeEdgeDim( const int newDim );
  //! Save tensors into given folder.
  /*!
   *  \param[in] folder The folder where tensor is stored.
   */
  void saveTensors( const std::string &folder ) const;
  //! Load tensors from given folder.
  /*!
   *  \param[in] folder The folder where tensor is load.
   */
  bool loadTensors( const std::string &folder );
  //corner bond: clockwise
  //edge bond: clockwise, but put inside bond the last
  //! boundary bond dimension
  int edgeDim;
  //UniTensor<T> initCorner( const UniTensor<T>& core, int idir, int edgeDim );
  //UniTensor<T> initEdge( const UniTensor<T>& core, int idir, int edgeDim );
};

//!Base class for CTM. 
/*!This class handle the convergence of corner transfer matrix.
 * In principal 2D TNS with translation invariance can be represented by repeated unit cell consist of reduced tensors.
 * For each reduced tensor within the unit cell, a nineTens is defined to stored the reduced tensor and its environment.
 * I implemented the CTM method to iteratively find the environment which can approximate the effect of the 2D infinite TNS.
 * One can define a derived class for computing environment for specific TNS ansatz or lattice.
 * For example, enIpeps is defined for ipeps on square lattice and enPess is defined for pess on kagome lattice.
 * enRedPess is defined for dimension-reduced environment for pess on kagome lattice.
 * Since the measurement procedures does vary for different ansatz or lattice, one should defined the measurement method 
 * and ctm iteration function in the derived class, as which is done in enIpeps, enPess and enRedPess.
*/
template<typename T>
class ctmBase{
  public:
    //! Constructor.
    /*!
     *  \param[in] nrowIn Number of rows of unit cell.
     *  \param[in] ncolIn Number of columns of unit cell.
     *  \param[in] edgeDimIn Boundary bond dimension of CTM. 
    */
    ctmBase( const int nrowIn, const int ncolIn, const int edgeDimIn );
    //! Constructor.
    /*!
     *  \param[in] nrowIn Number of rows of unit cell.
     *  \param[in] ncolIn Number of columns of unit cell.
     *  \param[in] unit Tensors consist of the unit cell. 
     *  The tensors stored in unit cell should be specified form left to right and up to down in the unit cell.
     *  \param[in] edgeDimIn Boundary bond dimension of CTM. 
    */
    ctmBase( const int nrowIn, const int ncolIn, const std::vector<uni10::UniTensor<T>>& unit , const int edgeDimIn );
    //virtual void ctmIteration( const int niter, const double errTol, uni10::UniTensor<T> &ham ) = 0;
    //! Change the boundary bond dimension of the CTM.
    /*!
     *  If newDim is smaller than the original boundary bond dimension, the excessive element will be truncated. 
     *  If newDim is larger than the original one, the empty elements will be filled by zeros.
     *  \param[in] newDim new boundary bond dimension.
    */
    void changeEdgeDim( const int newDim );
    //! Save tensors into a given path.
    /*!
     *  \param[in] folder Path to save the tensors.
    */
    void saveTensors( const std::string &folder ) const;
    //! Load tensors from a given path.
    /*!
     *  \param[in] folder Path to load the tensors.
     *  \return If the tensors are being saved successfully.
    */
    bool loadTensors( const std::string &folder );
    //! Perform a step of RG, i.e. left move, right move, up move, and down move successively.
    /*!
     *  \param[in] ifsvd If use svd to find isometries. If it is set to false, then QR deomposition is employee to find isometries.
    */
    void ctmrgAstep( bool ifsvd );
    //! Set the tensors within the unit cell and it's corresponding environment tensor (stored as nineTens).
    /*!
     *  The environment tensors will be reinitialized.
     *  \param[in] unit Tensors within unit cell.
    */
    void setGroups( const std::vector<uni10::UniTensor<T>> &unit );
    //! Set the tensors within the unit cell and it's corresponding environment tensor (stored as nineTens).
    /*!
     *  The environment tensors will not be reinitialized.
     *  \param[in] unit Tensors within unit cell.
    */
    void updateUnit( const std::vector<uni10::UniTensor<T>> &unit );
    //! Print the tensor diagrams for all the tensors.
    void printAll();
  private:
    const int nrow;
    const int ncol;
    std::vector<uni10::UniTensor<T>> findIsometrySVD( std::vector<int> order, const int idir ); 
    std::vector<uni10::UniTensor<T>> findIsometryQR ( std::vector<int> order, const int idir ); 
    void move( const int idir, bool ifsvd );
  protected:
    //! Boundary bond dimension.
    int ctmDim;
    //! Tensors consist of unit cell and its environemnts (4 C's and 4 T's). Each stored as a nineTens.
    std::vector<nineTens<T>> groups;//tensor from left to right, top to down
};

#endif
