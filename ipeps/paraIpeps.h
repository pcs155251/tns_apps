#ifndef __PARAIPEPS_INCLUDED__
#define __PARAIPEPS_INCLUDED__

#include <cmath>
#include <sys/stat.h> //for mkdir
#include "../parser/Parser.h"

//! input parameters for iPEPS
/*!
 *  Example of input files: *.rc
 *  Variables without documentation are auxiliary variables cannot be specified and used in input *.rc files.
 */
struct paraIpeps
{
  //! Constructor
  paraIpeps( int argc, char **argv );
  //! Whether load wave function or not.
  int ifloadWv;
  //! Whether perform simple update or not.
  int ifsu;
  //! Whether perform full update or not.
  int iffu;
  //! Whether save environment or not.
  int ifsaveEn;
  //! Virtural bond dimension.
  int virDim;
  //! Boundary bond dimension of loaded environment
  int edgeDimL;
  //! Boundary bond dimension of current environment 
  int edgeDimU;
  //! Maximum iteration step of simple update
  int suIter;
  //! Maximum iteration step of full update
  int fuIter;
  //! Maximum iteration step of CTM
  int ctmIter;
  //! Convergence criterion of simple update
  double suStop;
  //! Convergence criterion of full update
  double fuStop;
  //! Convergence criterion of CTM iteration
  double ctmStop;
  //! Sequence of time steps for imaginary time evolution
  std::vector<double> taus;
  //! Transverse fields wish to simulate
  std::vector<double> fields;
  //! Root folder
  std::string rootFolder
  //! Folder where wave function is loaded
  std::string wvLoadFolder;
  //! Folder where environment tensors are loaded
  std::string enLoadFolder;
  //! Folder where data and tensors are saved.
  std::string saveFolder;

  std::string tensorFolder, dataFolder;
  double fieldStart, fieldEnd, fieldDelta, fieldLoad;
};

#endif
