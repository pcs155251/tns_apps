#ifndef __PARAPESS_H_INCLUDED__
#define __PARAPESS_H_INCLUDED__

#include <sys/stat.h> //for mkdir
#include "../nsy_opt/operator.h"
#include "../parser/Parser.h"
#include <vector>
#include <string>

//! input parameters for PESS
/*
 *  Example of input files: *.rc
 *  Variables without documentation are auxiliary variables cannot be specified and used in input *.rc files.
 */
struct paraPess{
  //! Concstructor
  paraPess( int argc, char* argv [] );
  //! Root folder
  std::string rootFolder
  //! Folder where wave function is loaded.
  std::string loadFolder;
  //! Folder where environment tensors are loaded.
  std::string enviFolder;
  //! Folder where data and tensors are saved.
  std::string saveFolder;
  //! Whether find ground state or not.
  int ifGs;
  //! Whether find environment or not.
  int ifEn;
  //! Maximum iteration step of CTM.
  int ctmStep;
  //! Number of boundary bond dimension. The boundary bond dimension is increase as (start+slope*i)*dimD*dimD, where i in range (0,edgeDimN)
  int edgeDimN;
  double start;
  double slope;
  //! Whether use svd to find isometries in CTM or not.
  int ifsvd;
  //! Maximum iteration step of simple update.
  int gsIter;
  //! DM strength
  double field;
  //! Criterion of convergence with respect to boundary bond dimension.
  double edgeErr;
  //! Criterion of convergence of imaginary time evolution.
  double gsErr;
  //! Virtual bond dimensions being simulated
  std::vector<int> dimDs;
  //! A sequence of time step.
  std::vector<double> taus;
  int dimDload;
  std::string tensorFolder, dataFolder, fieldSubFolder;
};

#endif
