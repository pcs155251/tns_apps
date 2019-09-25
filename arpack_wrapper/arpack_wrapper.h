#ifndef __ARPACK_WRAPPER_H_INCLUDED__
#define __ARPACK_WRAPPER_H_INCLUDED__

#include "uni10.hpp"

extern "C" {
void dsaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);

void dseupd_(int *rvec, char *howmny, int *select, double *d, double *z, int *ldz, double *sigma, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);

void dnaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);

void dneupd_(int *rvec, char *howmny, int *select, double *dr, double *di, double *z, int *ldz, double *sigmar, double *sigmai, double *workev, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);

void znaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, std::complex<double> *resid, int *ncv, std::complex<double> *v, int *ldv, int *iparam, int *ipntr, std::complex<double> *workd, std::complex<double> *workl, int *lworkl, double *rwork, int *info);

//void zneupd_(int *rvec, char *All, int *select, std::complex<double> *d, std::complex<double> *z, int *ldz, std::complex<double> *sigma, std::complex<double> *workev, char *bmat, int *n, char *which, int *nev, double *tol, std::complex<double> *resid, int *ncv, std::complex<double> *v, int *ldv, int *iparam, int *ipntr, std::complex<double> *workd, std::complex<double> *workl, int *lworkl, double *rwork, int *info);

void zneupd_(int *rvec, char *howmny, int *select, std::complex<double> *d, std::complex<double> *z, int *ldz, std::complex<double> *sigma, std::complex<double> *workev, char *bmat, int *n, char *which, int *nev, double *tol, std::complex<double> *resid, int *ncv, std::complex<double> *v, int *ldv, int *iparam, int *ipntr, std::complex<double> *workd, std::complex<double> *workl, int *lworkl, double *rwork, int *info);
}

int arpackEignRight( const uni10::Matrix<double> &A_m, std::complex<double>& EigVal, uni10::Matrix<std::complex<double>> &EigVec_m, const int max_iter, const double err_tol );

int arpackEignLeft( const uni10::Matrix<double> &A_m, std::complex<double>& EigVal, uni10::Matrix<std::complex<double>> &EigVec_m, const int max_iter, const double err_tol );

int arpackEignRight( const uni10::Matrix<std::complex<double>> &A_m, std::complex<double>& EigVal, uni10::Matrix<std::complex<double>> &EigVec_m, const int max_iter, const double err_tol );

int arpackEignLeft( const uni10::Matrix<std::complex<double>> &A_m, std::complex<double>& EigVal, uni10::Matrix<std::complex<double>> &EigVec_m, const int max_iter, const double err_tol );


int arRealEignRight( const uni10::Matrix<double> &A_m, double& EigVal, uni10::Matrix<double> &EigVec_m, const int max_iter, const double err_tol, const double imagLimit );

int arRealEignLeft( const uni10::Matrix<double> &A_m, double& EigVal, uni10::Matrix<double> &EigVec_m, const int max_iter, const double err_tol, const double imagLimit );

int arCompEignRight( const uni10::Matrix<std::complex<double>> &A_m, std::complex<double>& EigVal, uni10::Matrix<std::complex<double>> &EigVec_m, const int max_iter, double err_tol, const bool ifinit );

int arCompEignLeft( const uni10::Matrix<std::complex<double>> &A_m, std::complex<double>& EigVal, uni10::Matrix<std::complex<double>> &EigVec_m, const int max_iter, double err_tol, const bool ifinit );

void arpack_for_corner ( const uni10::UniTensor<double> &target, uni10::Matrix<double> &corner, unsigned int & max_iter, double err_tol );

/*
int arRealEignRightBMPS( double& EigVal, uni10::UniTensor<double> &EigVec_t, const int max_iter, double err_tol, const double imagLimit, uni10::Network &transfer_net, const uni10::UniTensor<double> &mpsUpL, const uni10::UniTensor<double> &mpsUpR, const uni10::UniTensor<double> &mpsDnR, const uni10::UniTensor<double> &mpsDnL, const uni10::UniTensor<double> &evolUp );
*/

int arEignRightBMPS( double& EigVal, uni10::UniTensor<double> &EigVec_t, const int max_iter, double err_tol, const double imagLimit, uni10::Network &transfer_net, const uni10::UniTensor<double> &mpsUpL, const uni10::UniTensor<double> &mpsUpR, const uni10::UniTensor<double> &mpsDnR, const uni10::UniTensor<double> &mpsDnL, const uni10::UniTensor<double> &evolUp );

int arEignRightBMPS( std::complex<double>& EigVal, uni10::UniTensor<std::complex<double>> &EigVec_t, const int max_iter, double err_tol, const double imagLimit, uni10::Network &transfer_net, const uni10::UniTensor<std::complex<double>> &mpsUpL, const uni10::UniTensor<std::complex<double>> &mpsUpR, const uni10::UniTensor<std::complex<double>> &mpsDnR, const uni10::UniTensor<std::complex<double>> &mpsDnL, const uni10::UniTensor<std::complex<double>> &evolUp );

int arRealEignRightBMPSCompress( double& EigVal, uni10::UniTensor<double> &EigVec_t, const int max_iter, double err_tol, const double imagLimit, uni10::Network &transfer_net, const uni10::UniTensor<double> &up0, const uni10::UniTensor<double> &up1, const uni10::UniTensor<double> &triangleUp, const uni10::UniTensor<double> &dn0, const uni10::UniTensor<double> &dn1, const uni10::UniTensor<double> &triangleDn );

int arCompEignRightBMPSCompress( std::complex<double> &EigVal, uni10::UniTensor<double> &EigVec_t, const int max_iter, double err_tol, const double imagLimit, uni10::Network &transfer_net, const uni10::UniTensor<double> &up0, const uni10::UniTensor<double> &up1, const uni10::UniTensor<double> &triangleUp, const uni10::UniTensor<double> &dn0, const uni10::UniTensor<double> &dn1, const uni10::UniTensor<double> &triangleDn );

int arRealEignRightFivePess( double& EigVal, uni10::UniTensor<double> &EigVec_t, const int max_iter, double err_tol, const double imagLimit, uni10::Network &transfer_net, const uni10::UniTensor<double> &mpsUpL, const uni10::UniTensor<double> &mpsUpR, const uni10::UniTensor<double> &mpsDnR, const uni10::UniTensor<double> &mpsDnL, const std::vector<uni10::UniTensor<double>> &evolUp );

void arRealSymEignRight( const uni10::Matrix<double> &amat, const int nEigVal, uni10::Matrix<double> &EigVal, uni10::Matrix<double> &EigVec, const int maxIter, const double errTol);

void arRealSVD( const uni10::Matrix<double> &amat, const int nEigVal, uni10::Matrix<double> &EigVal, uni10::Matrix<double> &EigVec, const int maxIter, const double errTol);


template<typename ... Args>
int arEignArgs( double& EigVal, uni10::UniTensor<double> &EigVec_t, int max_iter, double err_tol, double imagLimit, uni10::Network &transfer_net, const Args&... args){
  ///this routine find the max magnitude Eigenvalue and its Eigenvector with dnaupd and dneupd for a nonsymmetric matrix
  ///if one want to find other kind of or more Eigen pairs, tune which and nev.
  ///since the matrix is nonsymmetric, the Eigenpairs can be std::complex, so the input EigVal and EigVec_t should be std::complex.
  ///the std::complex Eigenvalue will appear with it's std::complex conjugate
  int dim = EigVec_t.ElemNum();
  EigVec_t.Randomize();
  int nev = 1;

  int ido = 0; ///reverse communication parameter, must be zero before iteration
  char bmat = 'I'; ///'I': standard Eigenproblem, 'G': generalized Eigenproblem 
  char which[] = {'L','M'}; ///type of asked Eigenvalues 
  int info = 0; double *resid = new double[dim];
  ///If info = 0, a randomly initial residual vector is used.
  ///If info = 1, resid contains the initial guess vector provided by user or from previous run
  ///on output, resid contains the final residual vector
  int ncv = 5; ///the number of Ritz vector, nev+2 <= ncv <= dim
  if (ncv>dim){
    ncv = dim;
  }
  else {}
  int ldv = dim; ///leading dimension of v
  double *v = new double[ldv*ncv];
  int iparam[11];
  iparam[0] = 1;        ///Specifies the shift strategy (1->exact)
  iparam[2] = max_iter; ///Maximum number of iterations
  iparam[6] = 1;        ///Sets the mode of dsaupd.
  int ipntr[14];
  double *workd = new double[3*dim];
  int lworkl = 3*ncv*(ncv+2); ///LWORKL must be at least 3*NCV**2 + 6*NCV .
  double *workl = new double[lworkl];
  ///Parameters for dgemv
  double alpha = 1.0e0; double beta = 0.0e0; int inc = 1;

  ///start iteration
  dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  unsigned int count = 1;
  while( ido!=99 ){ ///for mode 1 (iparam[0]==1), ido only  can be 1 or 99
    //dgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    double *EigVec_p = EigVec_t.GetElem();
    memcpy(EigVec_p, workd+ipntr[0]-1, dim*sizeof(double));
    uni10::UniTensor<double> tempt;
    ContractArgs( tempt, transfer_net, args..., EigVec_t );
    double *tempt_p = tempt.GetElem();
    memcpy(workd+ipntr[1]-1, tempt_p, dim*sizeof(double));
    dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    count++;
  }
  if( info!=0 ){
      std::cerr<<"error with dnaupd, info = "<<info<<std::endl;
  }
  else{}
  ///calculate Eigenvalue and Eigenvector by dneupd
  int rvec = 1; ///0: only Eigenvalue, 1: also Eigenvector
  char howmny = 'A'; ///how many Eigenvectors to calculate: 'A' => nev Eigenvectors
  int *select = new int[ncv];/// when howmny == 'A', this is used as workspace to reorder the Eigenvectors
  double *dr = new double[nev+1];
  double *di = new double[nev+1];
  double *workev = new double[3*ncv];
  double *z = new double[dim*(nev+1)];
  double sigmar, sigmai;
  dneupd_(&rvec, &howmny, select, dr, di, z, &ldv, &sigmar, &sigmai, workev, &bmat, &dim, which, &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  if ( info != 0 ){
      std::cerr<<"Error with dneupd, info = "<<info<<std::endl;
  }
  else {}

  ///extract results and informations
  EigVal = dr[0];
  //assert( fabs( di[0] ) < imagLimit );
  EigVec_t.SetElem( z );
  //for (int i=0; i!=dim; i++){
  //  EigVec_t[i] = z[i];
  //}
  ///free memory
  delete [] resid;
  delete [] v;
  delete [] workd;
  delete [] workl;
  delete [] select;
  delete [] workev;
  delete [] dr;
  delete [] di;
  delete [] z;
  return count;
}

template<typename ... Args>
int arEignArgs( double& EigVal, uni10::UniTensor<std::complex<double>> &EigVec_t, int max_iter, double err_tol, double imagLimit, uni10::Network &transfer_net, const Args&... args){
  ///this routine find the max magnitude Eigenvalue and its Eigenvector with dnaupd and dneupd for a nonsymmetric matrix
  ///if one want to find other kind of or more Eigen pairs, tune which and nev.
  ///since the matrix is nonsymmetric, the Eigenpairs can be std::complex, so the input EigVal and EigVec_t should be std::complex.
  ///the std::complex Eigenvalue will appear with it's std::complex conjugate
  int dim = EigVec_t.ElemNum();
  //EigVec_t.Randomize();
  int nev = 1;

  int ido = 0; ///reverse communication parameter, must be zero before iteration
  char bmat = 'I'; ///'I': standard Eigenproblem, 'G': generalized Eigenproblem 
  char which[] = {'L','M'}; ///type of asked Eigenvalues 
  int info = 0; double *resid = new double[dim];
  ///If info = 0, a randomly initial residual vector is used.
  ///If info = 1, resid contains the initial guess vector provided by user or from previous run
  ///on output, resid contains the final residual vector
  int ncv = 5; ///the number of Ritz vector, nev+2 <= ncv <= dim
  if (ncv>dim){
    ncv = dim;
  }
  else {}
  int ldv = dim; ///leading dimension of v
  double *v = new double[ldv*ncv];
  int iparam[11];
  iparam[0] = 1;        ///Specifies the shift strategy (1->exact)
  iparam[2] = max_iter; ///Maximum number of iterations
  iparam[6] = 1;        ///Sets the mode of dsaupd.
  int ipntr[14];
  double *workd = new double[3*dim];
  int lworkl = 3*ncv*(ncv+2); ///LWORKL must be at least 3*NCV**2 + 6*NCV .
  double *workl = new double[lworkl];
  ///Parameters for dgemv
  double alpha = 1.0e0; double beta = 0.0e0; int inc = 1;

  ///start iteration
  dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  unsigned int count = 1;
  while( ido!=99 ){ ///for mode 1 (iparam[0]==1), ido only  can be 1 or 99
    //dgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    uni10::UniTensor<double> tempt(EigVec_t.bond());
    double *temp_p = tempt.GetElem();
    memcpy(temp_p, workd+ipntr[0]-1, dim*sizeof(double));
    uni10::UniTensor<double> tempd;
    tempt.SetElem( workd+ipntr[0]-1 );
    ContractArgs( tempd, transfer_net, args..., tempt );
    double *tempd_p = tempd.GetElem();
    memcpy(workd+ipntr[1]-1, tempd_p, dim*sizeof(double));
    dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    count++;
    /*
    double *EigVec_p = EigVec_t.GetElem();
    memcpy(EigVec_p, workd+ipntr[0]-1, dim*sizeof(double));
    uni10::UniTensor<double> tempt;
    ContractArgs( tempt, transfer_net, args..., EigVec_t );
    double *tempt_p = tempt.GetElem();
    memcpy(workd+ipntr[1]-1, tempt_p, dim*sizeof(double));
    dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    count++;
    */
  }
  if( info!=0 ){
      std::cerr<<"error with dnaupd, info = "<<info<<std::endl;
  }
  else{}
  ///calculate Eigenvalue and Eigenvector by dneupd
  int rvec = 1; ///0: only Eigenvalue, 1: also Eigenvector
  char howmny = 'A'; ///how many Eigenvectors to calculate: 'A' => nev Eigenvectors
  int *select = new int[ncv];/// when howmny == 'A', this is used as workspace to reorder the Eigenvectors
  double *dr = new double[nev+1];
  double *di = new double[nev+1];
  double *workev = new double[3*ncv];
  double *z = new double[dim*(nev+1)];
  double sigmar, sigmai;
  dneupd_(&rvec, &howmny, select, dr, di, z, &ldv, &sigmar, &sigmai, workev, &bmat, &dim, which, &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  if ( info != 0 ){
      std::cerr<<"Error with dneupd, info = "<<info<<std::endl;
  }
  else {}

  ///extract results and informations
  EigVal = dr[0];
  //assert( fabs( di[0] ) < imagLimit );
  std::complex<double> *celem = new std::complex<double>[dim];
  for (int i=0; i!=dim; i++){
    celem[i].real(z[i]);
    celem[i].imag(z[dim+i]);
  }
  EigVec_t.SetElem( celem );
  delete [] celem;
  //free memory
  delete [] resid;
  delete [] v;
  delete [] workd;
  delete [] workl;
  delete [] select;
  delete [] workev;
  delete [] dr;
  delete [] di;
  delete [] z;
  return count;
}

template<typename ... Args>
int arEignArgs( std::complex<double>& EigVal, uni10::UniTensor<std::complex<double>> &EigVec_t, int max_iter, double err_tol, double imagLimit, uni10::Network &transfer_net, const Args&... args){
//int arCompEignRight( const Matrix<std::complex<double>> &A_m, std::complex<double>& EigVal, Matrix<std::complex<double>> &EigVec_m, const int max_iter, double err_tol, const bool ifinit ){
  //this routine find the max magnitude Eigenvalue and its Eigenvector with znaupd and zneupd for a nonhermitian matrix
  //if one want to find other kind of or more Eigen pairs, tune which and nev.
  int dim = EigVec_t.ElemNum();
  EigVec_t.Randomize();
  int nev = 1;

  int ido = 0; ///reverse communication parameter, must be zero before iteration
  char bmat = 'I'; ///'I': standard Eigenproblem, 'G': generalized Eigenproblem 
  char which[] = {'L','M'}; ///type of asked Eigenvalues 
  int info=0; std::complex<double> *resid = new std::complex<double> [dim];
  ///If info = 0, a randomly initial residual vector is used.
  ///If info = 1, resid contains the initial guess vector provided by user or from previous run
  ///on output, resid contains the final residual vector
  int ncv = 5; ///the number of Ritz vector, nev+2 <= ncv <= dim
  if (ncv>dim){
    ncv = dim;
  } else {}
  int ldv = dim; ///leading dimension of v
  std::complex<double> *v = new std::complex<double>[dim*ncv];
  int iparam[11];
  iparam[0] = 1;        ///Specifies the shift strategy (1->exact)
  iparam[2] = max_iter; ///Maximum number of iterations
  iparam[3] = 1;        //current code only for 1
  iparam[6] = 1;        ///Sets the mode of znaupd.
  int ipntr[14];
  std::complex<double> *workd = new std::complex<double>[3*dim];
  int lworkl = 3*ncv*(ncv+2); ///LWORKL must be at least 3*NCV**2 + 6*NCV .
  std::complex<double> *workl = new std::complex<double>[lworkl];
  double *rwork = new double[ncv];

  ///start iteration
  //dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  znaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
  unsigned int count = 1;
  ///Parameters for zgemv
  std::complex<double> alpha( 1.0, 0); std::complex<double> beta( 0, 0); int inc = 1;
  while( ido!=99 ){ ///for mode 1 (iparam[0]==1), ido only  can be 1 or 99
    //zgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    std::complex<double> *EigVec_p = EigVec_t.GetElem();
    memcpy(EigVec_p, workd+ipntr[0]-1, dim*sizeof(std::complex<double>));
    uni10::UniTensor<std::complex<double>> tempt;
    ContractArgs( tempt, transfer_net, args..., EigVec_t );
    std::complex<double> *tempt_p = tempt.GetElem();
    memcpy(workd+ipntr[1]-1, tempt_p, dim*sizeof(std::complex<double>));
    znaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
    count++;
  }
  if( info!=0 ){
      std::cerr<<"error with znaupd, info = "<<info<<std::endl;
  } else{}
  ///calculate Eigenvalue and Eigenvector by dneupd
  int rvec = 1; ///0: only Eigenvalue, 1: also Eigenvector
  char howmny = 'A'; ///how many Eigenvectors to calculate: 'A' => nev Eigenvectors
  int *select = new int[ncv];/// when howmny == 'A', this is used as workspace to reorder the Eigenvectors
  std::complex<double> *d = new std::complex<double>[nev+1];
  std::complex<double> *z = new std::complex<double>[dim*nev];
  std::complex<double> sigma;
  std::complex<double> *workev = new std::complex<double>[2*ncv];
  //dneupd_(&rvec, &howmny, select, dr, di, z, &ldv, &sigmar, &sigmai, workev, &bmat, &dim, which, &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  zneupd_(&rvec, &howmny, select, d, z, &ldv, &sigma, workev, &bmat, &dim, which, &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
  if ( info != 0 ){
      std::cerr<<"error with zneupd, info = "<<info<<std::endl;
  }
  else {}
  ///extract results and informations
  EigVal = d[0];
  EigVec_t.SetElem( z );
  ///free memory
  delete [] resid;
  delete [] v;
  delete [] workd;
  delete [] workl;
  delete [] rwork;
  delete [] select;
  delete [] d;
  delete [] z;
  delete [] workev;
  return count;
}

/*
//void arRealSymEignRight( const Matrix<double> &amat, const int nEigVal, Matrix<double> &EigVal, Matrix<double> &EigVec, const int maxIter, const double errTol){
template<typename ... Args>
void arRealSymEignArgs( const int dim, const int nEigVal, Matrix<double> &EigVal, Matrix<double> &EigVec, const int maxIter, const double errTol, uni10::Network &transfer_net, const Args&... args){
  ///should specify eigenvector dimension and initialize before calling 
  ///For symmetric matrix this function finds given number of Eigenvalue and its Eigenvector with dsaupd and dseupd.
  int ido = 0; 
  char bmat = 'I';
  char which[] = {'L','M'}; 
  int nev = nEigVal;
  double tol = errTol;
  double *resid = new double[dim];
  int ncv = dim; 
  double *v = new double[dim*ncv];
  int ldv = dim; 

  int iparam[11];
  iparam[0] = 1;
  iparam[2] = maxIter; 
  iparam[6] = 1;

  int ipntr[11];
  double *workd = new double[3*dim];
  int lworkl = 2*ncv*(ncv+8);
  double *workl = new double[lworkl];
  int info = 0;

  ///Parameters for dgemv
  double alpha = 1.0e0; double beta = 0.0e0; int inc = 1;

  ///start iteration
  dsaupd_(&ido, &bmat, &dim, &which[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  while( ido!=99 ){ ///for mode 1 (iparam[0]==1), ido only  can be 1 or 99
    ///
    //dgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    std::complex<double> *EigVec_p = EigVec_t.GetElem();
    memcpy(EigVec_p, workd+ipntr[0]-1, dim*sizeof(std::complex<double>));
    uni10::UniTensor<std::complex<double>> tempt;
    ContractArgs( tempt, transfer_net, args..., EigVec_t );
    std::complex<double> *tempt_p = tempt.GetElem();
    memcpy(workd+ipntr[1]-1, tempt_p, dim*sizeof(std::complex<double>));
    ///
    dsaupd_(&ido, &bmat, &dim, &which[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  }
  if( info!=0 ){
      cerr<<"error with dsaupd, info = "<<info<<endl;
  }
  else{}

  ///calculate Eigenvalue and Eigenvector by dneupd
  int rvec = 1; ///0: only Eigenvalue, 1: also Eigenvector
  char howmny = 'A'; ///how many Eigenvectors to calculate: 'A' => nev Eigenvectors
  int *select = new int[ncv];/// when howmny == 'A', this is used as workspace to reorder the Eigenvectors
  double *d = new double[nev];
  double *z = new double[dim*nev];
  double sigma;
  int ldz = dim;
  //memcpy( z, v, dim*nev*sizeof(double) );
  dseupd_(&rvec, &howmny, select, d, z, &ldz, &sigma, &bmat, &dim, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  if ( info != 0 ){
      cerr<<"Error with dseupd, info = "<<info<<endl;
  }
  else {}
  delete [] resid;
  delete [] workd;
  delete [] v;
  delete [] workl;
  delete [] select;
  delete [] d;
  delete [] z;

  ///extract results and informations
  double *finald = new double[nev];
  EigVal = Matrix<double> (nEigVal, nEigVal, true);
  delete [] finald;

  double *finalz = new double[dim*nev];
  for (int i=0; i!=nev; i++){
    finald[i] = d[nev-i-1];
  }
  for (int i=0; i!=dim; i++){
    for (int j=0; j!=nev; j++){
      finalz[i*nev+j] = z[(nev-j-1)*dim+i];
    }
  }
  EigVec = Matrix<double> (    dim, nEigVal, false);
  EigVec.SetElem(finalz);
  delete [] finalz;
}
*/

std::vector<uni10::Matrix<double>> arRealSvd( const uni10::Matrix<double> &amat, const int nRankIn, const int maxIter, const double errTol );

#endif
