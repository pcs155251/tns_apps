#include "arpack_wrapper.h"

using namespace std;
using namespace uni10;

int arpackEignRight( const Matrix<complex<double>> &A_m, complex<double>& EigVal, Matrix<complex<double>> &EigVec_m, const int max_iter, double err_tol ){
  const bool ifinit=false;
  //this routine find the max magnitude Eigenvalue and its Eigenvector with znaupd and zneupd for a nonhermitian matrix
  //if one want to find other kind of or more eigen pairs, tune which and nev.
  assert( A_m.col() == A_m.row() );
  complex<double> *A = A_m.GetElem();
  int nev = 1;

  int dim = A_m.row(); ///dimension of Eigenproblem
  int ido = 0; ///reverse communication parameter, must be zero before iteration
  char bmat = 'I'; ///'I': standard Eigenproblem, 'G': generalized Eigenproblem 
  char which[] = {'L','M'}; ///type of asked Eigenvalues 
  int info; complex<double> *resid = new complex<double> [dim];
  if (ifinit){
    info = 1;
    memcpy( resid, EigVec_m.GetElem(), dim*sizeof( complex<double> ) );
    EigVec_m = Matrix<complex<double>> ( A_m.row(), 1);
  } 
  else {
    info = 0;
    EigVec_m = Matrix<complex<double>> ( A_m.row(), 1);
  }
  ///If info = 0, a randomly initial residual vector is used.
  ///If info = 1, resid contains the initial guess vector provided by user or from previous run
  ///on output, resid contains the final residual vector
  int ncv = dim; ///the number of Ritz vector, nev+2 <= ncv <= dim
  int ldv = dim; ///leading dimension of v
  complex<double> *v = new complex<double>[dim*ncv];
  int iparam[11];
  iparam[0] = 1;        ///Specifies the shift strategy (1->exact)
  iparam[2] = max_iter; ///Maximum number of iterations
  iparam[3] = 1;        //current code only for 1
  iparam[6] = 1;        ///Sets the mode of znaupd.
  int ipntr[14];
  complex<double> *workd = new complex<double>[3*dim];
  int lworkl = 3*ncv*(ncv+2); ///LWORKL must be at least 3*NCV**2 + 6*NCV .
  complex<double> *workl = new complex<double>[lworkl];
  double *rwork = new double[ncv];

  ///start iteration
  //dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  znaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
  unsigned int count = 1;
  ///Parameters for zgemv
  complex<double> alpha( 1.0, 0); complex<double> beta( 0, 0); int inc = 1;
  while( ido!=99 ){ ///for mode 1 (iparam[0]==1), ido only  can be 1 or 99
    //dgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    //dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    zgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    znaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
    count++;
  }
  if( info!=0 ){
      cerr<<"error with znaupd, info = "<<info<<endl;
  } else{}
  ///calculate Eigenvalue and Eigenvector by dneupd
  int rvec = 1; ///0: only Eigenvalue, 1: also Eigenvector
  char howmny = 'A'; ///how many Eigenvectors to calculate: 'A' => nev Eigenvectors
  int *select = new int[ncv];/// when howmny == 'A', this is used as workspace to reorder the Eigenvectors
  complex<double> *d = new complex<double>[nev+1];
  complex<double> *z = new complex<double>[dim*nev];
  complex<double> sigma;
  complex<double> *workev = new complex<double>[2*ncv];
  //dneupd_(&rvec, &howmny, select, dr, di, z, &ldv, &sigmar, &sigmai, workev, &bmat, &dim, which, &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  zneupd_(&rvec, &howmny, select, d, z, &ldv, &sigma, workev, &bmat, &dim, which, &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
  if ( info != 0 ){
      cerr<<"error with zneupd, info = "<<info<<endl;
  }
  else {}
  ///extract results and informations
  EigVal = d[0];
  EigVec_m.SetElem( z );
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

int arpackEignRight( const Matrix<double> &A_m, complex<double>& EigVal, Matrix<complex<double>> &EigVec_m, const int max_iter, double err_tol ){
  const double imagLimit = 1.0e-13;
  ///this routine find the max magnitude Eigenvalue and its Eigenvector with dnaupd and dneupd for a nonsymmetric matrix
  ///if one want to find other kind of or more Eigen pairs, tune which and nev.
  ///since the matrix is nonsymmetric, the Eigenpairs can be complex, so the input EigVal and EigVec_m should be complex.
  ///the complex Eigenvalue will appear with it's complex conjugate
  assert( A_m.col() == A_m.row() );
  double *A = A_m.GetElem();
  EigVec_m = Matrix<double> ( A_m.row(), 1);
  int nev = 1;

  int dim = A_m.row(); ///dimension of Eigenproblem
  int ido = 0; ///reverse communication parameter, must be zero before iteration
  char bmat = 'I'; ///'I': standard Eigenproblem, 'G': generalized Eigenproblem 
  char which[] = {'L','R'}; ///type of asked Eigenvalues 
  int info = 0; double *resid = new double[dim];
  ///If info = 0, a randomly initial residual vector is used.
  ///If info = 1, resid contains the initial guess vector provided by user or from previous run
  ///on output, resid contains the final residual vector
  int ncv = dim; ///the number of Ritz vector, nev+2 <= ncv <= dim
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
    dgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    count++;
  }
  if( info!=0 ){
      cerr<<"error with dnaupd, info = "<<info<<endl;
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
      cerr<<"Error with dneupd, info = "<<info<<endl;
  }
  else {}

  ///extract results and informations
  EigVal = complex<double> (dr[0],di[0]);
  if (di[0]<1.0e-14){
    for (int i=0; i!=dim; i++){
      EigVec_m[i] = z[i];
    }
  }
  else {
    for (int i=0; i!=dim; i++){
      EigVec_m[i] = complex<double> (z[i],z[dim+i]);
    }
  }
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

int arpackEignLeft( const Matrix<complex<double>> &A_m, complex<double>& EigVal, Matrix<complex<double>> &EigVec_m, const int max_iter, double err_tol ){
  Matrix<complex<double>> at_m = Transpose( A_m );
  int count = arpackEignRight( at_m, EigVal, EigVec_m, max_iter, err_tol );
  EigVec_m = Transpose( EigVec_m );
  return count;
}

int arpackEignLeft( const Matrix<double> &A_m, complex<double>& EigVal, Matrix<complex<double>> &EigVec_m, const int max_iter, double err_tol ){
  Matrix<double> at_m = Transpose( A_m );
  int count = arpackEignRight( at_m, EigVal, EigVec_m, max_iter, err_tol );
  EigVec_m = Transpose( EigVec_m );
  return count;
}

int arCompEignRight( const Matrix<complex<double>> &A_m, complex<double>& EigVal, Matrix<complex<double>> &EigVec_m, const int max_iter, double err_tol, const bool ifinit ){
  //this routine find the max magnitude Eigenvalue and its Eigenvector with znaupd and zneupd for a nonhermitian matrix
  //if one want to find other kind of or more Eigen pairs, tune which and nev.
  assert( A_m.col() == A_m.row() );
  complex<double> *A = A_m.GetElem();
  int nev = 1;

  int dim = A_m.row(); ///dimension of Eigenproblem
  int ido = 0; ///reverse communication parameter, must be zero before iteration
  char bmat = 'I'; ///'I': standard Eigenproblem, 'G': generalized Eigenproblem 
  char which[] = {'L','M'}; ///type of asked Eigenvalues 
  int info; complex<double> *resid = new complex<double> [dim];
  if (ifinit){
    info = 1;
    memcpy( resid, EigVec_m.GetElem(), dim*sizeof( complex<double> ) );
    EigVec_m = Matrix<complex<double>> ( A_m.row(), 1);
  } 
  else {
    info = 0;
    EigVec_m = Matrix<complex<double>> ( A_m.row(), 1);
  }
  ///If info = 0, a randomly initial residual vector is used.
  ///If info = 1, resid contains the initial guess vector provided by user or from previous run
  ///on output, resid contains the final residual vector
  int ncv = int( sqrt( dim ) ); ///the number of Ritz vector, nev+2 <= ncv <= dim
  if (ncv>dim){
    ncv = dim;
  } else {}
  int ldv = dim; ///leading dimension of v
  complex<double> *v = new complex<double>[dim*ncv];
  int iparam[11];
  iparam[0] = 1;        ///Specifies the shift strategy (1->exact)
  iparam[2] = max_iter; ///Maximum number of iterations
  iparam[3] = 1;        //current code only for 1
  iparam[6] = 1;        ///Sets the mode of znaupd.
  int ipntr[14];
  complex<double> *workd = new complex<double>[3*dim];
  int lworkl = 3*ncv*(ncv+2); ///LWORKL must be at least 3*NCV**2 + 6*NCV .
  complex<double> *workl = new complex<double>[lworkl];
  double *rwork = new double[ncv];

  ///start iteration
  //dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  znaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
  unsigned int count = 1;
  ///Parameters for zgemv
  complex<double> alpha( 1.0, 0); complex<double> beta( 0, 0); int inc = 1;
  while( ido!=99 ){ ///for mode 1 (iparam[0]==1), ido only  can be 1 or 99
    //dgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    //dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    zgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    znaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
    count++;
  }
  if( info!=0 ){
      cerr<<"error with znaupd, info = "<<info<<endl;
  } else{}
  ///calculate Eigenvalue and Eigenvector by dneupd
  int rvec = 1; ///0: only Eigenvalue, 1: also Eigenvector
  char howmny = 'A'; ///how many Eigenvectors to calculate: 'A' => nev Eigenvectors
  int *select = new int[ncv];/// when howmny == 'A', this is used as workspace to reorder the Eigenvectors
  complex<double> *d = new complex<double>[nev+1];
  complex<double> *z = new complex<double>[dim*nev];
  complex<double> sigma;
  complex<double> *workev = new complex<double>[2*ncv];
  //dneupd_(&rvec, &howmny, select, dr, di, z, &ldv, &sigmar, &sigmai, workev, &bmat, &dim, which, &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  zneupd_(&rvec, &howmny, select, d, z, &ldv, &sigma, workev, &bmat, &dim, which, &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
  if ( info != 0 ){
      cerr<<"error with zneupd, info = "<<info<<endl;
  }
  else {}
  ///extract results and informations
  EigVal = d[0];
  EigVec_m.SetElem( z );
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

int arRealEignRight( const Matrix<double> &A_m, double& EigVal, Matrix<double> &EigVec_m, const int max_iter, double err_tol, const double imagLimit ){
  ///this routine find the max magnitude Eigenvalue and its Eigenvector with dnaupd and dneupd for a nonsymmetric matrix
  ///if one want to find other kind of or more Eigen pairs, tune which and nev.
  ///since the matrix is nonsymmetric, the Eigenpairs can be complex, so the input EigVal and EigVec_m should be complex.
  ///the complex Eigenvalue will appear with it's complex conjugate
  assert( A_m.col() == A_m.row() );
  double *A = A_m.GetElem();
  EigVec_m = Matrix<double> ( A_m.row(), 1);
  int nev = 1;

  int dim = A_m.row(); ///dimension of Eigenproblem
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
    dgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    count++;
  }
  if( info!=0 ){
      cerr<<"error with dnaupd, info = "<<info<<endl;
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
      cerr<<"Error with dneupd, info = "<<info<<endl;
  }
  else {}

  ///extract results and informations
  EigVal = dr[0];
  assert( fabs( di[0] ) < imagLimit );
  for (int i=0; i!=dim; i++){
    EigVec_m[i] = z[i];
  }
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

int arCompEignLeft( const Matrix<complex<double>> &A_m, complex<double>& EigVal, Matrix<complex<double>> &EigVec_m, const int max_iter, double err_tol, const bool ifinit ){
  Matrix<complex<double>> at_m = Transpose( A_m );
  int count = arCompEignRight( at_m, EigVal, EigVec_m, max_iter, err_tol, ifinit );
  EigVec_m = Transpose( EigVec_m );
  return count;
}

int arRealEignLeft( const Matrix<double> &A_m, double& EigVal, Matrix<double> &EigVec_m, const int max_iter, double err_tol, const double imagLimit ){
  Matrix<double> at_m = Transpose( A_m );
  int count = arRealEignRight( at_m, EigVal, EigVec_m, max_iter, err_tol, imagLimit );
  EigVec_m = Transpose( EigVec_m );
  return count;
}

void arpack_for_corner ( const UniTensor<double> &target, Matrix<double> &corner, unsigned int & max_iter, double err_tol ){
  
  Network Eigen_net("./Networks/operation_for_corner.net");

  int dim = target.bond(0).dim();
  int n = dim*dim*dim*dim;
  
  vector<Bond> corner_bds = { Bond( BD_IN, dim ), Bond( BD_OUT, dim) };
  UniTensor<double> cornerTs (corner_bds);

  int ido = 0;
  char bmat = 'I';
  char which[] = {'L','M'};
  int info = 0; 
  double *resid = new double[n];
  
  int nev = 1;
  int ncv = 5;
  if (ncv > n) ncv = n;
  
  int ldv = n;  
  double *v = new double[n*ncv];
  int iparam[11];
  iparam[0] = 1;
  iparam[2] = max_iter;
  iparam[6] = 1;
  int ipntr[14];
  double *workd = new double[3*n];
  int lworkl = 3*ncv*(ncv+2);
  double *workl = new double[lworkl];

  double alpha = 1.0e0; double beta = 0.0e0; int inc = 1;

  
  dnaupd_(&ido, &bmat, &n, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  unsigned int count = 1;
  while (ido!=99){

    cornerTs.SetElem( workd + ipntr[0] - 1);
    ContractArgs(cornerTs, Eigen_net, target, cornerTs);
    double *cornerTs_p = cornerTs.GetElem();
    memcpy(workd+ipntr[1]-1, cornerTs_p, n*sizeof(double)); 
    dnaupd_(&ido, &bmat, &n, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    count++;
  }
  if ( info!=0 ) cerr << "Error:dnaupd, info = " << info << endl;
  
  int rvec = 1;
  char howmny = 'A';
  int *select = new int[ncv];
  double *dr = new double[nev+1];
  double *di = new double[nev+1];
  double *workev = new double[3*ncv];
  double *z = new double[n*(nev+1)];
  double sigmar, sigmai;
  dneupd_(&rvec, &howmny, select, dr, di, z, &ldv, &sigmar, &sigmai, workev, &bmat, &n, which, &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  if ( info!=0 ) cerr << "Error:dneupd, info = " << info << endl;
  complex<double> Eigval;
  Eigval.real(dr[0]);
  Eigval.imag(di[0]);
  if (di[0]!=0) cerr << "Error:Eigenvalue is complex" << endl; 
  cornerTs.SetElem(z);
  corner = cornerTs.GetBlock();
  ContractArgs(cornerTs, Eigen_net, target, cornerTs);
  cout << corner; 
  delete [] resid;
  delete [] v;
  delete [] workd;
  delete [] workl;
  delete [] select;
  delete [] workev;
  delete [] dr;
  delete [] di;
  delete [] z;
}

int arEignRightBMPS( double& EigVal, UniTensor<double> &EigVec_t, const int max_iter, double err_tol, const double imagLimit, Network &transfer_net, const UniTensor<double> &mpsUpL, const UniTensor<double> &mpsUpR, const UniTensor<double> &mpsDnL, const UniTensor<double> &mpsDnR, const UniTensor<double> &evolUp ){
  ///this routine find the max magnitude Eigenvalue and its Eigenvector with dnaupd and dneupd for a nonsymmetric matrix
  ///if one want to find other kind of or more Eigen pairs, tune which and nev.
  ///since the matrix is nonsymmetric, the Eigenpairs can be complex, so the input EigVal and EigVec_t should be complex.
  ///the complex Eigenvalue will appear with it's complex conjugate
  const int dimMPS = mpsUpL.bond(1).dim();
  //double *A = A_m.GetElem();
  int nev = 1;

  int dim = dimMPS*dimMPS; ///dimension of Eigenproblem
  EigVec_t = UniTensor<double> ( vector<Bond> (2, Bond( BD_IN, dimMPS )) );
  //EigVec_t.Randomize();
  vector<double> initElem( dimMPS*dimMPS, 1.0 );
  EigVec_t.SetElem( initElem );
  int ido = 0; ///reverse communication parameter, must be zero before iteration
  char bmat = 'I'; ///'I': standard Eigenproblem, 'G': generalized Eigenproblem 
  char which[] = {'L','M'}; ///type of asked Eigenvalues 
  int info = 1; double *resid = new double[dim];
  fill_n( resid, dim, 1.0 );
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
    UniTensor<double> tempt;
    ContractArgs( tempt, transfer_net, mpsUpL, mpsUpR, mpsDnL, mpsDnR, evolUp, EigVec_t );
    EigVec_t = tempt; 
    double *EigVec_p = EigVec_t.GetElem();
    memcpy(workd+ipntr[1]-1, EigVec_p, dim*sizeof(double));
    dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    count++;
  }
  if( info!=0 ){
      cerr<<"error with dnaupd, info = "<<info<<endl;
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
      cerr<<"Error with dneupd, info = "<<info<<endl;
  }
  else {}

  ///extract results and informations
  EigVal = dr[0];
  assert( fabs( di[0] ) < imagLimit );
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

int arEignRightBMPS( complex<double>& EigVal, UniTensor<complex<double>> &EigVec_t, const int max_iter, double err_tol, const double imagLimit, uni10::Network &transfer_net, const UniTensor<complex<double>> &mpsUpL, const UniTensor<complex<double>> &mpsUpR, const UniTensor<complex<double>> &mpsDnR, const UniTensor<complex<double>> &mpsDnL, const UniTensor<complex<double>> &evolUp ){
//int arCompEignRight( const Matrix<complex<double>> &A_m, complex<double>& EigVal, Matrix<complex<double>> &EigVec_m, const int max_iter, double err_tol, const bool ifinit ){
  //this routine find the max magnitude Eigenvalue and its Eigenvector with znaupd and zneupd for a nonhermitian matrix
  //if one want to find other kind of or more Eigen pairs, tune which and nev.
  const int dimMPS = mpsUpL.bond(1).dim();
  int nev = 1;

  int dim = dimMPS*dimMPS;
  EigVec_t = UniTensor<complex<double>> ( vector<Bond> (2, Bond( BD_IN, dimMPS )) );
  vector<complex<double>> initElem( dimMPS*dimMPS, 1.0 );
  EigVec_t.SetElem( initElem );
  int ido = 0; ///reverse communication parameter, must be zero before iteration
  char bmat = 'I'; ///'I': standard Eigenproblem, 'G': generalized Eigenproblem 
  char which[] = {'L','M'}; ///type of asked Eigenvalues 
  int info=1; complex<double> *resid = new complex<double> [dim];
  fill_n( resid, dim, 1.0 );
  ///If info = 0, a randomly initial residual vector is used.
  ///If info = 1, resid contains the initial guess vector provided by user or from previous run
  ///on output, resid contains the final residual vector
  int ncv = 5; ///the number of Ritz vector, nev+2 <= ncv <= dim
  if (ncv>dim){
    ncv = dim;
  } else {}
  int ldv = dim; ///leading dimension of v
  complex<double> *v = new complex<double>[dim*ncv];
  int iparam[11];
  iparam[0] = 1;        ///Specifies the shift strategy (1->exact)
  iparam[2] = max_iter; ///Maximum number of iterations
  iparam[3] = 1;        //current code only for 1
  iparam[6] = 1;        ///Sets the mode of znaupd.
  int ipntr[14];
  complex<double> *workd = new complex<double>[3*dim];
  int lworkl = 3*ncv*(ncv+2); ///LWORKL must be at least 3*NCV**2 + 6*NCV .
  complex<double> *workl = new complex<double>[lworkl];
  double *rwork = new double[ncv];

  ///start iteration
  //dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  znaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
  unsigned int count = 1;
  ///Parameters for zgemv
  complex<double> alpha( 1.0, 0); complex<double> beta( 0, 0); int inc = 1;
  while( ido!=99 ){ ///for mode 1 (iparam[0]==1), ido only  can be 1 or 99
    //dgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    //dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    //zgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
    UniTensor<complex<double>> tempt;
    ContractArgs( tempt, transfer_net, mpsUpL, mpsUpR, mpsDnL, mpsDnR, evolUp, EigVec_t );
    EigVec_t = tempt; 
    znaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
    count++;
  }
  if( info!=0 ){
      cerr<<"error with znaupd, info = "<<info<<endl;
  } else{}
  ///calculate Eigenvalue and Eigenvector by dneupd
  int rvec = 1; ///0: only Eigenvalue, 1: also Eigenvector
  char howmny = 'A'; ///how many Eigenvectors to calculate: 'A' => nev Eigenvectors
  int *select = new int[ncv];/// when howmny == 'A', this is used as workspace to reorder the Eigenvectors
  complex<double> *d = new complex<double>[nev+1];
  complex<double> *z = new complex<double>[dim*nev];
  complex<double> sigma;
  complex<double> *workev = new complex<double>[2*ncv];
  //dneupd_(&rvec, &howmny, select, dr, di, z, &ldv, &sigmar, &sigmai, workev, &bmat, &dim, which, &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  zneupd_(&rvec, &howmny, select, d, z, &ldv, &sigma, workev, &bmat, &dim, which, &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
  if ( info != 0 ){
      cerr<<"error with zneupd, info = "<<info<<endl;
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

int arRealEignRightBMPSCompress( double& EigVal, UniTensor<double> &EigVec_t, const int max_iter, double err_tol, const double imagLimit, uni10::Network &transfer_net, const UniTensor<double> &up0, const UniTensor<double> &up1, const UniTensor<double> &triangleUp, const UniTensor<double> &dn0, const UniTensor<double> &dn1, const UniTensor<double> &triangleDn ){
  ///this routine find the max magnitude Eigenvalue and its Eigenvector with dnaupd and dneupd for a nonsymmetric matrix
  ///if one want to find other kind of or more Eigen pairs, tune which and nev.
  ///since the matrix is nonsymmetric, the Eigenpairs can be complex, so the input EigVal and EigVec_t should be complex.
  ///the complex Eigenvalue will appear with it's complex conjugate
  const int dimVir = up0.bond(0).dim();
  const int dimMPS = up0.bond(1).dim();
  //double *A = A_m.GetElem();
  int nev = 1;

  int dim = dimMPS*dimMPS*dimVir*dimVir; ///dimension of Eigenproblem
  vector<Bond> vectBds = {Bond(BD_IN, dimMPS), Bond( BD_IN, 1), Bond( BD_IN, dimVir), Bond(BD_IN, dimVir), Bond(BD_IN, dimMPS) };
  EigVec_t = UniTensor<double> (vectBds);
  vector<double> initElem( dim, 1.0 );
  EigVec_t.SetElem( initElem );
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
    UniTensor<double> tempt;
    //ContractArgs( tempt, transfer_net, mpsUpL, mpsUpR, mpsDnL, mpsDnR, evolUp, EigVec_t );
    ContractArgs( tempt, transfer_net, up0, up1, triangleUp, triangleUp, dn0, dn1, triangleDn, triangleDn, EigVec_t );
    EigVec_t = tempt; 
    double *EigVec_p = EigVec_t.GetElem();
    memcpy(workd+ipntr[1]-1, EigVec_p, dim*sizeof(double));
    dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    count++;
  }
  if( info!=0 ){
      cerr<<"error with dnaupd, info = "<<info<<endl;
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
      cerr<<"Error with dneupd, info = "<<info<<endl;
  }
  else {}

  ///extract results and informations
  EigVal = dr[0];
  assert( fabs( di[0] ) < imagLimit );
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

int arCompEignRightBMPSCompress( complex<double> &EigVal, UniTensor<double> &EigVec_t, const int max_iter, double err_tol, const double imagLimit, uni10::Network &transfer_net, const UniTensor<double> &up0, const UniTensor<double> &up1, const UniTensor<double> &triangleUp, const UniTensor<double> &dn0, const UniTensor<double> &dn1, const UniTensor<double> &triangleDn ){
  ///this routine find the max magnitude Eigenvalue and its Eigenvector with dnaupd and dneupd for a nonsymmetric matrix
  ///if one want to find other kind of or more Eigen pairs, tune which and nev.
  ///since the matrix is nonsymmetric, the Eigenpairs can be complex, so the input EigVal and EigVec_t should be complex.
  ///the complex Eigenvalue will appear with it's complex conjugate
  const int dimVir = up0.bond(0).dim();
  const int dimMPS = up0.bond(1).dim();
  //double *A = A_m.GetElem();
  int nev = 1;

  int dim = dimMPS*dimMPS*dimVir*dimVir; ///dimension of Eigenproblem
  vector<Bond> vectBds = {Bond(BD_IN, dimMPS), Bond( BD_IN, 1), Bond( BD_IN, dimVir), Bond(BD_IN, dimVir), Bond(BD_IN, dimMPS) };
  EigVec_t = UniTensor<double> (vectBds);
  UniTensor<double> vectTen(vectBds);
  vector<double> initElem( dim, 1.0/sqrt(double(dim)) );
  vectTen.SetElem( initElem );
  int ido = 0; ///reverse communication parameter, must be zero before iteration
  char bmat = 'I'; ///'I': standard Eigenproblem, 'G': generalized Eigenproblem 
  char which[] = {'L','M'}; ///type of asked Eigenvalues 
  int info = 1; double *resid = new double[dim];
  fill_n( resid, dim, 1.0 );
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
    UniTensor<double> tempt;
    //ContractArgs( tempt, transfer_net, mpsUpL, mpsUpR, mpsDnL, mpsDnR, evolUp, EigVec_t );
    ContractArgs( tempt, transfer_net, up0, up1, triangleUp, triangleUp, dn0, dn1, triangleDn, triangleDn, vectTen );
    vectTen = tempt; 
    double *EigVec_p = vectTen.GetElem();
    memcpy(workd+ipntr[1]-1, EigVec_p, dim*sizeof(double));
    dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    count++;
  }
  if( info!=0 ){
      cerr<<"error with dnaupd, info = "<<info<<endl;
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
      cerr<<"Error with dneupd, info = "<<info<<endl;
  }
  else {}

  ///extract results and informations
  EigVal = complex<double> (dr[0],di[0]);
  //assert( fabs( di[0] ) < imagLimit );
  //EigVec_t.SetElem( z );
  double *values = new double [dim];
  for (int i=0; i!=dim; i++){
    values[i] = z[i];
    //values[i] = complex<double> (z[i], z[dim+i]);
  }
  //EigVec_t.SetElem( values );
  EigVec_t.SetElem( z );
  delete [] values;
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

int arRealEignRightFivePess( double& EigVal, UniTensor<double> &EigVec_t, const int max_iter, double err_tol, const double imagLimit, Network &transfer_net, const UniTensor<double> &mpsUpL, const UniTensor<double> &mpsUpR, const UniTensor<double> &mpsDnL, const UniTensor<double> &mpsDnR, const vector<UniTensor<double>> &evolUp ){
  ///this routine find the max magnitude Eigenvalue and its Eigenvector with dnaupd and dneupd for a nonsymmetric matrix
  ///if one want to find other kind of or more Eigen pairs, tune which and nev.
  ///since the matrix is nonsymmetric, the Eigenpairs can be complex, so the input EigVal and EigVec_t should be complex.
  ///the complex Eigenvalue will appear with it's complex conjugate
  const int dimMPS = mpsUpL.bond(1).dim();
  const int dimEvol = evolUp.at(0).bond(0).dim();
  //double *A = A_m.GetElem();
  int nev = 1;

  int dim = dimMPS*dimMPS; ///dimension of Eigenproblem
  EigVec_t = UniTensor<double> ( vector<Bond> { Bond( BD_IN, dimMPS ), Bond(BD_IN, dimEvol), Bond(BD_IN, dimMPS) } );
  EigVec_t.Randomize();
  //vector<double> initElem( dimEvol*dimMPS*dimMPS, 1.0 );
  //EigVec_t.SetElem( initElem );
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
    UniTensor<double> tempt;
    ContractArgs( tempt, transfer_net, mpsUpL, mpsUpR, mpsDnL, mpsDnR, evolUp.at(0), evolUp.at(1), EigVec_t );
    EigVec_t = tempt; 
    double *EigVec_p = EigVec_t.GetElem();
    memcpy(workd+ipntr[1]-1, EigVec_p, dim*sizeof(double));
    dnaupd_(&ido, &bmat, &dim, &which[0], &nev, &err_tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    count++;
  }
  if( info!=0 ){
      cerr<<"error with dnaupd, info = "<<info<<endl;
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
      cerr<<"Error with dneupd, info = "<<info<<endl;
  }
  else {}

  ///extract results and informations
  EigVal = dr[0];
  assert( fabs( di[0] ) < imagLimit );
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

void arRealSymEignRight( const Matrix<double> &amat, const int nEigVal, Matrix<double> &EigVal, Matrix<double> &EigVec, const int maxIter, const double errTol){
  ///For symmetric matrix this function finds given number of Eigenvalue and its Eigenvector with dsaupd and dseupd.

  int ido = 0; 
  char bmat = 'I';
  int dim = amat.col();
  char which[] = {'L','M'}; 
  double *A = amat.GetElem();
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
    dgemv((char*)"T", &dim, &dim, &alpha, A, &dim, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1, &inc); ///multiplication
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

  ///extract results and informations
  double *finald = new double[nev];
  double *finalz = new double[dim*nev];
  for (int i=0; i!=nev; i++){
    finald[i] = d[nev-i-1];
  }
  for (int i=0; i!=dim; i++){
    for (int j=0; j!=nev; j++){
      finalz[i*nev+j] = z[(nev-j-1)*dim+i];
    }
  }
  EigVal = Matrix<double> (nEigVal, nEigVal, true);
  EigVec = Matrix<double> (    dim, nEigVal, false);
  EigVal.SetElem(finald);
  EigVec.SetElem(finalz);
  ///free memory
  delete [] resid;
  delete [] v;
  delete [] workd;
  delete [] workl;
  delete [] select;
  delete [] d;
  delete [] z;
  delete [] finald;
  delete [] finalz;
}

void arRealSVD( const Matrix<double> &amat, const int nEigVal, Matrix<double> &EigVal, Matrix<double> &EigVec, const int maxIter, const double errTol){
  ///For symmetric matrix this function finds given number of Eigenvalue and its Eigenvector with dsaupd and dseupd.
  Matrix<double> aT = Transpose( amat );

  int ido = 0; 
  char bmat = 'I';
  int smallD = amat.col();
  int dim = 2*smallD;
  char which[] = {'L','A'}; 
  double *A = amat.GetElem();
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
    dgemv((char*)"T", &smallD, &smallD, &alpha, A, &smallD, workd+ipntr[0]-1+smallD, &inc, &beta, workd+ipntr[1]-1, &inc);
    dgemv((char*)"N", &smallD, &smallD, &alpha, A, &smallD, workd+ipntr[0]-1, &inc, &beta, workd+ipntr[1]-1+smallD, &inc);
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

  ///extract results and informations
  double *finald = new double[nev];
  double *finalz = new double[dim*nev];
  for (int i=0; i!=nev; i++){
    finald[i] = d[nev-i-1];
  }
  for (int i=0; i!=dim; i++){
    for (int j=0; j!=nev; j++){
      finalz[i*nev+j] = z[(nev-j-1)*dim+i];
    }
  }
  EigVal = Matrix<double> (nEigVal, nEigVal, true);
  EigVec = Matrix<double> (    dim, nEigVal, false);
  EigVal.SetElem(finald);
  EigVec.SetElem(finalz);
  ///free memory
  delete [] resid;
  delete [] v;
  delete [] workd;
  delete [] workl;
  delete [] select;
  delete [] d;
  delete [] z;
  delete [] finald;
  delete [] finalz;
}

std::vector<uni10::Matrix<double>> arRealSvd( const uni10::Matrix<double> &amat, const int nRankIn, const int maxIter, const double errTol ){
  int nRow = amat.row();
  int nCol = amat.col();
  //assert( nRow==nCol );
  int nEig = std::min(nRankIn, nCol);

  uni10::Matrix<double> aTa = Dot(Transpose(amat), amat);
  
  uni10::Matrix<double> EigVal, EigVec;
  arRealSymEignRight( aTa, nEig, EigVal, EigVec, maxIter, errTol );
  uni10::Matrix<double> smat( nEig, nEig, true );
  uni10::Matrix<double> sInv( nEig, nEig, true );
  uni10::Matrix<double> vTmat = Transpose( EigVec );
  for (int i=0; i!=EigVal.ElemNum(); i++){
    smat[i] = sqrt(EigVal.At(i,i));
    sInv[i] = (smat.At(i,i)<1.0e-14)? 0:1.0/smat.At(i,i);
  }
  uni10::Matrix<double> umat = Dot( amat, Dot(EigVec, sInv) );
  std::vector<uni10::Matrix<double>> result = {umat, smat, vTmat};
  return result;
}

