#ifndef _HAMILTONIAN_H_
#define _HAMILTONIAN_H_

#include "../nsy_opt/operator.h"

/*! 
 * \defgroup aux Auxiliary functions
*/


//!  Two site XXZ hamiltonian \ingroup aux 
uni10::UniTensor<double> XXZ(float Jx, float Jz, float spin = 0.5);

//!  Two site Heisenberg hamiltonian \ingroup aux
uni10::UniTensor<double> Heisenberg(float spin=0.5, double J=1.0);

//!  Transverse Ising hamiltonian \ingroup aux 
uni10::UniTensor<double> transverseIsing(float spin, float h, bool isAnti);

//! Pij = 1/4( Si dot Sj ) \ingroup aux 
uni10::UniTensor<double> Pij(float spin);

//! Four site J-Q model \ingroup aux
uni10::UniTensor<double> JQfourSites(double J, double Q);

//! Four site J1-J2 model \ingroup aux 
uni10::UniTensor<double> horizontalJ_verticalJp(double J, double Jp);

//! Create a product hamiltonian by two site hamiltonian with PBC. \ingroup aux
uni10::UniTensor<double> periodicHamiltonian(int N, const uni10::UniTensor<double>& H0);

uni10::UniTensor<std::complex<double>> periodicHamiltonian(int N, const uni10::UniTensor<std::complex<double>>& H0);

//! DM interatcion zhat dot (Si curl Sj ) \ingroup aux
uni10::UniTensor<double> DMinteraction();

uni10::UniTensor<std::complex<double>> xxz_DM( const double Jxy, const double Dmz );//set Jz=1.0

uni10::UniTensor<std::complex<double>> DM_upTriangle_inPlane( const double field );

uni10::Bond spin_bond(float spin, uni10::bondType btype, bool U1=false);
/*

std::vector<uni10::UniTensor<std::complex<double> > > KitaevModel( const double Jx=1.0, const double Jy=1.0, const double Jz=1.0 );

std::vector<uni10::UniTensor<std::complex<double> > > GammaModel( const double Jx=1.0, const double Jy=1.0, const double Jz=1.0 );

std::vector<uni10::UniTensor<double>> GammaModelTRS( const double Jx=1.0, const double Jy=1.0, const double Jz=1.0 );

uni10::UniTensor<std::complex<double> > KitaevModelSixSite( const double Jx=1.0, const double Jy=1.0, const double Jz=1.0 );

uni10::UniTensor<double> fiveSiteDM( const double field );

uni10::UniTensor<double> Heisenberg_stagger_field(float spin=0.5, double field=0);

uni10::UniTensor<std::complex<double>> Heisenberg_Sy(float spin=0.5, double field=0);

uni10::UniTensor<double> Heisenberg_U1(float spin=0.5, double J=1.0);

uni10::UniTensor<double> theModel(float spin, int i, double delta, double Dz, double hz, double dx);

uni10::UniTensor<double> directSum(size_t s1, size_t s2, const uni10::UniTensor<double>& H2, const uni10::UniTensor<double>& H_all);

uni10::UniTensor<double> J1J2model(double J1, double J2);

*/

#endif
