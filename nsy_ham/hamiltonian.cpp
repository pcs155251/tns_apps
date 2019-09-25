#include "hamiltonian.h"
using namespace std;

uni10::Bond spin_bond(float spin, uni10::bondType btype, bool U1){

  spin_check(spin);
  int dim = spin * 2 + 1;

  if(U1){
    std::vector<uni10::Qnum> qnums(dim);
    int halfint = true;
    if(spin == floor(spin))
      halfint = false;
    for(int i = 0; i < dim; i++){
      int s = spin - i;
      if(halfint){
        s = spin + 0.5 - i;
        if(s <= 0)
          s--;
      }
      qnums[i] = uni10::Qnum(s);
    }
    return uni10::Bond(btype, qnums);
  }
  else{

    return uni10::Bond(btype, dim);

  }

}

uni10::UniTensor<double> XXZ(float Jx, float Jz, float spin){

  uni10::Matrix<double> sp = matSp(spin);
  uni10::Matrix<double> sm = matSm(spin);
  uni10::Matrix<double> sz = matSz(spin);
  uni10::Matrix<double> ham = (double)Jz*uni10::Otimes(sz, sz);
  ham += Jx* 0.5 * (uni10::Otimes(sp, sm) + uni10::Otimes(sm, sp));
  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  std::vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor<double> H(bonds, "Heisenberg");
  H.PutBlock(ham);
  
  return H;

}

uni10::UniTensor<double> Heisenberg(float spin, double J){

  uni10::Matrix<double> sp = matSp(spin);
  uni10::Matrix<double> sm = matSm(spin);
  uni10::Matrix<double> sz = matSz(spin);
  uni10::Matrix<double> ham = uni10::Otimes(sz, sz);
  ham += 0.5 * (uni10::Otimes(sp, sm) + uni10::Otimes(sm, sp));
  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  std::vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor<double> H(bonds, "Heisenberg");
  H.PutBlock(ham);
  return J * H;

}

uni10::UniTensor<double> Pij(float spin){
  uni10::UniTensor<double> H = Heisenberg(spin,1.0);
  uni10::UniTensor<double> I = H;
  I.Identity();
  H = 0.25*I + (-1.0*H);
  return H;
}

uni10::UniTensor<double> JQfourSites(double J, double Q){
  //J term
  uni10::UniTensor<double> PtwoSites = Pij(0.5);
  uni10::UniTensor<double> Jterm = periodicHamiltonian( 4, PtwoSites );
  //Q term
  uni10::UniTensor<double> Q1 = Otimes( PtwoSites, PtwoSites );
  Q1.SetLabel( vector<int> {0,1,2,3,4,5,6,7} );
  uni10::UniTensor<double> Q2 = Permute( Q1, vector<int> {1,3,0,2,5,7,4,6}, 4 );
  uni10::UniTensor<double> Qterm = Q1 + Q2;
  Qterm.SetLabel( vector<int> {0,1,2,3,4,5,6,7} );
  Qterm = Permute( Qterm, vector<int> {0,2,3,1,4,6,7,5}, 4 );
  //H
  uni10::UniTensor<double> H = -J*Jterm - Q*Qterm;
  return H;
}

uni10::UniTensor<double> horizontalJ_verticalJp(double J, double Jp)
{
  uni10::UniTensor<double> PtwoSites = Pij(0.5);
  uni10::Matrix<double> mat = PtwoSites.GetBlock();
  uni10::Matrix<double> Jmat = Otimes( mat, mat );
  uni10::UniTensor<double> Jten( 
                                 vector<uni10::Bond> 
                                 {
                                   uni10::Bond( uni10::BD_IN, 2),
                                   uni10::Bond( uni10::BD_IN, 2),
                                   uni10::Bond( uni10::BD_IN, 2),
                                   uni10::Bond( uni10::BD_IN, 2),
                                   uni10::Bond( uni10::BD_OUT, 2),
                                   uni10::Bond( uni10::BD_OUT, 2),
                                   uni10::Bond( uni10::BD_OUT, 2),
                                   uni10::Bond( uni10::BD_OUT, 2)
                                  }
                                );
  Jten.SetLabel( vector<int> {1,2,3,4,-1,-2,-3,-4 } );
  uni10::UniTensor<double> Jpten = Jten;

  Jten.PutBlock(J*Jmat);
  Jpten.PutBlock(Jp*Jmat);
  Jpten = Permute( Jpten, vector<int> {2,3,4,1,-2,-3,-4,-1 }, 4 );
  Jten = Jten + Jpten;
  return Jten;
}

uni10::UniTensor<double> Heisenberg_stagger_field(float spin, double field){

  uni10::Matrix<double> id(2,2,true);
  id.Identity();
  uni10::Matrix<double> sp = matSp(spin);
  uni10::Matrix<double> sm = matSm(spin);
  uni10::Matrix<double> sz = matSz(spin);
  uni10::Matrix<double> ham = uni10::Otimes(sz, sz);
  ham += 0.5 * (uni10::Otimes(sp, sm) + uni10::Otimes(sm, sp));
  ham += 0.5*field* (uni10::Otimes(sz, id) + (-1.0)*uni10::Otimes(id, sz));
  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  std::vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor<double> H(bonds, "Heisenberg_stagger_field");
  H.PutBlock(ham);
  return H;

}

uni10::UniTensor<std::complex<double>> Heisenberg_Sy(float spin, double field){

  uni10::Matrix<double> id(2,2,true);
  id.Identity();
  uni10::Matrix<double> sp = matSp(spin);
  uni10::Matrix<double> sm = matSm(spin);
  uni10::Matrix<double> sz = matSz(spin);
  uni10::Matrix<std::complex<double>> sy = matSy(spin);
  uni10::Matrix<std::complex<double>> ham = uni10::Otimes(sz, sz);
  ham += 0.5 * (uni10::Otimes(sp, sm) + uni10::Otimes(sm, sp));
  ham += (1.0/3.0)*0.5*field* (uni10::Otimes(sy, id) + (-1.0)*uni10::Otimes(id, sy));
  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  std::vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor<std::complex<double>> H(bonds, "Heisenberg_Sy");
  H.PutBlock(ham);
  return H;

}

uni10::UniTensor<double> Heisenberg_U1(float spin, double J){

}

uni10::UniTensor<double> transverseIsing(float spin, float h, bool isAnti){

  uni10::Matrix<double> sx = matSx(spin);
  uni10::Matrix<double> sz = matSz(spin);
  uni10::Matrix<double> I(sx.row(), sx.col(), true);
  I.Identity();
  uni10::Matrix<double> ham = (isAnti) ? uni10::Otimes((double)2*sz, (double)2*sz) : ((double)-1.) * uni10::Otimes((double)2*sz, (double)2*sz); // Otimes(sigma_z, sizga_z);
  uni10::Matrix<double> sxl = uni10::Otimes((h/(double)2) * (double)2*sx, I);
  uni10::Matrix<double> sxr = uni10::Otimes(I, (h/(double)2) * (double)2*sx);
  ham = ham + 0.5 * (sxl + sxr) ;
  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  std::vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor<double> H(bonds, "transverseIsing");
  H.PutBlock(ham);
  return H;

}

uni10::UniTensor<double> theModel(float spin, int i, double delta, double Dz, double hz, double dx){

}

uni10::UniTensor<double> DMinteraction(){
  uni10::Matrix<double> sx=matSx(0.5);
  uni10::Matrix<double> sz=matSz(0.5);
  uni10::Matrix<double> outMatrix = uni10::Otimes( sz, sx ) + (-1.0)*uni10::Otimes( sx, sz );
  std::vector<uni10::Bond> twoSiteBds = {uni10::Bond( uni10::BD_IN, 2), uni10::Bond( uni10::BD_IN, 2), uni10::Bond( uni10::BD_OUT, 2), uni10::Bond( uni10::BD_OUT, 2) };
  uni10::UniTensor<double> DMhamiltonian( twoSiteBds );
  DMhamiltonian.PutBlock(outMatrix);
  return DMhamiltonian;
}

uni10::UniTensor<std::complex<double>> xxz_DM( const double Jxy, const double Dmz){
  uni10::Matrix<double> sx=matSx(0.5);
  uni10::Matrix<double> sz=matSz(0.5);
  uni10::Matrix<complex<double>> sy=matSy(0.5);
  uni10::Matrix<complex<double>> xxz_term = Jxy*( uni10::Otimes(sx,sx) + uni10::Otimes(sy,sy) ) + 1.0*( uni10::Otimes(sz,sz) );
  uni10::Matrix<complex<double>> DM_term = Dmz*( uni10::Otimes(sx,sy) + (-1.0)*uni10::Otimes( sy, sx ) );
  cout<<DM_term<<endl;

  std::vector<uni10::Bond> twoSiteBds = {uni10::Bond( uni10::BD_IN, 2), uni10::Bond( uni10::BD_IN, 2), uni10::Bond( uni10::BD_OUT, 2), uni10::Bond( uni10::BD_OUT, 2) };
  uni10::UniTensor<complex<double>> hamiltonian( twoSiteBds );
  hamiltonian.PutBlock(xxz_term+DM_term);
  return hamiltonian;
}

uni10::UniTensor<std::complex<double>> DM_upTriangle_inPlane( const double field ){
  uni10::Matrix<double> sx=matSx(0.5);
  uni10::Matrix<double> sz=matSx(0.5);
  uni10::Matrix<std::complex<double>> sy=matSy();
  // x component of Si curl Sj 
  uni10::Matrix<std::complex<double>> XcurlSiSj = uni10::Otimes( sy, sz ) + (-1.0)*uni10::Otimes( sz, sy );
  // z component of Si curl Sj 
  uni10::Matrix<std::complex<double>> ZcurlSiSj = uni10::Otimes( sx, sy ) + (-1.0)*uni10::Otimes( sy, sx );
  // Heisenberg term
  uni10::Matrix<double> heisenberg = Heisenberg(0.5,1.0).GetBlock();

  // local hamiltonian between sites
  uni10::Matrix<std::complex<double>> H01 = heisenberg +                0*field*XcurlSiSj +        field*ZcurlSiSj;
  uni10::Matrix<std::complex<double>> H12 = heisenberg + (-0.5*sqrt(3.0))*field*XcurlSiSj + (-0.5)*field*ZcurlSiSj;
  uni10::Matrix<std::complex<double>> H20 = heisenberg +    0.5*sqrt(3.0)*field*XcurlSiSj + (-0.5)*field*ZcurlSiSj;
  
  // Bonds of three site Hamiltonian
  std::vector<uni10::Bond> threeSiteBds(3,uni10::Bond(uni10::BD_IN,2));
  std::vector<uni10::Bond> threeSiteOutBds(3,uni10::Bond(uni10::BD_OUT,2));
  threeSiteBds.insert(threeSiteBds.end(), threeSiteOutBds.begin(), threeSiteOutBds.end());

  // three site hamiltonian
  uni10::Matrix<double> Id(2,2,true);
  Id.Identity();
  H01 = Otimes(H01, Id);
  uni10::UniTensor<std::complex<double>> Ham01(threeSiteBds);
  std::cout<<H01;
  std::cout<<Ham01;
  Ham01.PutBlock(H01);

  H12 = Otimes( Id,H12);
  uni10::UniTensor<std::complex<double>> Ham12(threeSiteBds);
  Ham12.PutBlock(H12);

  H20 = Otimes( Id,H20);
  uni10::UniTensor<std::complex<double>> Ham20(threeSiteBds);
  Ham20.PutBlock(H20);
  Ham20.SetLabel({2,3,1,-2,-3,-1});
  Ham20 = Permute(Ham20,{1,2,3,-1,-2,-3},3);

  return Ham01+Ham12+Ham20;
}

uni10::UniTensor<double> fiveSiteDM( const double field ){
  uni10::UniTensor<double> twoSiteHam = Heisenberg( 0.5, 1.0 );
  uni10::Matrix<double> twoSiteMat = twoSiteHam.GetBlock() + field*DMinteraction().GetBlock();
  twoSiteHam.PutBlock( twoSiteMat );
  uni10::UniTensor<double> threeSiteHam = periodicHamiltonian( 3, twoSiteHam );
  uni10::Matrix<double> twoSiteId( 4, 4, true );
  twoSiteId.Identity();
  uni10::Matrix<double> fiveSiteMat = Otimes( threeSiteHam.GetBlock(), twoSiteId );
  std::vector<uni10::Bond> allBds( 5, uni10::Bond( uni10::BD_IN, 2 ) );
  std::vector<uni10::Bond> outBds( 5, uni10::Bond( uni10::BD_OUT, 2 ) );
  allBds.insert( allBds.end(), outBds.begin(), outBds.end() );
  uni10::UniTensor<double> hamiltonian( allBds );
  hamiltonian.PutBlock( fiveSiteMat );
  uni10::UniTensor<double> hamiltonianCopy = Permute( hamiltonian, {3, 4, 2, 0, 1, 8, 9, 7, 5, 6}, 5 );
  hamiltonian = hamiltonian + hamiltonianCopy;
  return hamiltonian;
}

uni10::UniTensor<double> periodicHamiltonian(int N, const uni10::UniTensor<double>& H0){
  if(H0.BondNum() != 4){
    std::ostringstream err;
    err<<"Function periodicHamiltonain(N, H0) only for 2-site H0\n";
    throw std::runtime_error(err.str());
  }
  std::vector<uni10::Bond> bondI;
  bondI.push_back(H0.bond(0)), bondI.push_back(H0.bond(2));
  uni10::UniTensor<double> I(bondI);
  I.Identity();
  uni10::UniTensor<double> nH = H0;
  for(int i=0; i!=N-2; i++){
    nH = Otimes(nH, I);
  }
  uni10::UniTensor<double> perH = nH;
  std::vector<int> labels = perH.label();
  std::vector<int> per_labels(labels.size());
  int BondNum = per_labels.size();
  for(int l=0; l<BondNum/2; l++){
    per_labels[l] = labels[(l+1) % (BondNum/2)];
    per_labels[(BondNum/2) + l] = labels[BondNum/2 + ((l+1)%(BondNum/2))];
  }
  for(int h=0; h <N-1; h++){
    perH = Permute( perH, per_labels, BondNum/2);
    nH.PutBlock( nH.GetBlock() + perH.GetBlock() );
    perH.SetLabel(labels);
  }
  return nH;
}

uni10::UniTensor<complex<double>> periodicHamiltonian(int N, const uni10::UniTensor<complex<double>>& H0){
  if(H0.BondNum() != 4){ std::ostringstream err;
    err<<"Function periodicHamiltonain(N, H0) only for 2-site H0\n";
    throw std::runtime_error(err.str());
  }
  std::vector<uni10::Bond> bondI;
  bondI.push_back(H0.bond(0)), bondI.push_back(H0.bond(2));
  uni10::UniTensor<complex<double>> I(bondI);
  I.Identity();
  uni10::UniTensor<complex<double>> nH = H0;
  for(int i=0; i!=N-2; i++){
    nH = Otimes(nH, I);
  }
  uni10::UniTensor<complex<double>> perH = nH;
  std::vector<int> labels = perH.label();
  std::vector<int> per_labels(labels.size());
  int BondNum = per_labels.size();
  for(int l=0; l<BondNum/2; l++){
    per_labels[l] = labels[(l+1) % (BondNum/2)];
    per_labels[(BondNum/2) + l] = labels[BondNum/2 + ((l+1)%(BondNum/2))];
  }
  for(int h=0; h <N-1; h++){
    perH = Permute( perH, per_labels, BondNum/2);
    nH.PutBlock( nH.GetBlock() + perH.GetBlock() );
    perH.SetLabel(labels);
  }
  return nH;
}

std::vector<uni10::UniTensor<std::complex<double> > > KitaevModel( const double Jx, const double Jy, const double Jz ){
  uni10::Matrix<std::complex<double>> xx=Otimes(matSx(),matSx());
  uni10::Matrix<std::complex<double>> yy=Otimes(matSy(),matSy());
  uni10::Matrix<std::complex<double>> zz=Otimes(matSz(),matSz());
  std::vector<uni10::Bond> twoSiteHamBds = { uni10::Bond( uni10::BD_IN, 2),uni10::Bond( uni10::BD_IN, 2),uni10::Bond( uni10::BD_OUT, 2),uni10::Bond( uni10::BD_OUT, 2) };
  uni10::UniTensor<std::complex<double>> xxHam(twoSiteHamBds);
  uni10::UniTensor<std::complex<double>> yyHam(twoSiteHamBds);
  uni10::UniTensor<std::complex<double>> zzHam(twoSiteHamBds);
  xxHam.PutBlock(Jx*xx);
  yyHam.PutBlock(Jy*yy);
  zzHam.PutBlock(Jz*zz);
  return std::vector<uni10::UniTensor<std::complex<double> > > {xxHam, yyHam, zzHam};
}

std::vector<uni10::UniTensor<std::complex<double> > > GammaModel( const double Jx, const double Jy, const double Jz ){
  uni10::Matrix<std::complex<double>> xx=Otimes(matSy(),matSz())+Otimes(matSz(),matSy());
  uni10::Matrix<std::complex<double>> yy=Otimes(matSx(),matSz())+Otimes(matSz(),matSx());;
  uni10::Matrix<std::complex<double>> zz=Otimes(matSx(),matSy())+Otimes(matSy(),matSx());;
  std::vector<uni10::Bond> twoSiteHamBds = { uni10::Bond( uni10::BD_IN, 2),uni10::Bond( uni10::BD_IN, 2),uni10::Bond( uni10::BD_OUT, 2),uni10::Bond( uni10::BD_OUT, 2) };
  uni10::UniTensor<std::complex<double>> xxHam(twoSiteHamBds);
  uni10::UniTensor<std::complex<double>> yyHam(twoSiteHamBds);
  uni10::UniTensor<std::complex<double>> zzHam(twoSiteHamBds);
  xxHam.PutBlock(Jx*xx);
  yyHam.PutBlock(Jy*yy);
  zzHam.PutBlock(Jz*zz);
  return std::vector<uni10::UniTensor<std::complex<double> > > {xxHam, yyHam, zzHam};
}

std::vector<uni10::UniTensor<double>> GammaModelTRS( const double Jx, const double Jy, const double Jz ){
  uni10::Matrix<std::complex<double>> xxC=Otimes(matSy(),matSz())+Otimes(matSz(),matSy());
  uni10::Matrix<std::complex<double>> yyC=Otimes(matSx(),matSz())+Otimes(matSz(),matSx());
  uni10::Matrix<std::complex<double>> zzC=Otimes(matSx(),matSy())+Otimes(matSy(),matSx());
  std::vector<uni10::Bond> twoSiteHamBds = { uni10::Bond( uni10::BD_IN, 2),uni10::Bond( uni10::BD_IN, 2),uni10::Bond( uni10::BD_OUT, 2),uni10::Bond( uni10::BD_OUT, 2) };
  uni10::Matrix<std::complex<double>> umat = twoSiteTRStransform();
  uni10::Matrix<std::complex<double>> uinv = Dagger(umat);
  xxC = Dot(uinv,Dot(xxC,umat));
  yyC = Dot(uinv,Dot(yyC,umat));
  zzC = Dot(uinv,Dot(zzC,umat));
  uni10::Matrix<double> xx = xxC;
  uni10::Matrix<double> yy = yyC;
  uni10::Matrix<double> zz = zzC;
  uni10::UniTensor<double> xxHam(twoSiteHamBds);
  uni10::UniTensor<double> yyHam(twoSiteHamBds);
  uni10::UniTensor<double> zzHam(twoSiteHamBds);
  xxHam.PutBlock(Jx*xx);
  yyHam.PutBlock(Jy*yy);
  zzHam.PutBlock(Jz*zz);
  return std::vector<uni10::UniTensor<double > > {xxHam, yyHam, zzHam};
}

uni10::UniTensor<std::complex<double> > KitaevModelSixSite( const double Jx, const double Jy, const double Jz ){
  //the factor 4 is modify from spin to sigma
  uni10::Matrix<std::complex<double>> xx=4*Jx*Otimes(matSx(),matSx());
  uni10::Matrix<std::complex<double>> yy=4*Jy*Otimes(matSy(),matSy());
  uni10::Matrix<std::complex<double>> zz=4*Jz*Otimes(matSz(),matSz());
  uni10::Matrix<std::complex<double>> foureSiteId(16,16,true);
  foureSiteId.Identity();
  uni10::Matrix<std::complex<double>> sixSitexx = Otimes(xx,foureSiteId);
  uni10::Matrix<std::complex<double>> sixSiteyy = Otimes(yy,foureSiteId);
  uni10::Matrix<std::complex<double>> sixSitezz = Otimes(zz,foureSiteId);

  std::vector<uni10::Bond> sixSiteHamBds( 6, uni10::Bond( uni10::BD_IN, 2) );
  std::vector<uni10::Bond> sixSiteHamOutBds( 6, uni10::Bond( uni10::BD_OUT, 2) );
  sixSiteHamBds.insert( sixSiteHamBds.end(), sixSiteHamOutBds.begin(), sixSiteHamOutBds.end() );
  uni10::UniTensor<std::complex<double>> sixSiteHam(sixSiteHamBds);
  sixSiteHam.SetLabel({1,2,3,4,5,6,-1,-2,-3,-4,-5,-6});
  sixSiteHam.PutBlock(sixSitexx); 
  sixSiteHam = Permute( sixSiteHam, {2,3,4,5,6,1,-2,-3,-4,-5,-6,-1},6 );
  sixSiteHam.PutBlock(sixSiteHam.GetBlock()+sixSiteyy); 
  sixSiteHam = Permute( sixSiteHam, {3,4,5,6,1,2,-3,-4,-5,-6,-1,-2},6 );
  sixSiteHam.PutBlock(sixSiteHam.GetBlock()+sixSitezz); 
  sixSiteHam = Permute( sixSiteHam, {1,2,3,4,5,6,-1,-2,-3,-4,-5,-6},6 );

  uni10::Matrix<std::complex<double>> xxyyzz = sixSiteHam.GetBlock();
  sixSiteHam = Permute( sixSiteHam, {4,5,6,1,2,3,-4,-5,-6,-1,-2,-3},6 );
  sixSiteHam.PutBlock(sixSiteHam.GetBlock()+xxyyzz); 
  sixSiteHam.SetLabel({1,2,3,4,5,6,-1,-2,-3,-4,-5,-6});
  return sixSiteHam;
}

/*
uni10::UniTensor<double> J1J2model(double J1, double J2){

}

uni10::UniTensor<double> directSum(size_t s1, size_t s2, const uni10::UniTensor<double>& H2, const uni10::UniTensor<double>& H_all){

}
*/
