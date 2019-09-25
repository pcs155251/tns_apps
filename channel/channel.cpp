#include "channel.h"
using namespace std;
using namespace uni10;

channel::channel( const int dimMPSIn, const UniTensor<complex<double>> &ansazIn, const UniTensor<double> &hamiltonianIn ):dimPhy{ansazIn.bond(0).dim()}, dimMPS(dimMPSIn), ansaz( ansazIn ), hamiltonian( hamiltonianIn ), edges( vector<oneSiteiMPS<complex<double>>> (4) ), corners( vector<UniTensor<complex<double>>> (4) ), caps( vector<UniTensor<complex<double>>> (4) ), sups( vector<UniTensor<complex<double>>> (4) ), bluecaps( vector<UniTensor<complex<double>>> (4)), bluesI( vector<UniTensor<complex<double>>> (4)), bluesII( vector<UniTensor<complex<double>>> (4)), reds( vector<UniTensor<complex<double>>> (4)) {  
  dimD = ansaz.bond(1).dim();

  enum{ iLeft, iTop, iRight, iDown };
  for (int iDir=iLeft; iDir<=iDown; iDir++){
    edges.at(iDir) = oneSiteiMPS<complex<double>>( dimMPS, ansaz );
    vector<int> ansLab = ansaz.label();
    ansLab.push_back( ansLab[1] );
    ansLab.erase( ansLab.begin()+1 );
    ansaz = Permute( ansaz, ansLab, 1 );
  }
  
  for (int iDir=iLeft; iDir<=iDown; iDir++){
    findOneBlueI(iDir);
    findOneBlueII(iDir);
    findOneRed(iDir);
  }

}

void channel::findOneBlueI( const int iDir ){
  Network blueI_net("./Networks/blueI.net");
  UniTensor<complex<double>> conjAnsaz = Conj(edges[iDir].getAnsaz());
  ContractArgs( bluesI[iDir], blueI_net, edges[iDir].getAnsaz(), edges[iDir].getAnsaz(), hamiltonian, conjAnsaz, conjAnsaz );

  vector<int> oldLab = bluesI[iDir].label();
  vector<int> combineLab(2);
  for (int i=0; i!=6; i++){
    combineLab[0] = oldLab[i];
    combineLab[1] = oldLab[i+6];
    bluesI[iDir].CombineBond(combineLab);
  }
}

void channel::findOneBlueII( const int iDir ){
  vector<Bond> idBds(4,Bond(BD_IN,dimD));
  idBds.insert( idBds.begin(), Bond(BD_IN,dimPhy) );
  vector<Bond> idOutBds(4,Bond(BD_OUT,dimD));
  idOutBds.insert( idOutBds.begin(), Bond(BD_OUT,dimPhy) );
  idBds.insert( idBds.end(), idOutBds.begin(), idOutBds.end() );
  UniTensor<double> iden(idBds);
  iden.Identity();
  Network blueII_net("./Networks/blueII.net");
  UniTensor<complex<double>> conjAnsaz = Conj(edges[iDir].getAnsaz());
  ContractArgs( bluesII[iDir], blueII_net, edges[iDir].getAnsaz(), edges[iDir].getAnsaz(), hamiltonian, iden, conjAnsaz );

  vector<int> oldLab = bluesII[iDir].label();
  vector<int> combineLab(2);
  for (int i=0; i!=6; i++){
    combineLab[0] = oldLab[i];
    combineLab[1] = oldLab[i+6];
    bluesII[iDir].CombineBond(combineLab);
  }
}

void channel::findOneRed( const int iDir ){
  vector<Bond> idBds(4,Bond(BD_IN,dimD));
  idBds.insert( idBds.begin(), Bond(BD_IN,dimPhy) );
  vector<Bond> idOutBds(4,Bond(BD_OUT,dimD));
  idOutBds.insert( idOutBds.begin(), Bond(BD_OUT,dimPhy) );
  idBds.insert( idBds.end(), idOutBds.begin(), idOutBds.end() );
  UniTensor<double> iden(idBds);
  iden.Identity();

  Network red_net("./Networks/red.net");
  ContractArgs( reds[iDir], red_net, edges[iDir].getAnsaz(), iden );

  vector<int> oldLab = reds[iDir].label();
  vector<int> combineLab(2);
  for (int i=0; i!=4; i++){
    combineLab[0] = oldLab[i];
    combineLab[1] = oldLab[i+4];
    reds[iDir].CombineBond(combineLab);
  }
}

void channel::findAllEnv( const int iterMax, const double canonicalErr, const double arpErr ){
  //find all mps
  enum{ iLeft, iTop, iRight, iBot };
  for ( int iDir=iLeft; iDir<=iBot; iDir++){
    edges.at(iDir).fixedPoint( iterMax, canonicalErr, arpErr );
  }
  //find all corner
  enum{ iLT, iRT, iRB, iLB };
  for ( int iDir=iLT; iDir<=iLB; iDir++ ){
    findOneCorner( iDir, iterMax, arpErr );
  }
  //find all caps
  vector<complex<double>> EignValues(4);
  printf( "EignValues of different channel\n");
  for ( int iDir=iLeft; iDir<=iBot; iDir++ ){
    findOneCap( iDir, iterMax, arpErr, EignValues[iDir] );
    printf( "%4d%22.14e + %22.14e i\n", iDir, EignValues.at(iDir).real(), EignValues.at(iDir).imag() );
  }
  assert( fabs( (EignValues.at(iLeft) - EignValues.at(iRight))/EignValues.at(iLeft) )<1.0e-13 );
  assert( fabs( (EignValues.at(iTop) - EignValues.at(iBot))/EignValues.at(iTop) )<1.0e-13 );
  edges.at(iLeft) *= 1.0/EignValues.at(iTop);
  edges.at(iTop) *= 1.0/EignValues.at(iLeft);
  complex<double> allNorm = findNorm( );
  corners.at(0) *= 1.0/allNorm;
}

UniTensor<complex<double>> channel::findGradient(){
  //find some auxiliary tensor for gradient
  double expectEnergy = measureEnergy(false).real();
  UniTensor<double> subtract = hamiltonian;
  subtract.Identity();
  hamiltonian += -1.0*expectEnergy*subtract;

  for ( int iDir=0; iDir!=4; iDir++ ){
    findOneRed( iDir );
    findOneBlueI( iDir );
    findOneBlueII( iDir );
    findOneSup( iDir, 0 );
    findOneBlueCap( iDir );
  }

  //collect all terms of gradient
  vector<UniTensor<complex<double>>> gradients(4,edges[0].getAnsaz());
  for ( int iDir=0; iDir!=4; iDir++ ){
    gradients[iDir].Zeros();
    gradients[iDir] += grad0(iDir);
    gradients[iDir] += grad1(iDir);
    gradients[iDir] += grad2(iDir);
    gradients[iDir] += grad3(iDir);
    gradients[iDir] += grad4(iDir);
    gradients[iDir] += grad5(iDir);
    vector<int> gradientLab(5);
    gradientLab[0] = 0;
    gradientLab[1] = (iDir)%4+1;
    gradientLab[2] = (iDir+1)%4+1;
    gradientLab[3] = (iDir+2)%4+1;
    gradientLab[4] = (iDir+3)%4+1;
    gradients[iDir].SetLabel(gradientLab);
    /*
    cout<<grad0(iDir)<<endl;
    cout<<grad1(iDir)<<endl;
    cout<<grad2(iDir)<<endl;
    cout<<grad3(iDir)<<endl;
    cout<<grad4(iDir)<<endl;
    cout<<grad5(iDir)<<endl;
    exit(0);
    */
  }
  UniTensor<complex<double>> totalGradient( gradients[0].bond() );
  totalGradient.Zeros();
  for (int iDir=0; iDir!=4; iDir++){
    totalGradient += Permute( gradients[iDir], {0,1,2,3,4}, 1 );
  }

  hamiltonian += 1.0*expectEnergy*subtract;
  return totalGradient;
}

complex<double> channel::measureEnergy( const bool ifShow ) const {
  vector<complex<double>> energy(4);
  complex<double> average( 0, 0);
  cout<<"inside"<<endl;
  for (int i=0; i!=4; i++){
    energy.at(i) = twoSiteOpExpec( hamiltonian, i );
    average += energy.at(i);
  }
  cout<<"inside"<<endl;
  average *= 0.25;
  if (ifShow){
    printf( "\n" );
    for (int i=0; i!=4; i++){
      printf( "%12.8f+%8.2ei", energy.at(i).real(), energy.at(i).imag() );
    }
    printf( "%12.8f+%8.2ei\n", average.real(), average.imag() );

    //write result
    FILE *result = fopen( "result.dat", "a");
    for (int i=0; i!=4; i++){
      fprintf( result, "%12.8f+%8.2ei", energy.at(i).real(), energy.at(i).imag() );
    }
    fprintf( result, "%12.8f+%8.2ei\n", average.real(), average.imag() );
    fflush( result );
    fclose( result );

  } else {}
  return average;
}

//private function  
void channel::findOneCorner( const int iDir, unsigned int max_iter, double err_tol ){
  const oneSiteiMPS<complex<double>> &relTopEdge = edges[(iDir+1)%4];
  const oneSiteiMPS<complex<double>> &relLeftEdge = edges[iDir];

  UniTensor<complex<double>> target;
  complex<double> EigVal;
  
  UniTensor<complex<double>> ConjaL = Conj( relTopEdge.getAL() );
  UniTensor<complex<double>> ConjaR = Conj( relLeftEdge.getAR() );
  Network target_net("./Networks/matrix_for_corner.net");
  ContractArgs( target, target_net, relTopEdge.getIsometryL(), relLeftEdge.getIsometryR(), relLeftEdge.getEvolOperator(), ConjaL, ConjaR); 
  Matrix<complex<double>> cornerMat;
  arCompEignRight( target.GetBlock(), EigVal, cornerMat, max_iter, err_tol, 0);
  vector<Bond> cornerBds = { Bond( BD_IN, dimMPS ), Bond( BD_IN, dimMPS ) };
  corners[iDir] = UniTensor<complex<double>> (cornerBds);
  corners[iDir].PutBlock( cornerMat );
}

void channel::findOneCap( const int iDir, unsigned int max_iter, double err_tol, complex<double> &EignValue ){
  const oneSiteiMPS<complex<double>>& Ledge = edges[(iDir+3)%4];
  const oneSiteiMPS<complex<double>>& Redge = edges[(iDir+1)%4];
  const oneSiteiMPS<complex<double>>& Cedge = edges[iDir];
  UniTensor<complex<double>> transfer;
  complex<double> EigVal;
  Network transfer_net("./Networks/transferForCap.net");
  //UniTensor<complex<double>> evo = Permute( Cedge.getEvolOperator(), 0);
  ContractArgs( transfer, transfer_net, Ledge.getAR(), Redge.getAL(),Cedge.getEvolOperator());
  Matrix<complex<double>> capMat;
  arCompEignRight( transfer.GetBlock(), EigVal, capMat, max_iter, err_tol, 0);
  EignValue = EigVal;
  vector<Bond> capBds = { Bond( BD_IN, dimMPS ), Bond( BD_IN, dimD*dimD ), Bond( BD_IN, dimMPS ) };
  caps[iDir] = UniTensor<complex<double>> (capBds);
  caps[iDir].PutBlock( capMat );
}

void channel::findOneSup( const int iDir, const double q_y ){
  const oneSiteiMPS<complex<double>> &Ledge = edges[(iDir+3)%4];
  const oneSiteiMPS<complex<double>> &Redge = edges[(iDir+1)%4];
  const oneSiteiMPS<complex<double>> &Cedge = edges[iDir];
  UniTensor<complex<double>> transfer;
  Network transfer_net("./Networks/transferForCap.net");
  ContractArgs( transfer, transfer_net, Ledge.getAR(), Redge.getAL(), Cedge.getEvolOperator());
  transfer = Transpose( transfer );

  Matrix<complex<double>> leftVect;
  complex<double> EigVal;
  int max_iter = 10000;
  double err_tol = 1.0e-14;
  arCompEignRight( transfer.GetBlock(), EigVal, leftVect, max_iter, err_tol, 0);
  
  UniTensor<complex<double>> cap_dg( caps[iDir].bond() );
  cap_dg.PutBlock( leftVect );
  cap_dg.SetLabel( {-1,-2,-3} );
  caps[iDir].SetLabel( {-1,-2,-3});
  UniTensor<complex<double>> norm = Contract( caps[iDir], cap_dg );
  cap_dg *= 1.0/norm[0];
  caps[iDir].SetLabel( {1,2,3});
  UniTensor<complex<double>> capcap = Contract( caps[iDir], cap_dg);
  capcap = Permute( capcap, {-1,-2,-3,1,2,3}, 3 );
  
  complex<double> ci(0,1);
  Matrix<complex<double>> tmp, id;
  id = capcap.GetBlock();
  id.Identity();

  if (q_y==0){
    //tmp = Inverse( I-transfer.GetBlock()+capcap.GetBlock() );
    tmp = Inverse( id-transfer.GetBlock()+capcap.GetBlock() ) -capcap.GetBlock();
    sups[iDir].Assign( transfer.bond() );
    sups[iDir].PutBlock( tmp );
  }
  else {
    // need to check
    tmp = Inverse( id - exp(ci*q_y)*(transfer - capcap).GetBlock() );
    sups[iDir].Assign( transfer.bond() );
    sups[iDir].PutBlock( tmp );
  }

  // test spectral radius
/*
  complex<double> EigVal;
  Matrix<complex<double>> EigVec;
  unsigned int max_iter = 10000;
  double err_tol = 1e-12;
  arCompEignRight( (transfer-capcap).GetBlock(), EigVal, EigVec, max_iter, err_tol, 0);
  cout << EigVal << endl;
*/
}

void channel::findOneBlueCap( const int iDir ){ 
  const oneSiteiMPS<complex<double>> &Ledge = edges[(iDir+3)%4];
  const oneSiteiMPS<complex<double>> &Redge = edges[(iDir+1)%4];
  Network bluecap_net("./Networks/bluecap.net");
  ContractArgs( bluecaps[iDir], bluecap_net, caps[iDir], Ledge.getAR(), Ledge.getAR(), Redge.getAL(), Redge.getAL(), bluesI[iDir], sups[iDir]);
}

UniTensor<complex<double>> channel::grad0( const int iDir ){
  int s0 = iDir;
  int s1 = (iDir+1)%4;
  int s2 = (iDir+2)%4;
  int s3 = (iDir+3)%4;
  UniTensor<complex<double>> grad0;
  Network grad0_net("./Networks/grad0.net");
  ContractArgs( grad0, grad0_net, caps.at(s0), caps.at(s1), caps.at(s2), caps.at(s3), corners.at(s0), corners.at(s1), corners.at(s2), corners.at(s3), edges[s3].getAR(), edges[s1].getAL(), bluesII.at(s0) );
  return grad0;
}

UniTensor<complex<double>> channel::grad1( const int iDir ){
  int s0 = iDir;
  int s1 = (iDir+1)%4;
  int s2 = (iDir+2)%4;
  int s3 = (iDir+3)%4;
  UniTensor<complex<double>> grad1;
  Network grad1_net("./Networks/grad1.net");
  ContractArgs( grad1, grad1_net, bluecaps.at(s0), caps.at(s1), caps.at(s2), caps.at(s3), corners.at(s0), corners.at(s1), corners.at(s2), corners.at(s3), reds[s0]);
  return grad1;
} 

UniTensor<complex<double>> channel::grad2( const int iDir ){
  int s0 = iDir;
  int s1 = (iDir+1)%4;
  int s2 = (iDir+2)%4;
  int s3 = (iDir+3)%4;
  UniTensor<complex<double>> grad2;
  Network grad2_net("./Networks/grad2.net");
  ContractArgs( grad2, grad2_net, edges[s0].getEvolOperator(), bluecaps.at(s3), caps.at(s0), caps.at(s1), caps.at(s2), corners.at(s0), corners.at(s1), corners.at(s2), corners.at(s3), edges[s3].getAR(), edges[s1].getAL(), sups.at(s0), reds[s0]);
  return grad2;
}

UniTensor<complex<double>> channel::grad3( const int iDir ){
  int s0 = iDir;
  int s1 = (iDir+1)%4;
  int s2 = (iDir+2)%4;
  int s3 = (iDir+3)%4;
  UniTensor<complex<double>> grad3;
  Network grad3_net("./Networks/grad2.net");
  ContractArgs( grad3, grad3_net, edges[s0].getEvolOperator(), caps.at(s3), caps.at(s0), bluecaps.at(s1), caps.at(s2), corners.at(s0), corners.at(s1), corners.at(s2), corners.at(s3), edges[s3].getAR(), edges[s1].getAL(), sups.at(s0), reds[s0]);
  return grad3;
}

UniTensor<complex<double>> channel::grad4( const int iDir ){
  int s0 = iDir;
  int s1 = (iDir+1)%4;
  int s2 = (iDir+2)%4;
  int s3 = (iDir+3)%4;
  UniTensor<complex<double>> grad4;
  Network grad4_net("./Networks/grad4.net");
  ContractArgs( grad4, grad4_net, caps.at(s0), caps.at(s1), caps.at(s2), caps.at(s3), corners.at(s0), corners.at(s1), corners.at(s2), corners.at(s3), edges[s3].getAR(), edges[s1].getAL(), edges[s0].getAL(), edges[s2].getAR(), sups.at(s0), bluesI[s3], reds[s0]);
  return grad4; 
}

UniTensor<complex<double>> channel::grad5( const int iDir ){
  int s0 = iDir;
  int s1 = (iDir+1)%4;
  int s2 = (iDir+2)%4;
  int s3 = (iDir+3)%4;
  UniTensor<complex<double>> grad5;
  Network grad5_net("./Networks/grad5.net");
  ContractArgs( grad5, grad5_net, caps.at(s0), caps.at(s1), caps.at(s2), caps.at(s3), corners.at(s0), corners.at(s1), corners.at(s2), corners.at(s3), edges[s3].getAR(), edges[s1].getAL(), edges[s0].getAR(), edges[s2].getAL(), sups.at(s0), bluesI.at(s1), reds[s0]);
  return grad5;
}

complex<double> channel::findNorm( ){
  UniTensor<complex<double>> Norm;
  Network Norm_net("./Networks/normForChannel.net");
  ContractArgs( Norm, Norm_net, corners.at(0), corners.at(1), corners.at(2), corners.at(3), caps.at(0), caps.at(1), caps.at(2), caps.at(3), edges.at(0).getEvolOperator() );
  return Norm[0];
}

UniTensor<double> channel::twoSiteCluster( UniTensor<double> twoSiteOp, const int iDir ) const {
  assert( (iDir>=0)&&(iDir<=3) );
  UniTensor<double> tempAnsaz = ansaz;
  UniTensor<double> cluster;
  Network twoSiteCluster_net("./Networks/twoSiteOpCluster.net");

  vector<int> oldLab = tempAnsaz.label();
  for ( int i=0; i!=iDir; i++){
    oldLab.push_back( oldLab.at(1) );
    oldLab.erase( oldLab.begin()+1 );
  }
  tempAnsaz = Permute( tempAnsaz, oldLab, 1 );

  ContractArgs( cluster, twoSiteCluster_net, tempAnsaz, tempAnsaz, twoSiteOp, tempAnsaz, ansaz );
  combineTwoLayer( cluster );
  return cluster;
}

complex<double> channel::twoSiteOpExpec( const UniTensor<double> &twoSiteOp, const int iDir ) const {
  UniTensor<complex<double>> expectation;
  Network twoSiteOpExpec_net("./Networks/twoSiteOpExpec.net");
  UniTensor<double> cluster = twoSiteCluster( twoSiteOp, iDir );
  UniTensor<complex<double>> clusterC( cluster.bond() );
  clusterC.PutBlock( cMatrix( cluster.GetBlock() ) );
  //ContractArgs( expectation, twoSiteOpExpec_net, corners.at(0), corners.at(1), corners.at(2), corners.at(3), caps.at(0), caps.at(1), caps.at(2), caps.at(3), edges.at(1).getAL(), edges.at(3).getAR(), clusterC );
  ContractArgs( expectation, twoSiteOpExpec_net, corners.at(iDir), corners.at((iDir+1)%4), corners.at((iDir+2)%4), corners.at((iDir+3)%4), caps.at(iDir), caps.at((iDir+1)%4), caps.at((iDir+2)%4), caps.at((iDir+3)%4), edges.at((iDir+1)%4).getAL(), edges.at((iDir+3)%4).getAR(), clusterC );
  //ContractArgs( expectation, twoSiteOpExpec_net, cornerLT, cornerRT, cornerRB, cornerLB, capLeft, capTop, capRight, capBot, TopEdge.getAL(), BotEdge.getAR(), clusterC );
  return expectation[0];
}

void channel::SaveEnvTensors( const string &rootFolder ){
  mkdir( rootFolder.c_str(), 0755 );
  for (int i=0; i!=4; i++){
    char temp[4];
    sprintf( temp, "%1d/", i );
    string subFolder = rootFolder+"edge"+string(temp);
    mkdir( subFolder.c_str(), 0755 );
    edges.at(i).SaveTensors( subFolder );
  }
  SaveVecOfTensors( caps, rootFolder, "caps");
  SaveVecOfTensors( corners, rootFolder, "corners");
}

void channel::LoadEnvTensors( const string &rootFolder ){
  for (int i=0; i!=4; i++){
    char temp[4];
    sprintf( temp, "%1d/", i );
    string subFolder = rootFolder+"edge"+string(temp);
    edges.at(i).LoadTensors( subFolder );
  }
  LoadVecOfTensors( caps, rootFolder, "caps" );
  LoadVecOfTensors( corners, rootFolder, "corners" );
}

UniTensor<double> channel::numericalGradientReal( const double delta, const int iterMax, const double canonicalErr, const double arpErr  ){
  complex<double> energyOld = measureEnergy( false );
  int number = ansaz.ElemNum();

  vector<double> gradElem(number);
  UniTensor<complex<double>> originalAnsaz = ansaz;
  for (int i=0; i!=number; i++){
    ansaz = originalAnsaz;
    vector<double> elem( number, 0);
    elem[i] = delta;
    UniTensor<double> addOn( ansaz.bond() );
    addOn.SetRawElem( elem );
    addTensor( addOn );
    findAllEnv( iterMax, canonicalErr, arpErr );
    complex<double> energyNew = measureEnergy( false );
    gradElem[i] = (energyNew-energyOld).real()/delta;
  }

  ansaz = originalAnsaz;
  findAllEnv( iterMax, canonicalErr, arpErr );

  UniTensor<double> gradient( ansaz.bond() );
  gradient.SetRawElem( gradElem );
  return gradient;
} 


UniTensor<complex<double>> channel::numericalGradientComplex( const double delta, const int iterMax, const double canonicalErr, const double arpErr  ){
  complex<double> energyOld = measureEnergy( false );
  int number = ansaz.ElemNum();

  vector<complex<double>> gradElem(number);
  UniTensor<complex<double>> originalAnsaz = ansaz;
  for (int i=0; i!=number; i++){
    ansaz = originalAnsaz;
    vector<double> elem( number, 0);
    elem[i] = delta;
    UniTensor<double> addOn( ansaz.bond() );
    addOn.SetRawElem( elem );
    addTensor( addOn );
    findAllEnv( iterMax, canonicalErr, arpErr );
    complex<double> energyNew = measureEnergy( false );
    gradElem[i] = (energyNew-energyOld)/delta;

    ansaz = originalAnsaz;
    vector<complex<double>> elemc( number, complex<double> (0,0));
    elemc[i] = complex<double> (0,delta);
    UniTensor<complex<double>> addC( ansaz.bond() );
    addC.SetRawElem( elemc );
    addTensor( addC );
    findAllEnv( iterMax, canonicalErr, arpErr );
    energyNew = measureEnergy( false );
    gradElem[i].imag( (energyNew-energyOld).real()/delta );
  }

  findAllEnv( iterMax, canonicalErr, arpErr );
  UniTensor<complex<double>> gradient( ansaz.bond() );
  gradient.SetRawElem( gradElem );
  return gradient;
} 

void channel::normalize( ){
  ansaz *= 1.0/Norm( ansaz.GetBlock() );
}

void channel::addTensor( const UniTensor<complex<double>> &addOn ){
  ansaz = ansaz+addOn;
 // ansaz *= 1.0/maxAbsElem( ansaz );

  enum{ iLeft, iTop, iRight, iDown };
  for (int iDir=iLeft; iDir<=iDown; iDir++){
    edges.at(iDir) = oneSiteiMPS<complex<double>>( dimMPS, ansaz );
    vector<int> ansLab = ansaz.label();
    ansLab.push_back( ansLab[1] );
    ansLab.erase( ansLab.begin()+1 );
    ansaz = Permute( ansaz, ansLab, 1 );
  }
  
  /*
  for (int iDir=iLeft; iDir<=iDown; iDir++){
    findOneBlueI(iDir);
    findOneBlueII(iDir);
    findOneRed(iDir);
  }
  */

}

/*
void channel::testNormalization(){
  cout<<"norm: "<<setprecision(8)<<findNorm();
}
*/

void channel::test(){
  //find some auxiliary tensor for gradient
  double expectEnergy = measureEnergy(false).real();
  //complex<double> expectEnergy = measureEnergy(false);
  UniTensor<double> subtract = hamiltonian;
  subtract.Identity();
  hamiltonian += -1.0*expectEnergy*subtract;

  enum{ iLeft, iTop, iRight, iBot };
  for ( int iDir=iLeft; iDir<=iBot; iDir++ ){
    findOneRed( iDir );
    findOneBlueI( iDir );
    findOneBlueII( iDir );
    findOneSup( iDir, 0 );
    findOneBlueCap( iDir );
  }

  /*
  //collect all terms of gradient
  vector<UniTensor<complex<double>>> gradients(4,edges[0].getAnsaz());
  for ( int iDir=iLeft; iDir<=iBot; iDir++ ){
    gradients[iDir].Zeros();
    gradients[iDir] += grad0(iDir);
    gradients[iDir] += grad1(iDir);
    gradients[iDir] += grad2(iDir);
    gradients[iDir] += grad3(iDir);
    gradients[iDir] += grad4(iDir);
    gradients[iDir] += grad5(iDir);
    vector<int> gradientLab(5);
    gradientLab[0] = 0;
    gradientLab[1] = (iDir)%4+1;
    gradientLab[2] = (iDir+1)%4+1;
    gradientLab[3] = (iDir+2)%4+1;
    gradientLab[4] = (iDir+3)%4+1;
    gradients[iDir].SetLabel(gradientLab);
  }
  UniTensor<complex<double>> totalGradient =  gradients[0];
  totalGradient.Zeros();
  for (int iDir=iLeft; iDir<=iBot; iDir++){
    totalGradient += Permute( gradients[iDir], {0,1,2,3,4}, 1 );
  }
  */
  hamiltonian += 1.0*expectEnergy*subtract;
}
