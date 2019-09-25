#include "oneSiteiMPS.h"
using namespace std;
using namespace uni10;

/*
oneSiteiMPS<T>::oneSiteiMPS( const string &rootFolder, const string &subFolder ): gamma( UniTensor<double> () ), lambda( Matrix<double> () ), evolOperator( UniTensor<double> () ), isometryLeft( UniTensor<double> () ), isometryRight( UniTensor<double> () ),\
gammaLab( vector<int> {-2, 1, 2}), isometryLLab( vector<int> { 3, 1, -1} ), isometryRLab( vector<int> { 4, 2, -3 } ), evolOperLab( vector<int> { -1, -2, -3, -4 } ){
  LoadTensors( rootFolder, subFolder );
}
*/

template<typename T> 
oneSiteiMPS<T>::oneSiteiMPS( ): gammaLab( vector<int> {-2, 1, 2} ){
  //default constructor
}

template<typename T> 
oneSiteiMPS<T>::oneSiteiMPS( const oneSiteiMPS &input ){
  //copy constructor
  //this->evolOperLab = input.evolOperLab;
  this->gammaLab = input.gammaLab;
  this->aC = input.aC;
  this->aL = input.aL;
  this->aR = input.aR;
  this->isometryLeft = input.isometryLeft;
  this->isometryRight = input.isometryRight;
  this->cmat = input.cmat;
  this->evolComplex = input.evolComplex;
}

template<typename T> 
oneSiteiMPS<T>::oneSiteiMPS( const int horizDimIn, const UniTensor<complex<double>> &ansazIn ): horizDim{horizDimIn}, gammaLab( vector<int> {-2, 1, 2}) {
//oneSiteiMPS<T>::oneSiteiMPS( const int horizDimIn, const UniTensor<double> &evolOperatorIn ): horizDim{horizDimIn}, gammaLab( vector<int> {-2, 1, 2}) {
  assert( ansazIn.BondNum() == 5 );

  //ansaz = UniTensor<complex<double>> ( ansazIn.bond() );
  //ansaz.PutBlock( cMatrix( ansazIn.GetBlock() ) );
  ansaz = ansazIn;

  //evolComplex = UniTensor<complex<double>> ( evolOperatorIn.bond() );
  //evolComplex.PutBlock( cMatrix( evolOperatorIn.GetBlock() ) );
  vector<int> tracing_posits = { 0 };
  evolComplex = ContractDagger( ansaz, tracing_posits );
  combineTwoLayer( evolComplex );
  vertiDim = evolComplex.bond(0).dim();
  evolComplex.SetLabel({0,1,2,3});

  vector<Bond> edgeBond{ Bond(BD_IN, vertiDim), Bond(BD_OUT, horizDimIn), Bond(BD_OUT, horizDimIn) };
  aC = UniTensor<complex<double>> ( edgeBond );
  aC.Randomize();
  aC.SetLabel( gammaLab );

  vector<Bond> aLbds = { Bond( BD_IN, vertiDim ), Bond( BD_IN, horizDim ), Bond( BD_OUT, horizDim ) };
  vector<Bond> aRbds = { Bond( BD_IN, horizDim ), Bond( BD_OUT, vertiDim ), Bond( BD_OUT, horizDim ) };
  aL = UniTensor<complex<double>> (aLbds);
  aR = UniTensor<complex<double>> (aRbds);

  vector<Bond> isoLefBds = { Bond( BD_IN, horizDim ), Bond( BD_OUT, vertiDim ), Bond( BD_OUT, horizDim ) };
  vector<Bond> isoRigBds = { Bond( BD_IN, horizDim ), Bond( BD_IN, vertiDim ), Bond( BD_OUT, horizDim ) };
  isometryLeft = UniTensor<complex<double>> (isoLefBds);
  isometryRight = UniTensor<complex<double>> (isoRigBds);
  isometryLeft.Randomize();
  isometryRight.Randomize();
}

template<typename T> 
void oneSiteiMPS<T>::SaveTensors( const string &rootFolder ){
  mkdir( rootFolder.c_str(), 0755 );
  aC.Save( rootFolder + "aC.ten" );
  aL.Save( rootFolder + "aL.ten" );
  aR.Save( rootFolder + "aR.ten" );
  isometryLeft.Save( rootFolder + "isoL.ten" );
  isometryRight.Save( rootFolder + "isoR.ten" );
  cmat.Save( rootFolder + "cmat.ten" );
  evolComplex.Save( rootFolder + "evolComplex.ten" );
}

template<typename T> 
void oneSiteiMPS<T>::LoadTensors( const string &rootFolder ){
  string folder = rootFolder;

  aC.Load( rootFolder + "aC.ten" );
  aL.Load( rootFolder + "aL.ten" );
  aR.Load( rootFolder + "aR.ten" );
  isometryLeft.Load( rootFolder + "isoL.ten" );
  isometryRight.Load( rootFolder + "isoR.ten" );
  cmat.Load( rootFolder + "cmat.ten" );
  evolComplex.Load( rootFolder + "evolComplex.ten" );
  horizDim = aC.bond(1).dim();
  vertiDim = aC.bond(0).dim();
}

template<typename T> 
void oneSiteiMPS<T>::setEvolOperator( const UniTensor<double> &evolOperatorIn ){
  //evolOperator = evolOperatorIn;
  vertiDim = evolOperatorIn.bond(0).dim();

  evolComplex = UniTensor<complex<double>> ( evolOperatorIn.bond() );
  evolComplex.PutBlock( cMatrix( evolOperatorIn.GetBlock() ) );
}

template<typename T> 
void oneSiteiMPS<T>::canonical( const int maxIter, const double errTor ){
  bool if_canonical_converge = false;
  //to be implement: Arnoldi
  aC = Permute( aC, { 1, -2, 2}, 1 );
  vector<Matrix<complex<double>>> lqMatrix = Lq( aC.GetBlock() );
  Matrix<complex<double>> aRmat = lqMatrix.at(1);
  Matrix<complex<double>> cold = lqMatrix.at(0);
  cold *= 1.0/Norm( cold );

  const int ElemNum = cold.col()*cold.row();
  Matrix<complex<double>> aLmat;

  printf( "finding canonical form\n" );
  for ( int i=0; i!=maxIter; i++){
    Matrix<complex<double>> lform = Permute( aC, { -2, 1, 2}, 2 ).GetBlock();
    vector<Matrix<complex<double>>> qrMatrix = qrPositveR( lform );
    cmat = qrMatrix.at(1);
    aLmat = qrMatrix.at(0);
    cmat *= 1.0/Norm( cmat );
    Matrix<complex<double>> cform = Dot( cmat, aRmat );
    aC.PutBlock( cform );
    double diff = Norm( cmat + (complex<double>(-1.0))*cold );
    printf( "%6s%8d%6s%12.4e\n", "step", i, "diff", diff );
    if ( diff > errTor ){
      //keep iterate
      cold = cmat;
    } 
    else {
      if_canonical_converge = true;
      break;
    }
  }
  printf( "\n" );

  assert( if_canonical_converge == true );
  aC = Permute( aC, gammaLab, 1 );
  aL.PutBlock(aLmat);
  aR.PutBlock(aRmat);
}

template<typename T> 
void oneSiteiMPS<T>::printTensors( ){
}

template<typename T> 
void oneSiteiMPS<T>::testCanonical( ){
  aC = Permute( aC, { -2, 1, 2}, 2 );
  Matrix<complex<double>> aLmat = aL.GetBlock();
  Matrix<complex<double>> aRmat = aR.GetBlock();

  cout<<"check c*aR =? aL*c"<<endl;
  Matrix<complex<double>> right = aC.GetBlock();
  Matrix<complex<double>> left = Dot( aLmat, cmat );
  printf( "%12.4e\n\n", Norm( right + (complex<double>(-1.0))*left ) );

  cout<<"check left canonical";
  cout<<Dot( Dagger(aLmat), aLmat )<<endl;
  
  cout<<"check right canonical";
  cout<<Dot( aRmat, Dagger(aRmat) )<<endl;

  //recover
  aC = Permute( aC, gammaLab, 1 );
}

template<typename T> 
void oneSiteiMPS<T>::isometries( const double arpackErr ){
  UniTensor<complex<double>> tranMatL, tranMatR;

  UniTensor<complex<double>> DaggerAL = Dagger( aL );
  UniTensor<complex<double>> DaggerAR = Dagger( aR );
  //finding transfer matrix and dominant Eigenvector
  Network tranMatL_net("./Networks/tranMatL.net");
  ContractArgs( tranMatL, tranMatL_net, aL, evolComplex, DaggerAL );
  Network tranMatR_net("./Networks/tranMatR.net");
  ContractArgs( tranMatR, tranMatR_net, aR, evolComplex, DaggerAR );

  const int arpackIter = 10000;
  Matrix<complex<double>> leftVector, righVector;

  //use current as initial guess
  UniTensor<complex<double>> temp = isometryLeft;
  temp.SetLabel({3,1,2});
  temp = Permute( temp, {1,2,3}, 0);
  leftVector = temp.GetBlock();
  temp = isometryRight;
  temp.SetLabel({1,2,3});
  temp = Permute( temp, {1,2,3}, 3);
  righVector = temp.GetBlock();

  arCompEignRight( tranMatR.GetBlock(), lambdaR, righVector, arpackIter, arpackErr, 1 );
  arCompEignLeft(  tranMatL.GetBlock(), lambdaL, leftVector, arpackIter, arpackErr, 1 );
  vector<Bond> fLbonds = { Bond( BD_OUT, horizDim ), Bond( BD_OUT, vertiDim ), Bond( BD_OUT, horizDim ) };
  vector<Bond> fRbonds = { Bond( BD_IN, horizDim ), Bond( BD_IN, vertiDim ), Bond( BD_IN, horizDim ) };
  UniTensor<complex<double>> fLten( fLbonds );
  UniTensor<complex<double>> fRten( fRbonds );
  fLten.PutBlock( leftVector );
  fRten.PutBlock( righVector );
  fLten.SetLabel( {1, 2, 3} );
  fRten.SetLabel( {1, 2, 3} );
  isometryLeft = Permute( fLten, { 3, 1, 2}, 1 );
  isometryRight = Permute( fRten, { 1, 2, 3}, 2 );
}

template<typename T> 
void oneSiteiMPS<T>::Normalize( ){
  UniTensor<complex<double>> cten( cmat );
  UniTensor<complex<double>> NormTensor;
  Network NormTensor_net("./Networks/normTensor.net");
  Matrix<complex<double>> cdag = Dagger( cmat );
  UniTensor<complex<double>> cdagten( cdag );
  ContractArgs( NormTensor, NormTensor_net, cten, isometryLeft, isometryRight, cdagten );
  isometryLeft *= 1.0/Norm( NormTensor.GetBlock() );
}

template<typename T> 
void oneSiteiMPS<T>::improveAcC( const int arpackIter, const double arpackErr ){
  Matrix<complex<double>> vectAC, vectC; 
  UniTensor<complex<double>> updateAC;

  //use initial guess
  UniTensor<complex<double>> temp = Permute( aC, 3);
  vectAC = temp.GetBlock();
  temp = UniTensor<complex<double>> (cmat);
  temp = Permute( temp, 2);
  vectC = temp.GetBlock();

  Network updateAC_net("./Networks/updateAC.net");
  ContractArgs( updateAC, updateAC_net, isometryLeft, isometryRight, evolComplex );
  arCompEignRight( updateAC.GetBlock(), lambdaAC, vectAC, arpackIter, arpackErr, 1 );
  complex<double> factor = lambdaAC/lambdaL;
  lambdaAC = lambdaL;
  isometryLeft *= 1.0/factor;

  UniTensor<complex<double>> updateC;
  vector<int> ilLab = isometryLeft.label();
  vector<int> irLab = isometryRight.label();
  isometryLeft.SetLabel( { 4, 1, 3} );
  isometryRight.SetLabel( { 2, 3, 5} );
  updateC = Contract( isometryLeft, isometryRight );
  updateC = Permute( updateC, {4, 5, 1, 2}, 2);
  isometryLeft.SetLabel( ilLab );
  isometryLeft.SetLabel( irLab );
  arCompEignRight( updateC.GetBlock(), lambdaC, vectC, arpackIter, arpackErr, 1 );

  //put into aC and cmat
  aC = Permute( aC, 3);
  aC.PutBlock(vectAC);
  aC = Permute( aC, 1);

  UniTensor<complex<double>> cten( cmat );
  cten = Permute( cten, 2);
  cten.PutBlock(vectC);
  cten = Permute( cten, 1);
  cmat = cten.GetBlock();
}

template<typename T> 
void oneSiteiMPS<T>::fixedPoint( const int maxIter, const double errTor, const double arpackErr ){
  bool if_fixedPoint_converge = false;
  canonical( maxIter, errTor );
  printf( "finding fixed point\n" );
  double delta = 0.1;
  for (int i=0; i!=maxIter; i++ ){
    isometries( delta );

    //Normalize( );

    const int arpackIter = 10000; //can be tune
    improveAcC( arpackIter, delta );

    //update AL AR
    aC = Permute( aC, 2);
    vector<Matrix<complex<double>>> qArA = qrPositveR( aC.GetBlock() );
    vector<Matrix<complex<double>>> qCrC = qrPositveR( cmat );
    aL.PutBlock( Dot( qArA.at(0), Dagger(qCrC.at(0)) ) );
    double rDiff = Norm( qArA.at(1) + (complex<double>)(-1.0)*qCrC.at(1) );

    aC = Permute( aC, {1, -2, 2}, 1);
    vector<Matrix<complex<double>>> lAqA = lqPositveL( aC.GetBlock() );
    vector<Matrix<complex<double>>> lCqC = lqPositveL( cmat );
    aR.PutBlock( Dot( Dagger( lCqC.at(1)), lAqA.at(1) ) );
    double lDiff = Norm( lAqA.at(0) + (complex<double>)(-1.0)*lCqC.at(0) );

    delta = max( lDiff, rDiff );
    
    //recover and determine convergence
    aC = Permute( aC, {-2, 1, 2}, 1 );
    printf( "%8d%12.4e%6s%12.4e\n", i, rDiff, "", lDiff );
    if (delta > arpackErr){ 
      //keep iterate
    } 
    else {
      if_fixedPoint_converge = true;
      break;
    }
  }
  assert( if_fixedPoint_converge == true );
  isometries( arpackErr );

  complex<double> lambdaAC;
  Matrix<complex<double>> vectAC;
  UniTensor<complex<double>> updateAC;
  Network updateAC_net("./Networks/updateAC.net");
  ContractArgs( updateAC, updateAC_net, isometryLeft, isometryRight, evolComplex );
  arCompEignRight( updateAC.GetBlock(), lambdaAC, vectAC, 10000, arpackErr, 0 );
  complex<double> factor = lambdaAC/lambdaL;
  lambdaAC = lambdaL;
  isometryLeft *= 1.0/factor;

  printf( "\n" );
}

template<typename T> 
void oneSiteiMPS<T>::testFixedPoint( ){
  cout<<"test isometries"<<endl;
  UniTensor<complex<double>> tranMatL, tranMatR;
  UniTensor<complex<double>> DaggerAL = Dagger( aL );
  UniTensor<complex<double>> DaggerAR = Dagger( aR );

  Network tranMatL_net("./Networks/tranMatL.net");
  ContractArgs( tranMatL, tranMatL_net, aL, evolComplex, DaggerAL );
  UniTensor<complex<double>> temp = Permute( isometryLeft, {1,2,3}, 0 );
  cout<<"left vector diff: "<<Norm( lambdaL*temp.GetBlock() - Dot( temp.GetBlock(), tranMatL.GetBlock()) )<<endl;

  Network tranMatR_net("./Networks/tranMatR.net");
  ContractArgs( tranMatR, tranMatR_net, aR, evolComplex, DaggerAR );
  temp = Permute( isometryRight, 3);
  cout<<"right vector diff: "<<Norm( lambdaR*temp.GetBlock() - Dot(tranMatR.GetBlock(),temp.GetBlock()) )<<endl;

  cout<<"test Ac"<<endl;
  Matrix<complex<double>> vectAC, vectC; 
  UniTensor<complex<double>> updateAC;
  Network updateAC_net("./Networks/updateAC.net");
  ContractArgs( updateAC, updateAC_net, isometryLeft, isometryRight, evolComplex );
  temp = Permute( aC, 3);
  cout<<"Ac diff: "<<Norm( lambdaAC*temp.GetBlock() - Dot(updateAC.GetBlock(),temp.GetBlock()) )<<endl;

  cout<<"test C"<<endl;
  UniTensor<complex<double>> updateC;
  vector<int> ilLab = isometryLeft.label();
  vector<int> irLab = isometryRight.label();
  isometryLeft.SetLabel( { 4, 1, 3} );
  isometryRight.SetLabel( { 2, 3, 5} );
  updateC = Contract( isometryLeft, isometryRight );
  updateC = Permute( updateC, {4, 5, 1, 2}, 2);
  UniTensor<complex<double>> cten( cmat );
  cten = Permute( cten, 2);
  cout<<"C diff: "<<Norm( lambdaC*cten.GetBlock() - Dot(updateC.GetBlock(),cten.GetBlock()) )<<endl;

  cout<<"L "<<lambdaL<<"R "<<lambdaR<<"C "<<lambdaC<<"Ac "<<lambdaAC<<endl;
}

template class oneSiteiMPS<double>;
template class oneSiteiMPS<complex<double>>;
