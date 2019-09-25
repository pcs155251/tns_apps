#include "ctmBase.h"

using namespace std;
using namespace uni10;

#define IF_DEBUG false

template<typename T>
nineTens<T>::nineTens( const UniTensor<T> &coreIn, const int edgeDimIn ): 
core( coreIn )
,edgeDim{edgeDimIn}
,corners( vector<UniTensor<T>> (4, UniTensor<T> ( vector<Bond> (2,Bond(BD_IN,edgeDimIn))) ) )
,edges( vector<UniTensor<T>> {
                             UniTensor<T> ( vector<Bond> { Bond(BD_IN,edgeDimIn), Bond(BD_IN,edgeDimIn), coreIn.bond(1) } ),
                             UniTensor<T> ( vector<Bond> { Bond(BD_IN,edgeDimIn), Bond(BD_IN,edgeDimIn), coreIn.bond(2) } ),
                             UniTensor<T> ( vector<Bond> { Bond(BD_IN,edgeDimIn), Bond(BD_IN,edgeDimIn), coreIn.bond(3) } ),
                             UniTensor<T> ( vector<Bond> { Bond(BD_IN,edgeDimIn), Bond(BD_IN,edgeDimIn), coreIn.bond(0) } ) } )
{
  assert( core.BondNum()==4 );
  for (int i=0; i!=4; i++){
    corners[i].Randomize();
    edges  [i].Randomize();
  }
}

/*
template<typename T>
UniTensor<T> nineTens<T>::initCorner( const UniTensor<T>& core, int idir, int edgeDim )
{
  assert( idir<=3 );
  assert( idir>=0 );
  
  //todo more sophisicated init
  UniTensor<T> coreTmp = core;
  coreTmp.SetLabel( vector<int> {1,2,3,4} );
  //coreTmp.SetLabel( vector<int> {1,-1,2,-2,3,-3,4,-4} );
  int ib1 = (idir  )%4+1;
  int ib2 = (idir+1)%4+1;
  int ib3 = (idir+2)%4+1;
  int ib4 = (idir+3)%4+1;
  coreTmp = Permute( coreTmp, vector<int> { ib1, -ib1, ib2, -ib2, ib3, -ib3, ib4, -ib4}, 8 );
  int dim1 = coreTmp.bond( 0 ).dim();
  int dim2 = coreTmp.bond( 2 ).dim();
  Matrix<T> mat1( dim1, dim1, true );
  Matrix<T> mat2( dim2, dim2, true );
  mat1.Identity();
  mat2.Identity();
  UniTensor<T> ten1( mat1 );
  UniTensor<T> ten2( mat2 );
  ten1.SetLabel( vector<int> {ib1,-ib1} );
  ten2.SetLabel( vector<int> {ib2,-ib2} );
  UniTensor<T> corner = Contract( coreTmp, ten1, false ); 
  corner = Contract( corner, ten2, false ); 
  corner.CombineBond( vector<int> {ib3, -ib3} );
  corner.CombineBond( vector<int> {ib4, -ib4} );
  changeOneBond( corner, 0, edgeDim );
  changeOneBond( corner, 1, edgeDim );
  UniTensor<T> corner; 
  return corner;
}
*/

/*
template<typename T>
UniTensor<T> nineTens<T>::initEdge( const UniTensor<T>& core, int idir, int edgeDim )
{
  assert( idir<=3 );
  assert( idir>=0 );
  
  // todo more sophisicated init
  UniTensor<T> coreTmp = core;
  coreTmp.SetLabel( vector<int> {1,-1,2,-2,3,-3,4,-4} );
  int ib1 = (idir+1)%4+1;
  int ib2 = (idir+2)%4+1;
  int ib3 = (idir+3)%4+1;
  int ib4 = (idir  )%4+1;
  coreTmp = Permute( coreTmp, vector<int> { ib1, -ib1, ib2, -ib2, ib3, -ib3, ib4, -ib4}, 8 );
  int dim1 = coreTmp.bond( 0 ).dim();
  Matrix<T> mat1( dim1, dim1, true );
  mat1.Identity();
  UniTensor<T> ten1( mat1 );
  ten1.SetLabel( vector<int> {ib1,-ib1} );
  UniTensor<T> edge = Contract( coreTmp, ten1, false ); 
  edge.CombineBond( vector<int> {ib2, -ib2} );
  edge.CombineBond( vector<int> {ib4, -ib4} );
  changeOneBond( edge, 0, edgeDim );
  changeOneBond( edge, 3, edgeDim );
  edge = Permute( edge, vector<int> { ib2, ib4, ib3, -ib3 }, 4 );
  edge.CombineBond( vector<int> {ib3, -ib3} );
  UniTensor<T> edge;
  return edge;
}
*/

template<typename T>
void nineTens<T>::normalize( ){
  for (int i=0; i!=4; i++){
    edges.at(i) *= 1.0/maxAbsElem(edges.at(i));
    corners.at(i) *= 1.0/maxAbsElem(corners.at(i));
  }
}

template<typename T>
void nineTens<T>::changeEdgeDim( const int newDim ){
  edgeDim = newDim;
  for (int i=0; i!=4; i++){
    changeOneBond( edges.at(i), 0, newDim );
    changeOneBond( edges.at(i), 1, newDim );
    changeOneBond( corners.at(i), 0, newDim );
    changeOneBond( corners.at(i), 1, newDim );
  }
}

template<typename T>
void nineTens<T>::saveTensors( const std::string &folder )const {
  mkdir( folder.c_str(), 0755 );
  SaveVecOfTensors( edges, folder, "T");
  SaveVecOfTensors( corners, folder, "C");
}

template<typename T>
bool nineTens<T>::loadTensors( const std::string &folder )
{
  bool loadT = LoadVecOfTensors( edges, folder, "T");
  bool loadC = LoadVecOfTensors( corners, folder, "C");
  edgeDim = corners.at(0).bond(0).dim();
  return loadT&&loadC;
}

template<typename T>
void nineTens<T>::printAll( ) const
{
  core.PrintDiagram();
  for (int i=0; i!=corners.size(); ++i)
  {
    cout<<"C"<<i<<endl;
    corners[i].PrintDiagram();
  }
  for (int i=0; i!=edges.size(); ++i)
  {
    cout<<"T"<<i<<endl;
    edges[i].PrintDiagram();
  }
}

template<typename T>
ctmBase<T>::ctmBase( const int nrowIn, const int ncolIn, const int ctmDimIn ): 
nrow{nrowIn},
ncol{ncolIn},
ctmDim{ctmDimIn}
{
}

template<typename T>
void ctmBase<T>::setGroups( const vector<UniTensor<T>> &unit )
{
  assert( unit.size()==(nrow*ncol) );
  groups.clear();
  groups.reserve(unit.size());
  for (int i=0; i!=unit.size(); ++i )
  {
    groups.push_back( nineTens<T> ( unit[i],ctmDim ) );
  }
}


template<typename T>
ctmBase<T>::ctmBase( const int nrowIn, const int ncolIn, const vector<UniTensor<T>>& unit, const int ctmDimIn ): 
nrow{nrowIn},
ncol{ncolIn},
ctmDim{ctmDimIn}
{
  assert( unit.size()==(nrow*ncol) );
  groups.clear();
  groups.reserve(unit.size());
  ///TODO
  ///generalize the initialization to aritrary unit cell
  for (int i=0; i!=unit.size(); ++i )
  {
    groups.push_back( nineTens<T> ( unit[i],ctmDim ) );
  }
  /*
  groups.push_back( nineTens<T> ( unit[1],ctmDim ) );
  groups.push_back( nineTens<T> ( unit[0],ctmDim ) );
  groups.push_back( nineTens<T> ( unit[0],ctmDim ) );
  groups.push_back( nineTens<T> ( unit[1],ctmDim ) );
  groups[0].core = unit[0];
  groups[1].core = unit[1];
  groups[2].core = unit[2];
  groups[3].core = unit[3];
  */
}

template<typename T>
vector<UniTensor<T>> ctmBase<T>::findIsometryQR( vector<int> order, const int idir ){
  int s0 = idir;
  int s1 = (idir+1)%4;
  int s2 = (idir+2)%4;
  int s3 = (idir+3)%4;

  char buffer[16];
  sprintf( buffer, "%i", idir );
  string netString = string("../ctmBase/ctmBaseNets/isoUp") + string(buffer) + string(".net");
  Network isoUp_net( netString );
  netString = string("../ctmBase/ctmBaseNets/isoDn") + string(buffer) + string(".net");
  Network isoDn_net( netString );

  Matrix<T> pinv, pmat;
  int indim = 0;
  {
    Matrix<T> rUp;
    {
      if (IF_DEBUG){
        cout<<"contract upper half"<<flush<<endl;
      } else {}
      UniTensor<T> all;
      ContractArgs( all, isoUp_net, 
                    groups.at(order.at(0)).core, groups.at(order.at(0)).corners.at(s0), 
                    groups.at(order.at(0)).edges.at(s0), groups.at(order.at(0)).edges.at(s3),
                    groups.at(order.at(1)).core, groups.at(order.at(1)).corners.at(s1), 
                    groups.at(order.at(1)).edges.at(s1), groups.at(order.at(1)).edges.at(s0));
      indim = all.bond(3).dim();
      if (IF_DEBUG){
        cout<<"upper qr"<<flush<<endl;
      } else {}
      {
        vector<Matrix<T>> qrUp = Qr_nm( all.GetBlock() );
        rUp = qrUp[1];
      }
    }

    Matrix<T> rDnTran;
    {
      if (IF_DEBUG){
        cout<<"contract lower half"<<flush<<endl;
      } else {}
      UniTensor<T> all;
      ContractArgs( all, isoDn_net, 
                    groups.at(order.at(3)).core, groups.at(order.at(3)).corners.at(s2), 
                    groups.at(order.at(3)).edges.at(s2), groups.at(order.at(3)).edges.at(s1),
                    groups.at(order.at(2)).core, groups.at(order.at(2)).corners.at(s3), 
                    groups.at(order.at(2)).edges.at(s3), groups.at(order.at(2)).edges.at(s2));
      if (IF_DEBUG){
        cout<<"lower svd"<<flush<<endl;
      } else {}
      {
        vector<Matrix<T>> qrDn = Qr_nm( all.GetBlock() );
        rDnTran = Transpose( qrDn[1] );
      }
    }

    //using full svd
    if (IF_DEBUG){
      cout<<"final svd"<<flush<<endl;
    } else {}
    vector<Matrix<T>> usv = Sdd( Dot(rUp,rDnTran) );
    Resize( usv[0], usv[0].row(),      ctmDim, INPLACE);
    Resize( usv[1],      ctmDim,       ctmDim, INPLACE);
    Resize( usv[2],      ctmDim, usv[2].col(), INPLACE);

    Matrix<T> s_sqrt_inv = usv[1];
    for (int i=0; i!=ctmDim; i++){
      s_sqrt_inv[i] = 1.0/sqrt(usv[1].At(i,i));
    }

    pmat = Dot( s_sqrt_inv, Dot(Dagger(usv[0]),rUp) );
    pinv = Dot( rDnTran, Dot(Dagger(usv[2]),s_sqrt_inv) );
  }

  if (IF_DEBUG){
    cout<<"allocating isometries"<<flush<<endl;
  } else {}
  vector<Bond> bds = {Bond(BD_IN,ctmDim), Bond(BD_OUT,ctmDim), Bond(BD_OUT,indim) };
  UniTensor<T> pten(bds);
  pten.PutBlock(pmat);
  UniTensor<T> pinvTen = Transpose(pten);
  pinvTen.PutBlock(pinv);

  return vector<UniTensor<T>> {pten,pinvTen};
}

template<typename T>
vector<UniTensor<T>> ctmBase<T>::findIsometrySVD( vector<int> order, const int idir ){
  int s0 = idir;
  int s1 = (idir+1)%4;
  int s2 = (idir+2)%4;
  int s3 = (idir+3)%4;

  char buffer[16];
  sprintf( buffer, "%i", idir );
  string netString = string("../ctmBase/ctmBaseNets/isoUp") + string(buffer) + string(".net");
  Network isoUp_net( netString );
  netString = string("../ctmBase/ctmBaseNets/isoDn") + string(buffer) + string(".net");
  Network isoDn_net( netString );

  Matrix<T> pinv, pmat;
  int indim = 0;
  {
    Matrix<T> rUp;
    {
      if (IF_DEBUG){
        cout<<"contract upper half"<<flush<<endl;
      } else {}
      UniTensor<T> all;
      ContractArgs( all, isoUp_net, 
                    groups.at(order.at(0)).core, groups.at(order.at(0)).corners.at(s0), 
                    groups.at(order.at(0)).edges.at(s0), groups.at(order.at(0)).edges.at(s3),
                    groups.at(order.at(1)).core, groups.at(order.at(1)).corners.at(s1), 
                    groups.at(order.at(1)).edges.at(s1), groups.at(order.at(1)).edges.at(s0));
      indim = all.bond(3).dim();
      all *= 1.0/maxAbsElem( all );
      if (IF_DEBUG){
        cout<<"upper svd"<<flush<<endl;
      } else {}
      {
        Matrix<T> uten, sten, vTten;
        Svd( all.ConstGetBlock(), uten, sten, vTten, INPLACE );
        rUp = Dot(sten,vTten);
      }
    }

    Matrix<T> rDnTran;
    {
      if (IF_DEBUG){
        cout<<"contract lower half"<<flush<<endl;
      } else {}
      UniTensor<T> all;
      ContractArgs( all, isoDn_net, 
                    groups.at(order.at(3)).core, groups.at(order.at(3)).corners.at(s2), 
                    groups.at(order.at(3)).edges.at(s2), groups.at(order.at(3)).edges.at(s1),
                    groups.at(order.at(2)).core, groups.at(order.at(2)).corners.at(s3), 
                    groups.at(order.at(2)).edges.at(s3), groups.at(order.at(2)).edges.at(s2));
      all *= 1.0/maxAbsElem( all );
      if (IF_DEBUG){
        cout<<"lower svd"<<flush<<endl;
      } else {}
      {
        Matrix<T> uten, sten, vTten;
        Svd( all.ConstGetBlock(), uten, sten, vTten, INPLACE );
        rDnTran = Transpose( Dot(sten,vTten) );
      }
    }

    //using full svd
    if (IF_DEBUG){
      cout<<"final svd"<<flush<<endl;
    } else {}
    Matrix<T> uten, sten, vTten;
    Sdd( Dot( rUp,rDnTran), uten, sten, vTten, INPLACE );

    double factor = 1.0/sten[0];
    sten *= factor;
    rUp *= sqrt(factor);
    rDnTran *= sqrt(factor);

    Resize( uten , uten.row(),      ctmDim, INPLACE);
    Resize( sten ,     ctmDim,      ctmDim, INPLACE);
    Resize( vTten,     ctmDim, vTten.col(), INPLACE);

    Matrix<T> s_sqrt_inv = sten;
    for (int i=0; i!=ctmDim; i++){
      if ( (sten.At(i,i)/sten.At(0,0))>1.0e-14 )
      {
        s_sqrt_inv[i] = 1.0/sqrt(sten.At(i,i));
      }
      else 
      {
        s_sqrt_inv[i] = 1.0/sqrt(sten.At(i,i));
        //s_sqrt_inv[i] = 0;
      }
    }

    pmat = Dot( s_sqrt_inv, Dot(Dagger(uten),rUp) );
    pinv = Dot( rDnTran, Dot(Dagger(vTten),s_sqrt_inv) );
  }

  if (IF_DEBUG){
    cout<<"allocating isometries"<<flush<<endl;
  } else {}
  vector<Bond> bds = {Bond(BD_IN,ctmDim), Bond(BD_OUT,ctmDim), Bond(BD_OUT,indim) };
  UniTensor<T> pten(bds);
  pten.PutBlock(pmat);
  UniTensor<T> pinvTen = Transpose(pten);
  pinvTen.PutBlock(pinv);

  return vector<UniTensor<T>> {pten,pinvTen};
}

vector<int> buildTable( const int nrow, const int ncol, const int irow, const int icol, const int idir )
{
   //irow icol is the index of rotated matrix, so it should be specify in reverse limit when idir=1 or 3
   assert( idir>=0 && idir<=3 );
   const int nsite = nrow*ncol;
   vector<int> table( nsite );
   int irowp1, icolp1;
   int nrowTmp=nrow;
   int ncolTmp=ncol;
   if (idir%2)
   { 
      nrowTmp = ncol;
      ncolTmp = nrow;
   } else {}
   assert( irow<nrowTmp && icol<ncolTmp );
   irowp1 = (irow+1)%nrowTmp;
   icolp1 = (icol+1)%ncolTmp;
   switch(idir)
   {
     case 0:
       for ( int i=0; i!=nsite; ++i )
       {
          table[i] = i;
       }
       break;
     case 1:
       for ( int ir=0; ir!=nrowTmp; ++ir )
       {
          for ( int ic=0; ic!=ncolTmp; ++ic )
          {
            int index = nrow*ir + ic;
            table[index] = ncol-1-ir+ncol*ic;
          }
       }
       break;
     case 2:
       for ( int i=0; i!=nsite; ++i )
       {
          table[nsite-i-1] = i;
       }
       break;
     case 3:
       for ( int ir=0; ir!=nrowTmp; ++ir )
       {
          for ( int ic=0; ic!=ncolTmp; ++ic )
          {
            int index = nrow*ir + ic;
            table[index] = (nrow-1)*ncol+ir-ncol*ic;
          }
       }
       /*
       for (int ii=0; ii!=nsite; ++ii)
       {
          cout<<table[ii]<<" ";
       }
       exit(0);
       */
       break;
   }
   vector<int> out = { 
                       table[irow*ncolTmp+icol],   table[irow*ncolTmp+icolp1],
                       table[irowp1*ncolTmp+icol], table[irowp1*ncolTmp+icolp1]
                     };
   return out;
}

template<typename T>
void ctmBase<T>::move( const int idir, bool ifsvd )
{
  assert( idir<=3&&idir>=0 );
  int nrowTmp=nrow;
  int ncolTmp=ncol;
  if ( idir%2 )
  {  nrowTmp = ncol;
     ncolTmp = nrow;
  } else {}
  for (int icol=0; icol!=ncolTmp; ++icol )
  {
    //networks
    int s0 = idir;
    int s1 = (idir+1)%4;
    int s2 = (idir+2)%4;
    int s3 = (idir+3)%4;
    Network C0_net( "../ctmBase/ctmBaseNets/C0.net" );
    char buffer[8];
    sprintf( buffer, "%i", idir );
    string netString = string("../ctmBase/ctmBaseNets/T3_") + string(buffer) + string(".net");
    Network T3_net( netString );
    Network C3_net( "../ctmBase/ctmBaseNets/C3.net" );

    //isometries 
    vector<vector<UniTensor<T>>> isometries(nrowTmp);
    vector<vector<int>> tables(nrowTmp);
    for (int irow=0; irow!=nrowTmp; ++irow )
    {
      tables[irow] = buildTable( nrow, ncol, irow, icol, idir ); 
      if (ifsvd)
      {
        isometries[irow] = findIsometrySVD( tables[irow], idir );
      }
      else
      {
        isometries[irow] = findIsometryQR( tables[irow], idir );
      }
    }

    //merge tensors
    for (int irow=0; irow!=nrowTmp; ++irow )
    {
      int oldSite = tables[irow][0];
      int newSite = tables[irow][1];
      int irowp1 = (irow+1)%nrowTmp;
      int irowm1 = (irow+nrowTmp-1)%nrowTmp;
      ContractArgs( groups[newSite].corners[s0], C0_net, groups[oldSite].corners[s0], groups[oldSite].edges[s0], isometries[irowm1][1] );
      ContractArgs( groups[newSite].edges[s3],   T3_net, isometries[irowm1][0], groups[oldSite].edges[s3], groups[oldSite].core, isometries[irow][1] );
      ContractArgs( groups[newSite].corners[s3], C3_net, isometries[irow][0], groups[oldSite].corners[s3], groups[oldSite].edges[s2] );
    }
  }
}

template<typename T>
void ctmBase<T>::ctmrgAstep( bool ifsvd ){
  move( 0, ifsvd );
  move( 2, ifsvd );
  move( 1, ifsvd );
  move( 3, ifsvd );
  for (int i=0; i!=nrow*ncol; i++){
    groups[i].normalize();
  }
}

template<typename T>
void ctmBase<T>::changeEdgeDim( const int newDim ){
  cout<<"change edge dim"<<flush<<endl;
  ctmDim = newDim;
  for (int i=0; i!=groups.size(); ++i){
    groups.at(i).changeEdgeDim( newDim );
  }
}

template<typename T>
void ctmBase<T>::saveTensors( const std::string &rootFolder ) const
{

  char buffer[16];
  sprintf( buffer, "chi%04i/", ctmDim );
  string chiSubFolder = rootFolder + string( buffer ); 
  mkdir( chiSubFolder.c_str(), 0755 );
  for (int i=0; i!=groups.size(); i++){
    sprintf( buffer, "groups%i/", i );
    string subFolder = chiSubFolder + string( buffer );
    mkdir( subFolder.c_str(), 0755 );
    groups.at(i).saveTensors( subFolder );
  }
}

template<typename T>
bool ctmBase<T>::loadTensors( const std::string &rootFolder )
{
  bool if_load = true;

  char buffer[16];
  sprintf( buffer, "chi%04i/", ctmDim );
  string chiSubFolder = rootFolder + string( buffer ); 
  for (int i=0; i!=groups.size(); i++){
    sprintf( buffer, "groups%i/", i );
    string subFolder = chiSubFolder + string( buffer );
    if_load = if_load&&groups.at(i).loadTensors( subFolder );
  }
  return if_load;
}

template<typename T>
void ctmBase<T>::updateUnit( const vector<UniTensor<T>>& unit )
{
  cout<<groups.size()<<" "<<unit.size()<<endl;
  assert( groups.size() == unit.size() );
  for (int i=0; i!=groups.size(); ++i )
  {
    groups[i].core = unit[i];
  }
}

template<typename T>
void ctmBase<T>::printAll( )
{
  for (int i=0; i!=groups.size(); ++i )
  {
    cout<<"group "<<i<<endl;
    groups[i].printAll();
  }
}

template struct nineTens<double>;
template class ctmBase<double>;
