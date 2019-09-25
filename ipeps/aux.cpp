using namespace std;
using namespace uni10;

template<typename T>
void meaSimpTimQuants( ipeps<T>& tim, UniTensor<T>& ham, double field, string &dataPath )
{
  char buffer[16];
  sprintf( buffer, "h%04i/", int(round(1000*field)) );
  string path = dataPath + string(buffer);
  mkdir( path.c_str(), 0755 );

  sprintf( buffer, "D%02i/", tim.getVirDim() );
  path = path + string(buffer);
  mkdir( path.c_str(), 0755 );

  int isite = 1;
  Matrix<double> sigmax = 2.0*matSx();
  Matrix<double> sigmaz = 2.0*matSz();
  UniTensor<double> tenX( sigmax );
  UniTensor<double> tenZ( sigmaz );

  T eAver = 0;
  vector<T> ebond(4);
  for (int i=0; i!=4; ++i )
  {
    ebond[i] = tim.meaSimpTwoSiteOp( ham, i );
    eAver += ebond[i];
  }
  eAver *= 0.25;

  T mx = tim.meaSimpOneSiteOp( tenX, isite );
  T mz = tim.meaSimpOneSiteOp( tenZ, isite );

  ifstream check0( path+"all.dat" );
  if ( check0.peek() == ifstream::traits_type::eof() )
  {
    FILE *file = fopen( (path+"all.dat").c_str(), "w" );
    fprintf( file, "%6s%14s%14s%14s%14s%14s%14s%14s\n", "field", "eAver", "mx", "mz", "e0", "e1", "e2", "e3" );
    fflush( stdout );
    fclose( file );
  } else {}

  FILE *file = fopen( (path+"all.dat").c_str(), "a" );
  fprintf( file, "%6.2f%14.8f%14.8f%14.8f%14.8f%14.8f%14.8f%14.8f\n", field, eAver, mx, mz, ebond[0], ebond[1], ebond[2], ebond[3] );
  fflush( stdout );
  fclose( file );

  printf( "%6s%14s%14s%14s%14s%14s%14s%14s\n", "field", "eAver", "mx", "mz", "e0", "e1", "e2", "e3" );
  printf( "%6.2f%14.8f%14.8f%14.8f%14.8f%14.8f%14.8f%14.8f\n", field, eAver, mx, mz, ebond[0], ebond[1], ebond[2], ebond[3] );
  fflush( stdout );
}

template<typename T>
void meaFullTimQuants( ipeps<T>& tim, UniTensor<T>& ham, double field, string &dataPath )
{
  char buffer[16];
  sprintf( buffer, "h%04i/", int(1000*field) );
  string path = dataPath + string(buffer);
  mkdir( path.c_str(), 0755 );

  sprintf( buffer, "D%02i/", tim.getVirDim() );
  path = path + string(buffer);
  mkdir( path.c_str(), 0755 );

  int isite = 1;
  Matrix<double> sigmax = 2.0*matSx();
  Matrix<double> sigmaz = 2.0*matSz();
  UniTensor<double> tenX( sigmax );
  UniTensor<double> tenZ( sigmaz );

  T eAver = 0;
  vector<T> ebond(4);
  for (int ibond=0; ibond!=4; ++ibond )
  {
    double normVal = tim.meaTwoSiteNorm( tim.getGammas(), ibond );
    ebond[ibond] = tim.meaTwoSiteExp( tim.getGammas(), ham, ibond )/normVal;
    eAver += 0.25*ebond[ibond];
  }

  T normVal = tim.meaOneSiteNorm( tim.getGammas(), isite );
  T mx = tim.meaOneSiteExp( tim.getGammas(), tenX, isite )/normVal;
  T mz = tim.meaOneSiteExp( tim.getGammas(), tenZ, isite )/normVal;

  ifstream check0( path+"all.dat" );
  if ( check0.peek() == ifstream::traits_type::eof() )
  {
    FILE *file = fopen( (path+"all.dat").c_str(), "w" );
    fprintf( file, "%6s%14s%14s%14s%14s%14s%14s%14s\n", "field", "eAver", "mx", "mz", "e0", "e1", "e2", "e3" );
    fflush( stdout );
    fclose( file );
  } else {}

  FILE *file = fopen( (path+"all.dat").c_str(), "a" );
  fprintf( file, "%6.2f%14.8f%14.8f%14.8f%14.8f%14.8f%14.8f%14.8f\n", field, eAver, mx, mz, ebond[0], ebond[1], ebond[2], ebond[3] );
  fflush( stdout );
  fclose( file );
}

