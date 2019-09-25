#include "quant.h"

using namespace std;

void quantities::printAll( const string &dataFolder, const double field, const int dimD, const int dimCTM )
{
  calculate();
  printEnergy ( dataFolder, field, dimD, dimCTM );
  printMag    ( dataFolder, field, dimD, dimCTM );
}

void quantities::calculate()
{
  magAvg = 0;
  for (int i=0; i!=3; i++){
    magEach[i] = sqrt( spinComp[i][0]*spinComp[i][0] + spinComp[i][1]*spinComp[i][1] );
    magAvg += magEach[i];
  }
  magAvg /= 3.0;
  for (int i=0; i!=3; i++){
    int j = (i+1)%3;
    spinAngleCos[i] = (spinComp[i][0]*spinComp[j][0]+spinComp[i][1]*spinComp[j][1])/(magEach[i]*magEach[j]);
  }
  double upper = max( energyUp, energyDn );
  energyDn = min( energyUp, energyDn );
  energyUp = upper;
  energyDiff = fabs( energyDn - energyUp );
  ePerSite = (energyUp+energyDn)/3.0;
}

void quantities::writeColName( const string &dataFolder, const double field, const int dimD ){
  char buffer[32];
  sprintf( buffer, "J%04i/", int(round(1000*field)) );
  string fieldSubFolder = dataFolder + string(buffer);
  mkdir( fieldSubFolder.c_str(), 0755);

  sprintf( buffer, "D%03i/", dimD );
  string dimSubFolder = fieldSubFolder + string( buffer );
  mkdir( dimSubFolder.c_str(), 0755);

  FILE *fenergy = fopen( (dimSubFolder+"e.dat").c_str(), "a");
  fprintf( fenergy, "%4s%6s%12s%12s%12s%14s%12s%12s%12s%10s\n", "chi", "ctm", "energyUp", "energyDn", "energyDiff", "e", "trunErrEvol", "trunErrMPS", "M", "time" );
  fflush( fenergy );
  fclose( fenergy );

  FILE *fmagAll = fopen( (dimSubFolder+"magAll.dat").c_str(), "a");
  fprintf( fmagAll, "%4s%6s", "chi", "ctm" );
  for (int iSite=0; iSite!=3; iSite++){
    fprintf( fmagAll, "%9dmag", iSite );
  }
  fprintf( fmagAll, "%12s\n", "M");
  fflush( fmagAll );
  fclose( fmagAll );
}

void quantities::printEnergy( const string &dataFolder, const double field, const int dimD, const int dimCTM ){
  char buffer[32];
  sprintf( buffer, "J%04i/", int(round(1000*field)) );
  string fieldSubFolder = dataFolder + string(buffer);
  mkdir( fieldSubFolder.c_str(), 0755);

  sprintf( buffer, "D%03i/", dimD );
  string dimSubFolder = fieldSubFolder + string( buffer );
  mkdir( dimSubFolder.c_str(), 0755);

  double trunErrEvol=0, trunErrMPS=0;

  FILE *fenergy = fopen( (dimSubFolder+"e.dat").c_str(), "a");
  printf( "%4s%6s%12s%12s%12s%14s%12s%12s%12s%10s\n", "chi", "ctm", "energyUp", "energyDn", "energyDiff", "e", "trunErrEvol", "trunErrMPS", "M", "time" );
  printf( "%4d%6d%12.7f%12.7f%12.4e%14.9f%12.4e%12.4e%12.7f%10.2e\n", dimD, dimCTM, energyUp, energyDn, energyDiff, ePerSite, trunErrEvol, trunErrMPS, magAvg, elapsedT );
  fflush(stdout);
  fprintf( fenergy, "%4d%6d%12.7f%12.7f%12.4e%14.9f%12.4e%12.4e%12.7f%10.2e", dimD, dimCTM, energyUp, energyDn, energyDiff, ePerSite, trunErrEvol, trunErrMPS, magAvg, elapsedT );
  fprintf( fenergy, "\n" );
  fflush( fenergy );
  fclose( fenergy );
}

void quantities::printMag( const string &dataFolder, const double field, const int dimD, const int dimCTM ){
  char buffer[32];
  sprintf( buffer, "J%04i/", int(round(1000*field)) );
  string fieldSubFolder = dataFolder + string(buffer);
  mkdir( fieldSubFolder.c_str(), 0755);

  sprintf( buffer, "D%03i/", dimD );
  string dimSubFolder = fieldSubFolder + string( buffer );
  mkdir( dimSubFolder.c_str(), 0755);

  printf( "%4s%6s", "chi", "ctm" );
  for (int iSite=0; iSite!=3; iSite++){
    printf( "%11dx%11dz%9dmag", iSite, iSite, iSite );
  }
  printf("%12s\n", "M");
  fflush(stdout);
  printf( "%4d%6d", dimD, dimCTM );
  for (int iSite=0; iSite!=3; iSite++){
    printf( "%12.7f%12.7f%12.7f", spinComp.at(iSite).at(0), spinComp.at(iSite).at(1), magEach.at(iSite) );
  }
  printf("%12.7f\n", magAvg);
  fflush(stdout);

  FILE *fmagAll = fopen( (dimSubFolder+"magAll.dat").c_str(), "a");
  fprintf( fmagAll, "%4d%6d", dimD, dimCTM );
  for (int iSite=0; iSite!=3; iSite++){
    fprintf( fmagAll,"%12.7f", magEach.at(iSite) );
  }
  fprintf( fmagAll,"%12.7f\n", magAvg);
  fflush( fmagAll );
  fclose( fmagAll );

}

