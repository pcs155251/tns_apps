#ifndef __QUANT_H_INCLUDED__
#define __QUANT_H_INCLUDED__

#include <sys/stat.h> //for mkdir
#include <vector>
#include <stdio.h>
#include <math.h>
#include <string>

class quantities{
  public:
    quantities( ): spinComp( 3, std::vector<double> (2,0) ), magEach( std::vector<double> (3,0) ), spinAngleCos( std::vector<double> (3,0) ) { }
    std::vector<std::vector<double>> spinComp;
    std::vector<double> magEach;
    std::vector<double> spinAngleCos;
    double ePerSite=0, energyUp=0, energyDn=0, energyDiff=0, magAvg=0, elapsedT=0;
    void printAll( const std::string &dataFolder, const double field, const int dimD, const int dimCTM );
    void writeColName( const std::string &dataFolder, const double field, const int dimD );
  private:
    void calculate();
    void printEnergy ( const std::string &dataFolder, const double field, const int dimD, const int dimCTM );
    void printMag    ( const std::string &dataFolder, const double field, const int dimD, const int dimCTM );
};

#endif
