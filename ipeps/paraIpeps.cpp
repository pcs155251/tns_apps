#include "paraIpeps.h"

using namespace std;

paraIpeps::paraIpeps( int argc, char **argv ){
  Parser pars;
  pars.Bind("rootFolder",rootFolder);
  pars.Bind("@wvLoadFolder",wvLoadFolder);
  pars.Bind("enLoadFolder",enLoadFolder);
  pars.Bind("@saveFolder",saveFolder);
  pars.Bind("ifloadWv",ifloadWv);
  pars.Bind("ifsu",ifsu);
  pars.Bind("iffu",iffu);
  pars.Bind("ifsaveEn",ifsaveEn);
  pars.Bind("virDim",virDim);
  pars.Bind("edgeDimL",edgeDimL);
  pars.Bind("edgeDimU",edgeDimU);
  pars.Bind("ctmIter",ctmIter);
  pars.Bind("suStop",suStop);
  pars.Bind("suIter",suIter);
  pars.Bind("fuIter",fuIter);
  pars.Bind("fuStop",fuStop);
  pars.Bind("ctmStop",ctmStop);
  pars.Bind("taus",taus);
  pars.Bind("@fieldLoad",fieldLoad);
  pars.Bind("fieldStart",fieldStart);
  pars.Bind("fieldEnd",fieldEnd);
  pars.Bind("fieldDelta",fieldDelta);
  pars.Parse(string(argv[1]));

  cout<<"loading input.rc"<<endl;
  pars.PrintVars();

  double fieldTmp = fieldStart;
  while ( fieldTmp<=fieldEnd )
  {
    fields.push_back( fieldTmp );
    fieldTmp += fieldDelta;
  }

  for ( int i=2; i!=argc; i++ ){
    if ( string(argv[i]) == "-fieldLoad" )
    {
      i++;
      fieldLoad = atof( argv[i] );
    }
    else if ( string(argv[i]) == "-wvLoadFolder" )
    {
      i++;
      wvLoadFolder = string( argv[i] );
    }
    else if ( string(argv[i]) == "-saveFolder" )
    {
      i++;
      saveFolder = string( argv[i] );
    }
    else {}
  }
  cout<<"loading argv"<<endl;
  pars.PrintVars();

  mkdir( rootFolder.c_str(), 0755 );
  mkdir( (rootFolder+"OutputT/").c_str(), 0755 );
  mkdir( (rootFolder+"data/").c_str(), 0755 );
  wvLoadFolder = rootFolder + "OutputT/" + wvLoadFolder;
  enLoadFolder = rootFolder + "OutputT/" + enLoadFolder;
  tensorFolder = rootFolder + "OutputT/" + saveFolder;
  dataFolder = rootFolder + "data/" + saveFolder;
  mkdir( tensorFolder.c_str(), 0755 );
  mkdir( dataFolder.c_str(), 0755 );
}

