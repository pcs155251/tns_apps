#include "paraPess.h"

using namespace std;

paraPess::paraPess( int argc, char* argv[] ){
  Parser pars;
  pars.Bind("rootFolder",rootFolder);
  pars.Bind("@saveFolder",saveFolder);
  pars.Bind("loadFolder",loadFolder);
  pars.Bind("@enviFolder",enviFolder);
  pars.Bind("@ctmStep",ctmStep);
  pars.Bind("@edgeDimN",edgeDimN);
  pars.Bind("@edgeErr",edgeErr);
  pars.Bind("@start",start);
  pars.Bind("@slope",slope);
  pars.Bind("@field",field);
  pars.Bind("@dimDload",dimDload);
  pars.Bind("@dimDs",dimDs);
  pars.Bind("@ifsvd",ifsvd);
  pars.Bind("ifEn",ifEn);
  pars.Bind("ifGs",ifGs);
  pars.Bind("gsIter",gsIter);
  pars.Bind("gsErr",gsErr);
  pars.Bind("taus",taus);
  pars.Parse(string(argv[1]));

  cout<<"loading input.rc"<<endl;
  pars.PrintVars();

  int dimDstart=0, dimDend=0;
  for ( int i=2; i!=argc; i++ ){
    if ( string(argv[i]) == "-ctmStep" ){
      i++;
      ctmStep = atoi( argv[i] );
    }
    else if ( string(argv[i]) == "-edgeDimN" ){
      i++;
      edgeDimN = atoi( argv[i] );
    }
    else if ( string(argv[i]) == "-edgeErr" ){
      i++;
      edgeErr = atof( argv[i] );
    }
    else if ( string(argv[i]) == "-start" ){
      i++;
      start = atof( argv[i] );
    }
    else if ( string(argv[i]) == "-slope" ){
      i++;
      slope = atof( argv[i] );
    }
    else if ( string(argv[i]) == "-field" ){
      i++;
      field = atof( argv[i] );
    }
    else if ( string(argv[i]) == "-dimDload" ){
      i++;
      dimDload = atoi( argv[i] );
    }
    else if ( string(argv[i]) == "-dimDstart" ){
      i++;
      dimDstart = atoi( argv[i] );
    }
    else if ( string(argv[i]) == "-dimDend" ){
      i++;
      dimDend = atoi( argv[i] );
    }
    else if ( string(argv[i]) == "-saveFolder" ){
      i++;
      saveFolder = argv[i];
    }
    else if ( string(argv[i]) == "-enviFolder" ){
      i++;
      enviFolder = argv[i];
    }
    else if ( string(argv[i]) == "-ifsvd" ){
      i++;
      ifsvd = atoi( argv[i] );
    }
    else {}
  }
  

  if ((dimDstart!=0)&&(dimDstart!=0)){
    dimDs.clear();
    for (int i=dimDstart; i<dimDend; i++){
      dimDs.push_back( i );
    }
  } else {}

  cout<<"loading argv"<<endl;
  pars.PrintVars();

  //create folders
  mkdir( rootFolder.c_str(), 0755 );
  mkdir( (rootFolder+"OutputT/").c_str(), 0755 );
  mkdir( (rootFolder+"data/").c_str(), 0755 );

  loadFolder = rootFolder + "OutputT/" + loadFolder;
  enviFolder = rootFolder + "OutputT/" + enviFolder;
  tensorFolder = rootFolder + "OutputT/" + saveFolder;
  dataFolder = rootFolder + "data/" + saveFolder;
  mkdir( tensorFolder.c_str(), 0755 );
  mkdir( dataFolder.c_str(), 0755 );
  char buffer[16];
  sprintf( buffer, "J%04i/", int( 1000*field) );
  mkdir( dataFolder.c_str(), 0755);
  mkdir( tensorFolder.c_str(), 0755);
}

