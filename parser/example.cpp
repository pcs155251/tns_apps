#include "Parser.hpp"
using namespace std;

int main(int argc, char* argv[]){
  Parser pars;
  //The values to be parse:
  int MPSnumber, maxIter;
  string modelName;
  //string fieldSubFolder;
  double field, criterion, evoStop;
  vector<int> chis;
  vector<double> taus;
  //vector<int> dimDS;

  //Bind the value with keyword:
  pars.Bind("Model",modelName);
  pars.Bind("MPSnumber",MPSnumber);
  pars.Bind("maxIter",maxIter);
  pars.Bind("field",field);
  pars.Bind("criterion",evoStop);
  pars.Bind("chi",chis);
  pars.Bind("tau",taus);

  pars.PrintBinds();
  pars.PrintVars();
  //Start parse:
  pars.Parse(string(argv[1]));
  //Layout the vals:
  pars.PrintBinds();
  pars.PrintVars();
  return 0 ;
}
