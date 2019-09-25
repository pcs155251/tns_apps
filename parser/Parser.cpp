#include "Parser.h"
using namespace std;

string RemoveBlnk(string istr){

    if(istr.empty()){ 
      istr.clear(); 
      return istr;
    } else {}

    while(istr.at(0)==' '){
      istr.erase(istr.begin());
      if(istr.empty()){ 
        istr.clear(); 
        return istr;
      } else {}
    }

    while(istr.at(istr.size()-1)==' '){
      istr.erase(istr.end()-1);
      if(istr.empty()){ 
        istr.clear(); 
        return istr;
      } else {}
    }

    /*
    if(istr.find(" ") != std::string::npos){ 
      istr.clear(); 
      return istr;
    } else {}
    */

    return istr;

}

int Parser::TryRead(const string &rawkey,const std::string &readin_VarString){

  map<string,unsigned int>::iterator it = Matcher.find(RemoveBlnk(rawkey));
  if(it == Matcher.end()) {
    return 0;
  } else {}
  unsigned int idx = it->second;

  string rhs = RemoveBlnk(readin_VarString);
  if(rhs.empty()){
    return 0;
  } else {}
    
  //checkType:
  if(Keys[idx].pType == BOOL_type){
    bool* tmp = (bool*)Keys[idx].pVal;
    if (rhs=="true"){
      *tmp = true;
    }
    else if (rhs=="false"){
      *tmp = false;
    }
    else {
      cout<< "invalid input of BOOL_type: "<<rhs<<endl;
    }
  }
  else if(Keys[idx].pType == INT_type){
    int* tmp = (int*)Keys[idx].pVal;
    *tmp = atol(rhs.c_str());
  }
  else if (Keys[idx].pType == UINT_type){
    unsigned int* tmp = (unsigned int*)Keys[idx].pVal;
    *tmp = atol(rhs.c_str());
  }
  else if(Keys[idx].pType == LNG_type){
    long* tmp = (long*)Keys[idx].pVal;
    *tmp = atol(rhs.c_str());
  } 
  else if(Keys[idx].pType == ULNG_type){
    unsigned long* tmp = (unsigned long*)Keys[idx].pVal;
    *tmp = atol(rhs.c_str());
  }
  else if(Keys[idx].pType == FLT_type){
    float* tmp = (float*)Keys[idx].pVal;
    *tmp = atof(rhs.c_str());
  }
  else if(Keys[idx].pType == DBL_type){
    double* tmp = (double*)Keys[idx].pVal;
    *tmp = atof(rhs.c_str());
  }
  else if(Keys[idx].pType == STR_type){
    string* tmp = (string*)Keys[idx].pVal;
    *tmp = rhs;
  }
  else if(Keys[idx].pType == LLNG_type){
    long long* tmp = (long long*)Keys[idx].pVal;
    *tmp = atoll(rhs.c_str());
  }
  else if(Keys[idx].pType == VECT_INT_type){
    vector<int>* tmp = (vector<int> *)Keys[idx].pVal;
    istringstream iss(rhs);
    string token;
    while (getline(iss, token, ' ')){
      tmp->push_back(atol(token.c_str()));
    }
  }
  else if(Keys[idx].pType == VECT_DBL_type){
    vector<double>* tmp = (vector<double> *)Keys[idx].pVal;
    istringstream iss(rhs);
    string token;
    while (getline(iss, token, ' ')){
      tmp->push_back(atof(token.c_str()));
    }
  }
  else{
    cout << "[ERROR][Parser::TryRead] invalid Type." << endl;
    exit(1);
  }
  Keys[idx].isRead = 1;
  return 1;
}

void Parser::Parse(const string &rc_fname){
  FILE * pFile;
  char line[256];
  char *tmp;
  string lhs,rhs;
  pFile = fopen(rc_fname.c_str() , "rb");

  //check exist
  if(pFile != NULL){
    //loop line
    while(1){
      //chk EOF & read line
      if(fgets(line,sizeof(line),pFile)==NULL) {
        break;
      } else {}
      if(line[0]=='#') {
        continue;
      } else {}
      //get lhs(key)
      tmp = strtok(line,"=");
      lhs = string(tmp);
      if(tmp != NULL){
        //get rhs(val)
        tmp = strtok(NULL,"=");
        if(tmp ==NULL ) {
          continue;
        } else {}
        tmp = strtok(tmp,"\n");
        if(tmp == NULL) {
          continue;
        } else {}
        rhs = string(tmp);
        TryRead(lhs,rhs);
      } else {}
    }
    fclose(pFile);
  }
  else{
    cout << "[ERROR][Parser::Parse] invalid .rc file" << endl;
    exit(1);
  }
}

void Parser::PrintBinds(){
  cout<<endl<<"number of binded element: "<<Keys.size()<<endl;
  cout<<left<<setw(16)<<"parameter name"<<setw(8)<<"isRead"<<setw(16)<<"type"<<endl;;
  for(int i=0;i<Keys.size();i++){
  	cout<<setw(16)<<Keys[i].key<<setw(8)<<Keys[i].isRead<<setw(16)<<DictType[Keys[i].pType]<<endl;
  }
  cout<<endl;
}

void Parser::PrintVars(){
  cout<<endl<<"parameters:"<<endl<<"number of binded element: "<< Keys.size() << endl;
  string stat;

  for(int idx=0;idx<Keys.size();idx++){
    if(Keys[idx].isRead== 0){
    cout << setw(16) << left << Keys[idx].key << " [** No Parsed **] " << endl;
    continue;
    }
    //checkType:
    if(Keys[idx].pType ==  BOOL_type){
      bool* tmp =  (bool*)Keys[idx].pVal;
      cout << setw(16) << left << Keys[idx].key << "= " <<  *tmp << endl;
    }
    else if(Keys[idx].pType ==  INT_type){
      int* tmp =  (int*)Keys[idx].pVal;
      cout << setw(16) << left << Keys[idx].key << "= " <<  *tmp << endl;
    }
    else if(Keys[idx].pType ==  UINT_type){
      unsigned int* tmp =  (unsigned int*)Keys[idx].pVal;
      cout << setw(16) << left << Keys[idx].key << "= " <<  *tmp << endl;
    }
    else if(Keys[idx].pType ==  LNG_type){
      long* tmp =  (long*)Keys[idx].pVal;
      cout << setw(16) << left << Keys[idx].key << "= " <<  *tmp << endl;
    }
    else if(Keys[idx].pType ==  ULNG_type){
      unsigned long* tmp =  (unsigned long*)Keys[idx].pVal;
      cout << setw(16) << left << Keys[idx].key << "= " <<  *tmp << endl;
    }
    else if(Keys[idx].pType ==  FLT_type){
      float* tmp =  (float*)Keys[idx].pVal;
      cout << setw(16) << left << Keys[idx].key << "= " <<  *tmp << endl;
    }
    else if(Keys[idx].pType ==   DBL_type){
      double* tmp =  (double*)Keys[idx].pVal;
      cout << setw(16) << left << Keys[idx].key << "= " <<  *tmp << endl;
    }
    else if(Keys[idx].pType ==  STR_type){
      string* tmp =  (string*)Keys[idx].pVal;
      cout << setw(16) << left << Keys[idx].key << "= " <<  *tmp << endl;
    }
    else if(Keys[idx].pType ==  LLNG_type){
      long long* tmp =  (long long*)Keys[idx].pVal;
      cout << setw(16) << left << Keys[idx].key << "= " <<  *tmp << endl;
    }
    else if(Keys[idx].pType ==  VECT_INT_type){
      vector<int>* tmp =  (vector<int> *)Keys[idx].pVal;
      cout << setw(16) << left << Keys[idx].key << "= ";
      for (vector<int>::iterator it=tmp->begin(); it!=tmp->end(); it++){
        cout<<setw(6)<<*it;
      }
      cout<<endl;
    }
    else if(Keys[idx].pType ==  VECT_DBL_type){
      vector<double>* tmp =  (vector<double> *)Keys[idx].pVal;
      cout << setw(16) << left << Keys[idx].key << "= ";
      for (vector<double>::iterator it=tmp->begin(); it!=tmp->end(); it++){
        cout<<setw(6)<<*it;
      }
      cout<<endl;
    }
    else{
      cout << "[ERROR][Parser::PrintVars] invalid Type." << endl<<endl;;
      exit(1);
    }
  }
  cout<<endl;
}

void Parser::Remove(const std::string &key){
  map<std::string,unsigned int>::iterator it;
  it = Matcher.find(key);
  if(it!=Matcher.end()){
    unsigned int idx = it->second;
    Keys.erase(Keys.begin()+idx);
    Matcher.erase(it);
  } else {}
}

void Parser::Check_All(){
  int ff =0;
  for(int i=0;i<Keys.size();i++){
    if(Keys[i].isRead==0){
        cout << "[ERROR][Parser::Check_All] < "  << Keys[i].key << " > No read." << endl;
        ff = 1;
    } else {}
  }
  if(ff) {
    exit(1);
  } else {}
}

void Parser::Bind(const std::string &key, bool &Var){
  Matcher[key] = Keys.size();
  rElem tmp = {NULL,"",BOOL_type,0 } ;
  tmp.pVal = (void*)&Var;
  tmp.key = key;
  Keys.push_back(tmp);
 
}

void Parser::Bind(const std::string &key, int &Var){
  Matcher[key] = Keys.size();
  rElem tmp = {NULL,"",INT_type,0 } ;
  tmp.pVal = (void*)&Var;
  tmp.key = key;
  Keys.push_back(tmp);
 
}

void Parser::Bind(const std::string &key, float &Var){
  Matcher[key] = Keys.size();
  rElem tmp = {NULL,"",FLT_type,0 } ;
  tmp.pVal = (void*)&Var;
  tmp.key = key;
  Keys.push_back(tmp);
}

void Parser::Bind(const std::string &key, unsigned int &Var){
  Matcher[key] = Keys.size();
  rElem tmp = {NULL,"",UINT_type,0 } ;
  tmp.pVal = (void*)&Var;
  tmp.key = key;
  Keys.push_back(tmp);
}

void Parser::Bind(const std::string &key, unsigned long &Var){
  Matcher[key] = Keys.size();
  rElem tmp = {NULL,"",ULNG_type,0 } ;
  tmp.pVal = (void*)&Var;
  tmp.key = key;
  Keys.push_back(tmp);
}

void Parser::Bind(const std::string &key, long &Var){
  Matcher[key] = Keys.size();
  rElem tmp = {NULL,"",LNG_type,0 } ;
  tmp.pVal = (void*)&Var;
  tmp.key = key;
  Keys.push_back(tmp);
}

void Parser::Bind(const std::string &key, double &Var){
  Matcher[key] = Keys.size();
  rElem tmp = {NULL,"",DBL_type,0 } ;
  tmp.pVal = (void*)&Var;
  tmp.key = key;
  Keys.push_back(tmp);
}

void Parser::Bind(const std::string &key, std::string &Var){
  Matcher[key] = Keys.size();
  rElem tmp = {NULL,"",STR_type,0 } ;
  tmp.pVal = (void*)&Var;
  tmp.key = key;
  Keys.push_back(tmp);
}

void Parser::Bind(const std::string &key, long long &Var){
  Matcher[key] = Keys.size();
  rElem tmp = {NULL,"",LLNG_type,0 } ;
  tmp.pVal = (void*)&Var;
  tmp.key = key;
  Keys.push_back(tmp);
}

void Parser::Bind(const std::string &key, vector<int> &Var){
  Matcher[key] = Keys.size();
  rElem tmp = {NULL,"",VECT_INT_type,0 } ;
  tmp.pVal = (void*)&Var;
  tmp.key = key;
  Keys.push_back(tmp);
}

void Parser::Bind(const std::string &key, vector<double> &Var){
  Matcher[key] = Keys.size();
  rElem tmp = {NULL,"",VECT_DBL_type,0 } ;
  tmp.pVal = (void*)&Var;
  tmp.key = key;
  Keys.push_back(tmp);
}
