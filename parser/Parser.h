#ifndef PARSER_HPP_INCLUDED
#define PARSER_HPP_INCLUDED

#include <map>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>

//! A class to parse input parmeters.
/*!
 *  One can specify the input parameters in inputFile and Parser will handle the parsing automatically.
 *  For example,
 *  \code{.cpp}
 *  Parser pars;
 *  int dimD;
 *  pars.bind( "virtual bond dimension", dimD );
 *  pars.Parse( "inputFile" );
 *  \endcode
 *  And in inputFile:
 *  \code{.py}
 *  virtual bond dimension = 3
 *  \endcode
 *  will do the job.
 */
class Parser{
  private:
    enum {
      BOOL_type,
      INT_type,
      FLT_type,
      UINT_type,
      ULNG_type,
      LNG_type,
      DBL_type,
      STR_type,
      LLNG_type,
      VECT_INT_type,
      VECT_DBL_type
    };
    std::map<unsigned int,std::string> DictType;
  public :
    Parser(){	Init(); }
    void Init(){
      DictType[BOOL_type] = "bool";
      DictType[INT_type] = "int";
      DictType[FLT_type] = "float";
      DictType[UINT_type] = "unsigned int";
      DictType[ULNG_type] = "unsigned long";
      DictType[LNG_type] = "long";
      DictType[DBL_type] = "double";
      DictType[STR_type] = "string";
      DictType[LLNG_type] = "long long";
      DictType[VECT_INT_type] = "vector int";
      DictType[VECT_DBL_type] = "vector double";
    }
    struct rElem{
      void* pVal;
      std::string key;
      int pType;
      bool isRead;
    };
    std::map<std::string,unsigned int > Matcher;
    std::vector< rElem > Keys;
    //! Bind a key with a bool variable.
    /*!
     *  \param[in] Key being binded.
     *  \param[in/out] Variable being bind.
     */
    void Bind(const std::string &key, bool &Var);
    //! \overload
    void Bind(const std::string &key, int &Var);
    //! \overload
    void Bind(const std::string &key, float &Var);
    //! \overload
    void Bind(const std::string &key, unsigned int &Var);
    //! \overload
    void Bind(const std::string &key, unsigned long &Var);
    //! \overload
    void Bind(const std::string &key, long &Var);
    //! \overload
    void Bind(const std::string &key, double &Var);
    //! \overload
    void Bind(const std::string &key, std::string &Var);
    //! \overload
    void Bind(const std::string &key, long long &Var);
    //! \overload
    void Bind(const std::string &key, std::vector<int> &Var);
    //! \overload
    void Bind(const std::string &key, std::vector<double> &Var);
    int TryRead(const std::string &rawkey,const std::string &readin_VarString);
    //! Parse parameters from given file.
    /*!
     *  \param[in] rc_name Input files for parameters.
     *
     */
    void Parse(const std::string &rc_fname);
    //! Print the stirngs being binded.
    void PrintBinds();
    //! Print the stirngs being and its binding variables.
    /*!
     *  Example of input file:
     *  \code{.py}
     *  bindedString = 3
     *  \endcode
     *  Where bindedString is binded with an int variable in the program and is set to 3 after parsing.
     */
    void PrintVars();
    void Remove(const std::string &keys);
    void Check_All();
};

#endif // SSE_HPP_INCLUDED
