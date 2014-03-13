import os,sys
from UpdateMakefile import UM
from datetime import date

def generateHeader(name):
    headFile = "../../TbAlgorithms/src/Tb"+name+".h"
    print "generating "+headFile
    d = date.today()
    
    h = open(headFile,'w')
    s = """
#ifndef TB_%s_H 
#define TB_%s_H 1

/* Header for %s
 *
 * %s
 *
 */

#include "../../TbKernel/src/TbBaseClass.h"

class Tb%s : public TbBaseClass {

public: 

  /// Constructor
  Tb%s(const std::string& name);
  /// Destructor
  virtual ~Tb%s(){};

  bool configuration();
  bool initialize(AlgVec);    ///< Algorithm initialization
  bool execute(AlgVec );       ///< Algorithm execution
  bool end_event();
  bool finalize();      ///< Algorithm finalization

 private:

};
#endif
    """ % (name.upper(), name.upper(), name, d, name, name, name ) 

    h.write(s)
    h.close()

    print headFile+" generated!"

def generateCpp(name):
    cppFile = "../../TbAlgorithms/src/Tb"+name+".cpp"
    print "generating "+cppFile
    d = date.today()
    
    h = open(cppFile,'w')
    s = """
#include "Tb%s.h"


/** @file Tb%s.cpp
 *
 *  Implementation of class : Tb%s
 *
 */


//=============================================================================
/// Standard constructor
//=============================================================================
Tb%s::Tb%s(const std::string& name)  
  {

}


bool Tb%s::configuration(){
  return true;

}

//=============================================================================
/// Initialization
//=============================================================================
bool Tb%s::initialize(AlgVec  algos) {

  return true;

}

//=============================================================================
/// Main execution
//=============================================================================
bool Tb%s::execute(AlgVec algos) {

  return true;

}

//=============================================================================
// End of Event
//=============================================================================
bool Tb%s::end_event(){
 
  return true;
}

//=============================================================================
/// Finalize
//=============================================================================
bool Tb%s::finalize() {

  return true;
}

    """ % (name, name, name, name, name, name, name, name, name, name)

    h.write(s)
    h.close()

    print cppFile+" generated!"

def main():
    if len(sys.argv) != 2:
        print "usage: python generateClass.py <ClassName>"
        sys.exit()
    generateHeader(sys.argv[1])
    generateCpp(sys.argv[1])
    UM()

main()
