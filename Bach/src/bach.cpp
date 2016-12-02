/* bach.cpp
 *   Christoph Hombach
 *   main program to run BACH
 *
 *   usage: bin/bach    <xml-ConfigFile>
 *
 *   Program reads in algorithms and corresponding constants and executes those
 */

#define ATTR_SET ".<xmlattr>"

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <stdio.h>

//Add header file to Algorithm here

#include "TbGeometrySvc.h"
#include "TbBaseClass.h"
#include "TbDecoder.h"
#include "TbClustering.h"
#include "TbToyData.h"
#include "TbPatternRecognition.h"
#include "TbPlotTool.h"
#include "TbAlignment.h"

using namespace std;
using namespace ROOT::Math;
using namespace boost;
using namespace boost::property_tree;

const ptree& empty_ptree(){
    static ptree t;
    return t;
}

int n_evts;

vector< pair<string, map< string, vector< string> > > > Algorithms;

void Logo() {
  cout << std::setw(10) << " " << endl;
  cout << std::setw(10) << "    o                o       "<<endl;
  cout << std::setw(10) << "    o                o       "<<endl;
  cout << std::setw(10) << "    ooo     oo   ooo o oo    "<<endl;
  cout << std::setw(10) << "    o  o  o   o o    oo  o   "<<endl;
  cout << std::setw(10) << "    o  o  o   o o    o   o   "<<endl;
  cout << std::setw(10) << "    o  o  o   o o    o   o   "<<endl;
  cout << std::setw(10) << "    ooo    oo o  ooo o   o   "<<endl;
  cout << std::setw(10) << "-----------------------------"<<endl;
}

// Read in .xml-file and store elements in Algorithm_Container

void parse_xml(char* xmlfile) {

  ptree tree;
  read_xml(xmlfile, tree);
  const ptree & algorithms = tree.get_child("Algorithms");
  n_evts =  algorithms.get< int >("<xmlattr>.NoOfEvt");
  BOOST_FOREACH(const ptree::value_type & a1, algorithms) {
    string al = a1.first;
    if (al == "<xmlattr>" || al == "<xmlcomment>") continue;

    const ptree & algo = algorithms.get_child(al);
    map<string, vector< string> > Const_Map;
    BOOST_FOREACH(const ptree::value_type &a2, algo) {
      if(!std::strcmp(a2.first.c_str(), "constant")) {

        string name  = a2.second.get< std::string >("<xmlattr>.name");

        string value = a2.second.get< std::string >("<xmlattr>.value");

        string type  = a2.second.get< std::string >("<xmlattr>.type");
        vector< string> Values;
        Values.push_back(value);
        Values.push_back(type);

        Const_Map[name]=Values;
      }
    }
    Algorithms.push_back(make_pair(al,Const_Map));
  }
}

//Main execution

int main(int argc, char *argv[]) {
  Logo();
  if ( argc != 2 ){
    cout << "Usage: bin/bach <xml.-Configfile>" << endl;
    return 0;
  }

  //Read .xml-file
  char *xmlfile = argv[1];
  AlgVec Algorithm_Container;
  parse_xml(xmlfile);

  // Initialize Algorithms
  for ( vector< pair<string, map< string, vector< string> > > >::iterator it1 = Algorithms.begin(); it1 != Algorithms.end(); ++it1) {

    //Add your algorithm here
    if ( (*it1).first == "TbGeometrySvc" ) {
      TbGeometrySvc *tbgeo = new TbGeometrySvc((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbgeo));
    } else if ( (*it1).first == "TbDecoder" ) {
      TbDecoder *tbdec = new TbDecoder((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbdec));
    } else if ( (*it1).first == "TbClustering" ) {
      TbClustering *tbdec = new TbClustering((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbdec));
    } else if ( (*it1).first == "TbToyData" ) {
      TbToyData *tbtoy = new TbToyData((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbtoy));
    } else if ( (*it1).first == "TbPatternRecognition" ) {
      TbPatternRecognition *tbpr = new TbPatternRecognition((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbpr));
    } else if ( (*it1).first == "TbPlotTool" ) {
      TbPlotTool *tbpt = new TbPlotTool((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbpt));
    } else if ( (*it1).first == "TbAlignment" ) {
      TbAlignment *tbagn = new TbAlignment((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbagn));
    } else {
      cout << "Algorithm " << (*it1).first << " not known!" << std::endl;
      exit(1);
    }

    (Algorithm_Container.back()).second->configuration();

    //Read in constants from .xml file (Overrides default-values)

    for ( map< string, vector< string> >::iterator it2 = (*it1).second.begin(); it2 != (*it1).second.end(); ++it2) {
      string constant = (*it2).first;
      string type = (*it2).second[1];
      string value = (*it2).second[0];

      if (type == "D") {
        double con_value = atof(value.c_str());
        (Algorithm_Container.back()).second->Const_D(constant, con_value);
      }
      else if (type == "I") {
        int con_value = atoi(value.c_str());
        (Algorithm_Container.back()).second->Const_I(constant, con_value);
      }
      else if (type == "S") {
        string con_value = value;
        (Algorithm_Container.back()).second->Const_S(constant.c_str(), con_value);
      }
    }

    //Initialise algorithms
    cout << "Initialising --> " << (Algorithm_Container.back()).first << endl;
    (Algorithm_Container.back()).second->initialize(Algorithm_Container);

  }

  //Run algorithms over nevts events
  for (int i = 0; i < n_evts; ++i) {
    cout << "Processing event nr " << i << std::endl;

    for (AlgVec::iterator iter = Algorithm_Container.begin(); iter != Algorithm_Container.end();++iter) {
      cout << "\tExecute " << (*iter).first << endl;
      (*iter).second->execute(Algorithm_Container);
    }

    for (AlgVec::iterator iter_e = Algorithm_Container.begin(); iter_e != Algorithm_Container.end();++iter_e) {
      (*iter_e).second->end_event();
    }
  }

  //Finalize algorithms
  for (AlgVec::iterator iter_f = Algorithm_Container.begin(); iter_f != Algorithm_Container.end();++iter_f) {
    (*iter_f).second->finalize();
  }

  Algorithm_Container.clear();
  cout << "You're done!" << endl;

  return 0;
}
