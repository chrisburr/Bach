

#define ATTR_SET ".<xmlattr>"

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <stdio.h>


#include "../../TbKernel/src/TbGeometrySvc.h"
#include "../../TbKernel/src/TbBaseClass.h"
#include "../../TbAlgorithms/src/TbDecoder.h"
#include "../../TbAlgorithms/src/TbClustering.h"
#include "../../TbAlgorithms/src/TbToyData.h"
#include "../../TbAlgorithms/src/TbPatternRecognition.h"
#include "../../TbAlgorithms/src/TbPlotTool.h"
#include "../../TbAlgorithms/src/TbAlignment.h"
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

void parse_xml(char* xmlfile) {

  ptree tree;
  read_xml(xmlfile, tree);
  const ptree & algorithms = tree.get_child("Algorithms");
  n_evts =  algorithms.get< int >("<xmlattr>.NoOfEvt") ;
  BOOST_FOREACH(const ptree::value_type & a1, algorithms){
    string al = a1.first;
    if (al == "<xmlattr>") continue;
    const ptree & algo = algorithms.get_child(al);
    map<string, vector< string> > Const_Map; 
    BOOST_FOREACH(const ptree::value_type &a2, algo){
      if(!std::strcmp(a2.first.c_str(), "constant")){
	
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


int main(int argc, char *argv[])
  
{
  if ( argc!=2){
    cout << "Usage: bin/bach <xml.-Configfile>" << endl;
    return 0;}
  //Read .xml-file
  char *xmlfile = argv[1];
  AlgVec  Algorithm_Container;
  parse_xml(xmlfile);
  
  // Initialize Algorithms
  for ( vector< pair<string, map< string, vector< string> > > >::iterator it1 = Algorithms.begin();
	it1 != Algorithms.end(); ++it1){

    if ( (*it1).first == "TbGeometrySvc" ){
        
      TbGeometrySvc *tbgeo = new TbGeometrySvc((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbgeo));
      
    }
    else if ( (*it1).first == "TbDecoder" ){
     
   
      TbDecoder *tbdec = new TbDecoder((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbdec));
      
    }
    else if ( (*it1).first == "TbClustering" ){
      
      TbClustering *tbdec = new TbClustering((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbdec));
      
    }
    else if ( (*it1).first == "TbToyData" ){
      
      TbToyData *tbtoy = new TbToyData((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbtoy));
      
    }
    else if ( (*it1).first == "TbPatternRecognition" ){
     
      TbPatternRecognition *tbpr = new TbPatternRecognition((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbpr));
      
    }
    else if ( (*it1).first == "TbPlotTool" ){
     
      TbPlotTool *tbpt = new TbPlotTool((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbpt));
      
    }
    else if ( (*it1).first == "TbAlignment" ){
     
      TbAlignment *tbagn = new TbAlignment((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbagn));
      
    }
    else {cout << "Algorithm " << (*it1).first << " not known!" << std::endl;
	exit(1);}
    (Algorithm_Container.back()).second->configuration();
    for ( map< string, vector< string> >::iterator it2 = (*it1).second.begin();
	  it2 != (*it1).second.end(); ++it2){
      string constant = (*it2).first;
      string type = (*it2).second[1];
      string value = (*it2).second[0];
      
      if (type == "D")       {
	  double con_value = atof(value.c_str());
	  (Algorithm_Container.back()).second->Const_D(constant, con_value);
	}
	else if (type == "I")  {
	  int    con_value = atoi(value.c_str());
	  (Algorithm_Container.back()).second->Const_I(constant, con_value);
	}
	else if (type == "S")  {
	  string con_value = value;
	  (Algorithm_Container.back()).second->Const_S(constant.c_str(), con_value);
	  }
      
      
      
    }
    cout << "Initialising --> " << (Algorithm_Container.back()).first << endl;
    (Algorithm_Container.back()).second->initialize(Algorithm_Container);
   
  }
  //Run algorithms over nevts events
  
  for (int i = 0; i < n_evts; ++i){
    cout << "Processing event nr " << i << std::endl;
    for (AlgVec::iterator iter = Algorithm_Container.begin();
	 iter != Algorithm_Container.end();++iter)
      {
	(*iter).second->execute(Algorithm_Container);
      }
    for (AlgVec::iterator iter_e = Algorithm_Container.begin();
	 iter_e != Algorithm_Container.end();++iter_e)
      {
	(*iter_e).second->end_event();
      }
  }
  //Finalize algorithms
  for (AlgVec::iterator iter_f = Algorithm_Container.begin();
       iter_f != Algorithm_Container.end();++iter_f)
    {
      (*iter_f).second->finalize();
    }
  Algorithm_Container.clear();
  cout << "You're done!" << endl;

  return 0;
}
