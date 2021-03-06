#define ATTR_SET ".<xmlattr>"

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <iostream>
#include <map>
#include <stdio.h>
#include <string>
#include <vector>

#include "AlignmentExampleObjects.h"
#include "DD4hep/Factories.h"
#include "DD4hep/Objects.h"
#include "DDAlign/GlobalAlignmentWriter.h"
#include "DDAlign/GlobalAlignmentCache.h"

#include "TbAlignment.h"
#include "TbBaseClass.h"
#include "TbClustering.h"
#include "TbDecoder.h"
#include "TbGeometrySvc.h"
#include "TbPatternRecognition.h"
#include "TbPlotTool.h"
#include "TbToyData.h"

using namespace std;
using namespace ROOT::Math;
using namespace boost;
using namespace boost::property_tree;

const ptree &empty_ptree() {
  static ptree t;
  return t;
}

int n_evts;

vector<pair<string, map<string, vector<string>>>> Algorithms;

static void Logo() {
  cout << std::setw(10) << " " << endl;
  cout << std::setw(10) << "    o                o       " << endl;
  cout << std::setw(10) << "    o                o       " << endl;
  cout << std::setw(10) << "    ooo     oo   ooo o oo    " << endl;
  cout << std::setw(10) << "    o  o  o   o o    oo  o   " << endl;
  cout << std::setw(10) << "    o  o  o   o o    o   o   " << endl;
  cout << std::setw(10) << "    o  o  o   o o    o   o   " << endl;
  cout << std::setw(10) << "    ooo    oo o  ooo o   o   " << endl;
  cout << std::setw(10) << "-----------------------------" << endl;
}

// Read in .xml-file and store elements in Algorithm_Container

static void parse_xml(string xmlfile) {
  ptree tree;
  read_xml(xmlfile, tree);
  const ptree &algorithms = tree.get_child("Algorithms");
  n_evts = algorithms.get<int>("<xmlattr>.NoOfEvt");
  BOOST_FOREACH (const ptree::value_type &a1, algorithms) {
    string al = a1.first;
    if (al == "<xmlattr>" || al == "<xmlcomment>")
      continue;

    const ptree &algo = algorithms.get_child(al);
    map<string, vector<string>> Const_Map;
    BOOST_FOREACH (const ptree::value_type &a2, algo) {
      if (!std::strcmp(a2.first.c_str(), "constant")) {
        string name = a2.second.get<std::string>("<xmlattr>.name");

        string value = a2.second.get<std::string>("<xmlattr>.value");

        string type = a2.second.get<std::string>("<xmlattr>.type");
        vector<string> Values;
        Values.push_back(value);
        Values.push_back(type);

        Const_Map[name] = Values;
      }
    }
    Algorithms.push_back(make_pair(al, Const_Map));
  }
}

dd4hep::DetElement element_from_path(dd4hep::Detector &lcdd, string path) {

  if (path.substr(0, 7) != "/world/") {
    cout << "Invalid path given: " << path << endl;
    return dd4hep::DetElement();
  }
  // Remove "/world/" from the path
  path = path.substr(6);
  // Iterate through the path to find the requested element
  dd4hep::DetElement elm = lcdd.world();
  for (size_t start = 0, pos = 0; start++ < path.size(); start = pos) {
    pos = path.find('/', start);
    elm = elm.child(path.substr(start, pos-1));
  }
  cout << elm.path() << endl;
  return elm;
}

// Main execution
static int run_bach(dd4hep::Detector &lcdd, int argc, char **argv) {
  Logo();

  // Parse the arguments
  string input_fn, delta_fn, config_fn;
  bool arg_error = false;
  for (int i = 0; i < argc && argv[i]; ++i) {
    if (0 == ::strncmp("-input", argv[i], 4))
      input_fn = argv[++i];
    else if (0 == ::strncmp("-deltas", argv[i], 5))
      delta_fn = argv[++i];
    else if (0 == ::strncmp("-config", argv[i], 5))
      config_fn = argv[++i];
    else
      arg_error = true;
  }
  if (arg_error || input_fn.empty() || delta_fn.empty() || config_fn.empty()) {
    /// Help printout describing the basic command line interface
    cout << "Usage: -plugin <name> -arg [-arg]                            \n"
            "\tname:   factory name     Bach_main                         \n"
            "\t-input   <string>        Geometry file                     \n"
            "\t-deltas  <string>        Alignment deltas (Conditions)     \n"
            "\t-config  <string>        Configuration xml file            \n";
    ::exit(EINVAL);
  }

  // First we load the geometry
  lcdd.fromXML(input_fn);
  dd4hep::AlignmentExamples::installManager(lcdd);

  dd4hep::cond::ConditionsManager condMgr = dd4hep::cond::ConditionsManager::from(lcdd);
  const void *delta_args[] = {delta_fn.c_str(), 0}; // Better zero-terminate

  lcdd.apply("DD4hep_ConditionsXMLRepositoryParser", 1, (char **)delta_args);
  // Now the deltas are stored in the conditions manager in the proper IOV pools
  const dd4hep::IOVType *iov_typ = condMgr.iovType("run");
  if (0 == iov_typ) {
    dd4hep::except("ConditionsPrepare", "++ Unknown IOV type supplied.");
  }
  dd4hep::IOV req_iov(iov_typ, 1500); // IOV goes from run 1000 ... 2000
  std::shared_ptr<dd4hep::cond::ConditionsContent> content(new dd4hep::cond::ConditionsContent());
  std::unique_ptr<dd4hep::cond::ConditionsSlice> slice =
      std::make_unique<dd4hep::cond::ConditionsSlice>(condMgr, content);
  auto cres = condMgr.prepare(req_iov, *slice);

  // ++++++++++++++++++++++++ Compute the transformation matrices
  dd4hep::align::AlignmentsCalculator calc;
  dd4hep::ConditionsHashMap alignments;
  std::map<dd4hep::DetElement, dd4hep::Delta> deltas;
  dd4hep::DetectorScanner(dd4hep::align::deltaCollector(*slice, deltas), lcdd.world());
  slice->pool->flags |= dd4hep::cond::UserPool::PRINT_INSERT;
  dd4hep::align::AlignmentsCalculator::Result  ares = calc.compute(deltas, *slice);

  // Read the configuration xml
  AlgVec Algorithm_Container;
  parse_xml(config_fn);

  // Initialize Algorithms
  for (vector<pair<string, map<string, vector<string>>>>::iterator it1 =
           Algorithms.begin();
       it1 != Algorithms.end(); ++it1) {
    // Add your algorithm here
    if ((*it1).first == "TbGeometrySvc") {
      TbGeometrySvc *tbgeo = new TbGeometrySvc((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbgeo));
    } else if ((*it1).first == "TbDecoder") {
      TbDecoder *tbdec = new TbDecoder((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbdec));
    } else if ((*it1).first == "TbClustering") {
      TbClustering *tbdec = new TbClustering((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbdec));
    } else if ((*it1).first == "TbToyData") {
      TbToyData *tbtoy = new TbToyData((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbtoy));
    } else if ((*it1).first == "TbPatternRecognition") {
      TbPatternRecognition *tbpr = new TbPatternRecognition((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbpr));
    } else if ((*it1).first == "TbPlotTool") {
      TbPlotTool *tbpt = new TbPlotTool((*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbpt));
    } else if ((*it1).first == "TbAlignment") {
      TbAlignment *tbagn = new TbAlignment(lcdd, (*it1).first);
      Algorithm_Container.push_back(make_pair((*it1).first, tbagn));
    } else {
      cout << "Algorithm " << (*it1).first << " not known!" << std::endl;
      exit(1);
    }

    (Algorithm_Container.back()).second->configuration();

    // Read in constants from .xml file (Overrides default-values)

    for (map<string, vector<string>>::iterator it2 = (*it1).second.begin();
         it2 != (*it1).second.end(); ++it2) {
      string constant = (*it2).first;
      string type = (*it2).second[1];
      string value = (*it2).second[0];

      if (type == "D") {
        double con_value = atof(value.c_str());
        (Algorithm_Container.back()).second->Const_D(constant, con_value);
      } else if (type == "I") {
        int con_value = atoi(value.c_str());
        (Algorithm_Container.back()).second->Const_I(constant, con_value);
      } else if (type == "S") {
        string con_value = value;
        (Algorithm_Container.back())
            .second->Const_S(constant.c_str(), con_value);
      }
    }

    // Initialise algorithms
    cout << "Initialising --> " << (Algorithm_Container.back()).first << endl;
    if (it1->first == "TbGeometrySvc") {
      auto tbgeo = dynamic_cast<TbGeometrySvc *>(Algorithm_Container.back().second);
      tbgeo->initialize(lcdd.world(), Algorithm_Container);
    } else {
      (Algorithm_Container.back()).second->initialize(Algorithm_Container);
    }
  }

  // Run algorithms over nevts events
  for (int i = 0; i < n_evts; ++i) {
    cout << "Processing event nr " << i << std::endl;

    for (AlgVec::iterator iter = Algorithm_Container.begin();
         iter != Algorithm_Container.end(); ++iter) {
      cout << "\tExecute " << (*iter).first << endl;
      (*iter).second->execute(*slice, Algorithm_Container);
    }

    for (AlgVec::iterator iter_e = Algorithm_Container.begin();
         iter_e != Algorithm_Container.end(); ++iter_e) {
      (*iter_e).second->end_event();
    }
  }

  // Finalize algorithms
  for (AlgVec::iterator iter_f = Algorithm_Container.begin();
       iter_f != Algorithm_Container.end(); ++iter_f) {
    (*iter_f).second->finalize(*slice);
  }

  Algorithm_Container.clear();
  cout << "You're done!" << endl;

  return 0;
}

// first argument is the type from the xml file
DECLARE_APPLY(Bach_main, run_bach)
