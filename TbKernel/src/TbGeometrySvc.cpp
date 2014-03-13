#define ATTR_SET ".<xmlattr>"

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <stdio.h>
#include <fstream>
#include <sstream>

#include "Math/Translation3D.h"
#include "Math/RotationZYX.h"
#include "Math/Transform3D.h"
#include "TMath.h"

#include "TbGeometrySvc.h"

using namespace boost;
using namespace boost::property_tree;

TbGeometrySvc::TbGeometrySvc(const std::string& name) {m_name = name;}

TbGeometrySvc::~TbGeometrySvc() {}

bool TbGeometrySvc::configuration(){  
  
  
  Const_D("PitchX",0.055);
  Const_D("PitchY",0.055);
  Const_I("NoOfPixelX",256);
  Const_I("NoOfPixelY",256);
  Const_D("Thick"     ,0.3);
  Const_S("GeometryFile","");
  
  
  return true;
}

bool TbGeometrySvc::initialize(AlgVec algos) {

  if (!readConditionsXML()) {
    
    std::cout << "ERROR: Cannot import alignment conditions" << std::endl;
    return false;
  }
  //  PrintAlgorithm();
  return true;

}

bool TbGeometrySvc::finalize() {
  

  std::cout << "DEBUG: finalize" << std::endl;
  return false;
  
}

XYZPoint TbGeometrySvc::localToGlobal(const XYZPoint& p,
				      const std::string& id) {
  XYZPoint pGlobal;
  if (Modules.count(id) > 0) {
    
    pGlobal = Modules[id]->Transform() * p;
    return pGlobal;
  } else { 
    std::cout << "ERROR: Transform for " << id << " not found" << std::endl;
  }

  return pGlobal;
 
}

XYZPoint TbGeometrySvc::globalToLocal(const XYZPoint& p,
				      const std::string& id) {

  XYZPoint pLocal;
  if (Modules.count(id) > 0) {
    XYZPoint pLocal = Modules[id]->Transform().Inverse() * p;
    return pLocal;
  } else {
    std::cout << "ERROR: Transform for " << id << " not found" << std::endl;
  }
  return pLocal;
 
}

bool TbGeometrySvc::readConditions() {
  
  std::ifstream infile;
  infile.open(Const_S("GeometryFile").c_str(), std::ios::in);
  if (!infile) {
    std::cerr << "Can't open input file " << Const_S("GeometryFile") << std::endl;
    exit(1);
  }
  else std::cout << "Open input file " << Const_S("GeometryFile") << std::endl;
  std::string line;
  while (!infile.eof()) {
    std::getline(infile, line);
    size_t p = line.find_first_not_of(' ');
    if (std::string::npos == p) {
      // Empty line.
      continue;
    }
    // Trim whitespace from beginning of line.
    line = line.substr(p);
    // Skip comments.
    if (line[0] == '#') continue;
    unsigned int column = 0;
    std::string chip;
    double dx = 0., dy = 0., dz = 0.;
    double rx = 0., ry = 0., rz = 0.;
    std::string token;
    std::stringstream ss(line);
    while (std::getline(ss, token, ' ')) {
      if (token.empty()) continue;
      switch (column) {
        case 0:
          chip = token;
	  Modules[chip] = new TbModule();
          break;
        case 1:
          dx = atof(token.c_str());
	  break;
        case 2:
          dy = atof(token.c_str());
	  break;
        case 3:
          dz = atof(token.c_str());
	  break;
        case 4:
          rx = atof(token.c_str());
	  break;
        case 5:
          ry = atof(token.c_str());
	  break;
        case 6:
          rz = atof(token.c_str());
	  break;
      }
      ++column;
    }
    
    Modules[chip]->SetAlignment(chip, dx, dy, dz, rx, ry, rz, 0., 0., 0., 0., 0., 0.);    
    
  }

  return true;

}

bool TbGeometrySvc::readConditionsXML() {
  
  const char* xmlfile = Const_S("GeometryFile").c_str();
  ptree tree;
  read_xml(xmlfile, tree);
  const ptree & geometry = tree.get_child("Geometry");
  BOOST_FOREACH(const ptree::value_type & a1, geometry){
   std::string chip = a1.second.get< std::string >("<xmlattr>.name");

   float x = a1.second.get< float >("<xmlattr>.X"); 
   float y = a1.second.get< float >("<xmlattr>.Y");
   float z = a1.second.get< float >("<xmlattr>.Z");
   float rx = a1.second.get< float >("<xmlattr>.RX");
   float ry = a1.second.get< float >("<xmlattr>.RY");
   float rz = a1.second.get< float >("<xmlattr>.RZ");
   float dx = a1.second.get< float >("<xmlattr>.dX"); 
   float dy = a1.second.get< float >("<xmlattr>.dY");
   float dz = a1.second.get< float >("<xmlattr>.dZ");
   float drx = a1.second.get< float >("<xmlattr>.dRX");
   float dry = a1.second.get< float >("<xmlattr>.dRY");
   float drz = a1.second.get< float >("<xmlattr>.dRZ");

   Modules[chip] = new TbModule(); 
   Modules[chip]->SetAlignment(chip, x, y, z, rx, ry, rz, dx, dy, dz, drx, dry, drz);
   }  
return true;
} 

void TbGeometrySvc::writeConditionsXML() {
  const char* xmlfile = Const_S("GeometryFile").c_str();
  std::cout << "Writing geometry file " << xmlfile << std::endl; 
  std::ofstream myfile;
  myfile.open(xmlfile);
  myfile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n\n<Geometry>\n\n";
  std::map<std::string,TbModule* >::iterator itm = Modules.begin();
  for (;itm != Modules.end(); ++itm) {
    myfile << "\t<Module name=\'"<<(*itm).first
	   <<"\' X=\'"<<(*itm).second->X()<<"\' Y=\'"<< (*itm).second->Y() 
	   <<"\' Z=\'"<<(*itm).second->Z()<<"\' RX=\'"<< (*itm).second->RotX() 
	   <<"\' RY=\'"<<(*itm).second->RotY()<<"\' RZ=\'"<< (*itm).second->RotZ() 
	   <<"\' dX=\'"<<(*itm).second->dX()<<"\' dY=\'"<< (*itm).second->dY() 
	   <<"\' dZ=\'"<<(*itm).second->dZ()<<"\' dRX=\'"<< (*itm).second->dRotX() 
	   <<"\' dRY=\'"<<(*itm).second->dRotY()<<"\' dRZ=\'"<< (*itm).second->dRotZ() 
	   <<"\' />" << std::endl;
  }

  myfile << "</Geometry>";
}
