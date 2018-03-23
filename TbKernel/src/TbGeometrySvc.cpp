#include <fstream>
#include <sstream>
#include <stdio.h>

#include "DD4hep/Factories.h"
#include "TbGeometrySvc.h"

TbGeometrySvc::TbGeometrySvc(const std::string &name) { m_name = name; }

TbGeometrySvc::~TbGeometrySvc() {}

bool TbGeometrySvc::configuration() {
  Const_D("PitchX", 0.055);
  Const_D("PitchY", 0.055);
  Const_I("NoOfPixelX", 256);
  Const_I("NoOfPixelY", 256);
  Const_D("Thick", 0.3);

  return true;
}

std::vector<DetElement> get_alignables(DetElement::Children children) {
  std::cout << " Getting alignables " << std::endl;
  std::vector<DetElement> aligned;
  for (auto it = children.begin(); it != children.end(); ++it) {
    std::cout << " |-> " << it->second.path() << std::endl;
    if (it->second.path() != "/world/Telescope") {
    // if (it->second.hasConditions()) {
      aligned.push_back(it->second);
    }

    auto grandchildren = it->second.children();
    std::cout << "    > has " << it->second.children().size() << " children" << std::endl;
    if (grandchildren.begin() != grandchildren.end()) {
      auto aligned2 = get_alignables(grandchildren);
      aligned.insert(aligned.end(), aligned2.begin(), aligned2.end());
    }
  }
  return aligned;
}

bool TbGeometrySvc::initialize(DetElement world, AlgVec algos) {
  Modules = get_alignables(world.children());
  return true;
}

bool TbGeometrySvc::finalize(dd4hep::cond::ConditionsSlice &slice) {
  std::cout << "DEBUG: finalize" << std::endl;
  return false;
}
