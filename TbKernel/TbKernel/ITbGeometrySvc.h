#ifndef TBKERNEL_ITBGEOMETRYSVC_H
#define TBKERNEL_ITBGEOMETRYSVC_H

#include "Math/Point3D.h"
#include "Math/Transform3D.h"

#include "../include/TbBaseClass.h"
/** @class ITbGeometrySvc ITbGeometrySvc.h
 *
 *
 */
using namespace ROOT::Math;

class ITbGeometrySvc : public TbBaseClass {

public:
  virtual XYZPoint localToGlobal(const XYZPoint &p, const std::string &det) = 0;
  virtual XYZPoint globalToLocal(const XYZPoint &p, const std::string &det) = 0;
};

#endif
