#ifndef __TBGEOMETRYSVC_H__
#define __TBGEOMETRYSVC_H__

#include <map>
#include <string>

#include "DD4hep/Factories.h"

#include "TbBaseClass.h"

/** @class TbGeometrySvc TbGeometrySvc.h
 *
 */

using namespace DD4hep::Geometry;

class TbGeometrySvc : public TbBaseClass {
public:
  TbGeometrySvc(const std::string &name);
  virtual ~TbGeometrySvc();
  virtual bool configuration();
  virtual bool initialize(DetElement, AlgVec);
  virtual bool finalize(DD4hep::Conditions::ConditionsSlice &);

  std::vector<DetElement> Modules;

private:
  std::string m_name;
};

#endif
