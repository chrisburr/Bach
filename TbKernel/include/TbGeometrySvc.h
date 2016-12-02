#ifndef __TBGEOMETRYSVC_H__
#define __TBGEOMETRYSVC_H__

#include <map>
#include <string>

#include "../TbKernel/ITbGeometrySvc.h"

#include "Math/Point3D.h"
#include "Math/Transform3D.h"

#include "TbModule.h"

/** @class TbGeometrySvc TbGeometrySvc.h
 *
 */

using namespace ROOT::Math;
class TbGeometrySvc : public ITbGeometrySvc {

public:
  TbGeometrySvc(const std::string& name);
  virtual ~TbGeometrySvc();
  virtual bool configuration();
  virtual bool initialize(AlgVec);
  virtual bool finalize();

  virtual XYZPoint localToGlobal(const XYZPoint& p,
                                        const std::string& id);
  virtual XYZPoint globalToLocal(const XYZPoint& p,
                                        const std::string& id);

  virtual void  writeConditionsXML(std::string);
  std::map<std::string,TbModule* > Modules;
  std::vector<TbModule*> Modules_sort;
private:



  std::string m_name;
  bool readConditions();
  bool readConditionsXML();



};



#endif
