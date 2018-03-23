
#ifndef TB_PLOTTOOL_H
#define TB_PLOTTOOL_H 1

#include "TbAlignment.h"
#include "TbGeometrySvc.h"
#include "TbPatternRecognition.h"
#include "TbROOT.h"
#include "TbTrackAlgorithms.h"
#include "DD4hep/Factories.h"
/* Header for PlotTool
 *
 * 2014-02-25
 *
 */

#include "TbBaseClass.h"

using namespace ROOT::Math;
using namespace dd4hep;

class TbPlotTool : public TbBaseClass {
 public:
  /// Constructor
  TbPlotTool(const std::string &name);
  /// Destructor
  virtual ~TbPlotTool(){};

  bool configuration();
  bool initialize(AlgVec);  ///< Algorithm initialization
  bool execute(dd4hep::cond::ConditionsSlice &, AlgVec);     ///< Algorithm execution
  bool end_event();
  bool finalize(dd4hep::cond::ConditionsSlice &);  ///< Algorithm finalization
  TbBaseClass *find(AlgVec vec, std::string name) {
    for (AlgVec::iterator it = vec.begin(); it != vec.end(); ++it) {
      if ((*it).first == name) return (*it).second;
    }
    std::cout << "Couldn't find " << name << std::endl;
    return NULL;
  }

 private:
  TbROOT *tbroot;
  mutable TbGeometrySvc *m_geomSvc;
  mutable TbPatternRecognition *m_patternrec;
  TbTrackAlgorithms *tral;
  AlgVec m_algos;
};
#endif
