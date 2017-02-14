
#ifndef TB_ALIGNMENT_H
#define TB_ALIGNMENT_H 1

/* Header for Alignment
 *
 * 2014-02-28
 *
 */

#include "DD4hep/AlignmentsPrinter.h"
#include "DD4hep/Conditions.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetAlign.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Factories.h"
#include "DD4hep/Memory.h"
#include "DDAlign/AlignmentsForward.h"
#include "DDAlign/AlignmentsRegister.h"
#include "DDAlign/DDAlignForwardCall.h"
#include "DDAlign/DDAlignUpdateCall.h"
#include "DDCond/ConditionsManager.h"
#include "DDAlign/AlignmentsManager.h"
#include "DDCond/ConditionsSlice.h"
#include "Millepede.h"
#include "TbBaseClass.h"
#include "TbGeometrySvc.h"
#include "TbPatternRecognition.h"
#include "TbTrackAlgorithms.h"

using namespace ROOT::Math;
using namespace DD4hep::Geometry;

class TbAlignment : public TbBaseClass {
public:
  /// Constructor
  TbAlignment(DD4hep::Alignments::AlignmentsManager, const std::string &name);
  /// Destructor
  virtual ~TbAlignment();

  bool configuration();
  bool initialize(AlgVec); ///< Algorithm initialization
  bool execute(DD4hep::Conditions::ConditionsSlice &, AlgVec);    ///< Algorithm execution
  bool end_event();
  bool finalize(DD4hep::Conditions::ConditionsSlice &); ///< Algorithm finalization
  TbBaseClass *find(AlgVec vec, std::string name) {
    for (AlgVec::iterator it = vec.begin(); it != vec.end(); ++it) {
      if ((*it).first == name)
        return (*it).second;
    }
    std::cout << "Couldn't find " << name << std::endl;
    return NULL;
  }
  bool PutTrack(TbTrack *, int, int, bool *);
  bool PutTrack2(TbTrack *, int, int, bool *);
  TbTracks *GetTracks() { return m_trackcontainer; }
  int detectoridentifier(std::string id) {
    auto it = find_if(m_modulestoalign.begin(), m_modulestoalign.end(),
                      [&id](const DetElement &e) { return e.path() == id; });
    int detnr = it - m_modulestoalign.begin();
    return detnr;
  }

  TbGeometrySvc *GetGeom() { return m_geomSvc; }

private:
  mutable Millepede *m_millepede;
  mutable TbPatternRecognition *m_patternrec;
  DD4hep::Alignments::AlignmentsManager m_alignMgr;
  TbGeometrySvc *m_geomSvc;
  TbTrackAlgorithms *tral;
  bool m_debug;
  TbTracks *m_trackcontainer;
  std::map<std::string, int> m_detNum;
  std::vector<DetElement> m_modulestoalign;

  std::vector<double> ftx;
  std::vector<double> fty;
  std::vector<double> ftz;
  std::vector<double> frotx;
  std::vector<double> froty;
  std::vector<double> frotz;
  std::vector<double> fscaz;
  std::vector<double> shearx;
  std::vector<double> sheary;
};
#endif
