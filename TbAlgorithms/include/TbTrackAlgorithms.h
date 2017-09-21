#ifndef TB_TRACKALGORITHMS_H
#define TB_TRACKALGORITHMS_H 1

/* Header for TrackAlgorithms
 *
 * 2014-02-07
 *
 */

#include "TbGeometrySvc.h"
#include "TbBaseClass.h"
#include "TbTrack.h"
#include "DD4hep/Factories.h"

using namespace DD4hep::Geometry;
using DD4hep::Conditions::ConditionsSlice;

class TbTrackAlgorithms : public TbBaseClass {
 public:
  /// Constructor
  TbTrackAlgorithms(const std::string &name);
  /// Destructor
  virtual ~TbTrackAlgorithms(){};

  bool configuration();
  bool initialize(AlgVec);  ///< Algorithm initialization
  bool execute(ConditionsSlice &, AlgVec);     ///< Algorithm execution
  bool finalize(DD4hep::Conditions::ConditionsSlice &);          ///< Algorithm finalization

  void FitTrack(TbTrack *);

  Position getInterceptGlobal(TbTrack *, std::string, ConditionsSlice &);
  Position getInterceptGlobal(TbTrack *, DetElement, ConditionsSlice &);
  Position getInterceptGlobal_out(TbTrack *, std::string, ConditionsSlice &);
  Position getInterceptGlobal_out(TbTrack *, DetElement, ConditionsSlice &);
  Position getInterceptGlobal_in(TbTrack *, std::string, ConditionsSlice &);
  Position getInterceptGlobal_in(TbTrack *, DetElement, ConditionsSlice &);
  Position getInterceptLocal(TbTrack *, std::string, ConditionsSlice &);
  Position getInterceptLocal(TbTrack *, DetElement, ConditionsSlice &);
  /// Access geometry service on-demand
  TbGeometrySvc *Geom() const { return m_geomSvc; }
  void setGeom(TbGeometrySvc *geo) { m_geomSvc = geo; }

  private:
  mutable TbGeometrySvc *m_geomSvc;
};
#endif
