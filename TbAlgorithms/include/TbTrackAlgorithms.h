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

class TbTrackAlgorithms : public TbBaseClass {
 public:
  /// Constructor
  TbTrackAlgorithms(const std::string &name);
  /// Destructor
  virtual ~TbTrackAlgorithms(){};

  bool configuration();
  bool initialize(AlgVec);  ///< Algorithm initialization
  bool execute(AlgVec);     ///< Algorithm execution
  bool finalize();          ///< Algorithm finalization

  void FitTrack(TbTrack *);

  Position getInterceptGlobal(TbTrack *, std::string);
  Position getInterceptGlobal(TbTrack *, DetElement);
  Position getInterceptGlobal_out(TbTrack *, std::string);
  Position getInterceptGlobal_out(TbTrack *, DetElement);
  Position getInterceptGlobal_in(TbTrack *, std::string);
  Position getInterceptGlobal_in(TbTrack *, DetElement);
  Position getInterceptLocal(TbTrack *, std::string);
  Position getInterceptLocal(TbTrack *, DetElement);
  /// Access geometry service on-demand
  TbGeometrySvc *Geom() const { return m_geomSvc; }
  void setGeom(TbGeometrySvc *geo) { m_geomSvc = geo; }

  private:
  mutable TbGeometrySvc *m_geomSvc;
};
#endif
