
#ifndef TB_TRACKALGORITHMS_H 
#define TB_TRACKALGORITHMS_H 1

/* Header for TrackAlgorithms
 *
 * 2014-02-07
 *
 */

#include "../../TbKernel/src/TbBaseClass.h"
#include "../../TbKernel/src/TbTrack.h"
#include "../../TbKernel/src/TbGeometrySvc.h"


#include "Math/Point3D.h"
#include "Math/Vector3D.h"
using namespace ROOT::Math;
class TbTrackAlgorithms : public TbBaseClass {

public: 

  /// Constructor
  TbTrackAlgorithms(const std::string& name );
  /// Destructor
  virtual ~TbTrackAlgorithms(){};

  bool configuration();
  bool initialize(AlgVec);    ///< Algorithm initialization
  bool execute(AlgVec );       ///< Algorithm execution
  bool finalize();      ///< Algorithm finalization

 
  void FitTrack(TbTrack*);
  
  XYZPoint getInterceptGlobal(TbTrack* , const std::string);
  XYZPoint getInterceptGlobal_out(TbTrack* , const std::string);
  XYZPoint getInterceptGlobal_in(TbTrack* , const std::string);
  XYZPoint getInterceptLocal(TbTrack* , const std::string);
   /// Access geometry service on-demand
  TbGeometrySvc * Geom() const {
    return m_geomSvc;
  }
  void setGeom(TbGeometrySvc * geo) {m_geomSvc = geo;}

 private:
  mutable TbGeometrySvc * m_geomSvc;


};
#endif
    
