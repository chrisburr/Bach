
#ifndef TB_PATTERNRECOGNITION_H 
#define TB_PATTERNRECOGNITION_H 1

/* Header for PatternRecognition
 *
 * 2014-02-20
 *
 */

#include "../../TbKernel/src/TbBaseClass.h"
#include "../../TbKernel/src/TbTrack.h"
#include "../../TbKernel/src/TbGeometrySvc.h"
#include "TbTrackAlgorithms.h"
#include "TbClustering.h"
class TbPatternRecognition : public TbBaseClass {

public: 

  /// Constructor
  TbPatternRecognition(const std::string& name);
  /// Destructor
  virtual ~TbPatternRecognition();

  bool configuration();
  bool initialize(AlgVec);    ///< Algorithm initialization
  bool execute(AlgVec );       ///< Algorithm execution
  bool end_event();
  bool finalize();      ///< Algorithm finalization
  TbBaseClass* find( AlgVec vec , std::string name){
    for ( AlgVec::iterator it = vec.begin(); it!=vec.end(); ++it){
      if ((*it).first == name) return (*it).second;
    }
   }
   TbClustering * tbcluster() const {
      return m_tbcluster;
    }
    void tbcluster(TbClustering * tbc) {m_tbcluster = tbc;}
  TbTracks*  Tracks() {return m_tracks;}

  
 private:
    mutable TbClustering *  m_tbcluster;
    TbTracks  *m_tracks;
    TbGeometrySvc *m_geom;
    TbTrackAlgorithms *m_tral;
    std::string m_refmod;
    double m_dis;
    double m_chi2cut;
  };
#endif
    
