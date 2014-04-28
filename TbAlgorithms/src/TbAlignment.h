
#ifndef TB_ALIGNMENT_H 
#define TB_ALIGNMENT_H 1

/* Header for Alignment
 *
 * 2014-02-28
 *
 */

#include "../../TbKernel/src/TbBaseClass.h"
#include "../../TbKernel/src/Alignment/Millepede.h"
#include "../../TbKernel/src/TbGeometrySvc.h"
#include "TbPatternRecognition.h"
#include "TbTrackAlgorithms.h"
using namespace ROOT::Math;
class TbAlignment : public TbBaseClass {

public: 

  /// Constructor
  TbAlignment(const std::string& name);
  /// Destructor
  virtual ~TbAlignment();

  bool configuration();
  bool initialize(AlgVec);    ///< Algorithm initialization
  bool execute(AlgVec );       ///< Algorithm execution
  bool end_event();
  bool finalize();      ///< Algorithm finalization
  TbBaseClass* find( AlgVec vec , std::string name){
    for ( AlgVec::iterator it = vec.begin(); it!=vec.end(); ++it){
      if ((*it).first == name) return (*it).second;
    }
    std::cout << "Couldn't find " << name << std::endl;
    return NULL;
  }
  bool PutTrack( TbTrack* , int , int , bool* );
  bool PutTrack2( TbTrack* , int , int , bool* );
  TbTracks* GetTracks() {return m_trackcontainer;}
  int detectoridentifier(std::string id) {
    int detnr = -1;
    for (std::vector<TbModule*>::iterator itm = m_modulestoalign->begin();
	 itm != m_modulestoalign->end(); ++itm){
      if ((*itm)->id() == id) detnr = (*itm)->Nr();
    }
    return detnr;
  }
 
  TbGeometrySvc* GetGeom() {return m_geomSvc;}
 private:
  mutable Millepede * m_millepede;
  mutable TbPatternRecognition * m_patternrec;
  TbGeometrySvc * m_geomSvc;
  TbTrackAlgorithms *tral;
  bool m_debug;
  TbTracks *m_trackcontainer;
  std::map<int,double> TelescopeMap;
  std::map<std::string, int> m_detNum;
  std::vector<TbModule*> *m_modulestoalign;

  std::vector<double>   ftx;
  std::vector<double>   fty;
  std::vector<double>   ftz;
  std::vector<double>   frotx; 
  std::vector<double>   froty; 
  std::vector<double>   frotz;  
  std::vector<double>   fscaz;
  std::vector<double>   shearx; 
  std::vector<double>   sheary;
};
#endif
    
