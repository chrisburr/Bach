#ifndef TB_CLUSTER_H 
#define TB_CLUSTER_H 1
#include <string>
#include <vector>

#include "Math/Point3D.h"
using namespace ROOT::Math;
//#include "TMath.h"
#include "TbHit.h" 
class TbCluster {
 public:
  
  std::string id() {return m_id;}
  void id(std::string Id) {m_id = Id;}
  TbHits* HitContainer(){return m_hits;}
  void HitContainer(TbHits* hits) {m_hits = hits;}
  ROOT::Math::XYZPoint LocalPos() {return m_local_pos;}
  void LocalPos(ROOT::Math::XYZPoint lp) {m_local_pos = lp;}
  ROOT::Math::XYZPoint GlobalPos() {return m_global_pos;}
  void GlobalPos(ROOT::Math::XYZPoint gp) {m_global_pos = gp;}
  int XCount() {return m_xcount;}
  void XCount(int xc) {m_xcount = xc;}
  int YCount() {return m_ycount;}
  void YCount(int yc) {m_ycount = yc;}
 private:
  std::string m_id;
  int m_xcount;
  int m_ycount;
  TbHits* m_hits;
  ROOT::Math::XYZPoint m_local_pos;
  ROOT::Math::XYZPoint m_global_pos;
};


typedef std::vector<TbCluster*> TbClusters;
#endif
