#ifndef TB_CLUSTER_H
#define TB_CLUSTER_H 1
#include <string>
#include <vector>

#include "Math/Point3D.h"
#include "DD4hep/Factories.h"
using namespace ROOT::Math;
using namespace DD4hep::Geometry;

//#include "TMath.h"
#include "TbHit.h"
class TbCluster {
 public:
  std::string id() { return m_id; }
  void id(std::string Id) { m_id = Id; }
  TbHits *HitContainer() { return m_hits; }
  void HitContainer(TbHits *hits) { m_hits = hits; }
  Position LocalPos() { return m_local_pos; }
  void LocalPos(Position lp) { m_local_pos = lp; }
  Position GlobalPos() { return m_global_pos; }
  void GlobalPos(Position gp) { m_global_pos = gp; }
  int XCount() { return m_xcount; }
  void XCount(int xc) { m_xcount = xc; }
  int YCount() { return m_ycount; }
  void YCount(int yc) { m_ycount = yc; }

 private:
  std::string m_id;
  int m_xcount;
  int m_ycount;
  TbHits *m_hits;
  Position m_local_pos;
  Position m_global_pos;
};

typedef std::vector<TbCluster *> TbClusters;
#endif
