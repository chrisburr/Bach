#ifndef TB_TRACK_H
#define TB_TRACK_H 1
/* Header for Track
 *
 * 2014-02-07
 *
 */
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "TbCluster.h"
class TbTrack {

public:
  TbTrack() { m_clusters = new TbClusters; }
  void firstState(float x, float y, float z) {
    m_state.SetX(x);
    m_state.SetY(y);
    m_state.SetZ(z);
  }
  ROOT::Math::XYZPoint *firstState() { return &m_state; }

  void direction(float slope_xz, float slope_yz) {

    m_dir.SetX(slope_xz);
    m_dir.SetY(slope_yz);
    m_dir.SetZ(1.);
  }

  void direction(float x, float y, float z) {
    float r = sqrt(x * x + y * y + z * z);
    m_dir.SetX(x / r);
    m_dir.SetY(y / r);
    m_dir.SetZ(z / r);
  }
  ROOT::Math::XYZVector *direction() { return &m_dir; }

  float slopeXZ() { return m_dir.X() / m_dir.Z(); }
  float slopeYZ() { return m_dir.Y() / m_dir.Z(); }

  void chi2(float c2) { m_chi2 = c2; }
  float chi2() { return m_chi2; }

  void ndof(float ndf) { m_ndof = ndf; }
  float ndof() { return m_ndof; }

  TbClusters *Clusters() { return m_clusters; }

private:
  mutable TbClusters *m_clusters;
  // First measured state of the track
  ROOT::Math::XYZPoint m_state;
  // Direction of the track
  ROOT::Math::XYZVector m_dir;

  float m_chi2;
  float m_ndof;
};

typedef std::vector<TbTrack *> TbTracks;

#endif
