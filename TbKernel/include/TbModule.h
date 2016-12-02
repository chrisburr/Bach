#ifndef TB_MODULE_H
#define TB_MODULE_H 1

#include "Math/RotationZYX.h"
#include "Math/Transform3D.h"
#include "Math/Translation3D.h"

class TbModule {
 public:
  TbModule(){};
  virtual ~TbModule(){};

  void SetAlignment(std::string id, float x, float y, float z, float rotx,
                    float roty, float rotz, float dx, float dy, float dz,
                    float drotx, float droty, float drotz) {
    m_id = id;
    m_x = x;
    m_y = y;
    m_z = z;
    m_rotx = rotx;
    m_roty = roty;
    m_rotz = rotz;
    m_dx = dx;
    m_dy = dy;
    m_dz = dz;
    m_drotx = drotx;
    m_droty = droty;
    m_drotz = drotz;
    SetTransform();
  }
  void SetAlignment(float dx, float dy, float dz, float drotx, float droty,
                    float drotz) {
    m_dx = dx;
    m_dy = dy;
    m_dz = dz;
    m_drotx = drotx;
    m_droty = droty;
    m_drotz = drotz;
    SetTransform();
  }
  void SetTransform() {
    ROOT::Math::Translation3D shift(m_x + m_dx, m_y + m_dy, m_z + m_dz);
    ROOT::Math::RotationZYX rot(m_rotz + m_drotz, m_roty + m_droty,
                                m_rotx + m_drotx);
    ROOT::Math::Transform3D transform(rot, shift);
    m_transform = transform;
  }

  Transform3D Transform() { return m_transform; }

  std::string id() { return m_id; }
  float X() { return m_x; }
  float Y() { return m_y; }
  float Z() { return m_z; }
  float RotX() { return m_rotx; }
  float RotY() { return m_roty; }
  float RotZ() { return m_rotz; }
  float dX() { return m_dx; }
  float dY() { return m_dy; }
  float dZ() { return m_dz; }
  float dRotX() { return m_drotx; }
  float dRotY() { return m_droty; }
  float dRotZ() { return m_drotz; }

  void Nr(int nr) { m_nr = nr; }
  int Nr() { return m_nr; }

 private:
  std::string m_id;
  int m_nr;
  // Positions
  float m_x;
  float m_y;
  float m_z;
  float m_rotx;
  float m_roty;
  float m_rotz;
  // Alignment
  float m_dx;
  float m_dy;
  float m_dz;
  float m_drotx;
  float m_droty;
  float m_drotz;

  Transform3D m_transform;
};

#endif
