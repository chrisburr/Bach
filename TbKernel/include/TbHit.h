#ifndef TB_HIT_H
#define TB_HIT_H 1

#include <string>
#include <vector>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
using namespace ROOT::Math;
class TbHit {
public:
  TbHit() {
    m_mcx = 0.;
    m_mcy = 0.;
    m_incluster = false;
  }
  ~TbHit(){};
  std::string id() { return m_id; }
  int col() { return m_col; }
  int row() { return m_row; }
  float adc() { return m_adc; }

  int id_nr() { return m_idnr; }
  void id_nr(int in) { m_idnr = in; }
  bool incluster() { return m_incluster; }
  void setId(std::string id) { m_id = id; }
  void setCol(int col) { m_col = col; }
  void setRow(int row) { m_row = row; }
  void setAdc(float adc) { m_adc = adc; }
  void incluster(bool in) { m_incluster = in; }
  float mcx() { return m_mcx; }
  void mcx(float mx) { m_mcx = mx; }
  float mcy() { return m_mcy; }
  void mcy(float my) { m_mcy = my; }

  float mctrslxz() { return m_mctrslxz; }
  void mctrslxz(double mctsx) { m_mctrslxz = mctsx; }
  double mctrslyz() { return m_mctrslyz; }
  void mctrslyz(double mctsy) { m_mctrslyz = mctsy; }
  ROOT::Math::XYZPoint mctrintercept() { return m_mctrintercept; }
  void mctrintercept(ROOT::Math::XYZPoint mctri) { m_mctrintercept = mctri; }

private:
  std::string m_id;
  int m_idnr;
  int m_col;
  int m_row;
  float m_adc;
  float m_mcx;
  float m_mcy;
  bool m_incluster;
  double m_mctrslxz;
  double m_mctrslyz;
  ROOT::Math::XYZPoint m_mctrintercept;
};

typedef std::vector<TbHit *> TbHits;

#endif
