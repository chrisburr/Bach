
#include "TbTrackAlgorithms.h"
#include "DD4hep/Factories.h"

#include <math.h>

/** @file TbTrackAlgorithms.cpp
 *
 *  Implementation of class : TbTrackAlgorithms
 *
 *  Class for common Algorithms for tracks (interceptions, fit)
 *
 */

//=============================================================================
/// Standard constructor
//=============================================================================
TbTrackAlgorithms::TbTrackAlgorithms(const std::string &name) {}

bool TbTrackAlgorithms::configuration() { return true; }

//=============================================================================
/// Initialization
//=============================================================================
bool TbTrackAlgorithms::initialize(DD4hep::Geometry::LCDD &lcdd, AlgVec algos) { return true; }

//=============================================================================
/// Main execution
//=============================================================================
bool TbTrackAlgorithms::execute(AlgVec algos) { return true; }

//=============================================================================
/// Finalize
//=============================================================================
bool TbTrackAlgorithms::finalize() { return true; }

XYZPoint TbTrackAlgorithms::getInterceptGlobal(TbTrack *track,
                                               const std::string id) {
  XYZPoint planePointLocalCoords(0., 0., 0.);
  XYZPoint planePointGlobalCoords =
      Geom()->localToGlobal(planePointLocalCoords, id);

  XYZPoint planePointLocalCoords_2(0., 0., 1.);
  XYZPoint planePointGlobalCoords_2 =
      Geom()->localToGlobal(planePointLocalCoords_2, id);

  const double normal_x =
      planePointGlobalCoords_2.X() - planePointGlobalCoords.X();
  const double normal_y =
      planePointGlobalCoords_2.Y() - planePointGlobalCoords.Y();
  const double normal_z =
      planePointGlobalCoords_2.Z() - planePointGlobalCoords.Z();

  const double length =
      ((planePointGlobalCoords.X() - track->firstState()->X()) * normal_x +
       (planePointGlobalCoords.Y() - track->firstState()->Y()) * normal_y +
       (planePointGlobalCoords.Z() - track->firstState()->Z()) * normal_z) /
      (track->direction()->X() * normal_x + track->direction()->Y() * normal_y +
       track->direction()->Z() * normal_z);
  const float x_inter_global =
      track->firstState()->X() + length * track->direction()->X();
  const float y_inter_global =
      track->firstState()->Y() + length * track->direction()->Y();
  const float z_inter_global =
      track->firstState()->Z() + length * track->direction()->Z();

  PositionVector3D<Cartesian3D<double>> global(x_inter_global, y_inter_global,
                                               z_inter_global);
  return global;
}
XYZPoint TbTrackAlgorithms::getInterceptGlobal_out(TbTrack *track,
                                                   const std::string id) {
  XYZPoint planePointLocalCoords(0., 0., 0. + Geom()->Const_D("Thick") / 2.);
  XYZPoint planePointGlobalCoords =
      Geom()->localToGlobal(planePointLocalCoords, id);

  XYZPoint planePointLocalCoords_2(0., 0., 1. + Geom()->Const_D("Thick") / 2.);
  XYZPoint planePointGlobalCoords_2 =
      Geom()->localToGlobal(planePointLocalCoords_2, id);

  const double normal_x =
      planePointGlobalCoords_2.X() - planePointGlobalCoords.X();
  const double normal_y =
      planePointGlobalCoords_2.Y() - planePointGlobalCoords.Y();
  const double normal_z =
      planePointGlobalCoords_2.Z() - planePointGlobalCoords.Z();

  const double length =
      ((planePointGlobalCoords.X() - track->firstState()->X()) * normal_x +
       (planePointGlobalCoords.Y() - track->firstState()->Y()) * normal_y +
       (planePointGlobalCoords.Z() - track->firstState()->Z()) * normal_z) /
      (track->direction()->X() * normal_x + track->direction()->Y() * normal_y +
       track->direction()->Z() * normal_z);
  const float x_inter_global =
      track->firstState()->X() + length * track->direction()->X();
  const float y_inter_global =
      track->firstState()->Y() + length * track->direction()->Y();
  const float z_inter_global =
      track->firstState()->Z() + length * track->direction()->Z();
  PositionVector3D<Cartesian3D<double>> global(x_inter_global, y_inter_global,
                                               z_inter_global);
  return global;
}
XYZPoint TbTrackAlgorithms::getInterceptGlobal_in(TbTrack *track,
                                                  const std::string id) {
  XYZPoint planePointLocalCoords(0., 0., 0. - Geom()->Const_D("Thick") / 2.);
  XYZPoint planePointGlobalCoords =
      Geom()->localToGlobal(planePointLocalCoords, id);

  XYZPoint planePointLocalCoords_2(0., 0., 1. - Geom()->Const_D("Thick") / 2.);
  XYZPoint planePointGlobalCoords_2 =
      Geom()->localToGlobal(planePointLocalCoords_2, id);

  const double normal_x =
      planePointGlobalCoords_2.X() - planePointGlobalCoords.X();
  const double normal_y =
      planePointGlobalCoords_2.Y() - planePointGlobalCoords.Y();
  const double normal_z =
      planePointGlobalCoords_2.Z() - planePointGlobalCoords.Z();

  const double length =
      ((planePointGlobalCoords.X() - track->firstState()->X()) * normal_x +
       (planePointGlobalCoords.Y() - track->firstState()->Y()) * normal_y +
       (planePointGlobalCoords.Z() - track->firstState()->Z()) * normal_z) /
      (track->direction()->X() * normal_x + track->direction()->Y() * normal_y +
       track->direction()->Z() * normal_z);
  const float x_inter_global =
      track->firstState()->X() + length * track->direction()->X();
  const float y_inter_global =
      track->firstState()->Y() + length * track->direction()->Y();
  const float z_inter_global =
      track->firstState()->Z() + length * track->direction()->Z();
  PositionVector3D<Cartesian3D<double>> global(x_inter_global, y_inter_global,
                                               z_inter_global);
  return global;
}
XYZPoint TbTrackAlgorithms::getInterceptLocal(TbTrack *track,
                                              const std::string id) {
  XYZPoint global = getInterceptGlobal(track, id);

  XYZPoint local = Geom()->globalToLocal(global, id);

  return local;
}

void TbTrackAlgorithms::FitTrack(TbTrack *track) {
  double vecx[2] = {0., 0.};
  double vecy[2] = {0., 0.};
  double matx[2][2] = {{0., 0.}, {0., 0.}};
  double maty[2][2] = {{0., 0.}, {0., 0.}};

  TbClusters *clusters = track->Clusters();
  TbClusters::iterator it;
  int nd(clusters->size());
  const double error = 0.015;
  for (it = clusters->begin(); it != clusters->end(); ++it) {
    if ((*it) == 0) continue;

    const double x = (*it)->GlobalPos().X();
    const double y = (*it)->GlobalPos().Y();
    const double z = (*it)->GlobalPos().Z();

    const double wx = 1. / (error * error);
    const double wy = wx;
    vecx[0] += wx * x;
    vecx[1] += wx * x * z;
    vecy[0] += wy * y;
    vecy[1] += wy * y * z;
    matx[0][0] += wx;
    matx[1][0] += wx * z;
    matx[1][1] += wx * z * z;
    maty[0][0] += wy;
    maty[1][0] += wy * z;
    maty[1][1] += wy * z * z;
  }

  // Invert the matrix and compute track parameters.
  double detx = matx[0][0] * matx[1][1] - matx[1][0] * matx[1][0];
  double dety = maty[0][0] * maty[1][1] - maty[1][0] * maty[1][0];
  // Check for singularities.
  if (detx == 0. || dety == 0.) return;

  double slopex = (vecx[1] * matx[0][0] - vecx[0] * matx[1][0]) / detx;
  double slopey = (vecy[1] * maty[0][0] - vecy[0] * maty[1][0]) / dety;

  double interceptx = (vecx[0] * matx[1][1] - vecx[1] * matx[1][0]) / detx;
  double intercepty = (vecy[0] * maty[1][1] - vecy[1] * maty[1][0]) / dety;

  track->firstState(interceptx, intercepty, 0);
  track->direction(slopex, slopey, 1);

  // Compute chi2.
  float m_chi2 = 0.;
  for (it = clusters->begin(); it != clusters->end(); ++it) {
    const float x = (*it)->GlobalPos().X();
    const float y = (*it)->GlobalPos().Y();
    const float z = (*it)->GlobalPos().Z();

    const double wx = error;
    const double wy = error;

    const double dz = fabs(z - track->firstState()->Z());
    const float xfit = track->firstState()->X() + dz * track->slopeXZ();
    const float yfit = track->firstState()->Y() + dz * track->slopeYZ();
    const double dx = x - xfit;
    const double dy = y - yfit;
    m_chi2 += (dx * dx) / (wx * wx) + (dy * dy) / (wy * wy);
  }

  double m_ndof = 2 * nd - 4;
  track->chi2(m_chi2);
  track->ndof(m_ndof);
}
