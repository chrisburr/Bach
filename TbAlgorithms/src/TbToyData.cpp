
#include "TbToyData.h"
#include "TbTrackAlgorithms.h"

#include "TFile.h"
#include "TTree.h"

#include "TMath.h"

/** @file TbToyData.cpp
 *
 *  Implementation of class : TbToyData
 *
 *  Class to generate tracks, determine interception points and store
 * corresponding hits
 *
 */

//=============================================================================
/// Standard constructor
//=============================================================================
TbToyData::TbToyData(const std::string &name) {}
TbToyData::~TbToyData() { delete m_hits; }
bool TbToyData::configuration() {
  Const_D("MeanX", 0.0);
  Const_D("MeanY", 0.0);
  Const_D("SigmaX", 0.0001);
  Const_D("SigmaY", 0.0001);
  Const_B("Write_Txt", false);
  Const_S("TxtFile", "TxtFile.txt");
  Const_I("NoOfTracks", 10);
  Const_I("Seed", 0);
  Const_S("GeometryFile", "");
  return true;
}

//=============================================================================
/// Initialization
//=============================================================================
bool TbToyData::initialize(AlgVec algos) {
  TbGeometrySvc *geo = new TbGeometrySvc("misaligned");
  geomSvc(geo);
  geomSvc()->configuration();
  geomSvc()->Const_S("GeometryFile", Const_S("GeometryFile"));
  m_seed = Const_I("Seed");

  geomSvc()->initialize(algos);
  if (Const_B("Write_Txt")) {
    datei = fopen((char *)Const_S("TxtFile").c_str(), "w");
    if (datei == NULL) {
      std::cout << "Could not open " << Const_S("TxtFile") << std::endl;
      return false;
    }
  }
  m_hits = new TbHits();
  return true;
}

//=============================================================================
/// Main execution
//=============================================================================
bool TbToyData::execute(AlgVec algos) {
  // Get geometry parameter

  const double pitch_x = geomSvc()->Const_D("PitchX");
  const double pitch_y = geomSvc()->Const_D("PitchY");
  const double noofpix_x = geomSvc()->Const_I("NoOfPixelX");
  const double noofpix_y = geomSvc()->Const_I("NoOfPixelY");

  double x_min = -pitch_x * noofpix_x / 2.;
  double x_max = pitch_x * noofpix_x / 2.;
  double y_min = -pitch_y * noofpix_y / 2.;
  double y_max = pitch_y * noofpix_y / 2.;

  TbTrackAlgorithms *tral = new TbTrackAlgorithms("TrackAlgos");
  tral->setGeom(geomSvc());
  std::map<std::string, TbModule *> align = geomSvc()->Modules;
  typedef std::map<std::string, TbModule *>::iterator it_type;

  int n_event = 0;

  int hitidnr = 0;
  m_r.SetSeed(m_seed++);

  // Generate Tracks
  for (int j = 0; j < Const_I("NoOfTracks"); ++j) {
    double slopex = m_r.Gaus(Const_D("MeanX"), Const_D("SigmaX"));
    double slopey = m_r.Gaus(Const_D("MeanX"), Const_D("SigmaX"));

    double interx = m_r.Uniform(x_min, x_max);
    double intery = m_r.Uniform(y_min, y_max);

    TbTrack *m_tr = new TbTrack;
    m_tr->firstState(interx, intery, 0.);
    m_tr->direction(slopex, slopey);

    for (it_type iterator = align.begin(); iterator != align.end();
         iterator++) {
      // Get intercepts and corresponding pixel hits
      const std::string id = iterator->first;
      XYZPoint global_intercept = tral->getInterceptGlobal(m_tr, id);
      XYZPoint global_intercept_in = tral->getInterceptGlobal_in(m_tr, id);
      XYZPoint global_intercept_out = tral->getInterceptGlobal_out(m_tr, id);

      XYZPoint local_intercept = geomSvc()->globalToLocal(global_intercept, id);

      float pix_x = local_intercept.X() / pitch_x + noofpix_x / 2.;
      float pix_y = local_intercept.Y() / pitch_y + noofpix_y / 2.;

      XYZPoint local_intercept_in =
          geomSvc()->globalToLocal(global_intercept_in, id);
      float pix_x_in = local_intercept_in.X() / pitch_x + noofpix_x / 2.;
      float pix_y_in = local_intercept_in.Y() / pitch_y + noofpix_y / 2.;

      XYZPoint local_intercept_out =
          geomSvc()->globalToLocal(global_intercept_out, id);
      float pix_x_out = local_intercept_out.X() / pitch_x + noofpix_x / 2.;
      float pix_y_out = local_intercept_out.Y() / pitch_y + noofpix_y / 2.;

      // Delete hit, if it is outside the sensorregion
      if (pix_x_in < 0 || pix_y_in < 0 || pix_x_out < 0 || pix_y_out < 0 ||
          pix_x_in > noofpix_x || pix_y_in > noofpix_y ||
          pix_x_out > noofpix_x || pix_y_out > noofpix_y)
        continue;
      float length = TMath::Sqrt(
          TMath::Power(global_intercept_out.X() - global_intercept_in.X(), 2.) +
          TMath::Power(global_intercept_out.Y() - global_intercept_in.Y(), 2.) +
          TMath::Power(global_intercept_out.Z() - global_intercept_in.Z(), 2.));

      int adc = 100000 * length;

      // 1 hit cluster
      if (((int)pix_x_in == (int)pix_x_out) &&
          ((int)pix_y_in == (int)pix_y_out)) {
        TbHit *hit = new TbHit;
        hit->setId(id);
        hit->id_nr(hitidnr);
        hit->setCol((int)pix_y);
        hit->setRow((int)pix_x);
        hit->setAdc(adc);
        hit->mcx(global_intercept.X());
        hit->mcy(global_intercept.Y());
        hit->mctrslxz(slopex);
        hit->mctrslyz(slopey);
        hit->mctrintercept(*(m_tr->firstState()));
        m_hits->push_back(hit);
        ++hitidnr;
      }

      // 2 hit custer
      else if (((int)pix_x_in == (int)pix_x_out) &&
               ((int)pix_y_in != (int)pix_y_out)) {
        int ref_y = (int)pix_y_in + ((int)pix_y_out - (int)pix_y_in) / 2. + 1;
        float l1 = (ref_y - pix_y_in) / (pix_y_out - pix_y_in);
        float l2 = 1. - l1;
        TbHit *hita = new TbHit;
        hita->setId(id);
        hita->id_nr(hitidnr);
        hita->setCol((int)pix_y_in);
        hita->setRow((int)pix_x_in);
        hita->setAdc(l1 * adc);
        hita->mcx(global_intercept.X());
        hita->mcy(global_intercept.Y());
        hita->mctrslxz(slopex);
        hita->mctrslyz(slopey);
        hita->mctrintercept(*(m_tr->firstState()));
        m_hits->push_back(hita);
        ++hitidnr;
        TbHit *hitb = new TbHit;
        hitb->setId(id);
        hitb->id_nr(hitidnr);
        hitb->setCol((int)pix_y_out);
        hitb->setRow((int)pix_x_in);
        hitb->setAdc(l2 * adc);
        hitb->mcx(global_intercept.X());
        hitb->mcy(global_intercept.Y());
        hitb->mctrslxz(slopex);
        hitb->mctrslyz(slopey);
        hitb->mctrintercept(*(m_tr->firstState()));
        m_hits->push_back(hitb);
        ++hitidnr;
      } else if (((int)pix_y_in == (int)pix_y_out) &&
                 ((int)pix_x_in != (int)pix_x_out)) {
        int ref_x = (int)pix_x_in + ((int)pix_x_out - (int)pix_x_in) / 2. + 1;
        float l1 = (ref_x - pix_x_in) / (pix_x_out - pix_x_in);
        float l2 = 1. - l1;
        TbHit *hita = new TbHit;
        hita->setId(id);
        hita->id_nr(hitidnr);
        hita->setCol((int)pix_y_in);
        hita->setRow((int)pix_x_in);
        hita->setAdc(l1 * adc);
        hita->mcx(global_intercept.X());
        hita->mcy(global_intercept.Y());
        hita->mctrslxz(slopex);
        hita->mctrslyz(slopey);
        hita->mctrintercept(*(m_tr->firstState()));
        m_hits->push_back(hita);
        ++hitidnr;
        TbHit *hitb = new TbHit;
        hitb->setId(id);
        hitb->id_nr(hitidnr);
        hitb->setCol((int)pix_y_in);
        hitb->setRow((int)pix_x_out);
        hitb->setAdc(l2 * adc);
        hitb->mcx(global_intercept.X());
        hitb->mcy(global_intercept.Y());
        hitb->mctrslxz(slopex);
        hitb->mctrslyz(slopey);
        hitb->mctrintercept(*(m_tr->firstState()));
        m_hits->push_back(hitb);
        ++hitidnr;
      }

      // 3 hit cluster
      else {
        int ref_x = (int)pix_x_in + ((int)pix_x_out - (int)pix_x_in) / 2. + 1;
        int ref_y = (int)pix_y_in + ((int)pix_y_out - (int)pix_y_in) / 2. + 1;
        float l1 = (ref_x - pix_x_in) / (pix_x_out - pix_x_in);
        float l2 = (ref_y - pix_y_in) / (pix_y_out - pix_y_in);

        if (l1 > l2) l1 = 1. - l1;
        if (l2 > l1) l2 = 1. - l2;
        float l3 = 1. - l1 - l2;

        TbHit *hita = new TbHit;
        hita->setId(id);
        hita->id_nr(hitidnr);
        hita->setCol((int)pix_y_in);
        hita->setRow((int)pix_x_in);
        hita->setAdc(l1 * adc);
        hita->mcx(global_intercept.X());
        hita->mcy(global_intercept.Y());
        hita->mctrslxz(slopex);
        hita->mctrslyz(slopey);
        hita->mctrintercept(*(m_tr->firstState()));
        m_hits->push_back(hita);
        ++hitidnr;
        float mid_y = (pix_y_out - pix_y_in) / (pix_x_out - pix_x_in) *
                          (ref_x - pix_x_in) +
                      pix_y_in;
        float mid_x = (pix_x_out - pix_x_in) / (pix_y_out - pix_y_in) *
                          (ref_y - pix_y_in) +
                      pix_x_in;

        TbHit *hitb = new TbHit;
        hitb->setId(id);
        hitb->id_nr(hitidnr);
        hitb->setCol((int)mid_y);
        hitb->setRow((int)mid_x);
        hitb->setAdc(l3 * adc);
        hitb->mcx(global_intercept.X());
        hitb->mcy(global_intercept.Y());
        hitb->mctrslxz(slopex);
        hitb->mctrslyz(slopey);
        hitb->mctrintercept(*(m_tr->firstState()));
        m_hits->push_back(hitb);
        ++hitidnr;
        TbHit *hitc = new TbHit;
        hitc->setId(id);
        hitc->id_nr(hitidnr);
        hitc->setCol((int)pix_y_out);
        hitc->setRow((int)pix_x_out);
        hitc->setAdc(l2 * adc);
        hitc->mcx(global_intercept.X());
        hitc->mcy(global_intercept.Y());
        hitc->mctrslxz(slopex);
        hitc->mctrslyz(slopey);
        hitc->mctrintercept(*(m_tr->firstState()));
        m_hits->push_back(hitc);
        ++hitidnr;
      }
    }
    delete m_tr;
  }
  return true;
}

//=============================================================================
// End of Event
//=============================================================================
bool TbToyData::end_event() {
  m_hits->clear();
  return true;
}

//=============================================================================
/// Finalize
//=============================================================================
bool TbToyData::finalize() {
  std::cout << "TbToyData: finalize() " << std::endl;

  return true;
}
