#include <algorithm>
#include <map>

#include "TbGeometrySvc.h"
#include "TbROOT.h"
#include "DD4hep/Factories.h"

#include "TbDecoder.h"

// Local
#include "TbClustering.h"

/** @file TbClustering.cpp
 *
 *  Implementation of class : TbClustering
 *
 *  Class to generate cluster and their position from hitselection (either from
 * data or toy data)
 *
 */

//=============================================================================
/// Standard constructor
//=============================================================================
TbClustering::TbClustering(const std::string &name) : m_nEvents(0) {}
TbClustering::~TbClustering() {}
bool TbClustering::configuration() {
  Const_B("DoToyData", false);
  return true;
}

//=============================================================================
/// Initialization
//=============================================================================
bool TbClustering::initialize(AlgVec algos) {
  TbGeometrySvc *geo =
      dynamic_cast<TbGeometrySvc *>(find(algos, "TbGeometrySvc"));
  geomSvc(geo);

  for (auto elm : m_geomSvc->Modules) {
    m_clusters[elm.path()] = new TbClusters;
  }

  // Read hits from TbToyData
  if (Const_B("DoToyData")) {
    TbToyData *td = dynamic_cast<TbToyData *>(find(algos, "TbToyData"));
    toyData(td);
  }
  // Read hits from Data
  else {
    TbDecoder *dec = dynamic_cast<TbDecoder *>(find(algos, "TbDecoder"));

    decoder(dec);
  }

  return true;
}

//=============================================================================
/// Main execution
//=============================================================================
bool TbClustering::execute(dd4hep::cond::ConditionsSlice &slice, AlgVec algos) {
  TbHits *hits;

  if (!Const_B("DoToyData")) {
    hits = decoder()->getHits();
  }

  else {
    hits = toyData()->getHits();
  }

  if (!hits) {
    std::cout << "No hits in Decoder" << std::endl;
    return false;
  }

  // Get geometry parameter

  const double pitch_x = geomSvc()->Const_D("PitchX");
  const double pitch_y = geomSvc()->Const_D("PitchY");
  const double noofpix_x = geomSvc()->Const_I("NoOfPixelX");
  const double noofpix_y = geomSvc()->Const_I("NoOfPixelY");

  // Loop over hits
  TbHits::const_iterator ith;
  for (ith = hits->begin(); ith != hits->end(); ++ith) {
    TbCluster *cluster = new TbCluster;
    TbHits *hitcontainer = new TbHits;

    if ((*ith)->incluster() == true)
      continue;  // Check if hit is already in cluster
    hitcontainer->push_back((*ith));
    TbHits::const_iterator ith_2;

    // Find adjacent hits
    for (ith_2 = hits->begin(); ith_2 != hits->end(); ++ith_2) {
      if (ith == ith_2) continue;
      if ((*ith)->id() != (*ith_2)->id()) continue;
      if (((*ith)->col() == (*ith_2)->col() + 1 &&
           (*ith)->row() == (*ith_2)->row()) ||
          ((*ith)->col() == (*ith_2)->col() &&
           (*ith)->row() == (*ith_2)->row() + 1) ||
          ((*ith)->col() == (*ith_2)->col() - 1 &&
           (*ith)->row() == (*ith_2)->row()) ||
          ((*ith)->col() == (*ith_2)->col() &&
           (*ith)->row() == (*ith_2)->row() - 1) ||
          ((*ith)->col() == (*ith_2)->col() + 1 &&
           (*ith)->row() == (*ith_2)->row() + 1) ||
          ((*ith)->col() == (*ith_2)->col() - 1 &&
           (*ith)->row() == (*ith_2)->row() + 1) ||
          ((*ith)->col() == (*ith_2)->col() + 1 &&
           (*ith)->row() == (*ith_2)->row() - 1) ||
          ((*ith)->col() == (*ith_2)->col() - 1 &&
           (*ith)->row() == (*ith_2)->row() - 1)

              ) {
        (*ith_2)->incluster(true);

        hitcontainer->push_back((*ith_2));
      }
    }

    // Get local position of cluster (CoG-method)
    TbHits::const_iterator ihit;

    float num_row = 0.;
    float denum = 0.;
    float num_col = 0.;

    int size_x = 0;
    int size_y = 0;
    int temp_row = 0;
    int temp_col = 0;
    for (ihit = hitcontainer->begin(); ihit != hitcontainer->end(); ++ihit) {
      if ((*ihit)->row() != temp_row) {
        size_x += 1;
        temp_row = (*ihit)->row();
      }
      if ((*ihit)->col() != temp_col) {
        size_y += 1;
        temp_col = (*ihit)->col();
      }
      num_row += (*ihit)->adc() * ((*ihit)->row() + 0.5);
      num_col += (*ihit)->adc() * ((*ihit)->col() + 0.5);
      denum += (*ihit)->adc();
    }

    float xLocal = (num_row / denum - geomSvc()->Const_I("NoOfPixelX") / 2.) *
                   geomSvc()->Const_D("PitchX");
    float yLocal = (num_col / denum - geomSvc()->Const_I("NoOfPixelY") / 2.) *
                   geomSvc()->Const_D("PitchY");

    Position pLocal(xLocal, yLocal, 0.);
    Position pGlobal(0., 0., 0.);
    auto it = find_if(m_geomSvc->Modules.begin(), m_geomSvc->Modules.end(),
      [&ith] (const DetElement &e) {
        return e.path() == (*ith)->id();
      });

    Position planePointLocalCoords(0., 0., 0.);
    Position planePointGlobalCoords(0., 0., 0.);

    Alignment alignment = slice.get(*it, dd4hep::align::Keys::alignmentKey);
    alignment.data().localToWorld(pLocal, pGlobal);

    cluster->id((*ith)->id());
    cluster->LocalPos(pLocal);
    cluster->GlobalPos(pGlobal);
    cluster->HitContainer(hitcontainer);
    cluster->XCount(size_x);
    cluster->YCount(size_y);

    m_clusters[(*ith)->id()]->push_back(cluster);
  }

  return true;
}
//=============================================================================
// End of Event
//=============================================================================
bool TbClustering::end_event() {
  for (auto elm : m_geomSvc->Modules) {
    m_clusters[elm.path()]->clear();
  }
  return true;
}

//=============================================================================
/// Finalize
//=============================================================================
bool TbClustering::finalize(dd4hep::cond::ConditionsSlice &slice) { return true; }
