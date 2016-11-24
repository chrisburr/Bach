
#include "TbPatternRecognition.h"

#include "TMath.h"

/** @file TbPatternRecognition.cpp
 *
 *  Implementation of class : TbPatternRecognition
 *
 *  Class to combine cluster to a track and perform trackfit
 *
 */


//=============================================================================
/// Standard constructor
//=============================================================================
TbPatternRecognition::TbPatternRecognition(const std::string& name)
  {

}
TbPatternRecognition::~TbPatternRecognition()
  {
  delete m_tracks;
  }

bool TbPatternRecognition::configuration(){
  Const_S("ReferenceModule",m_refmod);
  Const_D("Distance"       ,m_dis   );
  Const_D("Chi2ndof-cut"   ,m_chi2cut);
  return true;

}

//=============================================================================
/// Initialization
//=============================================================================
bool TbPatternRecognition::initialize(AlgVec  algos) {
  m_tral = new TbTrackAlgorithms("TrackAlgos2");
  TbGeometrySvc *m_geom = dynamic_cast<TbGeometrySvc*> (find(algos,"TbGeometrySvc"));
  TbClustering *tbc = dynamic_cast<TbClustering*> (find(algos,"TbClustering"));
  tbcluster(tbc);
  m_tracks = new TbTracks;
  return true;

}

//=============================================================================
/// Main execution
//=============================================================================
bool TbPatternRecognition::execute(AlgVec algos) {
  std::map<std::string,TbClusters*> clustermap = tbcluster()->clusters();
  std::string reference_module = Const_S("ReferenceModule");  // Collect all cluster from reference module and search for close cluster on neighbor modules





  TbClusters::const_iterator itc;
  for (itc = clustermap[reference_module]->begin(); itc != clustermap[reference_module]->end(); ++itc) {

    float ref_x = (*itc)->GlobalPos().X();
    float ref_y = (*itc)->GlobalPos().Y();
    TbTrack* track = new TbTrack;

    track->Clusters()->push_back((*itc));
    for (std::map<std::string,TbClusters*>::iterator itm = clustermap.begin();
	 itm != clustermap.end(); ++ itm){
      if ((*itm).first == reference_module) continue;
      float  max_dist = Const_D("Distance");       //Radius in which we look for associated cluster
      TbCluster * closest_cl = new TbCluster;
      TbClusters::const_iterator itc2;
      for (itc2 = (*itm).second->begin(); itc2 != (*itm).second->end(); ++itc2) {
          float cl_x = (*itc2)->GlobalPos().X();
          float cl_y = (*itc2)->GlobalPos().Y();
          float distance = TMath::Sqrt( TMath::Power(ref_x - cl_x,2.) +  TMath::Power(ref_y - cl_y,2.) );

	  if (distance < max_dist){
	    max_dist = distance;
	    closest_cl = (*itc2);
	  }

      }
      if (closest_cl->GlobalPos().X() != 0 )
      track->Clusters()->push_back(closest_cl);
    }
    //Perform Track Fit
    if( track->Clusters()->size() == 8) {
      m_tral -> FitTrack( track);
      //Chi2/ndof-cut

      if (track->chi2()/track->ndof() < Const_D("Chi2ndof-cut")) {

	m_tracks->push_back(track);
      }
    }
  }

  return true;

}
//=============================================================================
// End of Event
//=============================================================================
bool TbPatternRecognition::end_event(){
  m_tracks->clear();
  return true;
}
//=============================================================================
/// Finalize
//=============================================================================
bool TbPatternRecognition::finalize() {
  delete m_tral;

  return true;
}


