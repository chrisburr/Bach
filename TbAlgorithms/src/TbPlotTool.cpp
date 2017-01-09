
#include "TbPlotTool.h"
#include "DD4hep/Factories.h"

/** @file TbPlotTool.cpp
 *
 *  Implementation of class : TbPlotTool
 *
 */

//=============================================================================
/// Standard constructor
//=============================================================================
TbPlotTool::TbPlotTool(const std::string &name) { tbroot = new TbROOT; }

bool TbPlotTool::configuration() {
  Const_S("OutputFile", "out/Histograms.root");
  return true;
}

//=============================================================================
/// Initialization
//=============================================================================
bool TbPlotTool::initialize(DD4hep::Geometry::LCDD &lcdd, AlgVec algos) {
  m_algos = algos;
  tral = new TbTrackAlgorithms("TrackAlgos");
  m_geomSvc = dynamic_cast<TbGeometrySvc *>(find(algos, "TbGeometrySvc"));
  tral->setGeom(m_geomSvc);
  // Define Histograms

  tbroot->File(Const_S("OutputFile").c_str());
  m_patternrec =
      dynamic_cast<TbPatternRecognition *>(find(algos, "TbPatternRecognition"));

  tbroot->AddHisto1D("track_slope_xz", "track_slope_xz", 20, -0.0002, 0.0002);
  tbroot->AddHisto1D("track_slope_yz", "track_slope_yz", 20, -0.0002, 0.0002);
  tbroot->AddHisto1D("track_slopemc_slope_xz", "track_slopemc_slope_xz", 75, 0.,
                     1.5);
  tbroot->AddHisto1D("track_slopemc_slope_yz", "track_slopemc_slope_yz", 75, 0.,
                     1.5);
  tbroot->AddHisto1D("track_chi2ndof", "track_chi2ndof", 50, 0., 5.);
  tbroot->AddHisto2D("track_firststate", "track_firststate", 2560, -7.04, 7.04,
                     2560, -7.04, 7.04);

  for (std::map<std::string, TbModule *>::iterator itr =
           m_geomSvc->Modules.begin();
       itr != m_geomSvc->Modules.end(); ++itr) {
    tbroot->AddHisto2D(("HitMap" + (*itr).first).c_str(),
                       ("HitMap" + (*itr).first).c_str(), 256, -7.04, 7.04, 256,
                       -7.04, 7.04);
    tbroot->AddHisto3D(("Residuals2D_x_" + (*itr).first).c_str(),
                       ("Residuals2D_x_" + (*itr).first).c_str(), 32, -7.04,
                       7.04, 32, -7.04, 7.04);
    tbroot->AddHisto3D(("Residuals2D_y_" + (*itr).first).c_str(),
                       ("Residuals2D_y_" + (*itr).first).c_str(), 32, -7.04,
                       7.04, 32, -7.04, 7.04);
    tbroot->AddHisto3D(("Residuals2D_x_afteralign" + (*itr).first).c_str(),
                       ("Residuals2D_x_afteralign" + (*itr).first).c_str(), 32,
                       -7.04, 7.04, 32, -7.04, 7.04);
    tbroot->AddHisto3D(("Residuals2D_y_afteralign" + (*itr).first).c_str(),
                       ("Residuals2D_y_afteralign" + (*itr).first).c_str(), 32,
                       -7.04, 7.04, 32, -7.04, 7.04);
    tbroot->AddHisto1D(("Residuals_x_" + (*itr).first).c_str(),
                       ("Residuals_x_" + (*itr).first).c_str(), 50, -0.1, 0.1);
    tbroot->AddHisto1D(("Residuals_y_" + (*itr).first).c_str(),
                       ("Residuals_y_" + (*itr).first).c_str(), 50, -0.1, 0.1);
    tbroot->AddHisto1D(("Residuals_x_afteralign" + (*itr).first).c_str(),
                       ("Residuals_x_afteralign" + (*itr).first).c_str(), 50,
                       -0.1, 0.1);
    tbroot->AddHisto1D(("Residuals_y_afteralign" + (*itr).first).c_str(),
                       ("Residuals_y_afteralign" + (*itr).first).c_str(), 50,
                       -0.1, 0.1);
    tbroot->AddHisto1D(("Clustersize_" + (*itr).first).c_str(),
                       ("Clustersize_" + (*itr).first).c_str(), 10, 0., 10.);
    tbroot->AddHisto1D(("Cluster_real-recox_" + (*itr).first).c_str(),
                       ("Cluster_real-recox_" + (*itr).first).c_str(), 120,
                       -1.2, 1.2);
    tbroot->AddHisto1D(("Cluster_real-recoy_" + (*itr).first).c_str(),
                       ("Cluster_real-recoy_" + (*itr).first).c_str(), 120,
                       -1.2, 1.2);
    tbroot->AddHisto1D(("Cluster_real-recox_1hit" + (*itr).first).c_str(),
                       ("Cluster_real-recox_1hit" + (*itr).first).c_str(), 120,
                       -1.2, 1.2);
    tbroot->AddHisto1D(("Cluster_real-recoy_1hit" + (*itr).first).c_str(),
                       ("Cluster_real-recoy_1hit" + (*itr).first).c_str(), 120,
                       -1.2, 1.2);
    tbroot->AddHisto1D(("Cluster_real-recox_2hit" + (*itr).first).c_str(),
                       ("Cluster_real-recox_2hit" + (*itr).first).c_str(), 120,
                       -1.2, 1.2);
    tbroot->AddHisto1D(("Cluster_real-recoy_2hit" + (*itr).first).c_str(),
                       ("Cluster_real-recoy_2hit" + (*itr).first).c_str(), 120,
                       -1.2, 1.2);
    tbroot->AddHisto1D(("Cluster_real-recox_3hit" + (*itr).first).c_str(),
                       ("Cluster_real-recox_3hit" + (*itr).first).c_str(), 120,
                       -1.2, 1.2);
    tbroot->AddHisto1D(("Cluster_real-recoy_3hit" + (*itr).first).c_str(),
                       ("Cluster_real-recoy_3hit" + (*itr).first).c_str(), 120,
                       -1.2, 1.2);
  }

  return true;
}

//=============================================================================
/// Main execution
//=============================================================================
bool TbPlotTool::execute(AlgVec algos) {
  // loop over tracks
  TbTracks *tracks = m_patternrec->Tracks();

  TbTracks::const_iterator it;
  for (it = tracks->begin(); it != tracks->end(); ++it) {
    tbroot->Histo1D("track_slope_xz")->Fill((*it)->slopeXZ());
    tbroot->Histo1D("track_slope_yz")->Fill((*it)->slopeYZ());
    tbroot->Histo1D("track_chi2ndof")->Fill((*it)->chi2() / (*it)->ndof());
    tbroot->Histo2D("track_firststate")
        ->Fill((*it)->firstState()->X(), (*it)->firstState()->Y());
    double true_slopexz =
        (*it)->Clusters()->at(0)->HitContainer()->at(0)->mctrslxz();
    double true_slopeyz =
        (*it)->Clusters()->at(0)->HitContainer()->at(0)->mctrslyz();

    tbroot->Histo1D("track_slopemc_slope_xz")
        ->Fill((*it)->slopeXZ() / true_slopexz);
    tbroot->Histo1D("track_slopemc_slope_yz")
        ->Fill((*it)->slopeYZ() / true_slopeyz);

    // loop over clusters
    TbClusters *clusters = (*it)->Clusters();

    TbClusters::const_iterator ic;
    for (ic = clusters->begin(); ic != clusters->end(); ++ic) {
      XYZPoint intercept = tral->getInterceptGlobal((*it), (*ic)->id());
      double res_x = (*ic)->GlobalPos().X() - intercept.X();
      double res_y = (*ic)->GlobalPos().Y() - intercept.Y();
      tbroot->Histo3D("Residuals2D_x_" + (*ic)->id())
          ->Fill((*ic)->GlobalPos().X(), (*ic)->GlobalPos().Y(), res_x);
      tbroot->Histo3D("Residuals2D_y_" + (*ic)->id())
          ->Fill((*ic)->GlobalPos().X(), (*ic)->GlobalPos().Y(), res_y);

      tbroot->Histo1D("Residuals_x_" + (*ic)->id())->Fill(res_x);
      tbroot->Histo1D("Residuals_y_" + (*ic)->id())->Fill(res_y);
      tbroot->Histo1D("Clustersize_" + (*ic)->id())
          ->Fill((*ic)->HitContainer()->size());

      double real_reco_x =
          ((*ic)->HitContainer()->at(0)->mcx() - (*ic)->GlobalPos().X()) /
          0.055;
      double real_reco_y =
          ((*ic)->HitContainer()->at(0)->mcy() - (*ic)->GlobalPos().Y()) /
          0.055;
      tbroot->Histo1D("Cluster_real-recox_" + (*ic)->id())->Fill(real_reco_x);
      tbroot->Histo1D("Cluster_real-recoy_" + (*ic)->id())->Fill(real_reco_y);

      if ((*ic)->XCount() == 1)
        tbroot->Histo1D("Cluster_real-recox_1hit" + (*ic)->id())
            ->Fill(real_reco_x);
      else if ((*ic)->XCount() == 2)
        tbroot->Histo1D("Cluster_real-recox_2hit" + (*ic)->id())
            ->Fill(real_reco_x);
      else if ((*ic)->XCount() == 3)
        tbroot->Histo1D("Cluster_real-recox_3hit" + (*ic)->id())
            ->Fill(real_reco_x);

      if ((*ic)->YCount() == 1)
        tbroot->Histo1D("Cluster_real-recoy_1hit" + (*ic)->id())
            ->Fill(real_reco_y);
      else if ((*ic)->YCount() == 2)
        tbroot->Histo1D("Cluster_real-recoy_2hit" + (*ic)->id())
            ->Fill(real_reco_y);
      else if ((*ic)->YCount() == 3)
        tbroot->Histo1D("Cluster_real-recoy_3hit" + (*ic)->id())
            ->Fill(real_reco_y);
    }
  }

  return true;
}
//=============================================================================
// End of Event
//=============================================================================
bool TbPlotTool::end_event() { return true; }
//=============================================================================
/// Finalize
//=============================================================================
bool TbPlotTool::finalize() {
  if (TbAlignment *align =
          dynamic_cast<TbAlignment *>(find(m_algos, "TbAlignment"))) {
    tral->setGeom(align->GetGeom());
    // loop over tracks
    TbTracks *tracks = align->GetTracks();

    TbTracks::const_iterator it;
    for (it = tracks->begin(); it != tracks->end(); ++it) {
      // loop over clusters
      TbClusters *clusters = (*it)->Clusters();
      TbClusters::const_iterator ic;
      for (ic = clusters->begin(); ic != clusters->end(); ++ic) {
        XYZPoint intercept = tral->getInterceptGlobal((*it), (*ic)->id());
        double res_x = (*ic)->GlobalPos().X() - intercept.X();
        double res_y = (*ic)->GlobalPos().Y() - intercept.Y();
        tbroot->Histo1D("Residuals_x_afteralign" + (*ic)->id())->Fill(res_x);
        tbroot->Histo1D("Residuals_y_afteralign" + (*ic)->id())->Fill(res_y);
        tbroot->Histo3D("Residuals2D_x_afteralign" + (*ic)->id())
            ->Fill((*ic)->GlobalPos().X(), (*ic)->GlobalPos().Y(), res_x);
        //	tbroot->Histo3D("Residuals2D_y_aftrtalign"+(*ic)->id())->Fill((*ic)->GlobalPos().X(),(*ic)->GlobalPos().Y(),res_y);
      }
    }
  }
  std::cout << "Writing root file " << Const_S("OutputFile") << std::endl;
  tbroot->File()->Write();
  tbroot->File()->Close();
  delete tral;
  return true;
}
