#include <map>
#include <algorithm>

#include "../../TbKernel/src/TbGeometrySvc.h"
#include "../../TbKernel/src/TbROOT.h"

#include "TbDecoder.h"

// Local
#include "TbClustering.h"

/** @file TbClustering.cpp
 *
 *  Implementation of class : TbClustering
 *
 *  Class to generate cluster and their position from hitselection (either from data or toy data)
 *
 */


//=============================================================================
/// Standard constructor
//=============================================================================
TbClustering::TbClustering(const std::string& name) : 
    m_nEvents(0) {

}
TbClustering::~TbClustering() {

}
bool TbClustering::configuration(){
  Const_B("DoToyData",true);
}

//=============================================================================
/// Initialization
//=============================================================================
bool TbClustering::initialize(AlgVec  algos) {
  TbGeometrySvc *geo = dynamic_cast<TbGeometrySvc*> (find(algos,"TbGeometrySvc")); 
  geomSvc(geo);
  
  for (std::map<std::string,TbModule*>::iterator itr = m_geomSvc->Modules.begin();
       itr!= m_geomSvc->Modules.end(); ++itr){
    
    m_clusters[(*itr).first] = new TbClusters;
  }

  //Read hits from TbToyData
  if ( Const_B("DoToyData"))  {
     TbToyData *td = dynamic_cast<TbToyData*> (find(algos,"TbToyData")); 
     toyData(td);
   }
  //Read hits from Data
  else {  
      TbDecoder *dec = dynamic_cast<TbDecoder*> (find(algos,"TbDecoder")); 

      decoder(dec);
  }

  return true;

}

//=============================================================================
/// Main execution
//=============================================================================
bool TbClustering::execute(AlgVec algos) {

  
  TbHits * hits;

  if (!Const_B("DoToyData")) {hits = decoder()->getHits();}

  else {hits = toyData()->getHits();}
  
  if (!hits) {
    std::cout << "No hits in Decoder" << std::endl;
    return false;
  }
  
  //Get geometry parameter

  const double pitch_x    = geomSvc()->Const_D("PitchX");
  const double pitch_y    = geomSvc()->Const_D("PitchY");
  const double noofpix_x = geomSvc()->Const_I("NoOfPixelX");
  const double noofpix_y = geomSvc()->Const_I("NoOfPixelY");  
  
  //Loop over hits
  TbHits::const_iterator ith;
  for (ith = hits->begin(); ith != hits->end(); ++ith) {
    TbCluster *cluster = new TbCluster;
    TbHits *hitcontainer = new TbHits;
    
    
    
    if ( (*ith)->incluster() == true ) continue; //Check if hit is already in cluster
    hitcontainer->push_back((*ith));
    TbHits::const_iterator ith_2;


    //Find adjacent hits
    for (ith_2 = hits->begin(); ith_2 != hits->end(); ++ith_2) {
      if (ith == ith_2) continue;
      if ((*ith)->id() != (*ith_2)->id()) continue;
      if ( ((*ith)->col() == (*ith_2)->col()+1 && (*ith)->row() == (*ith_2)->row() )      ||
         ( (*ith)->col() == (*ith_2)->col()    && (*ith)->row() == (*ith_2)->row()+1 )    ||
         ( (*ith)->col() == (*ith_2)->col()-1  && (*ith)->row() == (*ith_2)->row() )      ||
	 ( (*ith)->col() == (*ith_2)->col()    && (*ith)->row() == (*ith_2)->row()-1 )    ||
         ( (*ith)->col() == (*ith_2)->col()+1  && (*ith)->row() == (*ith_2)->row()+1 )    ||
         ( (*ith)->col() == (*ith_2)->col()-1  && (*ith)->row() == (*ith_2)->row()+1 )    ||
         ( (*ith)->col() == (*ith_2)->col()+1  && (*ith)->row() == (*ith_2)->row()-1 )    ||
         ( (*ith)->col() == (*ith_2)->col()-1  && (*ith)->row() == (*ith_2)->row()-1 )
	   
      ) {

        (*ith_2)->incluster(true);

	hitcontainer->push_back((*ith_2));
        
      }
      }
 
     //Get local position of cluster (CoG-method)
     TbHits::const_iterator ihit;

     float num_row = 0.;
     float denum = 0.;
     float num_col = 0.;
     
     
     for (ihit = hitcontainer->begin(); ihit != hitcontainer->end()  ;++ihit ){
       
       num_row += (*ihit)->adc()*((*ihit)->row()+0.5);
       num_col += (*ihit)->adc()*((*ihit)->col()+0.5);
       denum += (*ihit)->adc();
     }

     
     float xLocal = (num_row/denum - geomSvc()->Const_I("NoOfPixelX")/2.)*geomSvc()->Const_D("PitchX");
     float yLocal = (num_col/denum - geomSvc()->Const_I("NoOfPixelY")/2.)*geomSvc()->Const_D("PitchY");

     
     
     XYZPoint pLocal(xLocal, yLocal, 0.);
     XYZPoint pGlobal = geomSvc()->localToGlobal(pLocal, (*ith)->id());
     
     cluster->id((*ith)->id());
     cluster->LocalPos(pLocal);
     cluster->GlobalPos(pGlobal);
     cluster->HitContainer(hitcontainer);
     

    

    m_clusters[(*ith)->id()]->push_back(cluster);

  } 



return true;

}
//=============================================================================
// End of Event
//=============================================================================
bool TbClustering::end_event(){
  for (std::map<std::string,TbModule*>::iterator itr = geomSvc()->Modules.begin();
       itr!= geomSvc()->Modules.end(); ++itr){
    
    m_clusters[(*itr).first]->clear();
  }
  return true;
}

//=============================================================================
/// Finalize
//=============================================================================
bool TbClustering::finalize() {
  

 
  return true;

}


