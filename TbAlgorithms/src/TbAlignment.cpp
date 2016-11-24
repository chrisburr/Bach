
#include "TbAlignment.h"


/** @file TbAlignment.cpp
 *
 *  Implementation of class : TbAlignment
 *
 */


//=============================================================================
/// Standard constructor
//=============================================================================
TbAlignment::TbAlignment(const std::string& name)
  {

}
TbAlignment::~TbAlignment(){
  delete m_trackcontainer;
  delete m_millepede;
  delete m_patternrec;
  delete m_modulestoalign;
  delete tral;
}
bool TbAlignment::configuration(){

  Const_I("Iterations", 1);
  Const_S("FixedModule", "D09-W0108");
  return true;

}

//=============================================================================
/// Initialization
//=============================================================================
bool TbAlignment::initialize(AlgVec  algos) {

  tral = new TbTrackAlgorithms("TrackAlgos");
  m_geomSvc = dynamic_cast<TbGeometrySvc*> (find(algos,"TbGeometrySvc"));

  m_patternrec = dynamic_cast<TbPatternRecognition*> (find(algos,"TbPatternRecognition"));
  m_trackcontainer = new TbTracks;
  m_millepede = new Millepede;
  m_modulestoalign = new std::vector<TbModule*>;

  int nDetToAlign = 0;
  std::vector<double> z_sort;
  for (std::map<std::string,TbModule*>::iterator itr = m_geomSvc->Modules.begin();
       itr!= m_geomSvc->Modules.end(); ++itr){

      z_sort.push_back((*itr).second->Z());
  }
  sort(z_sort.begin(),z_sort.end());

  m_millepede->aligndut = false;
  int Modules[8] = {1,1,1,1,1,1,1,1};

  int counter=0;
  for (int i=0; i<z_sort.size(); ++i){

    if (Modules[i]!=0){
      for (std::map<std::string,TbModule*>::iterator itr2 = m_geomSvc->Modules.begin();
       itr2!= m_geomSvc->Modules.end(); ++itr2){
	if ((*itr2).second->Z() == z_sort[i]) {

	  (*itr2).second->Nr(nDetToAlign);

	  m_modulestoalign->push_back((*itr2).second);
	  ++nDetToAlign;
	}
	TelescopeMap.insert(std::make_pair(counter,z_sort[i]));


      }
    }
    else {
      TelescopeMap.insert(std::make_pair(counter,-1));
    }

    ++counter;
    }
  int fixed = m_geomSvc->Modules[Const_S("FixedModule")]->Nr();

  m_millepede->m_fixed = fixed;

  return true;

}

//=============================================================================
/// Main execution
//=============================================================================
bool TbAlignment::execute(AlgVec algos) {
  //Collect Tracks
  TbTracks *tracks = m_patternrec->Tracks();

  TbTracks::const_iterator it;
  for (it = tracks->begin(); it != tracks->end(); ++it) {

    m_trackcontainer->push_back( (*it) );

				}
  return true;

}

//=============================================================================
// End of Event
//=============================================================================
bool TbAlignment::end_event(){

  return true;
}

//=============================================================================
/// Finalize
//=============================================================================
bool TbAlignment::finalize() {

  std::cout << "Run Millepede!"<<std::endl;

  //Define Constraints, DOF, Sigma
  //int Cons[9]={0,0,0,0,0,0,0,0,0};
  int Cons[9]={0,1,0,1,0,0,1,1,1};  //activate constraint equations
  bool DOF[6] = {1,1,0,0,0,1};      //DoF to align
  double Sigm[6] = {0.1,0.1,0.1,0.08,0.08,0.08};  //Sigma for each DoF

  m_millepede->m_iteration = true;






  int nglo = m_modulestoalign->size();



  int nloc = 4;
  double startfact = 100;
  int nstd = 0;
  double res_cut = 0.06;
  double res_cut_init = 0.3;


  m_millepede->aligndut = false;
  m_millepede->dut = -1;

  for (int i = 0; i <  Const_I("Iterations"); ++i){
    std::cout << "Alignment ---> Iteration: " << i << std::endl;
    std::cout << "Number of grabbed tracks: " << m_trackcontainer->size() << std::endl;
    int n_fits = m_trackcontainer->size();
    m_millepede->InitMille(DOF,Sigm, nglo, nloc, startfact, nstd , res_cut, res_cut_init, n_fits);

    int Nstations  = nglo;   // Number of stations to be aligned (for VELO)
    int Nparams    = 6*nglo;   // Number of params to be aligned (for VELO)

    //
    // Parameters for constraint equations
    //

    double zmoy   = 0.0;
    double s_zmoy = 0.0;
    int nonzer    = 0;

    int CorrectTelescopeMap[TelescopeMap.size()];

    for (int i = 0; i<TelescopeMap.size(); ++i){
      CorrectTelescopeMap[i] = i;

      if(TelescopeMap[i] > 0 ){
	zmoy+=TelescopeMap[i];
	nonzer++;
      }

    }
    zmoy /= nonzer;
    for (int i = 0; i<TelescopeMap.size(); ++i){
      if(TelescopeMap[i] > 0){
	s_zmoy+=(TelescopeMap[i]-zmoy)*(TelescopeMap[i]-zmoy);
      }

    }
    s_zmoy /= nonzer;

    float m_slopex = 0.0;  //Average slope of all tracks for constraint
    float m_slopey = 0.0;
    float m_alpha = 0.0;   //Average z-rotation off all tracks for constraint

    for (int i=0;i<n_fits; ++i) {

      m_slopex+=m_trackcontainer->at(i)->slopeXZ();
      m_slopey+=m_trackcontainer->at(i)->slopeYZ();

      if ((m_trackcontainer->at(i)->firstState()->Y() > 0 && m_trackcontainer->at(i)->slopeXZ()<0) ||
      	  (m_trackcontainer->at(i)->firstState()->Y() < 0 && m_trackcontainer->at(i)->slopeXZ()>0))
      {
	  m_alpha+=sqrt(m_trackcontainer->at(i)->slopeXZ()*m_trackcontainer->at(i)->slopeXZ()
		       +m_trackcontainer->at(i)->slopeYZ()*m_trackcontainer->at(i)->slopeYZ()
		       /
		       (m_trackcontainer->at(i)->firstState()->X()*m_trackcontainer->at(i)->firstState()->X()
			+m_trackcontainer->at(i)->firstState()->Y()*m_trackcontainer->at(i)->firstState()->Y()));

	  }
	  else
	{m_alpha-=sqrt(m_trackcontainer->at(i)->slopeXZ()*m_trackcontainer->at(i)->slopeXZ()
		       +m_trackcontainer->at(i)->slopeYZ()*m_trackcontainer->at(i)->slopeYZ()
		       /
		       (m_trackcontainer->at(i)->firstState()->X()*m_trackcontainer->at(i)->firstState()->X()
			+m_trackcontainer->at(i)->firstState()->Y()*m_trackcontainer->at(i)->firstState()->Y()));



			}
	}
    m_slopex/=n_fits;
    m_slopey/=n_fits;
    m_alpha/=n_fits;


    std::cout << "m_alpha " << m_alpha << std::endl;

    //
    // Here we define the 9 constraints equations
    // according to the requested geometry
    //

    ftx.resize(Nparams);
    fty.resize(Nparams);
    ftz.resize(Nparams);
    frotx.resize(Nparams);
    froty.resize(Nparams);
    frotz.resize(Nparams);

    shearx.resize(Nparams);
    sheary.resize(Nparams);
    fscaz.resize(Nparams);
    for (int j=0; j<Nparams; j++)
      {
	ftx[j]    = 0.0;
	fty[j]    = 0.0;
	ftz[j]    = 0.0;
	frotx[j]  = 0.0;
	froty[j]  = 0.0;
	frotz[j]  = 0.0;
	shearx[j] = 0.0;
	sheary[j] = 0.0;
	fscaz[j]  = 0.0;

      }

    for (int j=0; j<TelescopeMap.size(); j++)
      {

	double z_station = TelescopeMap[j];



	if (z_station >= 0){


	  ftx[CorrectTelescopeMap[j]]                 = 1.0;
	  fty[CorrectTelescopeMap[j]+Nstations]       = 1.0;}

        ftz[CorrectTelescopeMap[j]+2*Nstations]     = 1.0;
        frotx[CorrectTelescopeMap[j]+3*Nstations]   = 1.0;//(z_station-zmoy)/s_zmoy;
        froty[CorrectTelescopeMap[j]+4*Nstations]   = 1.0;//(z_station-zmoy)/s_zmoy;
        frotz[CorrectTelescopeMap[j]+5*Nstations]   = (z_station-zmoy)/s_zmoy;
        shearx[CorrectTelescopeMap[j]]              = (z_station-zmoy)/s_zmoy;
        sheary[CorrectTelescopeMap[j]+Nstations]    = (z_station-zmoy)/s_zmoy;


        fscaz[CorrectTelescopeMap[j]+2*Nstations]   = (z_station-zmoy)/s_zmoy;



      }
    //  Here we put the constraints information in the basket


    if (Cons[0] && DOF[0])   m_millepede->ConstF(&ftx[0],     0.0);
    if (Cons[1] && DOF[0])   m_millepede->ConstF(&shearx[0],  -(nglo)*(m_slopex));//+5.21074e-05));
    if (Cons[2] && DOF[1])   m_millepede->ConstF(&fty[0],     0.0);
    if (Cons[3] && DOF[1])   m_millepede->ConstF(&sheary[0],  -(nglo)*(m_slopey));//-(nglo)+1.40022e-04));
    if (Cons[4] && DOF[2])   m_millepede->ConstF(&ftz[0],     0.0);
    if (Cons[5] && DOF[2])   m_millepede->ConstF(&fscaz[0],   0.0);
    if (Cons[6] && DOF[3])   m_millepede->ConstF(&frotx[0],   -(m_slopey));
    if (Cons[7] && DOF[4])   m_millepede->ConstF(&froty[0],   -(m_slopex));
    if (Cons[8] && DOF[5])   m_millepede->ConstF(&frotz[0],   -2.5*(m_alpha));//-2.5

    // That's it!

    //Feed Millepede with tracks

    for (int i=0;i<n_fits; ++i) {

      PutTrack2(m_trackcontainer->at(i), nglo, nloc, DOF);

    }
    m_millepede->mis_const->clear();
    m_millepede->mis_const->resize(6*nglo);
    m_millepede->mis_error->clear();
    m_millepede->mis_error->resize(6*nglo);
    m_millepede->mis_pull->clear();
    m_millepede->mis_pull->resize(6*nglo);


    m_millepede->MakeGlobalFit(m_millepede->mis_const,m_millepede->mis_error,m_millepede->mis_pull);

    std::cout << "Global Fit made!" << std::endl;
    // Update Geometry

    int index = 0;
    for (std::vector<TbModule*>::iterator itrm = m_modulestoalign->begin();
	 itrm!= m_modulestoalign->end(); ++itrm){
      double rotx;
      double roty;
      double rotz;
      double transx;
      double transy;
      double transz;

      transx = (*itrm)->dX()+m_millepede->m_par->at(index+0*nglo);
      transy = (*itrm)->dY()+m_millepede->m_par->at(index+1*nglo);
      transz = (*itrm)->dZ()+m_millepede->m_par->at(index+2*nglo);
      rotx   = (*itrm)->dRotX()+m_millepede->m_par->at(index+3*nglo);
      roty   = (*itrm)->dRotY()+m_millepede->m_par->at(index+4*nglo);
      rotz   = (*itrm)->dRotZ()+m_millepede->m_par->at(index+5*nglo);

      (*itrm)->SetAlignment(transx, transy, transz, rotx, roty, rotz);
      ++index;
    }

    //Update Tracks


    for (TbTracks::iterator itt = m_trackcontainer->begin(); itt != m_trackcontainer->end(); ++itt) {
      TbClusters* clusters = (*itt)->Clusters();
      TbClusters::iterator ic;
      for (ic = clusters->begin(); ic != clusters->end(); ++ic) {
	if ((*ic) == 0) continue;
	XYZPoint pLocal = (*ic)->LocalPos();
	XYZPoint pGlobal = m_geomSvc->localToGlobal(pLocal, (*ic)->id() );
	(*ic)->GlobalPos(pGlobal);
      }
      tral->setGeom(m_geomSvc);

      tral->FitTrack( (*itt) );


    }
  }
  m_geomSvc->writeConditionsXML("geom/Geom_output.xml");
  delete m_millepede;

return true;
}


bool  TbAlignment::PutTrack( TbTrack *track, int nglo, int nloc, bool m_DOF[])
{

  //Feed Millepede with track equations
  bool sc;


  int Nmodules = nglo;   // Number of modules to be aligned
  int Nlocal = nloc;
  int Nparams = 6*Nmodules;

  std::vector<double>   derLC;
  std::vector<double>   derGB;
  std::vector<double>   derNonLin;
  std::vector<double>   derNonLin_i;

  derGB.clear();
  derGB.resize(Nparams);          // Vector containing the global derivatives
  derNonLin.clear();
  derNonLin.resize(Nparams);      // Global derivatives non linearly related to residual
  derNonLin_i.clear();
  derNonLin_i.resize(Nparams);    // Global derivatives non linearly related to residual
  derLC.clear();
  derLC.resize(Nlocal);           // Vector containing the local derivatives



  for (int i=0; i<Nparams; i++) {derGB[i]=0.; derNonLin[i]=0.; derNonLin_i[i]=0.;}
  for (int i=0; i<Nlocal; i++)  {derLC[i]=0.;}


  double track_params[2*Nlocal+2];     // Vector containing the track parameters
  for (int i=0; i<2*Nlocal+2; i++) {track_params[i]=0.;}

  int n_station = 0;

  double x_cor  = 0.;
  double y_cor  = 0.;
  double z_cor  = 0.;
  double err_x  = 0.;
  double err_y  = 0.;


  // Now we iterate over each cluster on the track and sort those first
  TbClusters* cluster = track->Clusters();

  TbClusters::iterator clus = track->Clusters()->begin();
  std::vector< std::pair <float, TbCluster* > >* cluster_sort = new std::vector< std::pair <float, TbCluster* > >;

  for(; clus != track->Clusters()->end(); ++clus) {
    cluster_sort->push_back( std::make_pair ((*clus)->GlobalPos().Z() , (*clus) ) );

  }
  std::sort(cluster_sort->begin(), cluster_sort->end() );
   std::vector< std::pair <float, TbCluster* > >::iterator itclus = cluster_sort->begin();


  for(; itclus != cluster_sort->end(); ++itclus) {


    std::string dID = (*itclus).second->id();

    z_cor = (*itclus).second->GlobalPos().Z();
    x_cor = (*itclus).second->GlobalPos().X();
    y_cor = (*itclus).second->GlobalPos().Y();



    n_station = detectoridentifier(dID);

    err_x = 0.004;  //Expected error on x- and y-measurements
    err_y = 0.004;

    m_millepede->ZerLoc(&derGB[0],&derLC[0],&derNonLin[0],&derNonLin_i[0]);

    // LOCOAL 1st derivatives for the X equation

    derLC[0] = 1.0;
    derLC[1] = z_cor;
    derLC[2] = 0.0;
    derLC[3] = 0.0;

    // GLOBAL 1st derivatives (see LHCbnote-2005-101 for definition)

    if (m_DOF[0]) derGB[n_station]             = -1.0;    // dX
    if (m_DOF[1]) derGB[Nmodules+n_station]   =  0.0;    // dY
    if (m_DOF[2]) derGB[2*Nmodules+n_station] =  0.0;    // dZ
    if (m_DOF[3]) derGB[3*Nmodules+n_station] =  0.0;    // d_alpha
    if (m_DOF[4]) derGB[4*Nmodules+n_station] =  0.0;    // d_beta
    if (m_DOF[5]) derGB[5*Nmodules+n_station] =  y_cor;  // d_gamma

    if (m_DOF[0]) derNonLin[n_station]             =  0.0;      // dX
    if (m_DOF[1]) derNonLin[Nmodules+n_station]   =  0.0;      // dY
    if (m_DOF[2]) derNonLin[2*Nmodules+n_station] =  1.0;      // dZ
    if (m_DOF[3]) derNonLin[3*Nmodules+n_station] =  y_cor;    // d_alpha
    if (m_DOF[4]) derNonLin[4*Nmodules+n_station] =  x_cor;    // d_beta
    if (m_DOF[5]) derNonLin[5*Nmodules+n_station] =  0.0;      // d_gamma

    if (m_DOF[0]) derNonLin_i[n_station]             =  0.0;      // dX
    if (m_DOF[1]) derNonLin_i[Nmodules+n_station]   =  0.0;      // dY
    if (m_DOF[2]) derNonLin_i[2*Nmodules+n_station] =  0.0;      // dZ 1.0
    if (m_DOF[3]) derNonLin_i[3*Nmodules+n_station] =  1.0;      // d_alpha 1.0
    if (m_DOF[4]) derNonLin_i[4*Nmodules+n_station] =  1.0;      // d_beta 1.0
    if (m_DOF[5]) derNonLin_i[5*Nmodules+n_station] =  0.0;      // d_gamma
    sc = m_millepede->EquLoc(&derGB[0], &derLC[0], &derNonLin[0], &derNonLin_i[0],
    		     x_cor, err_x); // Store hits parameters
    if (! sc) {break;}
    m_millepede->ZerLoc(&derGB[0],&derLC[0],&derNonLin[0],&derNonLin_i[0]);

    // LOCAL 1st derivatives for the Y equation

    derLC[0] = 0.0;
    derLC[1] = 0.0;
    derLC[2] = 1.0;
    derLC[3] = z_cor;

    // GLOBAL 1st derivatives
    if (m_DOF[0]) derGB[n_station]             =  0.0;        // dX
    if (m_DOF[1]) derGB[Nmodules+n_station]   = -1.0;        // dY
    if (m_DOF[2]) derGB[2*Nmodules+n_station] = 0.0;         // dZ
    if (m_DOF[3]) derGB[3*Nmodules+n_station] = 0.0;         // d_alpha
    if (m_DOF[4]) derGB[4*Nmodules+n_station] = 0.0;         // d_beta
    if (m_DOF[5]) derGB[5*Nmodules+n_station] = -x_cor;      // d_gamma

    if (m_DOF[0]) derNonLin[n_station]             =  0.0;     // dX
    if (m_DOF[1]) derNonLin[Nmodules+n_station]   =  0.0;     // dY
    if (m_DOF[2]) derNonLin[2*Nmodules+n_station] =  1.0;     // dZ
    if (m_DOF[3]) derNonLin[3*Nmodules+n_station] =  y_cor;   // d_alpha
    if (m_DOF[4]) derNonLin[4*Nmodules+n_station] =  x_cor;   // d_beta
    if (m_DOF[5]) derNonLin[5*Nmodules+n_station] =  0.0;     // d_gamma

    if (m_DOF[0]) derNonLin_i[n_station]             =  0.0;      // dX
    if (m_DOF[1]) derNonLin_i[Nmodules+n_station]   =  0.0;      // dY
    if (m_DOF[2]) derNonLin_i[2*Nmodules+n_station] =  0.0;      // dZ 3.0
    if (m_DOF[3]) derNonLin_i[3*Nmodules+n_station] =  1.0;      // d_alpha 3.0
    if (m_DOF[4]) derNonLin_i[4*Nmodules+n_station] =  1.0;      // d_beta 3.0
    if (m_DOF[5]) derNonLin_i[5*Nmodules+n_station] =  0.0;      // d_gamma
    sc = m_millepede->EquLoc(&derGB[0], &derLC[0], &derNonLin[0], &derNonLin_i[0],
    y_cor, err_y); // Store hits parameters

    if (! sc) {break;}


  }


  sc = m_millepede->FitLoc(m_millepede->GetTrackNumber(),track_params,0);


  if(sc) {
    m_millepede->SetTrackNumber(m_millepede->GetTrackNumber()+1);

  }
  delete cluster_sort;

  return true;
}


bool  TbAlignment::PutTrack2( TbTrack *track, int nglo, int nloc, bool m_DOF[])
{

  //Feed Millepede with track equations
  bool sc;


  int Nmodules = nglo;   // Number of modules to be aligned
  int Nlocal = nloc;
  int Nparams = 6*Nmodules;

  std::vector<double>   derLC;
  std::vector<double>   derGB;
  std::vector<double>   derNonLin;
  std::vector<double>   derNonLin_i;

  derGB.clear();
  derGB.resize(Nparams);          // Vector containing the global derivatives
  derNonLin.clear();
  derNonLin.resize(Nparams);      // Global derivatives non linearly related to residual
  derNonLin_i.clear();
  derNonLin_i.resize(Nparams);    // Global derivatives non linearly related to residual
  derLC.clear();
  derLC.resize(Nlocal);           // Vector containing the local derivatives



  for (int i=0; i<Nparams; i++) {derGB[i]=0.; derNonLin[i]=0.; derNonLin_i[i]=0.;}
  for (int i=0; i<Nlocal; i++)  {derLC[i]=0.;}


  double track_params[2*Nlocal+2];     // Vector containing the track parameters
  for (int i=0; i<2*Nlocal+2; i++) {track_params[i]=0.;}

  int n_station = 0;

  double x_cor  = 0.;
  double y_cor  = 0.;
  double z_cor  = 0.;
  double err_x  = 0.;
  double err_y  = 0.;
  double err_z  = 0.;

  // Now we iterate over each cluster on the track and sort those first
  TbClusters* cluster = track->Clusters();

  TbClusters::iterator clus = track->Clusters()->begin();
  std::vector< std::pair <float, TbCluster* > >* cluster_sort = new std::vector< std::pair <float, TbCluster* > >;

  for(; clus != track->Clusters()->end(); ++clus) {
    cluster_sort->push_back( std::make_pair ((*clus)->GlobalPos().Z() , (*clus) ) );

  }
  std::sort(cluster_sort->begin(), cluster_sort->end() );
   std::vector< std::pair <float, TbCluster* > >::iterator itclus = cluster_sort->begin();


  for(; itclus != cluster_sort->end(); ++itclus) {


    std::string dID = (*itclus).second->id();
    float z_mod =  m_geomSvc->Modules[dID]->Z();

    z_cor = (*itclus).second->GlobalPos().Z();
    x_cor = (*itclus).second->GlobalPos().X();
    y_cor = (*itclus).second->GlobalPos().Y();

    float z_loc = z_mod-z_cor;
    //std::cout << z_loc << std::endl;
    n_station = detectoridentifier(dID);

    err_x = 0.004;  //Expected error on x- and y-measurements
    err_y = 0.004;
    err_z = 0.004;
    m_millepede->ZerLoc(&derGB[0],&derLC[0],&derNonLin[0],&derNonLin_i[0]);

    // LOCOAL 1st derivatives for the X equation

    derLC[0] = 1.0;
    derLC[1] = z_cor;
    derLC[2] = 0.0;
    derLC[3] = 0.0;
    //derLC[4] = 0.0;
    //derLC[5] = 0.0;
    // GLOBAL 1st derivatives (see LHCbnote-2005-101 for definition)

    if (m_DOF[0]) derGB[n_station]             = -1.0;    // dX
    if (m_DOF[1]) derGB[Nmodules+n_station]   =  0.0;    // dY
    if (m_DOF[2]) derGB[2*Nmodules+n_station] =  0.0;    // dZ
    if (m_DOF[3]) derGB[3*Nmodules+n_station] =  0.0;    // d_alpha
    if (m_DOF[4]) derGB[4*Nmodules+n_station] =  z_loc;    // d_beta
    if (m_DOF[5]) derGB[5*Nmodules+n_station] =  y_cor;  // d_gamma

    if (m_DOF[0]) derNonLin[n_station]            =  0.0;      // dX
    if (m_DOF[1]) derNonLin[Nmodules+n_station]   =  -5*Nmodules+n_station;      // dY
    if (m_DOF[2]) derNonLin[2*Nmodules+n_station] =  -4*Nmodules+n_station;      // dZ
    if (m_DOF[3]) derNonLin[3*Nmodules+n_station] =  0.0;      // d_alpha
    if (m_DOF[4]) derNonLin[4*Nmodules+n_station] =  -2*Nmodules+n_station;      // d_beta
    if (m_DOF[5]) derNonLin[5*Nmodules+n_station] =  -Nmodules+n_station;     // d_gamma

    if (m_DOF[0]) derNonLin_i[n_station]            =  0.0;      // dX
    if (m_DOF[1]) derNonLin_i[Nmodules+n_station]   =  0.0;      // dY
    if (m_DOF[2]) derNonLin_i[2*Nmodules+n_station] =  0.0;      // dZ 1.0
    if (m_DOF[3]) derNonLin_i[3*Nmodules+n_station] =  0.0;      // d_alpha 1.0
    if (m_DOF[4]) derNonLin_i[4*Nmodules+n_station] =  0.0;      // d_beta 1.0
    if (m_DOF[5]) derNonLin_i[5*Nmodules+n_station] =  0.0;      // d_gamma
    sc = m_millepede->EquLoc(&derGB[0], &derLC[0], &derNonLin[0], &derNonLin_i[0],
    		     x_cor, err_x); // Store hits parameters
    if (! sc) {break;}
    m_millepede->ZerLoc(&derGB[0],&derLC[0],&derNonLin[0],&derNonLin_i[0]);

    // LOCAL 1st derivatives for the Y equation

    derLC[0] = 0.0;
    derLC[1] = 0.0;
    derLC[2] = 1.0;
    derLC[3] = z_cor;
    //derLC[4] = 0.0;
    //derLC[5] = 0.0;
    // GLOBAL 1st derivatives
    if (m_DOF[0]) derGB[n_station]             =  0.0;        // dX
    if (m_DOF[1]) derGB[Nmodules+n_station]   = -1.0;        // dY
    if (m_DOF[2]) derGB[2*Nmodules+n_station] = 0.0;         // dZ
    if (m_DOF[3]) derGB[3*Nmodules+n_station] = z_loc;         // d_alpha
    if (m_DOF[4]) derGB[4*Nmodules+n_station] = 0.0;         // d_beta
    if (m_DOF[5]) derGB[5*Nmodules+n_station] = -x_cor;      // d_gamma

    if (m_DOF[0]) derNonLin[n_station]            =  5*Nmodules+n_station;      // dX
    if (m_DOF[1]) derNonLin[Nmodules+n_station]   =  0.0;      // dY
    if (m_DOF[2]) derNonLin[2*Nmodules+n_station] =  -3*Nmodules+n_station;      // dZ
    if (m_DOF[3]) derNonLin[3*Nmodules+n_station] =  -2*Nmodules+n_station;      // d_alpha
    if (m_DOF[4]) derNonLin[4*Nmodules+n_station] =  0.0;      // d_beta
    if (m_DOF[5]) derNonLin[5*Nmodules+n_station] =  n_station;     // d_gamma

    if (m_DOF[0]) derNonLin_i[n_station]            =  0.0;      // dX
    if (m_DOF[1]) derNonLin_i[Nmodules+n_station]   =  0.0;      // dY
    if (m_DOF[2]) derNonLin_i[2*Nmodules+n_station] =  0.0;      // dZ 1.0
    if (m_DOF[3]) derNonLin_i[3*Nmodules+n_station] =  0.0;      // d_alpha 1.0
    if (m_DOF[4]) derNonLin_i[4*Nmodules+n_station] =  0.0;      // d_beta 1.0
    if (m_DOF[5]) derNonLin_i[5*Nmodules+n_station] =  0.0;      // d_gamma
    sc = m_millepede->EquLoc(&derGB[0], &derLC[0], &derNonLin[0], &derNonLin_i[0],
    y_cor, err_y); // Store hits parameters

    if (! sc) {break;}
    m_millepede->ZerLoc(&derGB[0],&derLC[0],&derNonLin[0],&derNonLin_i[0]);

    // LOCAL 1st derivatives for the Z equation

    /*    derLC[0] = 0.0;
    derLC[1] = 0.0;
    derLC[2] = 0.0;
    derLC[3] = 0.0;
    derLC[4] = 0.0;
    derLC[5] = z_cor;
    // GLOBAL 1st derivatives
    if (m_DOF[0]) derGB[n_station]             =  0.0;        // dX
    if (m_DOF[1]) derGB[Nmodules+n_station]   = 0.0;        // dY
    if (m_DOF[2]) derGB[2*Nmodules+n_station] = -1.0;         // dZ
    if (m_DOF[3]) derGB[3*Nmodules+n_station] = -y_cor;         // d_alpha
    if (m_DOF[4]) derGB[4*Nmodules+n_station] = -x_cor;         // d_beta
    if (m_DOF[5]) derGB[5*Nmodules+n_station] = 0.0;      // d_gamma

    if (m_DOF[0]) derNonLin[n_station]            =  4*Nmodules+n_station;      // dX
    if (m_DOF[1]) derNonLin[Nmodules+n_station]   =  3*Nmodules+n_station;      // dY
    if (m_DOF[2]) derNonLin[2*Nmodules+n_station] =  0.0;      // dZ
    if (m_DOF[3]) derNonLin[3*Nmodules+n_station] =  Nmodules+n_station;      // d_alpha
    if (m_DOF[4]) derNonLin[4*Nmodules+n_station] =  n_station;      // d_beta
    if (m_DOF[5]) derNonLin[5*Nmodules+n_station] =  0.0;     // d_gamma

    if (m_DOF[0]) derNonLin_i[n_station]            =  0.0;      // dX
    if (m_DOF[1]) derNonLin_i[Nmodules+n_station]   =  0.0;      // dY
    if (m_DOF[2]) derNonLin_i[2*Nmodules+n_station] =  0.0;      // dZ 1.0
    if (m_DOF[3]) derNonLin_i[3*Nmodules+n_station] =  0.0;      // d_alpha 1.0
    if (m_DOF[4]) derNonLin_i[4*Nmodules+n_station] =  0.0;      // d_beta 1.0
    if (m_DOF[5]) derNonLin_i[5*Nmodules+n_station] =  0.0;      // d_gamma
    sc = m_millepede->EquLoc(&derGB[0], &derLC[0], &derNonLin[0], &derNonLin_i[0],
    z_cor, err_z); // Store hits parameters

    if (! sc) {break;}
    */
  }


  sc = m_millepede->FitLoc(m_millepede->GetTrackNumber(),track_params,0);


  if(sc) {
    m_millepede->SetTrackNumber(m_millepede->GetTrackNumber()+1);

  }
  delete cluster_sort;

  return true;
}
