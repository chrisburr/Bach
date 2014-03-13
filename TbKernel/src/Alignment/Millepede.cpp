// Include files 
#include <iostream>
#include <iomanip> 
#include <algorithm>
#include <fstream>

// local
#include "Millepede.h"
#include <math.h>
//-----------------------------------------------------------------------------
// Implementation file for class : Millepede
//
// 2012-06-19 : Christoph Hombach
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================

using namespace std;


 Millepede:: Millepede(  ) {




   m_debug = false;

  mis_const = new std::vector<double>;
  mis_error = new std::vector<double>;
  mis_pull = new std::vector<double>;
  
  indst = new std::vector<int>;
  arest = new std::vector<double>;
  arenl = new std::vector<double>;

  storeind = new std::vector<int> ;
  
  storeplace = new std::vector<int> ;
  storeare = new std::vector<double>;
  storenl = new std::vector<double>;
  psigm = new std::vector<double>;
  m_par = new std::vector<double>;
  
  clcmat = new std::vector<std::vector<double> >;
  corrm = new std::vector<std::vector<double> >;
  adercs = new std::vector<std::vector<double> >;
  clmat = new std::vector<std::vector<double> >;  
  cgmat = new std::vector<std::vector<double> >; 
  
  m_glmat = new std::vector<std::vector<double> >;
  m_glmatinv = new std::vector<std::vector<double> >;
  corrv = new std::vector<double>;
  pparm = new std::vector<double>;
  dparm = new std::vector<double>;

  scdiag = new std::vector<double>;
  blvec = new std::vector<double>;
  arhs = new std::vector<double>;
  diag = new std::vector<double>;
  bgvec = new std::vector<double>;

  indgb = new std::vector<int>;
  nlnpa = new std::vector<int>;
  indnz = new std::vector<int>;
  indbk = new std::vector<int>;
  indlc = new std::vector<int>;
  
  scflag = new std::vector<bool>;


}

//=============================================================================
// Destructor
//=============================================================================
 Millepede::~ Millepede() {
  

  delete mis_const;        mis_const = 0;
  delete mis_error;        mis_error  = 0;
  delete mis_pull;         mis_pull = 0;

  delete indst;            indst = 0;
  delete arest;            arest = 0;
  delete arenl;            arenl = 0;
  
  delete storeind;         storeind = 0;
  delete storeplace;       storeplace = 0;
  delete storeare;         storeare = 0;
  delete storenl;          storenl = 0;
  delete psigm;            psigm = 0;
  delete m_par;            m_par = 0;
    
  DeleteVector(clcmat);
  DeleteVector(corrm);          
  DeleteVector(adercs);         
  
  DeleteVector(clmat);
  DeleteVector(cgmat);
  
  DeleteVector( m_glmat);
  DeleteVector( m_glmatinv);
  delete corrv;            corrv = 0;
  delete pparm;            pparm = 0;
  delete dparm;            dparm = 0;

  delete scdiag;           scdiag = 0;
  delete blvec;            blvec = 0;
  delete arhs;             arhs = 0;
  delete diag;             diag = 0;
  delete bgvec;            bgvec = 0; 
  
  delete indgb;            indgb = 0;
  delete nlnpa;            nlnpa = 0;
  delete indnz;            indnz = 0;
  delete indbk;            indbk = 0;
  delete indlc;            indlc = 0;

  delete scflag;           scflag = 0;
  
  

} 

//=============================================================================
void  Millepede::DeleteVector(std::vector<std::vector<double> >* vector)
{
  
  if (vector){
    
    for (int i =0; i < vector->size(); i++){
      (vector->at(0)).clear();           
    }
    delete vector;                     vector =0;
  }
}


  

bool  Millepede::InitMille(bool DOF[], double Sigm[], int nglo
				, int nloc, double startfact, int nstd 
				, double res_cut, double res_cut_init, int n_fits)
{

  
  if(m_debug) cout << "" << endl;
  if(m_debug) cout << "----------------------------------------------------" << endl;
  if(m_debug) cout << "" << endl;
  if(m_debug) cout << "    Entering InitMille" << endl;
  if(m_debug) cout << "" << endl;
  if(m_debug) cout << "-----------------------------------------------------" << endl;
  if(m_debug) cout << "" << endl;
  int count = 0;
 

  
  ncs = 0;
  loctot  = 0;                        // Total number of local fits
  locrej  = 0;                        // Total number of local fits rejected
  cfactref  = 1.0;                    // Reference value for Chi^2/ndof cut

  SetTrackNumber(0);       // Number of local fits (starts at 0)

  m_residual_cut = res_cut;
  m_residual_cut_init = res_cut_init; 
 
  nagb	  = 6*nglo;    // Number of global derivatives
  nalc	  = nloc;       // Number of local derivatives
  nstdev  = nstd;     // Number of StDev for local fit chisquare cut

  m_par->clear();       // Vector containing the alignment constants
  m_par->resize(nagb);

  if(m_debug) cout << "Number of global parameters   : " << nagb << endl;
  if(m_debug) cout << "Number of local parameters    : " << nalc << endl;
  if(m_debug) cout << "Number of standard deviations : " << nstdev << endl;

  if (nagb>mglobl || nalc>mlocal)
 {
    if(m_debug) cout << "Two many parameters !!!!!" << endl;
    return false;
  }

  // All parameters initializations
 
 
  corrm->clear();
  corrm->resize(mglobl);
  clcmat->clear();
  clcmat->resize(mglobl);  
  adercs->clear();
  adercs->resize(mglobl);
  psigm->clear();
  psigm->resize(mglobl);
  cgmat->clear();
  cgmat->resize(mgl); 
  clmat->clear();
  clmat->resize(mlocal);  
  corrv->clear();
  corrv->resize(mglobl);    
  pparm->clear();
  pparm->resize(mglobl); 
  dparm->clear();
  dparm->resize(mglobl);
  scdiag->clear();
  scdiag->resize(mglobl);    
  scflag->clear();
  scflag->resize(mglobl);
  indgb->clear();
  indgb->resize(mglobl); 
  nlnpa->clear();
  nlnpa->resize(mglobl);
  indnz->clear();
  indnz->resize(mglobl);
  indbk->clear();
  indbk->resize(mglobl);
  diag->clear();    
  diag->resize(mgl);
  bgvec->clear(); 
  bgvec->resize(mgl);
  blvec->clear();
  blvec->resize(mlocal);
  indlc->clear();
  indlc->resize(mlocal);
  arhs->clear();
  arhs->resize(mcs);
 
  for (int i=0; i<mglobl; i++)
  {
    (corrm->at(i)).clear();
    (corrm->at(i)).resize(mglobl);
    (clcmat->at(i)).clear();   
    (clcmat->at(i)).resize(mglobl);
    (adercs->at(i)).clear();     
    (adercs->at(i)).resize(mglobl);    
    
  }
 
  
  for (int i=0; i<mgl; i++)
  {
    (cgmat->at(i)).clear();
    (cgmat->at(i)).resize(mgl);
   
  }

  for (int i=0; i<mlocal; i++)
  {
   
    (clmat->at(i)).clear();
    (clmat->at(i)).resize(mlocal);
    
  }

  
  for (int j=0; j<mcs;j++) arhs->at(j) = 0.;   

  // Then we fix all parameters...
  
  for (int j=0; j<nagb; j++)  {
   
    
    ParSig(j,0.0);
  }
  
  // ...and we allow them to move if requested

  for (int i=0; i<6; i++)
  {
    if(m_debug) cout << "GetDOF(" << i << ")= " << DOF[i] << endl;
    
    
                
    if (DOF[i]) {for (int j=i*nglo; j<(i+1)*nglo; j++) ParSig(j,Sigm[i]);}

  }
  
  //  ParSig(0,0.);
  //ParSig(nglo+0,0.);
  //ParSig(2*nglo,0.);
  //ParSig(3*nglo,0.);
  //ParSig(4*nglo,0.);
  //ParSig(5*nglo,0.);
  
  
  for (int i = 0;i<nglo;++i){
    if(aligndut && i != dut){
      //if(i == dut){
      ParSig(i,0.);
      ParSig(nglo+i,0.);
      ParSig(2*nglo+i,0.);
      ParSig(3*nglo+i,0.);
      ParSig(4*nglo+i,0.);
      ParSig(5*nglo+i,0.);
    }
    }
  
  if (dut >= 0 && dut < nglo &&  !aligndut )
  {
    if(m_debug) cout << "You are fixing module " << m_fixed << endl;

    ParSig(dut,0.);
    ParSig(nglo+dut,0.);
    ParSig(2*nglo+dut,0.);
    ParSig(3*nglo+dut,0.);
    ParSig(4*nglo+dut,0.);
    ParSig(5*nglo+dut,0.);
  }
    

  if (m_fixed >= 0 && m_fixed < nglo &&  !aligndut )
  {
    //if(m_debug) 
      cout << "You are fixing module " << m_fixed << endl;

    ParSig(m_fixed,0.);
    ParSig(nglo+m_fixed,0.);
    ParSig(2*nglo+m_fixed,0.);
    ParSig(3*nglo+m_fixed,0.);
    ParSig(4*nglo+m_fixed,0.);
    ParSig(5*nglo+m_fixed,0.);
  }

      //if(m_debug) cout << "Sigm(" << j << ")= " << psigm->at(j) << endl;

  // Activate iterations (if requested)

  itert   = 0;	// By default iterations are turned off
  cfactr  = startfact;
  if (m_iteration) InitUn(startfact);          

  arest->clear();  // Number of stored parameters when doing local fit
  //arest->resize(1000);
  arenl->clear(); // Linear or not
  //arenl->resize(1000);
  indst->clear(); 
  //indst->resize(1000);

  storeind->clear();
  storeare->clear();
  storenl->clear();
  storeplace->clear();

  // Memory allocation for the stores

  if(m_debug) cout << "Store size is " << n_fits*(nagb+nalc+3) << endl;
  
  
  storeind->reserve(2*n_fits*(nagb+nalc+3));
  storeare->reserve(2*n_fits*(nagb+nalc+3));
  storenl->reserve(2*n_fits*(nagb+nalc+3));
  storeplace->reserve(2*n_fits);

  if(m_debug) cout << "" << endl;
  if(m_debug) cout << "----------------------------------------------------" << endl;
  if(m_debug) cout << "" << endl;
  if(m_debug) cout << "    InitMille has been successfully called!" << endl;
  if(m_debug) cout << "" << endl;
  if(m_debug) cout << "-----------------------------------------------------" << endl;
  if(m_debug) cout << "" << endl;
		
  return true;
}


int     Millepede::GetTrackNumber()                      {return m_track_number;}
void    Millepede::SetTrackNumber(int value)             {m_track_number = value;}
/*
-----------------------------------------------------------
  PARGLO: initialization of global parameters
-----------------------------------------------------------

  index    = the index of the global parameter in the 
             result array (equivalent to dparm[]).

  param    = the starting value

-----------------------------------------------------------
*/

bool  Millepede::ParGlo(int index, double param)
{
 if (index<0 || index>=nagb)
   {return false;}
 else
 {pparm->at(index) = param;}

 return true;
}
    
/*
-----------------------------------------------------------
  CONSTF: define a constraint equation in Millepede
-----------------------------------------------------------

  dercs    = the row containing constraint equation 
             derivatives (put into the final matrix)

  rhs      = the lagrange multiplier value (sum of equation)	     

-----------------------------------------------------------
*/

bool  Millepede::ConstF(double dercs[], double rhs)
{  
  if (ncs>=mcs) // mcs is defined in Millepede.h
  {
    if(m_debug)cout << "Too many constraints !!!" << endl;
    return false;
  }
 	
  for (int i=0; i<nagb; i++) {(adercs->at(ncs))[i] = dercs[i];}
  
 	
  arhs->at(ncs) = rhs;
  ncs++ ;
  if(m_debug) cout << "Number of constraints increased to " << ncs << endl;
  return true;
}

/*
-----------------------------------------------------------
  PARSIG: define a constraint for a single global param
          param is 'encouraged' to vary within [-sigma;sigma] 
	  range
-----------------------------------------------------------

  index    = the index of the global parameter in the 
             result array (equivalent to dparm[]).

  sigma	   = value of the constraint (sigma <= 0. will 
             mean that parameter is FIXED !!!) 
 
-----------------------------------------------------------

*/

  
bool  Millepede::ParSig(int index, double sigma)
{


  
  if (index>=nagb) 
    {return false;}
  else
  {psigm->at(index) = sigma;}
  if(m_debug) cout << "--> Parsig " << index << " value " << psigm->at(index) << endl;
  return true;
}

/*
-----------------------------------------------------------
  INITUN: unit for iteration
-----------------------------------------------------------
  
  cutfac is used by Fitloc to define the Chi^2/ndof cut value

  A large cutfac value enables to take a wider range of tracks 
  for first iterations, which might be useful if misalignments
  are large.

  As soon as cutfac differs from 0 iteration are requested.
  cutfac is then reduced, from one iteration to the other,
  and iterations are stopped when it reaches the value 1.

  At least one more iteration is often needed in order to remove
  tracks containing outliers.
  
-----------------------------------------------------------
*/
 
bool  Millepede::InitUn(double cutfac)
{
  cfactr = std::max(1.0, cutfac);
  
  if(m_debug) cout << "Initial cut factor is  " << cfactr << endl;
  itert = 1; // Initializes the iteration process
  return true;
}


/*



/*
-----------------------------------------------------------
  ZERLOC: reset the derivative vectors
-----------------------------------------------------------

  dergb[1..nagb]		= global parameters derivatives
  dernl[1..nagb]		= global parameters 'non-linear' derivatives
  dernl_i[1..nagb]		= array linking 'non-linear' derivatives 
                                  and corresponding local param
  dergb[1..nalc]		= local parameters derivatives
 
-----------------------------------------------------------
*/
 
bool  Millepede::ZerLoc(double dergb[], double derlc[],double dernl[], double dernl_i[]) 
{
  
  for(int i=0; i<nalc; i++) {derlc[i] = 0.0;}
  for(int i=0; i<nagb; i++) {dergb[i] = 0.0;}
  for(int i=0; i<nagb; i++) {dernl[i] = 0.0;}
  for(int i=0; i<nagb; i++) {dernl_i[i] = 0.0;}

  return true;
}
/*
-----------------------------------------------------------
  EQULOC: write ONE equation in the matrices
-----------------------------------------------------------

  dergb[1..nagb]	= global parameters derivatives
  dernl[1..nagb]		= global parameters 'non-linear' derivatives
  dernl_i[1..nagb]		= array linking 'non-linear' derivatives 
                                  and corresponding local param
  derlc[1..nalc] 	= local parameters derivatives
  rmeas  		= measured value
  sigma 		= error on measured value (nothin to do with ParSig!!!)

-----------------------------------------------------------
*/

bool  Millepede::EquLoc(double dergb[], double derlc[], double dernl[], double dernl_i[], 
			     double rmeas, double sigma)
{	
  if (sigma<=0.0) // If parameter is fixed, then no equation
  {
    for (int i=0; i<nalc; i++)
    {
      derlc[i] = 0.0;
    }
    for (int i=0; i<nagb; i++)
    {
      dergb[i] = 0.0;
    }
    if(m_debug)cout << "Didn't use the track" << endl;

    return true;
  }
  
// Serious equation, initialize parameters
  	
  double wght =  1.0/(sigma*sigma);
  int nonzer  =  0;
  int ialc    = -1;
  int iblc    = -1;
  int iagb    = -1;
  int ibgb    = -1;
 
  for (int i=0; i<nalc; i++) // Retrieve local param interesting indices
  {
 
    
    if (derlc[i]!=0.0)
    {
      nonzer++;
      if (ialc == -1) ialc=i;	// first index
      iblc = i;       	     	// last index
    }
  }
	 
  if(m_debug)cout <<"derlc first and last index: " << ialc << " / " << iblc << endl;
	
  for (int i=0; i<nagb; i++)  // Idem for global parameters
  {
    //    cout << "dergb[" << i << "] "<<dergb[i] << endl;   
    if (dergb[i]!=0.0 || dernl[i]!=0.0)
    {
      nonzer++;
      if (iagb == -1) iagb=i;	// first index
      ibgb = i; 	     	// last index
    }
  }

  if(m_debug)cout <<"dergl first and last index: " << iagb << " / " << ibgb << endl;
  
  
  indst->push_back(-1);
  
  arest->push_back(rmeas);
  
  arenl->push_back(0.);

  
  if (ialc != -1)  // Just in case of constrained fit (PV for example)
  {
    for (int i=ialc; i<=iblc; i++)
    {
      if (derlc[i]!=0.0)
      {
        indst->push_back(i);
        arest->push_back(derlc[i]);
        arenl->push_back(0.0);
        derlc[i]   = 0.0;
      }
    }
  }  

  indst->push_back(-1);
  arest->push_back(wght);
  arenl->push_back(0.);

  
  if (iagb != -1)
  {
    for (int i=iagb; i<=ibgb; i++)
    {
      if (dergb[i]!=0.0 || dernl[i]!=0.0)
      {
        indst->push_back(i);
        arest->push_back(dergb[i]);
        
        
        if (dernl[i]!=0.0) 
        {
	  
          arenl->push_back(dernl[i]+(dernl_i[i]+0.5)*nonlin_param);
        }
        else
        {
          arenl->push_back(0.);
        }
        
        dergb[i]   = 0.0;
        dernl[i]   = 0.0;
      }
    }
  }	

  if(m_debug)cout << "Out Equloc --  NST = " << arest->size() << endl;
  return true; 	
}

/*
-----------------------------------------------------------
  FITLOC:  perform local params fit, once all the equations
           have been written by EquLoc
-----------------------------------------------------------

  n            = number of the fit, it is used to store 
                 fit parameters and then retrieve them 
		 for iterations (via STOREIND and STOREARE)

  track_params = contains the fitted track parameters and
                 related errors

  single_fit   = is an option, if it is set to 1, we don't 
                 perform the last loop. It is used to update 
		 the track parameters without modifying global
		 matrices

-----------------------------------------------------------
*/

bool  Millepede::FitLoc(int n, double track_params[], int single_fit)
{
// Few initializations
	
  int i, j, k, ik, ij, ist, nderlc, ndergl, ndf;
  int ja      = -1;
  int jb      = 0;
  int nagbn   = 0;
	
  double rmeas, wght, rms, cutval;

  double summ  = 0.0;
  int    nsum  = 0;
  nst   = 0; 
  nst   = arest->size();
 
  //for (i=0; i < arenl->size() ; ++i) cout << arenl->at(i) << endl;
  // Fill the track store at first pass
  
  
  if (itert < 2 && single_fit != 1)  // Do it only once //hier!! 
  {
    if(m_debug) cout  << "Store equation no: " << n << endl; 
       
    for (i=0; i<nst; i++)    // Store the track parameters
    {
      storeind->push_back(indst->at(i));
      storeare->push_back(arest->at(i));
      storenl->push_back(arenl->at(i));

      //if (arenl->at(i) != 0.) arest->at(i) = 0.0; // Reset global derivatives if non linear and first iteration 
    }

    arenl->clear();
    storeplace->push_back(storeind->size());

    if(m_debug) cout  << "StorePlace size = " << storeplace->at(n) << endl; 
    if(m_debug) cout  << "StoreInd size   = " << storeind->size() << endl; 
  }	

  blvec->clear();
  blvec->resize(nalc);
  
  for (i=0; i<nalc; i++) // reset local params
  {
  
    (clmat->at(i)).clear();
    (clmat->at(i)).resize(nalc);
    
    
  }
  indnz->clear();
  indnz->resize(nagb);
  
  for (i=0; i<nagb; i++) {indnz->at(i) = -1;} // reset mixed params


/*

  LOOPS : HOW DOES IT WORKS ?	

  Now start by reading the informations stored with EquLoc.
  Those informations are in vector INDST and AREST.
  Each -1 in INDST delimits the equation parameters:
  
  First -1  ---> rmeas in AREST 
  Then we have indices of local eq in INDST, and derivatives in AREST
  Second -1 ---> weight in AREST
  Then follows indices and derivatives of global eq.
  ....
  
  We took them and store them into matrices.
  
  As we want ONLY local params, we substract the part of the estimated value
  due to global params. Indeed we could have already an idea of these params,
  with previous alignment constants for example (set with PARGLO). Also if there
  are more than one iteration (FITLOC could be called by FITGLO)

*/

    
//
// FIRST LOOP : local track fit
//
	
  ist = 0;
  indst->push_back(-1);
	
  while (ist <= nst)
  {
    if (indst->at(ist) == -1)
    {
      if (ja == -1)     {ja = ist;}  // First  0 : rmeas
      else if (jb == 0) {jb = ist;}  // Second 0 : weight 
      else                           // Third  0 : end of equation  
      {
	
        rmeas	= arest->at(ja);
        
	wght 	= arest->at(jb);
        if(m_debug) cout  << "rmeas = " << rmeas << endl ;
        if(m_debug) cout  << "wght = " << wght << endl ;
        
        for (i=(jb+1); i<ist; i++)   // Now suppress the global part   
          // (only relevant with iterations)
        {

          
          j = indst->at(i);              // Global param indice
          if(m_debug) cout  << "dparm[" << j << "] = " << dparm->at(j) << endl;        
          if(m_debug) cout  << "Starting misalignment = " << pparm->at(j) << endl;        
          
          
          rmeas -= arest->at(i)*(pparm->at(j)+dparm->at(j));
        }

        if(m_debug) cout  << "rmeas after global stuff removal = " << rmeas << endl ;
				
        for (i=(ja+1); i<jb; i++)    // Finally fill local matrix and vector
        {
          j = indst->at(i);   // Local param indice (the matrix line) 
          blvec->at(j) += wght*rmeas*arest->at(i);  // See note Millepede for precisions b= sum w*d*z
         
          
          if(m_debug) cout  << "blvec[" << j << "] = " << blvec->at(j) << endl ;
					
          for (k=(ja+1); k<=i ; k++) // Symmetric matrix, don't bother k>j coeffs
          {
            ik = indst->at(k);						
            (clmat->at(j))[ik] += wght*arest->at(i)*arest->at(k);
	    
            if(m_debug) cout  << "clmat[" << j << "][" << ik << "] = " << (clmat->at(j))[ik] << endl;
          } 
        }  
        ja = -1;
        jb = 0;
        ist--;
      } // End of "end of equation" operations
    } // End of loop on equation
    ist++;
  } // End of loop on all equations used in the fit


//
// Local params matrix is completed, now invert to solve...
//
	
  nrank = 0;  // Rank is the number of nonzero diagonal elements 
  nrank = SpmInv(clmat, blvec, nalc, scdiag, scflag);
       	
  if(m_debug) cout  << "" << endl;
  if(m_debug) cout  << " __________________________________________________" << endl;
  if(m_debug) cout  << " Printout of local fit  (FITLOC)  with rank= "<< nrank << endl;
  if(m_debug) cout  << " Result of local fit :      (index/parameter/error)" << endl;
  
  for (i=0; i<nalc; i++)
  {
    if(m_debug) cout  << std::setprecision(4) << std::fixed;
    if(m_debug) cout  << std::setw(20) << i << "   /   " << std::setw(10) 
                      << blvec->at(i) << "   /   " << sqrt((clmat->at(i))[i]) << endl;	
  }


// Store the track params, errors

  for (i=0; i<nalc; i++)
  {
    track_params[2*i] = blvec->at(i);
    track_params[2*i+1] = sqrt(fabs((clmat->at(i))[i]));
  }

    
//
// SECOND LOOP : residual calculation
//
  
  ist = 0;
  ja = -1;
  jb = 0;

  while (ist <= nst)
  {
    if (indst->at(ist) == -1)
    {
      if (ja == -1)     {ja = ist;}  // First  0 : rmeas
      else if (jb == 0) {jb = ist;}  // Second 0 : weight 
      else                           // Third  0 : end of equation  
      {	
        rmeas	= arest->at(ja);
        wght 	= arest->at(jb);

        nderlc = jb-ja-1;    // Number of local derivatives involved
        ndergl = ist-jb-1;   // Number of global derivatives involved
        
        // Print all (for debugging purposes)
        
        if(m_debug) cout  << "" << endl;
        if(m_debug) cout  << std::setprecision(4) << std::fixed;
        if(m_debug) cout  << ". equation:  measured value " << std::setw(13) 
                          << rmeas << " +/- " << std::setw(13) << 1.0/sqrt(wght) << endl;
        if(m_debug) cout  << "Number of derivatives (global, local): " 
                          << ndergl << ", " << nderlc << endl;
        if(m_debug) cout  << "Global derivatives are: (index/derivative/parvalue) " << endl;
        
        for (i=(jb+1); i<ist; i++)
        {if(m_debug) cout  << indst->at(i) << " / " << arest->at(i) 
                           << " / " << pparm->at(indst->at(i)) << endl;} 
        if(m_debug) cout  << "Global Non-Lin derivatives are: (index/derivative/parvalue) " << endl;
        
        
        
        
        /*for (i=(jb+1); i<ist; i++)
        {    
          if(m_debug) cout  << indst->at(i) << " / " << arenl->at(i) 
          << " / " << pparm->at(indst->at(i)) << endl;} 
        */
        if(m_debug) cout  << "Local derivatives are: (index/derivative) " << endl;
        
        for (i=(ja+1); i<jb; i++) {if(m_debug) cout  << indst->at(i) << " / " << arest->at(i) << endl;}	  
        
        // Now suppress local and global parts to RMEAS;
        
        for (i=(ja+1); i<jb; i++) // First the local part 
        {
          j = indst->at(i);
          rmeas -= arest->at(i)*blvec->at(j);
        }
        
        for (i=(jb+1); i<ist; i++) // Then the global part
        {
          j = indst->at(i);
          rmeas -= arest->at(i)*(pparm->at(j)+dparm->at(j));
          
        }
        
        // rmeas contains now the residual value
        if(m_debug) cout  << "Residual value : "<< rmeas << endl;
        
        // reject the track if rmeas is too important (outlier)
        if (fabs(rmeas) >= m_residual_cut_init && itert <= 1)  //hier!!  
        {
          if(m_debug) cout  << "Rejecting track due to residual cut in iteration 0-1!!!!!" << endl;
          if (single_fit == 0) locrej++;      
          indst->clear(); // reset stores and go to the next track 
          //indst->resize(1000);
	  arest->clear();	  
	  //arest->resize(1000);
	  return false;
        }
        
        if (fabs(rmeas) >= m_residual_cut && itert > 1)  //hier!!   
        {
          if(m_debug) cout  << "Rejected track due to residual cut in iteration " << itert << "!!!!!" << endl;
          if (single_fit == 0) locrej++;      
          indst->clear(); // reset stores and go to the next track 
          //indst->resize(1000);
	  arest->clear();
	  //arest->resize(1000);
          return false;
        }
        
        summ += wght*rmeas*rmeas ; // total chi^2
        nsum++;                    // number of equations			
        ja = -1;
        jb = 0;
        ist--;
      } // End of "end of equation" operations
    }   // End of loop on equation
    ist++;
  } // End of loop on all equations used in the fit

  ndf = nsum-nrank;	
  rms = 0.0;
  
  if(m_debug) cout  << "Final chi square / degrees of freedom "<< summ << " / " << ndf << endl;
  
  if (ndf > 0) rms = summ/float(ndf);  // Chi^2/dof
  
  if (single_fit == 0) loctot++;	
  
  if (nstdev != 0 && ndf > 0 && single_fit != 1) // Chisquare cut
  {
    cutval = chindl(nstdev, ndf)*cfactr;
    
    if(m_debug) cout  << "Reject if Chisq/Ndf = " << rms << "  >  " << cutval << endl;
 
    if (rms > cutval) // Reject the track if too much...
    {
      if(m_debug) cout << "Rejected track because rms (" << rms << ") larger than " << cutval << " !!!!!" << endl;
      locrej++;      
      indst->clear(); // reset stores and go to the next track 
      //indst->resize(1000);
      arest->clear();
      //arest->resize(1000);
      return false;
    }
  }

  if (single_fit == 1) // Stop here if just updating the track parameters
  {
    indst->clear(); // Reset store for the next track 
    //indst->resize(1000);
    arest->clear();
    //arest->resize(1000);
    return true;
  }
  

  // Store the track number of DOFs (for the final chisquare estimation)

  track_params[2*nalc]   = float(ndf);
  track_params[2*nalc+1] = summ;

  //  
  // THIRD LOOP: local operations are finished, track is accepted 
  // We now update the global parameters (other matrices)
  //
  
  ist = 0;
  ja = -1;
  jb = 0;

  while (ist <= nst)
  {
    if (indst->at(ist) == -1)
    {
      if (ja == -1)     {ja = ist;}  // First  0 : rmeas
      else if (jb == 0) {jb = ist;}  // Second 0 : weight 
      else                           // Third  0 : end of equation  
      {	
        rmeas	= arest->at(ja);
        wght 	= arest->at(jb);
	

        for (i=(jb+1); i<ist; i++) // Now suppress the global part
        {
          j = indst->at(i);   // Global param indice
          rmeas -= arest->at(i)*(pparm->at(j)+dparm->at(j));
        }
        
        // First of all, the global/global terms (exactly like local matrix)
        
        for (i=(jb+1); i<ist; i++)  
        {
          j = indst->at(i);   // Global param indice (the matrix line)          
          
          bgvec->at(j) += wght*rmeas*arest->at(i);  // See note for precisions
          if(m_debug) cout  << "bgvec[" << j << "] = " << bgvec->at(j) << endl ;
          
          for (k=(jb+1); k<ist ; k++)
          {
            ik = indst->at(k);						
            (cgmat->at(j))[ik] += wght*arest->at(i)*arest->at(k);
            if(m_debug) cout  << "cgmat[" << j << "][" << ik << "] = " << (cgmat->at(j))[ik] << endl;
          } 
        }
        
        // Now we have also rectangular matrices containing global/local terms.
        
        for (i=(jb+1); i<ist; i++) 
        {
          j  = indst->at(i);  // Global param indice (the matrix line)
          ik = indnz->at(j);  // Index of index
          
          if (ik == -1)	  // New global variable
          {
            for (k=0; k<nalc; k++) {(clcmat->at(nagbn))[k] = 0.0;} // Initialize the row
            
            indnz->at(j) = nagbn;
            indbk->at(nagbn) = j;
            ik = nagbn;
            nagbn++;
          }
          
          for (k=(ja+1); k<jb ; k++) // Now fill the rectangular matrix
          {
            ij = indst->at(k);						
            (clcmat->at(ik))[ij] += wght*arest->at(i)*arest->at(k);
            if(m_debug) cout << "clcmat[" << ik << "][" << ij << "] = " << (clcmat->at(ik))[ij] << endl ;
          } 
        }
        ja = -1;
        jb = 0;
        ist--;
      } // End of "end of equation" operations
    }   // End of loop on equation
    ist++;
  } // End of loop on all equations used in the fit
	
  // Third loop is finished, now we update the correction matrices (see notes)

  SpAVAt(clmat, clcmat, corrm, nalc, nagbn);
  SpAX(clcmat, blvec, corrv, nalc, nagbn);

  for (i=0; i<nagbn; i++)
  {
    j = indbk->at(i);
    bgvec->at(j) -= corrv->at(i);
		
    for (k=0; k<nagbn; k++)
    {
      ik = indbk->at(k);
      (cgmat->at(j))[ik] -= (corrm->at(i))[k];
    }
  }
  //chombach
  /*cout << storeare->size() << "  " << arest->size() << endl;
  nst = arest->size();
  storeare->clear();
  for (i=0; i<nst; i++)    // Store the track parameters
    {storeare->push_back(arest->at(i));}
  cout << storeare->size() << "  " << arest->size() << endl;
  */
  indst->clear(); // Reset store for the next track 

  arest->clear();
  //arest->resize(1000);

  return true;
}

/*
-----------------------------------------------------------
  SPMINV:  obtain solution of a system of linear equations with symmetric matrix 
 	   and the inverse (using 'singular-value friendly' GAUSS pivot)
-----------------------------------------------------------

	Solve the equation :  V * X = B
	
	V is replaced by inverse matrix and B by X, the solution vector
-----------------------------------------------------------

s. V. Blobel Kap. 3
*/

int  Millepede::SpmInv(std::vector<std::vector<double> >* v,std::vector<double>* b, int n,
                              std::vector<double>* diag,std::vector<bool>*  flag)
{
  
  
  
  if (v->size()==mgl){

    
  
  int k;
  double vkk, *temp;
  double  *r, *c;
  double eps = 0.0000000000001;
  bool *used_param;

  r = new double[n];
  c = new double[n];

  temp = new double[n];
  used_param = new bool[n];

  for (int i=0; i<n; i++)
  {
    r[i] = 0.0;
    c[i] = 0.0;
    flag->at(i) = true;
    used_param[i] = true; 

    for (int j=0; j<=i; j++) {(v->at(j))[i] = (v->at(i))[j];}
  }
  
  // Small loop for matrix equilibration (gives a better conditioning) 

  for (int i=0; i<n; i++)
  {
    for (int j=0; j<n; j++)
    { 
      if (fabs((v->at(i))[j]) >= r[i]) r[i] = fabs((v->at(i))[j]); // Max elemt of row i
      if (fabs((v->at(j))[i]) >= c[i]) c[i] = fabs((v->at(j))[i]); // Max elemt of column i
    }
  }

  for (int i=0; i<n; i++)
  {
    if (0.0 != r[i]) r[i] = 1./r[i]; // Max elemt of row i
    if (0.0 != c[i]) c[i] = 1./c[i]; // Max elemt of column i

    if (eps >= r[i]) r[i] = 0.0; // Max elemt of row i not wihin requested precision
    if (eps >= c[i]) c[i] = 0.0; // Max elemt of column i not wihin requested precision
    
  }

  for (int i=0; i<n; i++) // Equilibrate the V matrix
  {
    for (int j=0; j<n; j++) {(v->at(i))[j] = sqrt(r[i])*(v->at(i))[j]*sqrt(c[j]);}
  }

  nrank = 0;

  // save diagonal elem absolute values 	
  for (int i=0; i<n; i++) 
  {
    diag->at(i) = fabs((v->at(i))[i]);

    if (r[i] == 0. && c[i] == 0.) // This part is empty (non-linear treatment with non constraints)
    {
      flag->at(i) = false;
      used_param[i] = false; 
    }
  }

  for (int i=0; i<n; i++)
  {
    vkk = 0.0;
    k = -1;
    
    for (int j=0; j<n; j++) // First look for the pivot, ie max unused diagonal element 
    {
      if (flag->at(j) && (fabs((v->at(j))[j])>std::max(fabs(vkk),eps)))
      {
       
        
        vkk = (v->at(j))[j];
        k = j;
      }
    }
	     
    if (k >= 0)    // pivot found
    {      
      if(m_debug) cout << "Pivot value :" << vkk << endl;
      nrank++;
      flag->at(k) = false; // This value is used
      vkk = 1.0/vkk;
      (v->at(k))[k] = -vkk; // Replace pivot by its inverse
     
      for (int j=0; j<n; j++)
      {
	for (int jj=0; jj<n; jj++)  
	{
	  if (j != k && jj != k && used_param[j] && used_param[jj] ) // Other elements (!!! do them first as you use old v[k][j]'s !!!)
	  {
	    (v->at(j))[jj] = (v->at(j))[jj] - vkk*(v->at(j))[k]*(v->at(k))[jj];
	  }					
	}					
      }

      for (int j=0; j<n; j++)
      {
	if (j != k && used_param[j]) // Pivot row or column elements 
	{
	  (v->at(j))[k] = ((v->at(j))[k])*vkk;	// Column
	  (v->at(k))[j] = ((v->at(k))[j])*vkk;	// Line
	}
      }	
    }
    else   // No more pivot value (clear those elements)
    {
      for (int j=0; j<n; j++)
      {
        if (flag->at(j))
	{
	  b->at(j) = 0.0;

	  for (int k=0; k<n; k++)
	  {
	    (v->at(j))[k] = 0.0;
	    (v->at(k))[j] = 0.0;
	  }
	}				
      }

      break;  // No more pivots anyway, stop here
    }
  }

  for (int i=0; i<n; i++) // Correct matrix V
  {
    for (int j=0; j<n; j++) {(v->at(i))[j] = sqrt(c[i])*(v->at(i))[j]*sqrt(r[j]);}
  }
	
  for (int j=0; j<n; j++)
  {
    temp[j] = 0.0;
    
    for (int jj=0; jj<n; jj++)  // Reverse matrix elements
    {
      (v->at(j))[jj] = -(v->at(j))[jj];
      temp[j] += (v->at(j))[jj]*b->at(jj);
    }					
  }

  for (int j=0; j<n; j++) {b->at(j) = temp[j];}	// The final result				

  delete []temp;
  delete []r;
  delete []c;
  delete []used_param;

  return nrank;
  }


//
// Same method but for local fit, so heavily simplified
//


  else if (v->size()==mlocal)
{

  
  int i, j, jj, k;
  double vkk, *temp;
  double eps = 0.0000000000001;

  temp = new double[n];

  
  for (i=0; i<n; i++)
  {
    flag->at(i) = true;
    diag->at(i) = fabs((v->at(i))[i]);     // save diagonal elem absolute values 	

    for (j=0; j<=i; j++)
    {
      (v->at(j))[i] = (v->at(i))[j] ;
    }
  }

  
  nrank = 0;

  for (i=0; i<n; i++)
  {
    vkk = 0.0;
    k = -1;
		
    for (j=0; j<n; j++) // First look for the pivot, ie max unused diagonal element 
    {
      if (flag->at(j) && (fabs((v->at(j))[j])>std::max(fabs(vkk),eps*diag->at(j))))
      {
        vkk = (v->at(j))[j];
        k = j;
      }
    }
		
    if (k >= 0)    // pivot found
    {
      nrank++;
      flag->at(k) = false;
      vkk = 1.0/vkk;
      (v->at(k))[k] = -vkk; // Replace pivot by its inverse
      
      for (j=0; j<n; j++)
      {
	for (jj=0; jj<n; jj++)  
	{
	  if (j != k && jj != k) // Other elements (!!! do them first as you use old v[k][j]'s !!!)
	  {
	    (v->at(j))[jj] = (v->at(j))[jj] - vkk*(v->at(j))[k]*(v->at(k))[jj];
	  }					
	}					
      }

      for (j=0; j<n; j++)
      {
	if (j != k) // Pivot row or column elements 
	{
	  (v->at(j))[k] = ((v->at(j))[k])*vkk; // Column
	  (v->at(k))[j] = ((v->at(k))[j])*vkk; // Line
	}
      }					
    }
    else  // No more pivot value (clear those elements)
    {
      for (j=0; j<n; j++)
      {
        if (flag->at(j))
	{
	  b->at(j) = 0.0;

	  for (k=0; k<n; k++)
	  {
	    (v->at(j))[k] = 0.0;
	  }
	}				
      }

      break;  // No more pivots anyway, stop here
    }
  }

  for (j=0; j<n; j++)
  {
    temp[j] = 0.0;
    
    for (jj=0; jj<n; jj++)  // Reverse matrix elements
    {
      (v->at(j))[jj] = -(v->at(j))[jj];
      temp[j] += (v->at(j))[jj]*b->at(jj);
    }					
  }

  for (j=0; j<n; j++)
  {	
    b->at(j) = temp[j];
  }					

  delete []temp;
  
  return nrank;
 }
  return 0;
  
}

/*
----------------------------------------------------------------
  CHINDL:  return the limit in chi^2/nd for n sigmas stdev authorized
----------------------------------------------------------------

  Only n=1, 2, and 3 are expected in input
----------------------------------------------------------------
*/

double  Millepede::chindl(int n, int nd)
{
  int m;
  double sn[3]        =	{0.47523, 1.690140, 2.782170};
  double table[3][30] = {{1.0000, 1.1479, 1.1753, 1.1798, 1.1775, 1.1730, 1.1680, 1.1630,
                          1.1581, 1.1536, 1.1493, 1.1454, 1.1417, 1.1383, 1.1351, 1.1321,
                          1.1293, 1.1266, 1.1242, 1.1218, 1.1196, 1.1175, 1.1155, 1.1136,
                          1.1119, 1.1101, 1.1085, 1.1070, 1.1055, 1.1040},
                         {4.0000, 3.0900, 2.6750, 2.4290, 2.2628, 2.1415, 2.0481, 1.9736,
                          1.9124, 1.8610, 1.8171, 1.7791, 1.7457, 1.7161, 1.6897, 1.6658,
                          1.6442, 1.6246, 1.6065, 1.5899, 1.5745, 1.5603, 1.5470, 1.5346,
                          1.5230, 1.5120, 1.5017, 1.4920, 1.4829, 1.4742},
                         {9.0000, 5.9146, 4.7184, 4.0628, 3.6410, 3.3436, 3.1209, 2.9468,
                          2.8063, 2.6902, 2.5922, 2.5082, 2.4352, 2.3711, 2.3143, 2.2635,
                          2.2178, 2.1764, 2.1386, 2.1040, 2.0722, 2.0428, 2.0155, 1.9901,
                          1.9665, 1.9443, 1.9235, 1.9040, 1.8855, 1.8681}};

  if (nd < 1)
  {
    return 0.0;
  }
  else
  {
    m = std::max(1,std::min(n,3));

    if (nd <= 30)
    {
      return table[m-1][nd-1];
    }
    else // approximation
    {
      return ((sn[m-1]+sqrt(float(2*nd-3)))*(sn[m-1]+sqrt(float(2*nd-3))))/float(2*nd-2);
    }
  }
}

/*
-----------------------------------------------------------
  SPAX
-----------------------------------------------------------

  multiply general M-by-N matrix A and N-vector X
 
  CALL  SPAX(A,X,Y,M,N)          Y =  A * X
                                 M   M*N  N
 
  where A = general M-by-N matrix (A11 A12 ... A1N  A21 A22 ...)
        X = N vector
        Y = M vector
-----------------------------------------------------------
*/
 
bool  Millepede::SpAX(std::vector<std::vector<double> >* a, std::vector<double>*x,std::vector<double>* y, int n, int m)
{
  int i,j; 
	
  for (i=0; i<m; i++)
  {
    y->at(i) = 0.0;	    // Reset final vector
			
    for (j=0; j<n; j++)
    {
      y->at(i) += (a->at(i))[j]*x->at(j);  // fill the vector
    }
  }
	
  return true;
}

/*
-----------------------------------------------------------
  SPAVAT
-----------------------------------------------------------

  multiply symmetric N-by-N matrix from the left with general M-by-N
  matrix and from the right with the transposed of the same  general
  matrix  to  form  symmetric  M-by-M   matrix.
  
                                                       T
  CALL SPAVAT(V,A,W,N,M)      W   =   A   *   V   *   A
   		           M*M     M*N     N*N     N*M
  
  where V = symmetric N-by-N matrix
        A = general N-by-M matrix
        W = symmetric M-by-M matrix
-----------------------------------------------------------
*/

  
bool  Millepede::SpAVAt(std::vector<std::vector<double> >* v,std::vector<std::vector<double> >* a,
                               std::vector<std::vector<double> >* w, int n, int m)
{
  for (int i=0; i<m; ++i)
  {
    for (int j=0; j<=i; ++j)   // Fill upper left triangle
    {
      (w->at(i))[j] = 0.0;      // Reset final matrix
			
      for (int k=0; k<n; ++k)
      {
	for (int l=0; l<n; ++l)
	{
	  (w->at(i))[j] += (a->at(i))[k]*(v->at(k))[l]*(a->at(j))[l];  // fill the matrix
	}
      }

      (w->at(j))[i] = (w->at(i))[j] ; // Fill the rest
    }
  }
	
  return true;
}

/*
-----------------------------------------------------------
  MAKEGLOBALFIT:  perform global params fit, once all the 'tracks'
                  have been fitted by FitLoc
-----------------------------------------------------------

  par[]        = array containing the computed global 
                 parameters (the misalignment constants)

  error[]      = array containing the error on global 
                 parameters (estimated by Millepede)

  pull[]        = array containing the corresponding pulls 

-----------------------------------------------------------
*/

bool  Millepede::MakeGlobalFit(std::vector<double>* par, std::vector<double>* error, std::vector<double>* pull)
{
  int i, j, nf, nvar;
  int itelim = 0;
  
  int nstillgood;
  
  double sum;
  
  double step[150];
  
  double trackpars[2*(mlocal+1)];
  
  int ntotal_start, ntotal;

  // Intended to compute the final global chisquare
  
  double final_cor = 0.0;
  double final_chi2 = 0.0;
  double final_ndof = 0.0;


  ntotal_start =  GetTrackNumber();
  if(m_debug) cout << "Doing global fit with " << ntotal_start << " tracks....." << endl;

  std::vector<double>* local_params = new std::vector<double>; // For non-linear treatment
  local_params->clear();
  local_params->resize(mlocal*ntotal_start);
  for (int i=0; i<mlocal*ntotal_start; i++) local_params->at(i) = 0.;
	
  if (itert <= 1) itelim=10;    // Max number of iterations

  for (i=0; i<nagb; i++)  if(m_debug) cout << "Psigm       = " << std::setprecision(5) << psigm->at(i) << endl;
  
  if(m_debug) cout << "itelim = " << itelim << endl;
  if(m_debug) cout << "itert = " << itert << endl;
  while (itert < itelim)  // Iteration for the final loop
  {
    if(m_debug) cout << "ITERATION --> " << itert << endl;
    
    ntotal = GetTrackNumber();
    if(m_debug) cout << "...using " << ntotal << " tracks..." << endl;
    
    final_cor = 0.0;
    final_chi2 = 0.0;
    final_ndof = 0.0;

    // Start by saving the diagonal elements
    
    for (i=0; i<nagb; i++) {diag->at(i) = (cgmat->at(i))[i];}

    //  Then we retrieve the different constraints: fixed parameter or global equation
    
    nf = 0; // First look at the fixed global params
    
    for (i=0; i<nagb; i++)
    {
      if (psigm->at(i) <= 0.0)   // fixed global param
      {
        nf++;

        for (j=0; j<nagb; j++)
        {
          (cgmat->at(i))[j] = 0.0;  // Reset row and column
          (cgmat->at(j))[i] = 0.0;
        }			
      }
      else if (psigm->at(i) > 0.0) {(cgmat->at(i))[i] += 1.0/(psigm->at(i)*psigm->at(i));}
    }
    
    
    nvar = nagb;  // Current number of equations	
    
    
    if(m_debug) cout << "Number of constraint equations : " << ncs << endl;
    
    if (ncs > 0) // Then the constraint equation
    {
      for (i=0; i<ncs; i++)
      {
        sum = arhs->at(i);
                
        for (j=0; j<nagb; j++)
        {
          (cgmat->at(nvar))[j] = float(ntotal)*(adercs->at(i))[j];
          (cgmat->at(j))[nvar] = float(ntotal)*(adercs->at(i))[j];          
          if(psigm->at(j) == 0.){
           if(m_debug) cout<<" --> Centi: set par to zero...."<<endl;
          
           (cgmat->at(nvar))[j] = 0.0;
           (cgmat->at(j))[nvar] = 0.0;
          }
          
          sum -= (adercs->at(i))[j]*(pparm->at(j)+dparm->at(j));
        }
        
        (cgmat->at(nvar))[nvar] = 0.0;
       
        bgvec->at(nvar) = float(ntotal)*sum;
        
        nvar++;
      }
    }
    
    
    // Intended to compute the final global chisquare
    
    
    if (itert > 1) 
    {
      for (j=0; j<nagb; j++)
      {
        
        if(m_debug) cout << "Psigm       = " << std::setprecision(5) << psigm->at(j) << endl;
        
        if(m_debug) cout << "diag. value : " << j << " = " << std::setprecision(5) << (cgmat->at(j))[j] << endl;
        
        for (i=0; i<nagb; i++)
        {
          if (psigm->at(i) > 0.0)
          {
            final_cor += step[j]*(cgmat->at(j))[i]*step[i]; 
            if (i == j) final_cor -= step[i]*step[i]/(psigm->at(i)*psigm->at(i));
          }
        }
      }
    }
    
    if(m_debug) cout << " Final coeff is " << final_cor << endl;		
    if(m_debug) cout << " Final NDOFs =  " << nagb << endl;
    
    //  The final matrix inversion
    
    
    nrank = SpmInv(cgmat, bgvec, nvar, scdiag, scflag);
    
    for (i=0; i<nagb; i++)
    {
      dparm->at(i) += bgvec->at(i);    // Update global parameters values (for iterations)
      
      if(m_debug) cout << "bgvec[" << i << "] = " << bgvec->at(i) << endl;
      if(m_debug) cout << "dparm[" << i << "] = " << dparm->at(i) << endl;
      
      if(m_debug) cout << "cgmat[" << i << "][" << i << "] = " << (cgmat->at(i))[i] << endl;
      if(m_debug) cout << "err = " << sqrt(fabs((cgmat->at(i))[i])) << endl;
      if(m_debug) cout << "cgmat * diag = " << std::setprecision(5) << (cgmat->at(i))[i]*diag->at(i) << endl;
      
      step[i] = bgvec->at(i);
      //      if (bgvec[i] < psigm[i]/10.)  psigm[i] = -1.;  // Fix parameter if variation becomes too small
      
      if (itert == 1 ) error->at(i) = (cgmat->at(i))[i]; // Unfitted error //hier!!
    }
    
    if(m_debug) cout << "" << endl;
    if(m_debug) cout << "The rank defect of the symmetric " << nvar << " by " << nvar 
                     << " matrix is " << nvar-nf-nrank << " (bad if non 0)" << endl;
    //    cout << itert << endl;
    if (itert == 0)  break;       

    itert++;

    if(m_debug) cout << "" << endl;
    if(m_debug) cout << "Total : " << loctot << " local fits, " 
                     << locrej << " rejected." << endl;
     cout << "Total : " << loctot << " local fits, " 
                     << locrej << " rejected." << endl;
    // Reinitialize parameters for iteration
    
    loctot = 0;
    locrej = 0;
    
    if (cfactr != cfactref && sqrt(cfactr) > 1.2*cfactref)
    {
      cfactr = sqrt(cfactr);
    }
    else
    {
      cfactr = cfactref;
      //      itert = itelim;
    }
    
    if (itert == itelim)  break;  // End of story         
    
    if(m_debug) cout << "Iteration " << itert << " with cut factor " << cfactr << endl;
    
    // Reset global variables
    
    for (i=0; i<nvar; i++)
    {
      bgvec->at(i) = 0.0;
      for (j=0; j<nvar; j++)
      {
        (cgmat->at(i))[j] = 0.0;
      }
    }
    
    //
    // We start a new iteration
    //
    
    // First we read the stores for retrieving the local params
    
    nstillgood = 0;	
    int k = 0;
    for (int s = 0; s < storeare->size(); ++s){
      if (storeind->at(s) == -1 && storeare->at(s)< 62500)
	//      std::cout << dparm->at(k) << "    " << storeare->at(s) << "   " <<storeind->at(s)<<"   " <<endl;
      //storeind->at(s)+=dparm->at(k);
      ++k;
      if (k == 48) k=0;
    }
    for (int i=0; i<ntotal_start; i++)
    {
      int rank_i = 0;
      int rank_f = 0;
      
      
      (i>0) ? rank_i = abs(storeplace->at(i-1)) : rank_i = 0;
      rank_f = storeplace->at(i);
      
      if(m_debug) cout << "Track " << i << " : " << endl;
      if(m_debug) cout  << "Starts at " << rank_i << endl;
      if(m_debug) cout  << "Ends at " << rank_f << endl;
      
      if (rank_f >= 0) // Fit is still OK
      {
        indst->clear();
	
        arest->clear();
        int k = 0;
        for (int j=rank_i; j<rank_f; j++)
        {
          indst->push_back(storeind->at(j));
          
	  //cout << storeare->at(j) << "    " << storenl->at(j) << endl;
	  ++k;
          if (storenl->at(j) == 0) arest->push_back(storeare->at(j));
	 
                    
          if (itert > 1 && storenl->at(j) != 0.) // Non-linear treatment (after two iterations)
          {  
            int local_index = int(storenl->at(j))/nonlin_param;
	    
            
	      
	    arest->push_back(storeare->at(j) + local_params->at(nalc*i+local_index)
                             *(storenl->at(j)-nonlin_param*(local_index+0.5)));
          }
        }
        //for (int g = 0 ; g < arest->size(); ++g) cout <<arest->at(g) << endl;
        for (int j=0; j<2*nalc; j++) {trackpars[j] = 0.;}	
        
        bool sc = FitLoc(i,trackpars,0); // Redo the fit
        
        for (int j=0; j<nalc; j++) {local_params->at(nalc*i+j) = trackpars[2*j];}
        
        if (sc) final_chi2 += trackpars[2*nalc+1]; 
        if (sc) final_ndof += trackpars[2*nalc]; 
        
        (sc) 
          ? nstillgood++
          : storeplace->at(i) = -rank_f; 
      }
    } // End of loop on fits
    
    if(m_debug) cout << " Final chisquare is " << final_chi2 << endl;		
    if(m_debug) cout << " Final chi/DOF =  " << final_chi2/(int(final_ndof)-nagb+ncs) << endl;
    
    SetTrackNumber(nstillgood);
    //PrtGlo();
  } // End of iteration loop
	
  PrtGlo(); // Print the final results
  
  for (j=0; j<nagb; j++)
  {
    par->at(j)   = dparm->at(j);
    m_par->at(j) = par->at(j);
    dparm->at(j) = 0.;
    pull->at(j)  = par->at(j)/sqrt(psigm->at(j)*psigm->at(j)-(cgmat->at(j))[j]);
    error->at(j) = sqrt(fabs((cgmat->at(j))[j]));
  }
  
  cout << std::setw(10) << " " << endl;
  cout << std::setw(10) << "            * o o                   o      " << endl;
  cout << std::setw(10) << "              o o                   o      " << endl;
  cout << std::setw(10) << "   o ooooo  o o o  oo  ooo   oo   ooo  oo  " << endl;
  cout << std::setw(10) << "    o  o  o o o o o  o o  o o  o o  o o  o " << endl;
  cout << std::setw(10) << "    o  o  o o o o oooo o  o oooo o  o oooo " << endl;
  cout << std::setw(10) << "    o  o  o o o o o    ooo  o    o  o o    " << endl;
  cout << std::setw(10) << "    o  o  o o o o  oo  o     oo   ooo  oo ++ ends." << endl;
  cout << std::setw(10) << "                       o                   " << endl;	  

  delete local_params;      local_params = 0;
  
  
  return true;
}



bool  Millepede::PrtGlo()
{
  double err, gcor;
	
  cout << "" << endl;
  cout << "   Result of fit for global parameters" << endl;
  cout << "   ===================================" << endl;
  cout << "    I       initial       final       differ   " 
         << "     lastcor        Error       gcor" << endl;
  cout << "-----------------------------------------" 
         << "------------------------------------------" << endl;
	
  for (int i=0; i<nagb; i++)
  {
    
    err = sqrt(fabs((cgmat->at(i))[i]));
    if ((cgmat->at(i))[i] < 0.0) err = -err;
    gcor = 0.0;

    if (i%(nagb/6) == 0)
    {
      cout << "-----------------------------------------" 
           << "------------------------------------------" << endl;
    }
    
 
    if (fabs((cgmat->at(i))[i]*diag->at(i)) > 0)
    {
      
      cout << std::setprecision(6) << std::fixed;
      gcor = sqrt(fabs(1.0-1.0/((cgmat->at(i))[i]*diag->at(i))));
      cout << std::setw(4) << i << "  / " << std::setw(10) << pparm->at(i) 
           << "  / " << std::setw(10) << pparm->at(i)+ dparm->at(i) 
           << "  / " << std::setw(13) << dparm->at(i) 
           << "  / " << std::setw(13) << bgvec->at(i) 
             << "  / " << std::setw(10) << std::setprecision(5) << err 
             << "  / " << std::setw(10) << gcor << endl;
    }
    else
    {
      cout << std::setw(4) << i << "  / " << std::setw(10) << "OFF" 
	     << "  / " << std::setw(10) << "OFF" 
	     << "  / " << std::setw(13) << "OFF" 
	     << "  / " << std::setw(13) << "OFF" 
	     << "  / " << std::setw(10) << "OFF" 
	     << "  / " << std::setw(10) << "OFF" << endl;
    }
  }


  for (int i=0; i<nagb; i++)
  {
    if(m_debug)cout << " i=" << i << "  sqrt(fabs(cgmat[i][i]))=" <<  sqrt(fabs((cgmat->at(i))[i]))
                    << " diag = " << diag->at(i) <<endl;
  }
  
  return true;
}

