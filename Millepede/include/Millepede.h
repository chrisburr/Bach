#ifndef MILLEPEDE_H
#define MILLEPEDE_H 1

// Include files
#include <map>
#include <stdlib.h>
#include <vector>
/** @class Millepede Millepede.h
 *
 *
 *  @author Christoph Hombach
 *  @date   2012-06-19
 */
class Millepede {
public:
  /// Standard constructor
  Millepede();

  virtual ~Millepede(); ///< Destructor

  bool InitMille(bool DOF[], double Sigm[], int nglo, int nloc,
                 double startfact, int nstd, double res_cut,
                 double res_cut_init, int n_fits);
  bool MakeGlobalFit(std::vector<double> *, std::vector<double> *,
                     std::vector<double> *);
  int GetTrackNumber();
  void SetTrackNumber(int value);
  void DeleteVector(std::vector<std::vector<double>> *);

  bool ZerLoc(double dergb[], double derlc[], double dernl[], double dernl_i[]);
  bool ParSig(int index, double sigma);
  bool InitUn(double cutfac);
  bool EquLoc(double dergb[], double derlc[], double dernl[], double dernl_i[],
              double rmeas, double sigma);
  bool FitLoc(int n, double track_params[], int single_fit);
  bool PrtGlo();
  bool ParGlo(int index, double param);
  bool ConstF(double dercs[], double rhs);

  int SpmInv(std::vector<std::vector<double>> *, std::vector<double> *, int n,
             std::vector<double> *, std::vector<bool> *);

  double chindl(int n, int nd);
  bool SpAX(std::vector<std::vector<double>> *, std::vector<double> *,
            std::vector<double> *, int n, int m);

  bool SpAVAt(std::vector<std::vector<double>> *,
              std::vector<std::vector<double>> *,
              std::vector<std::vector<double>> *, int n, int m);

  void DrawResiduals(int);

  void UpdateTracks();

  bool aligndut;
  int dut;

  int itert, nst, nfl, ncs, nstdev;
  int loctot, locrej, nagb, nalc, nrank;
  int store_row_size;

  int m_NofDOFs;
  int m_track_number;
  double m_residual_cut_init;
  double m_residual_cut;
  bool m_iteration;
  int m_fixed;
  double cfactr, cfactref;
  bool rotate;

  std::vector<double> *mis_const;
  std::vector<double> *mis_error;
  std::vector<double> *mis_pull;

  std::vector<int> *indst;
  std::vector<double> *arest;
  std::vector<double> *arenl;

  std::vector<int> *storeind;
  std::vector<int> *storeplace;
  std::vector<double> *storeare;

  std::vector<double> *storenl;
  std::vector<double> *psigm;
  std::vector<double> *m_par;

  std::vector<std::vector<double>> *clcmat;
  std::vector<std::vector<double>> *corrm;
  std::vector<std::vector<double>> *adercs;
  std::vector<std::vector<double>> *clmat;
  std::vector<std::vector<double>> *cgmat;

  bool m_debug;

  static const int mglobl = 2000; // Max. number of global parameters
  static const int mlocal = 50;   // Max. number of local parameters
  static const int mcs = 50;      // Max. number of constraint equations
  static const int mgl = 2050;    // mglobl+mlocal

  static const int nonlin_param = 1000000; // For non-linear terms treatment
  // See how it works in EquLoc() and MakeGlobalFit()
  // , Should be much larger than any of the derivatives

  static const int mobjects = 100000; // Max. number of objects to align

  std::vector<std::vector<double>> *m_glmat;
  std::vector<std::vector<double>> *m_glmatinv;
  std::vector<double> *corrv;
  std::vector<double> *pparm;
  std::vector<double> *dparm;

  std::vector<double> *scdiag;
  std::vector<double> *blvec;
  std::vector<double> *arhs;
  std::vector<double> *diag;
  std::vector<double> *bgvec;

  std::vector<int> *indgb;
  std::vector<int> *nlnpa;
  std::vector<int> *indnz;
  std::vector<int> *indbk;
  std::vector<int> *indlc;

  std::vector<bool> *scflag;

protected:
private:
};
#endif // TESTBEAMMILLEPEDE_H
