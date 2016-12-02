#ifndef TB_ROOT_H
#define TB_ROOT_H 1

#include <map>

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile2D.h"

#include "TFile.h"
class TbROOT {
 public:
  TbROOT(){};
  virtual ~TbROOT(){};
  void File(const char *name) { m_rootfile = new TFile(name, "recreate"); }
  TFile *File() { return m_rootfile; }

  void AddHisto2D(std::string Key, const char *name, Int_t nbinsx, Float_t xmin,
                  Float_t xmax, Int_t nbinsy, Float_t ymin, Float_t ymax) {
    m_histos2D[Key] =
        new TH2F(name, name, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  }
  TH2F *Histo2D(std::string Key) { return m_histos2D[Key]; }
  void AddHisto3D(std::string Key, const char *name, Int_t nbinsx, Float_t xmin,
                  Float_t xmax, Int_t nbinsy, Float_t ymin, Float_t ymax) {
    m_histos3D[Key] =
        new TProfile2D(name, name, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  }
  TProfile2D *Histo3D(std::string Key) { return m_histos3D[Key]; }

  void AddHisto1D(std::string Key, const char *name, Int_t nbinsx, Float_t xmin,
                  Float_t xmax) {
    m_histos1D[Key] = new TH1F(name, name, nbinsx, xmin, xmax);
  }
  TH1F *Histo1D(std::string Key) { return m_histos1D[Key]; }

 private:
  TFile *m_rootfile;
  std::map<std::string, TProfile2D *> m_histos3D;
  std::map<std::string, TH2F *> m_histos2D;
  std::map<std::string, TH1F *> m_histos1D;
};
#endif
