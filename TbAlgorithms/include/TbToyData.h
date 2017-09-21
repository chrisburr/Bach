
#ifndef TB_TOYDATA_H
#define TB_TOYDATA_H 1

/* Header for ToyData
 *
 * 2014-02-07
 *
 */
#include <iostream>
#include <fstream>

#include "TbBaseClass.h"
#include "TbGeometrySvc.h"
#include "TbHit.h"
#include "DD4hep/Factories.h"

#include <stdio.h>
#include <sstream>
#include "TRandom3.h"

using namespace ROOT::Math;
using namespace DD4hep::Geometry;

class TbToyData : public TbBaseClass {
 public:
  /// Constructor
  TbToyData(const std::string &name);
  /// Destructor
  virtual ~TbToyData();

  bool configuration();
  bool initialize(AlgVec);  ///< Algorithm initialization
  bool execute(DD4hep::Conditions::ConditionsSlice &, AlgVec);     ///< Algorithm execution
  bool end_event();
  bool finalize(DD4hep::Conditions::ConditionsSlice &);  ///< Algorithm finalization

  TbHits *getHits() { return m_hits; }

  TbBaseClass *find(AlgVec vec, std::string name) {
    for (AlgVec::iterator it = vec.begin(); it != vec.end(); ++it) {
      if ((*it).first == name) return (*it).second;
    }
    std::cout << "Couldn't find " << name << std::endl;
    return NULL;
  }

 private:
  mutable TbGeometrySvc *m_geomSvc;
  /// Access geometry service on-demand
  TbGeometrySvc *geomSvc() const { return m_geomSvc; }
  void geomSvc(TbGeometrySvc *geo) { m_geomSvc = geo; }

  TbHits *m_hits;
  int m_nevent;

  TRandom3 m_r;

  ULong_t m_seed;

  std::string m_geomfile;

  std::ofstream m_outfile;
};
#endif
