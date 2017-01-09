
#ifndef TB_TOYDATA_H
#define TB_TOYDATA_H 1

/* Header for ToyData
 *
 * 2014-02-07
 *
 */

#include "TbBaseClass.h"
#include "TbGeometrySvc.h"
#include "TbHit.h"
#include "DD4hep/Factories.h"

#include <stdio.h>
#include <sstream>
#include "TRandom3.h"
// typedef std::vector<TbHit*> TbHits;
using namespace ROOT::Math;
class TbToyData : public TbBaseClass {
 public:
  /// Constructor
  TbToyData(const std::string &name);
  /// Destructor
  virtual ~TbToyData();

  bool configuration();
  bool initialize(DD4hep::Geometry::LCDD&, AlgVec);  ///< Algorithm initialization
  bool execute(AlgVec);     ///< Algorithm execution
  bool end_event();
  bool finalize();  ///< Algorithm finalization

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

  TRandom3 m_r;

  ULong_t m_seed;

  std::string m_geomfile;

  FILE *datei;
};
#endif
