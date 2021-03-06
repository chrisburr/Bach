#ifndef TB_CLUSTERING_H
#define TB_CLUSTERING_H 1

//#include "TbDecoder.h"
#include "Math/Point3D.h"
#include "TbBaseClass.h"
#include "TbCluster.h"
#include "TbDecoder.h"
#include "TbGeometrySvc.h"
#include "TbROOT.h"
#include "TbToyData.h"
#include "DD4hep/Factories.h"
/** @class TbClustering TbClustering.h
 *
 */
using namespace ROOT::Math;
using namespace dd4hep;

class TbClustering : public TbBaseClass {
 public:
  /// Constructor
  TbClustering(const std::string &name);
  /// Destructor
  virtual ~TbClustering();
  virtual bool configuration();
  virtual bool initialize(AlgVec);  ///< Algorithm initialization
  virtual bool execute(dd4hep::cond::ConditionsSlice &, AlgVec);     ///< Algorithm execution
  bool end_event();
  virtual bool finalize(dd4hep::cond::ConditionsSlice &);  ///< Algorithm finalization

  TbBaseClass *find(AlgVec vec, std::string name) {
    for (AlgVec::iterator it = vec.begin(); it != vec.end(); ++it) {
      if ((*it).first == name) return (*it).second;
    }
    std::cout << "Couldn't find " << name << std::endl;
    return NULL;
  }

  /// Access geometry service on-demand
  TbGeometrySvc *geomSvc() const { return m_geomSvc; }
  void geomSvc(TbGeometrySvc *geo) { m_geomSvc = geo; }
  /// Access decoder on-demand
  TbDecoder *decoder() { return m_decoder; }
  void decoder(TbDecoder *dec) { m_decoder = dec; }
  /// Access clusters on-demand
  std::map<std::string, TbClusters *> clusters() { return m_clusters; }
  void clusters(std::map<std::string, TbClusters *> clus) { m_clusters = clus; }
  /// Access ToyData service on-demand
  TbToyData *toyData() const { return m_toydata; }
  void toyData(TbToyData *td) { m_toydata = td; }

 private:
  unsigned int m_nEvents;
  mutable TbGeometrySvc *m_geomSvc;
  mutable TbDecoder *m_decoder;
  mutable std::map<std::string, TbClusters *> m_clusters;
  mutable TbToyData *m_toydata;
};

#endif
