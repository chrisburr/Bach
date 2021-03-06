#ifndef TB_DECODER_H
#define TB_DECODER_H 1

#include <fstream>

// Tb
#include "TbBaseClass.h"
#include "TbHit.h"
#include "DD4hep/Factories.h"
/** @class TbDecoder TbDecoder.h
 *
 */

class TbDecoder : public TbBaseClass {
 public:
  /// Constructor
  TbDecoder(const std::string &name);
  /// Destructor
  virtual ~TbDecoder();

  bool configuration();
  bool initialize(AlgVec);  ///< Algorithm initialization
  bool execute(dd4hep::cond::ConditionsSlice &, AlgVec);     ///< Algorithm execution
  bool finalize(dd4hep::cond::ConditionsSlice &);          ///< Algorithm finalization

  TbHits *getHits() { return hits; }

 private:
  std::ifstream m_inputFile;
  std::string m_name;

  TbHits *hits;
  int m_lastEvent;
};
#endif
