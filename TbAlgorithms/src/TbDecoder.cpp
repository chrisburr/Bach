#include <sstream>

// Local
#include "TbDecoder.h"
#include "DD4hep/Factories.h"

//=============================================================================
/// Standard constructor
//=============================================================================
TbDecoder::TbDecoder(const std::string &name) { m_name = name; }

//=============================================================================
/// Destructor
//=============================================================================
TbDecoder::~TbDecoder() {}

//=============================================================================
/// Initialization
//=============================================================================

bool TbDecoder::configuration() {
  Const_S("Name", m_name);
  Const_S("InputFile", "../../Bach/cmt/hits.txt");
  return true;
}

bool TbDecoder::initialize(AlgVec algos) {
  // Open the specified file.
  m_inputFile.open(Const_S("InputFile").c_str(), std::ios::in);
  if (!m_inputFile.is_open()) {
    std::cout << "Cannot open file " << Const_S("InputFile") << std::endl;
    return false;
  }
  PrintAlgorithm();
  return true;
}

//=============================================================================
/// Main execution
//=============================================================================
bool TbDecoder::execute(DD4hep::Conditions::ConditionsSlice &slice, AlgVec algos) {
  if (!m_inputFile.is_open()) return false;
  // Container for storing hits.
  hits = new TbHits;

  std::string line;
  while (getline(m_inputFile, line)) {
    std::stringstream ss(line);
    int evt, col, row, adc;
    std::string id;
    ss >> evt >> id >> col >> row >> adc;

    if (m_lastEvent < 0) m_lastEvent = evt;

    TbHit *m_hit = new TbHit;
    m_hit->setId(id);
    m_hit->setCol(col);
    m_hit->setRow(row);
    m_hit->setAdc(adc);
    if (evt != m_lastEvent) {
      // Hit belongs to the next event.
      m_lastEvent = evt;
      break;
    }
    hits->push_back(m_hit);
  }

  if (m_inputFile.bad()) {
    std::cout << "Error reading input file " << Const_S("InputFile").c_str()
              << std::endl;
    return false;
  }
  return true;
}

//=============================================================================
/// Finalize
//=============================================================================
bool TbDecoder::finalize(DD4hep::Conditions::ConditionsSlice &slice) {
  if (m_inputFile.is_open()) m_inputFile.close();
  // delete hits;
  return true;
}
