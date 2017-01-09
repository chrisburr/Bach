#ifndef __TBBASECLASS_H__
#define __TBBASECLASS_H__

#include <stdlib.h>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "DD4hep/Factories.h"

class TbBaseClass {
 public:
  TbBaseClass() {}
  virtual ~TbBaseClass() {}

  virtual bool configuration() { return true; };
  virtual bool initialize(DD4hep::Geometry::LCDD&, std::vector<std::pair<std::string, TbBaseClass *>>) {
    return true;
  }
  virtual bool execute(std::vector<std::pair<std::string, TbBaseClass *>>) {
    return true;
  }
  virtual bool end_event() { return true; }
  virtual bool finalize() { return true; }

  void Const_S(std::string name, std::string value) {
    put(m_const_s, make_pair(name, value));
  }
  std::string Const_S(std::string name) {
    std::map<std::string, std::string>::const_iterator iter =
        m_const_s.find(name);
    if (iter == m_const_s.end()) {
      std::cout << "Couldn't find parameter " << name
                << "; \nCheck initialisation" << std::endl;
      exit(1);
    }

    return (*iter).second;
  }

  void Const_B(std::string name, bool value) {
    put(m_const_b, make_pair(name, value));
  }
  bool Const_B(std::string name) {
    std::map<std::string, bool>::const_iterator iter = m_const_b.find(name);
    if (iter == m_const_b.end()) {
      std::cout << "Couldn't find parameter " << name
                << "; \nCheck initialisation" << std::endl;
      exit(1);
    }
    return (*iter).second;
  }
  void Const_D(std::string name, double value) {
    put(m_const_d, make_pair(name, value));
  }
  double Const_D(std::string name) {
    std::map<std::string, double>::const_iterator iter = m_const_d.find(name);
    if (iter == m_const_d.end()) {
      std::cout << "Couldn't find parameter " << name
                << "; \nCheck initialisation" << std::endl;
      exit(1);
    }
    return (*iter).second;
  }
  void Const_I(std::string name, int value) {
    put(m_const_i, make_pair(name, value));
  }
  int Const_I(std::string name) {
    std::map<std::string, int>::const_iterator iter = m_const_i.find(name);
    if (iter == m_const_i.end()) {
      std::cout << "Couldn't find parameter " << name
                << "; \nCheck initialisation" << std::endl;
      exit(1);
    }
    return (*iter).second;
  }

  void PrintAlgorithm() {
    for (std::map<std::string, std::string>::const_iterator iter =
             m_const_s.begin();
         iter != m_const_s.end(); ++iter) {
      std::cout << (*iter).first << '\t' << (*iter).second << std::endl;
    }
    for (std::map<std::string, int>::const_iterator iter = m_const_i.begin();
         iter != m_const_i.end(); ++iter) {
      std::cout << (*iter).first << '\t' << (*iter).second << std::endl;
    }
    for (std::map<std::string, double>::const_iterator iter = m_const_d.begin();
         iter != m_const_d.end(); ++iter) {
      std::cout << (*iter).first << '\t' << (*iter).second << std::endl;
    }
  }

 protected:
 private:
  std::map<std::string, std::string> m_const_s;
  std::map<std::string, int> m_const_i;
  std::map<std::string, double> m_const_d;
  std::map<std::string, bool> m_const_b;
  std::string m_name;
  template <typename MapType>
  bool put(MapType &m, const typename MapType::value_type &v) {
    std::pair<typename MapType::iterator, bool> result = m.insert(v);
    if (!result.second) {
      result.first->second = v.second;
    }
    return result.second;
  }
};

typedef std::vector<std::pair<std::string, TbBaseClass *>> AlgVec;

#endif
