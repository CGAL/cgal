#ifndef _PFACET_H_
#define _PFACET_H_

#include "dpqueue.h"

template <class FT, class Facet_handle>
class CPFacet {
protected:
  Facet_handle m_facet;
  FT m_priority;
public:
  CPFacet() {
    m_facet = Facet_handle();
    m_priority = 0.0;
  }

  CPFacet(const Facet_handle &fh, const FT priority = 0.0) {
    m_facet = fh;
    m_priority = priority;
  }

  CPFacet(const CPFacet &pfacet) {
    m_facet = pfacet.facet();
    m_priority = pfacet.priority();
  }

  virtual ~CPFacet() {}

  CPFacet& operator = (const CPFacet& pfacet) {
    m_facet = pfacet.facet();
    m_priority = pfacet.priority();
    return *this;
  }

  bool operator == (const CPFacet& pfacet) const {
    return m_facet == pfacet.facet();
  }

  bool operator < (const CPFacet& pfacet) const {
    return m_facet < pfacet.facet();
  }

  bool operator >(const CPFacet& pfacet) const {
    return m_facet > pfacet.facet();
  }

  const Facet_handle facet() const { return m_facet; }
  const FT priority() const { return m_priority; }
};

#endif