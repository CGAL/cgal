#ifndef _PVERTEX_H_
#define _PVERTEX_H_

#include "dpqueue.h"

template <class FT, class Vertex_handle>
class CPVertex {
protected:
  Vertex_handle m_vertex;
  FT m_priority;
public:
  CPVertex() {
    m_vertex = Vertex_handle();
    m_priority = 0.0;
  }

  CPVertex(const Vertex_handle &he, const FT priority = 0.0) {
    m_vertex = he;
    m_priority = priority;
  }

  CPVertex(const CPVertex &pvertex) {
    m_vertex = pvertex.vertex();
    m_priority = pvertex.priority();
  }

  virtual ~CPVertex() {}

  CPVertex& operator = (const CPVertex& pvertex) {
    m_vertex = pvertex.vertex();
    m_priority = pvertex.priority();
    return *this;
  }

  bool operator == (const CPVertex& pvertex) const {
    return m_vertex == pvertex.vertex();
  }

  bool operator < (const CPVertex& pvertex) const {
    return m_vertex < pvertex.vertex();
  }

  bool operator >(const CPVertex& pvertex) const {
    return m_vertex > pvertex.vertex();
  }

  const Vertex_handle vertex() const { return m_vertex; }
  const FT priority() const { return m_priority; }
};

#endif