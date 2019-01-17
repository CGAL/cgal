#ifndef _PHALFEDGE_H_
#define _PHALFEDGE_H_

#include "dpqueue.h"

template <class FT, class Halfedge_handle>
class CPHalfedge {
protected:
  Halfedge_handle m_halfedge;
  FT              m_priority;
  //Point			      m_point;		
  FT m_x, m_y, m_z;   // for simulated edge collapse: the destination after relocation.

public:
  CPHalfedge() {
    m_halfedge = Halfedge_handle();
    m_priority = 0.0;
  }

  CPHalfedge(const Halfedge_handle& he, const FT priority = 0.0) {
    m_halfedge = he;
    m_priority = priority;
  }

  CPHalfedge(const Halfedge_handle& he, const FT priority, FT x, FT y, FT z) {
    m_halfedge = he;
    m_priority = priority;
    m_x = x;
    m_y = y;
    m_z = z;
  }

  CPHalfedge(const CPHalfedge& phedge) {
    m_halfedge = phedge.halfedge();
    m_priority = phedge.priority();
    //m_point = phedge.point();
    m_x = phedge.m_x;
    m_y = phedge.m_y;
    m_z = phedge.m_z;
  }

  virtual ~CPHalfedge() { }

  CPHalfedge& operator = (const CPHalfedge& phedge) {
    m_halfedge = phedge.halfedge();
    m_priority = phedge.priority();
    //m_point = phedge.point();
    m_x = phedge.m_x;
    m_y = phedge.m_y;
    m_z = phedge.m_z;
    return *this;
  }

  bool operator == (const CPHalfedge& phedge) const {
    return m_halfedge == phedge.halfedge();
  }

  bool operator < (const CPHalfedge& phedge) const {
    return m_halfedge < phedge.halfedge();
  }

  const Halfedge_handle halfedge() const { return m_halfedge; }
  const FT priority() const { return m_priority; }
  //const Point point() const { return m_point; }
  const FT x() const { return m_x; }
  const FT y() const { return m_y; }
  const FT z() const { return m_z; }
};

template <class T>
class CDPQueue_short : public DynamicPriorityQueue<T> {
public:
  CDPQueue_short() { }
  ~CDPQueue_short() { }
  bool compare(const T& a, const T& b) const  { return a.priority() < b.priority(); }
};

template <class T>
class CDPQueue_long : public DynamicPriorityQueue<T> {
public:
  CDPQueue_long() { }
  ~CDPQueue_long() { }
  bool compare(const T& a, const T& b) const  { return a.priority() > b.priority(); }
};

#endif
