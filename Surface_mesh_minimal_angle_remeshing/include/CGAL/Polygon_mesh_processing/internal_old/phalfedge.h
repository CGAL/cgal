// Copyright (c) 2019  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Kaimo Hu

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
