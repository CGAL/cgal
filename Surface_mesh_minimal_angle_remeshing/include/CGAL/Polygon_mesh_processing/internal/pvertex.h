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