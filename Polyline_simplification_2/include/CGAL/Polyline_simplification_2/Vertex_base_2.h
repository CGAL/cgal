// Copyright (c) 2012 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Philipp Moeller

#ifndef CGAL_POLYLINE_SIMPLIFICATION_VERTEX_BASE_2_H
#define CGAL_POLYLINE_SIMPLIFICATION_VERTEX_BASE_2_H

#include <CGAL/license/Polyline_simplification_2.h>


#include <CGAL/Triangulation_vertex_base_2.h>

namespace CGAL {

namespace Polyline_simplification_2 {

/// \ingroup PkgPolylineSimplification2Classes

/// A vertex base class with data members needed by the simplification algorithm.
/// \tparam Vb must be a model of the concept `TriangulationVertexBase_2`
/// \cgalModels `PolylineSimplificationVertexBase_2`.
  template<class K, class Vb = CGAL::Triangulation_vertex_base_2<K> >
class Vertex_base_2
  : public Vb
{
  typedef typename K::FT FT;
  typedef Vb                                            Base;
  typedef typename Base::Triangulation_data_structure   Tds;

  bool m_removable;
  FT m_cost;

#ifndef DOXYGEN_RUNNING
public:
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other      Vb2;
    typedef Vertex_base_2<K,Vb2>         Other;
  };

  Vertex_base_2()
    : Base(), m_removable(true), m_cost(-1.0)
  {}


  bool is_removable() const
  {
    return m_removable;
  }

  void set_removable(bool b)
  {
    m_removable = b;
  }

  FT cost() const
  {
    return m_cost;
  }

  void set_cost(const FT& ft)
  {
    m_cost = ft;
  }
#endif
};


} // Polyline_simplification_2
} // CGAL

#endif /* CGAL_POLYLINE_SIMPLIFICATION_VERTEX_BASE_2_H */
