// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_EDGE_LESS_H
#define CGAL_VORONOI_DIAGRAM_2_EDGE_LESS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

template<class Edge_t>
struct Edge_less
{
  typedef Edge_t        Edge;
  typedef bool          result_type;

  bool operator()(const Edge& e1, const Edge& e2) const {
    if ( e1.first != e2.first ) { return e1.first < e2.first; }
    return e1.second < e2.second;
  }
};


} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_EDGE_LESS_H
