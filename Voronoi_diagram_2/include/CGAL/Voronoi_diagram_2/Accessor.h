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

#ifndef CGAL_VORONOI_DIAGRAM_2_ACCESSOR_H
#define CGAL_VORONOI_DIAGRAM_2_ACCESSOR_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

template<class VDA>
class Accessor
{
public:
  // TYPES
  //------
  typedef typename VDA::Non_degenerate_faces_iterator
  Non_degenerate_faces_iterator;

  typedef typename VDA::Non_degenerate_edges_iterator
  Non_degenerate_edges_iterator;

  typedef typename VDA::Valid_edges_iterator
  Valid_edges_iterator;

  typedef typename VDA::Edge_rejector  Edge_rejector;
  typedef typename VDA::Face_rejector  Face_rejector;

  typedef typename VDA::Find_valid_vertex      Find_valid_vertex;

  // CONSTRUCTOR
  //------------
  Accessor(VDA* vda) : vda(vda) {}
  Accessor(const VDA* vda) : vda(const_cast<VDA*>(vda)) {}

  // METHODS
  //--------
  VDA* ptr() { return vda; }
  const VDA* ptr() const { return vda; }

private:
  VDA* vda;
};


} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL


#endif // CGAL_VORONOI_DIAGRAM_2_ACCESSOR_H
