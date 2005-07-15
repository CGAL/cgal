// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_ACCESSOR_H
#define CGAL_VORONOI_DIAGRAM_2_ACCESSOR_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

template<class VDA>
struct Accessor
{
  typedef typename VDA::Non_degenerate_faces_iterator
  Non_degenerate_faces_iterator;

  typedef typename VDA::Non_degenerate_edges_iterator
  Non_degenerate_edges_iterator;

  typedef typename VDA::Valid_edges_iterator
  Valid_edges_iterator;

  typedef typename VDA::Edge_degeneracy_tester Edge_degeneracy_tester;
  typedef typename VDA::Face_degeneracy_tester Face_degeneracy_tester;

#ifndef VDA_NO_CACHED_TESTERS
  typedef typename VDA::Cached_edge_degeneracy_tester
  Cached_edge_degeneracy_tester;

  typedef typename VDA::Cached_face_degeneracy_tester
  Cached_face_degeneracy_tester;
#endif

#ifndef VDA_NO_CACHED_TESTERS
  static const Cached_edge_degeneracy_tester& edge_tester(const VDA& vda) {
    return vda.cached_e_tester_;
  }

  static const Cached_face_degeneracy_tester& face_tester(const VDA& vda) {
    return vda.cached_f_tester_;
  }
#else
  static const Edge_degeneracy_tester& edge_tester(const VDA& vda) {
    return vda.tr_.edge_degeneracy_tester_object();
  }

  static const Face_degeneracy_tester& face_tester(const VDA& vda) {
    return vda.tr_.face_degeneracy_tester_object();
  }
#endif
};


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_ACCESSOR_H

