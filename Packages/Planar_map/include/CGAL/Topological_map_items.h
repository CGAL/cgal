// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.4-I-40 $
// release_date  : $CGAL_Date: 2001/12/28 $
//
// file          : include/CGAL/Topological_map_items.h
// package       : Planar_map (5.80)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Shai Hirsch <shaihi@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin halperin<@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_TOPOLOGICAL_MAP_ITEMS_H
#define CGAL_TOPOLOGICAL_MAP_ITEMS_H 1

#ifndef CGAL_HALFEDGEDS_VERTEX_BASE_H
#include <CGAL/HalfedgeDS_vertex_base.h>
#endif
#ifndef CGAL_HALFEDGEDS_HALFEDGE_BASE_H
#include <CGAL/HalfedgeDS_halfedge_base.h>
#endif
#ifndef CGAL_HALFEDGEDS_FACE_BASE_H
#include <CGAL/Topological_map_face_base.h>
#endif

CGAL_BEGIN_NAMESPACE

class Topological_map_items {
public:
  template < class Refs, class Traits >
  struct Vertex_wrapper {
    typedef typename Traits::Point_2 Point;
    typedef HalfedgeDS_vertex_base< Refs, Tag_true, Point > Vertex;
  };
  template < class Refs, class Traits>
  struct Halfedge_wrapper {
    typedef typename Traits::X_curve X_curve;
    typedef HalfedgeDS_halfedge_base< Refs >                Halfedge;
  };
  template < class Refs, class Triats>
  struct Face_wrapper {
    typedef Topological_map_face_list_base< Refs >          Face;
  };
};

CGAL_END_NAMESPACE
#endif // CGAL_TOPOLOGICAL_MAP_ITEMS_H
// EOF //
