// Copyright (c) 2001  Tel-Aviv University (Israel).
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
// Author(s)     : Shai Hirsch <shaihi@post.tau.ac.il>
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
