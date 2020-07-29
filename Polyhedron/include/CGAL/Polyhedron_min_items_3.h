// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>)

#ifndef CGAL_POLYHEDRON_MIN_ITEMS_3_H
#define CGAL_POLYHEDRON_MIN_ITEMS_3_H 1

#include <CGAL/license/Polyhedron.h>

#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_face_base.h>

namespace CGAL {

class Polyhedron_min_items_3 {
public:
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef typename Traits::Point_3 Point;
        typedef HalfedgeDS_vertex_base< Refs, Tag_false, Point> Vertex;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper {
        typedef HalfedgeDS_halfedge_base< Refs>                 Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        typedef HalfedgeDS_face_base< Refs, Tag_false>          Face;
    };
};

} //namespace CGAL
#endif // CGAL_POLYHEDRON_MIN_ITEMS_3_H //
// EOF //
