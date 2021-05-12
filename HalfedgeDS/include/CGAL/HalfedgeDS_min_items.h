// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_HALFEDGEDS_MIN_ITEMS_H
#define CGAL_HALFEDGEDS_MIN_ITEMS_H 1
#include <CGAL/HalfedgeDS_vertex_min_base.h>
#include <CGAL/HalfedgeDS_halfedge_min_base.h>
#include <CGAL/HalfedgeDS_face_min_base.h>

namespace CGAL {

class HalfedgeDS_min_items {
public:
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef HalfedgeDS_vertex_min_base< Refs>   Vertex;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper {
        typedef HalfedgeDS_halfedge_min_base< Refs> Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        typedef HalfedgeDS_face_min_base< Refs>     Face;
    };
};

} //namespace CGAL
#endif // CGAL_HALFEDGEDS_MIN_ITEMS_H //
// EOF //
