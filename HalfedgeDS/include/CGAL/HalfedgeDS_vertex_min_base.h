// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_HALFEDGEDS_VERTEX_MIN_BASE_H
#define CGAL_HALFEDGEDS_VERTEX_MIN_BASE_H 1

#include <CGAL/basic.h>

namespace CGAL {

template < class Refs>
class HalfedgeDS_vertex_min_base {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_vertex_min_base< Refs>    Base;
    typedef Tag_false                            Supports_vertex_halfedge;
    typedef Tag_false                            Supports_vertex_point;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Halfedge              Halfedge;
    typedef typename Refs::Face                  Face;
};

} //namespace CGAL

#endif // CGAL_HALFEDGEDS_VERTEX_MIN_BASE_H //
// EOF //
