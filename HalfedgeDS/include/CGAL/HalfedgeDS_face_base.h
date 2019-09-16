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

#ifndef CGAL_HALFEDGEDS_FACE_BASE_H
#define CGAL_HALFEDGEDS_FACE_BASE_H 1

#include <CGAL/basic.h>

namespace CGAL {

// We use Tag_false to indicate that no plane type is provided.

template < class Refs, class T = Tag_true, class Pln = Tag_false>
class HalfedgeDS_face_base;

template < class Refs >
class HalfedgeDS_face_base< Refs, Tag_false, Tag_false> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_face_base< Refs, Tag_false, Tag_false>  Base;
    typedef Tag_false                            Supports_face_halfedge;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Vertex                Vertex;
    typedef typename Refs::Halfedge              Halfedge;
    // Additional tags required by Polyhedron.
    typedef Tag_false                            Supports_face_plane;
    struct Plane_not_supported {};
    typedef Plane_not_supported                  Plane;
    // No longer required.
    // typedef Tag_false                            Supports_face_normal;
};

template < class Refs >
class HalfedgeDS_face_base< Refs, Tag_true, Tag_false> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_face_base< Refs, Tag_true, Tag_false>   Base;
    typedef Tag_true                             Supports_face_halfedge;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Vertex                Vertex;
    typedef typename Refs::Halfedge              Halfedge;
    // Additional tags required by Polyhedron.
    typedef Tag_false                            Supports_face_plane;
    struct Plane_not_supported {};
    typedef Plane_not_supported                  Plane;
    // No longer required.
    //typedef Tag_false                            Supports_face_normal;
private:
    Halfedge_handle hdg;
public:
    Halfedge_handle       halfedge()                        { return hdg; }
    Halfedge_const_handle halfedge() const                  { return hdg; }
    void                  set_halfedge( Halfedge_handle h)  { hdg = h; }
};

template < class Refs, class Pln >
class HalfedgeDS_face_base< Refs, Tag_false, Pln> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_face_base< Refs, Tag_false, Pln>     Base;
    typedef Tag_false                            Supports_face_halfedge;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Vertex                Vertex;
    typedef typename Refs::Halfedge              Halfedge;
    // Additional tags and types required by Polyhedron.
    typedef Tag_true                             Supports_face_plane;
    typedef Pln                                  Plane;
    // No longer required.
    //typedef Tag_true                             Supports_face_normal;
    //typedef Trts                                 Traits;
    //typedef typename Traits::Normal              Normal;
    //typedef typename Traits::Plane               Plane;
private:
    Plane  pln;
public:
    HalfedgeDS_face_base() {}
    HalfedgeDS_face_base( const Plane& g) : pln(g) {}
    Plane&                plane()                           { return pln; }
    const Plane&          plane() const                     { return pln; }
    // No longer required.
    // Normal              normal() const { return pln.orthogonal_vector();}
};

template < class Refs, class Pln >
class HalfedgeDS_face_base< Refs, Tag_true, Pln> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_face_base< Refs, Tag_true, Pln>      Base;
    typedef Tag_true                             Supports_face_halfedge;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Vertex                Vertex;
    typedef typename Refs::Halfedge              Halfedge;
    // Additional tags and types required by Polyhedron.
    typedef Tag_true                             Supports_face_plane;
    typedef Pln                                  Plane;
    // No longer required.
    //typedef Tag_true                             Supports_face_normal;
    //typedef Trts                                 Traits;
    //typedef typename Traits::Normal              Normal;
    //typedef typename Traits::Plane               Plane;
private:
    Halfedge_handle hdg;
    Plane           pln;
public:
    HalfedgeDS_face_base() {}
    HalfedgeDS_face_base( const Plane& g) : pln(g) {}
    Halfedge_handle       halfedge()                        { return hdg; }
    Halfedge_const_handle halfedge() const                  { return hdg; }
    void                  set_halfedge( Halfedge_handle h)  { hdg = h; }
    Plane&                plane()                           { return pln; }
    const Plane&          plane() const                     { return pln; }
    // No longer required.
    //Normal                normal() const { return pln.orthogonal_vector();}
};

} //namespace CGAL

#endif // CGAL_HALFEDGEDS_FACE_BASE_H //
// EOF //
