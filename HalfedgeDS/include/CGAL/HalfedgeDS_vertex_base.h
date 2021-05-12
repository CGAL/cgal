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

#ifndef CGAL_HALFEDGEDS_VERTEX_BASE_H
#define CGAL_HALFEDGEDS_VERTEX_BASE_H 1

#include <CGAL/basic.h>

namespace CGAL {

// We use Tag_false to indicate that no point type is provided.

template < class Refs, class T = Tag_true, class P = Tag_false>
class HalfedgeDS_vertex_base;

template < class Refs >
class HalfedgeDS_vertex_base< Refs, Tag_false, Tag_false> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_vertex_base< Refs, Tag_false, Tag_false>  Base;
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

template < class Refs>
class HalfedgeDS_vertex_base< Refs, Tag_true, Tag_false> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_vertex_base< Refs, Tag_true, Tag_false>   Base;
    typedef Tag_true                             Supports_vertex_halfedge;
    typedef Tag_false                            Supports_vertex_point;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Halfedge              Halfedge;
    typedef typename Refs::Face                  Face;
private:
    Halfedge_handle hdg;
public:
    Halfedge_handle       halfedge()                        { return hdg; }
    Halfedge_const_handle halfedge() const                  { return hdg; }
    void                  set_halfedge( Halfedge_handle h)  { hdg = h; }
};

template < class Refs, class P>
class HalfedgeDS_vertex_base< Refs, Tag_false, P> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_vertex_base< Refs, Tag_false, P>     Base;
    typedef Tag_false                            Supports_vertex_halfedge;
    typedef Tag_true                             Supports_vertex_point;
    typedef P                                    Point;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Halfedge              Halfedge;
    typedef typename Refs::Face                  Face;
private:
    Point   p;
public:
    HalfedgeDS_vertex_base() {}
    HalfedgeDS_vertex_base( const Point& pp) : p(pp) {}
    Point&                point()                           { return p; }
    const Point&          point() const                     { return p; }
};

template < class Refs, class P>
class HalfedgeDS_vertex_base< Refs, Tag_true, P> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_vertex_base< Refs, Tag_true, P>      Base;
    typedef Tag_true                             Supports_vertex_halfedge;
    typedef Tag_true                             Supports_vertex_point;
    typedef P                                    Point;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Halfedge              Halfedge;
    typedef typename Refs::Face                  Face;
private:
    Halfedge_handle hdg;
    Point           p;
public:
    HalfedgeDS_vertex_base() {}
    HalfedgeDS_vertex_base( const Point& pp) : p(pp) {}
    Halfedge_handle       halfedge()                        { return hdg; }
    Halfedge_const_handle halfedge() const                  { return hdg; }
    void                  set_halfedge( Halfedge_handle h)  { hdg = h; }
    Point&                point()                           { return p; }
    const Point&          point() const                     { return p; }
};

} //namespace CGAL

#endif // CGAL_HALFEDGEDS_VERTEX_BASE_H //
// EOF //
