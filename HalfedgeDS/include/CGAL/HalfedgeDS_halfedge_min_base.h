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

#ifndef CGAL_HALFEDGEDS_HALFEDGE_MIN_BASE_H
#define CGAL_HALFEDGEDS_HALFEDGE_MIN_BASE_H 1

#include <CGAL/basic.h>

namespace CGAL {

template < class Refs>
struct HalfedgeDS_halfedge_min_base_base {
    // Base_base will be used to access set_opposite(), which is
    // made private in the normal halfedge bases. Since halfedges
    // come always in pairs, managed by the HDS, the set_opposite()
    // member function is protected from the user.
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_halfedge_min_base_base<Refs> Base_base;
    typedef Tag_false                            Supports_halfedge_prev;
    typedef Tag_false                            Supports_halfedge_vertex;
    typedef Tag_false                            Supports_halfedge_face;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Vertex                Vertex;
    typedef typename Refs::Face                  Face;
private:
    Halfedge_handle  opp;
    Halfedge_handle  nxt;
public:
    Halfedge_handle       opposite()                        { return opp;}
    Halfedge_const_handle opposite() const                  { return opp;}
    void                  set_opposite( Halfedge_handle h)  { opp = h;}
    Halfedge_handle       next()                            { return nxt;}
    Halfedge_const_handle next() const                      { return nxt;}
    void                  set_next( Halfedge_handle h)      { nxt = h;}

    bool is_border() const { return false;}
        // is always false as long as faces are not supported.
};

template < class Refs>
class HalfedgeDS_halfedge_min_base
    : public HalfedgeDS_halfedge_min_base_base< Refs>
{
public:
    typedef typename Refs::Halfedge_handle Halfedge_handle;
    typedef HalfedgeDS_halfedge_min_base<Refs>       Base;
    typedef HalfedgeDS_halfedge_min_base_base<Refs>  Base_base;
private:
    void  set_opposite( Halfedge_handle h)  { Base_base::set_opposite(h);}
};

} //namespace CGAL

#endif // CGAL_HALFEDGEDS_HALFEDGE_MIN_BASE_H //
// EOF //
