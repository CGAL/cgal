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

#ifndef CGAL_HALFEDGEDS_ITEMS_DECORATOR_H
#define CGAL_HALFEDGEDS_ITEMS_DECORATOR_H 1

#include <CGAL/basic.h>
#include <cstddef>

namespace CGAL {

template < class p_HDS >
class HalfedgeDS_items_decorator {
public:

// TYPES
// ----------------------------------
    typedef p_HDS                                 HDS;
    typedef p_HDS                                 HalfedgeDS;
    typedef typename HDS::Traits                  Traits;
    typedef typename HDS::Vertex                  Vertex;
    typedef typename HDS::Halfedge                Halfedge;
    typedef typename HDS::Face                    Face;

    typedef typename HDS::Vertex_handle           Vertex_handle;
    typedef typename HDS::Vertex_const_handle     Vertex_const_handle;
    typedef typename HDS::Vertex_iterator         Vertex_iterator;
    typedef typename HDS::Vertex_const_iterator   Vertex_const_iterator;

    typedef typename HDS::Halfedge_handle         Halfedge_handle;
    typedef typename HDS::Halfedge_const_handle   Halfedge_const_handle;
    typedef typename HDS::Halfedge_iterator       Halfedge_iterator;
    typedef typename HDS::Halfedge_const_iterator Halfedge_const_iterator;

    typedef typename HDS::Face_handle             Face_handle;
    typedef typename HDS::Face_const_handle       Face_const_handle;
    typedef typename HDS::Face_iterator           Face_iterator;
    typedef typename HDS::Face_const_iterator     Face_const_iterator;

    typedef typename HDS::size_type               size_type;
    typedef typename HDS::difference_type         difference_type;
    typedef typename HDS::iterator_category       iterator_category;

// The following types are equal to either `Tag_true' or `Tag_false',
// dependent whether the named feature is supported or not.

    typedef typename HDS::Supports_vertex_halfedge
                                                  Supports_vertex_halfedge;
    typedef typename HDS::Supports_halfedge_prev  Supports_halfedge_prev;
    typedef typename HDS::Supports_halfedge_vertex
                                                  Supports_halfedge_vertex;
    typedef typename HDS::Supports_halfedge_face  Supports_halfedge_face;
    typedef typename HDS::Supports_face_halfedge  Supports_face_halfedge;

    typedef typename HDS::Supports_removal        Supports_removal;

protected:
    typedef typename Vertex::Base                 VBase;
    typedef typename Halfedge::Base               HBase;
    typedef typename Halfedge::Base_base          HBase_base;
    typedef typename Face::Base                   FBase;

public:

// CREATION
// ----------------------------------

    HalfedgeDS_items_decorator() {}

// Access Functions
// ----------------------------------

    Halfedge_handle get_vertex_halfedge( Vertex_handle v) const {
        // returns the incident halfedge of v if supported,
        // `Halfedge_handle()' otherwise.
        return get_vertex_halfedge( v, Supports_vertex_halfedge());
    }

    Vertex_handle get_vertex( Halfedge_handle h) const {
        // returns the incident vertex of h if supported,
        // `Vertex_handle()' otherwise.
        return get_vertex( h, Supports_halfedge_vertex());
    }

    Halfedge_handle get_prev( Halfedge_handle h) const {
        // returns the previous halfedge of h if supported,
        // `Halfedge_handle()' otherwise.
        return get_prev( h, Supports_halfedge_prev());
    }

    Halfedge_handle find_prev( Halfedge_handle h) const {
        // returns the previous halfedge of h. Uses the `prev()' method if
        // supported or performs a search around the face using `next()'.
        return find_prev( h, Supports_halfedge_prev());
    }

    Halfedge_handle find_prev_around_vertex( Halfedge_handle h) const {
        // returns the previous halfedge of h. Uses the `prev()' method if
        // supported or performs a search around the vertex using `next()'.
        return find_prev_around_vertex( h, Supports_halfedge_prev());
    }

    Face_handle get_face( Halfedge_handle h) const {
        // returns the incident face of h if supported,
        // `Face_handle()' otherwise.
        return get_face( h, Supports_halfedge_face());
    }

    Halfedge_handle get_face_halfedge( Face_handle f) const {
        // returns the incident halfedge of f if supported,
        // `Halfedge_handle()' otherwise.
        return get_face_halfedge( f, Supports_face_halfedge());
    }

// Const Access Functions
// ----------------------------------

    Halfedge_const_handle
    get_vertex_halfedge( Vertex_const_handle v) const {
        return get_vertex_halfedge( v, Supports_vertex_halfedge());
    }

    Vertex_const_handle get_vertex( Halfedge_const_handle h) const {
        return get_vertex( h, Supports_halfedge_vertex());
    }

    Halfedge_const_handle get_prev( Halfedge_const_handle h) const {
        return get_prev( h, Supports_halfedge_prev());
    }

    Halfedge_const_handle find_prev( Halfedge_const_handle h) const {
        return find_prev( h, Supports_halfedge_prev());
    }

    Halfedge_const_handle
    find_prev_around_vertex( Halfedge_const_handle h) const {
        return find_prev_around_vertex( h, Supports_halfedge_prev());
    }

    Face_const_handle get_face( Halfedge_const_handle h) const {
        return get_face( h, Supports_halfedge_face());
    }

    Halfedge_const_handle get_face_halfedge( Face_const_handle f) const {
        return get_face_halfedge( f, Supports_face_halfedge());
    }

// Modifying Functions (Primitives)
// ----------------------------------

    void set_vertex_halfedge( Vertex_handle v, Halfedge_handle g) const {
        // sets the incident halfedge of v to g.
        set_vertex_halfedge( v, g, Supports_vertex_halfedge());
    }

    void set_vertex_halfedge( Halfedge_handle h) const {
        // sets the incident halfedge of the vertex incident to h to h.
        set_vertex_halfedge( h, h, Supports_halfedge_vertex());
    }

    void set_vertex( Halfedge_handle h, Vertex_handle v) const {
        // sets the incident vertex of h to v.
        set_vertex(h, v, Supports_halfedge_vertex());
    }

    void set_prev( Halfedge_handle h, Halfedge_handle g) const {
        // sets the previous link of h to g.
        set_prev( h, g, Supports_halfedge_prev());
    }

    void set_face( Halfedge_handle h, Face_handle f) const {
        // sets the incident face of h to f.
        set_face(h, f, Supports_halfedge_face());
    }

    void set_face_halfedge( Face_handle f, Halfedge_handle g) const {
        // sets the incident halfedge of f to g.
        set_face_halfedge( f, g, Supports_face_halfedge());
    }

    void set_face_halfedge( Halfedge_handle h) const {
        // sets the incident halfedge of the face incident to h to h.
        set_face_halfedge( h, h, Supports_halfedge_face());
    }

// Modifying Functions (Composed)
// ----------------------------------

    void close_tip( Halfedge_handle h) const {
        // makes `h->opposite()' the successor of h.
        h->HBase::set_next( h->opposite());
        set_prev( h->opposite(), h);
    }

    void close_tip( Halfedge_handle h, Vertex_handle v) const {
        // makes `h->opposite()' the successor of h and sets the incident
        // vertex of h to v.
        h->HBase::set_next( h->opposite());
        set_prev( h->opposite(), h);
        set_vertex( h, v);
        set_vertex_halfedge( h);
    }

    void insert_tip( Halfedge_handle h, Halfedge_handle v) const {
        // inserts the tip of the edge h into the halfedges around the
        // vertex pointed to by v. Halfedge `h->opposite()' is the new
        // successor of v and `h->next()' will be set to `v->next()'. The
        // vertex of h will be set to the vertex v refers to if vertices
        // are supported.
        h->HBase::set_next( v->next());
        v->HBase::set_next( h->opposite());
        set_prev( h->next(), h);
        set_prev( h->opposite(), v);
        set_vertex( h, get_vertex( v));
    }

    void remove_tip( Halfedge_handle h) const {
        // removes the edge `h->next()->opposite()' from the halfedge
        // circle around the vertex referred to by h. The new successor
        // halfedge of h will be `h->next()->opposite()->next()'.
        h->HBase::set_next( h->next()->opposite()->next());
        set_prev( h->next(), h);
    }

    void insert_halfedge( Halfedge_handle h, Halfedge_handle f) const {
        // inserts the halfedge h between f and `f->next()'. The face of h
        // will be the one f refers to if faces are supported.
        h->HBase::set_next( f->next());
        f->HBase::set_next( h);
        set_prev( h->next(), h);
        set_prev( h, f);
        set_face( h, get_face( f));
    }

    void remove_halfedge( Halfedge_handle h) const {
        // removes edge `h->next()' from the halfedge circle around the
        // face referred to by h. The new successor of h will be
        // `h->next()->next()'.
        h->HBase::set_next( h->next()->next());
        set_prev( h->next(), h);
    }

    void set_vertex_in_vertex_loop( Halfedge_handle  , Vertex_handle  ,
                                    Tag_false) const {}
    void set_vertex_in_vertex_loop( Halfedge_handle h, Vertex_handle v,
                                    Tag_true) const {
        CGAL_assertion_code( std::size_t termination_count = 0;)
        Halfedge_handle end = h;
        do {
            CGAL_assertion( ++termination_count != 0);
            h->HBase::set_vertex( v);
            h = h->next()->opposite();
        } while ( h != end);
    }

    void
    set_vertex_in_vertex_loop( Halfedge_handle h, Vertex_handle v) const {
        // loops around the vertex incident to h and sets all vertex
        // pointers to v. Precondition: `h != Halfedge_handle()'.
        CGAL_precondition( h != Halfedge_handle());
        set_vertex_in_vertex_loop( h, v, Supports_halfedge_vertex());
    }

    void set_face_in_face_loop( Halfedge_handle  , Face_handle  ,
                                Tag_false) const {}
    void set_face_in_face_loop( Halfedge_handle h, Face_handle f,
                                Tag_true) const {
        CGAL_assertion_code( std::size_t termination_count = 0;)
        Halfedge_handle end = h;
        do {
            CGAL_assertion( ++termination_count != 0);
            h->HBase::set_face( f);
            h = h->next();
        } while ( h != end);
    }

    void set_face_in_face_loop( Halfedge_handle h, Face_handle f) const {
        // loops around the face incident to h and sets all face pointers
        // to f. Precondition: `h != Halfedge_handle()'.
        CGAL_precondition( h != Halfedge_handle());
        set_face_in_face_loop( h, f, Supports_halfedge_face());
    }

    Halfedge_handle flip_edge( Halfedge_handle h) const {
        // performs an edge flip. It returns h after
        // rotating the edge h one vertex in the direction of the face
        // orientation. Precondition: `h != Halfedge_handle()' and both
        // incident faces of h are triangles.
        CGAL_precondition( h != Halfedge_handle());
        CGAL_precondition( h == h->next()->next()->next());
        CGAL_precondition( h->opposite() ==
                    h->opposite()->next()->next()->next());
        Halfedge_handle hprev = h->next()->next();
        Halfedge_handle gprev = h->opposite()->next()->next();
        remove_tip( hprev);
        remove_tip( gprev);
        set_face_halfedge(  hprev);
        set_face_halfedge(  gprev);
        set_vertex_halfedge( hprev);
        set_vertex_halfedge( gprev);
        set_face( hprev->next(), hprev->face());
        set_face( gprev->next(), gprev->face());
        hprev = hprev->next();
        gprev = gprev->next();
        insert_tip( h, gprev);
        insert_tip( h->opposite(), hprev);
        CGAL_postcondition( h == h->next()->next()->next());
        CGAL_postcondition( h->opposite() ==
                     h->opposite()->next()->next()->next());
        return h;
    }

// Implementing These Functions.
// ====================================================
// Access Functions
// ----------------------------------

    Halfedge_handle
    get_vertex_halfedge( Vertex_handle  , Tag_false) const {
        return Halfedge_handle();
    }
    Halfedge_handle
    get_vertex_halfedge( Vertex_handle v, Tag_true) const {
        return v->halfedge();
    }

    Vertex_handle get_vertex( Halfedge_handle  , Tag_false) const {
        return Vertex_handle();
    }
    Vertex_handle get_vertex( Halfedge_handle h, Tag_true)  const {
        return h->vertex();
    }

    Halfedge_handle get_prev( Halfedge_handle  , Tag_false) const {
        return Halfedge_handle();
    }
    Halfedge_handle get_prev( Halfedge_handle h, Tag_true)  const {
        return h->HBase::prev();
    }

    Halfedge_handle find_prev( Halfedge_handle h, Tag_true) const {
        return h->HBase::prev();
    }
    Halfedge_handle find_prev( Halfedge_handle h, Tag_false) const {
        Halfedge_handle g = h;
        while ( g->next() != h)
            g = g->next();
        return g;
    }

    Halfedge_handle find_prev_around_vertex( Halfedge_handle h,
                                             Tag_true) const {
        return h->HBase::prev();
    }
    Halfedge_handle find_prev_around_vertex( Halfedge_handle h,
                                             Tag_false) const {
        Halfedge_handle g = h->opposite();
        while ( g->next() != h)
            g = g->next()->opposite();
        return g;
    }

    Face_handle get_face( Halfedge_handle  , Tag_false) const {
        return Face_handle();
    }
    Face_handle get_face( Halfedge_handle h, Tag_true)  const {
        return h->face();
    }

    Halfedge_handle get_face_halfedge( Face_handle  , Tag_false) const {
        return Halfedge_handle();
    }
    Halfedge_handle get_face_halfedge( Face_handle f, Tag_true) const {
        return f->halfedge();
    }

// Const Access Functions
// ----------------------------------

    Halfedge_const_handle
    get_vertex_halfedge( Vertex_const_handle  ,Tag_false) const {
        return Halfedge_const_handle();
    }
    Halfedge_const_handle
    get_vertex_halfedge( Vertex_const_handle v,Tag_true) const {
        return v->halfedge();
    }

    Vertex_const_handle
    get_vertex( Halfedge_const_handle  , Tag_false) const {
        return Vertex_const_handle();
    }
    Vertex_const_handle
    get_vertex( Halfedge_const_handle h, Tag_true)  const {
        return h->vertex();
    }

    Halfedge_const_handle
    get_prev( Halfedge_const_handle  , Tag_false) const {
        return Halfedge_const_handle();
    }
    Halfedge_const_handle
    get_prev( Halfedge_const_handle h, Tag_true)  const {
        return h->HBase::prev();
    }

    Halfedge_const_handle
    find_prev( Halfedge_const_handle h, Tag_true) const {
        return h->HBase::prev();
    }
    Halfedge_const_handle
    find_prev( Halfedge_const_handle h, Tag_false) const {
        Halfedge_const_handle g = h;
        while ( g->next() != h)
            g = g->next();
        return g;
    }

    Halfedge_const_handle
    find_prev_around_vertex( Halfedge_const_handle h, Tag_true) const {
        return h->HBase::prev();
    }
    Halfedge_const_handle
    find_prev_around_vertex( Halfedge_const_handle h, Tag_false) const {
        Halfedge_const_handle g = h->opposite();
        while ( g->next() != h)
            g = g->next()->opposite();
        return g;
    }

    Face_const_handle
    get_face( Halfedge_const_handle  , Tag_false) const {
        return Face_const_handle();
    }
    Face_const_handle
    get_face( Halfedge_const_handle h, Tag_true)  const {
        return h->face();
    }

    Halfedge_const_handle
    get_face_halfedge( Face_const_handle  , Tag_false) const {
        return Halfedge_const_handle();
    }
    Halfedge_const_handle
    get_face_halfedge( Face_const_handle f, Tag_true) const {
        return f->halfedge();
    }

// Modifying Function Primitives
// ----------------------------------

    void set_vertex_halfedge( Vertex_handle,
                              Halfedge_handle,
                              Tag_false) const {}
    void set_vertex_halfedge( Vertex_handle v,
                              Halfedge_handle g,
                              Tag_true)  const {
        v->VBase::set_halfedge(g);
    }

    void set_vertex_halfedge( Halfedge_handle,
                              Halfedge_handle,
                              Tag_false) const {}
    void set_vertex_halfedge( Halfedge_handle h,
                              Halfedge_handle g,
                              Tag_true) const {
        set_vertex_halfedge( h->vertex(), g);
    }

    void set_vertex( Halfedge_handle,   Vertex_handle,  Tag_false) const {}
    void set_vertex( Halfedge_handle h, Vertex_handle v, Tag_true) const {
        h->HBase::set_vertex(v);
    }

    void set_prev( Halfedge_handle,   Halfedge_handle,  Tag_false) const {}
    void set_prev( Halfedge_handle h, Halfedge_handle g, Tag_true) const {
        h->HBase::set_prev( g);
    }

    void set_face( Halfedge_handle,   Face_handle,  Tag_false) const {}
    void set_face( Halfedge_handle h, Face_handle f, Tag_true) const {
        h->HBase::set_face(f);
    }

    void set_face_halfedge( Face_handle,
                            Halfedge_handle,
                            Tag_false) const {}
    void set_face_halfedge( Face_handle f,
                            Halfedge_handle g,
                            Tag_true)  const {
        f->FBase::set_halfedge(g);
    }

    void set_face_halfedge( Halfedge_handle,
                            Halfedge_handle,
                            Tag_false) const {}
    void set_face_halfedge( Halfedge_handle h,
                            Halfedge_handle g,
                            Tag_true) const {
        set_face_halfedge( h->face(), g);
    }
};

} //namespace CGAL

#endif // CGAL_HALFEDGEDS_ITEMS_DECORATOR_H //
// EOF //
