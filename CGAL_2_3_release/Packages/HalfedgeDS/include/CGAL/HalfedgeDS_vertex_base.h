// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : HalfedgeDS_vertex_base.h
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: HalfedgeDS 3.3 (27 Sep 2000) $
// source        : hds_bases.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Halfedge Data Structure Base Class for Vertices.
// ============================================================================

#ifndef CGAL_HALFEDGEDS_VERTEX_BASE_H
#define CGAL_HALFEDGEDS_VERTEX_BASE_H 1

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

CGAL_BEGIN_NAMESPACE


#ifndef CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION

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

#else // CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION //

// Partial specialization doesn't work. We can factor out the
// Point parameter in a base class with full specialization
// on 'Tag_false', but we cannot get rid of the halfedge reference.
// So, we just waste the space and have it always.
//   Furthermore, it is likely to have a non-optimal memory
// price-tag for the base class as well if it is the empty base
// class for point type 'Tag_false', since empty structs probably
// consume at least a byte, probably a word.
//   See HalfedgeDS_face_min_base.h for an alternative.

// We use Tag_false to indicate that no point type is provided.

template <class Pt>
struct I_HalfedgeDS_vertex_base_point {
    Pt point;
    I_HalfedgeDS_vertex_base_point() {}
    I_HalfedgeDS_vertex_base_point( const Pt& pt) : point(pt) {}
    typedef Tag_true Supports_point;
    typedef Pt Point;
};
template <>
struct I_HalfedgeDS_vertex_base_point<Tag_false> {
    typedef Tag_false Supports_point;
    struct Point_not_supported {};
    typedef Point_not_supported Point;
};

template < class Refs, class T = Tag_true, class P = Tag_false>
class HalfedgeDS_vertex_base : public I_HalfedgeDS_vertex_base_point<P> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_vertex_base< Refs, T, P>  Base;
    typedef T                                    Supports_vertex_halfedge;
    typedef I_HalfedgeDS_vertex_base_point<P>    Point_base;
    typedef typename Point_base::Supports_point  Supports_vertex_point;
    typedef typename Point_base::Point           Point;
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
    HalfedgeDS_vertex_base() {}
    HalfedgeDS_vertex_base( const Point& pp)
        : I_HalfedgeDS_vertex_base_point<P>(pp) {}
    Halfedge_handle       halfedge()                      { return hdg; }
    Halfedge_const_handle halfedge() const                { return hdg; }
    void                  set_halfedge( Halfedge_handle h){ hdg = h; }
    Point&                point()           { return Point_base::point; }
    const Point&          point() const     { return Point_base::point; }
};

#endif // CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION //

CGAL_END_NAMESPACE

#endif // CGAL_HALFEDGEDS_VERTEX_BASE_H //
// EOF //
