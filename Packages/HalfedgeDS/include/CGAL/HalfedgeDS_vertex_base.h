// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_HALFEDGEDS_VERTEX_BASE_H
#define CGAL_HALFEDGEDS_VERTEX_BASE_H 1

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

CGAL_BEGIN_NAMESPACE


#if ! defined( CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION) \
 && ! defined( CGAL_CFG_SUNPRO_PARTIAL_SPEC_BUG)

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

#else //    ! defined( CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION)
      // && ! defined( CGAL_CFG_SUNPRO_PARTIAL_SPEC_BUG)

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

#endif //    ! defined( CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION)
       // && ! defined( CGAL_CFG_SUNPRO_PARTIAL_SPEC_BUG)

CGAL_END_NAMESPACE

#endif // CGAL_HALFEDGEDS_VERTEX_BASE_H //
// EOF //
