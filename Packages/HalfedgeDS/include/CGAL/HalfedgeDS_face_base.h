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
// file          : HalfedgeDS_face_base.h
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: HalfedgeDS 3.3 (27 Sep 2000) $
// source        : hds_bases.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Halfedge Data Structure Base Classes for Faces.
// ============================================================================

#ifndef CGAL_HALFEDGEDS_FACE_BASE_H
#define CGAL_HALFEDGEDS_FACE_BASE_H 1

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION

template < class Refs, class T = Tag_true, class Trts = void>
class HalfedgeDS_face_base;

template < class Refs >
class HalfedgeDS_face_base< Refs, Tag_false, void> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_face_base< Refs, Tag_false, void>  Base;
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
    typedef Tag_false                            Supports_face_normal;
};

template < class Refs >
class HalfedgeDS_face_base< Refs, Tag_true, void> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_face_base< Refs, Tag_true, void>   Base;
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
    typedef Tag_false                            Supports_face_normal;
private:
    Halfedge_handle hdg;
public:
    Halfedge_handle       halfedge()                        { return hdg; }
    Halfedge_const_handle halfedge() const                  { return hdg; }
    void                  set_halfedge( Halfedge_handle h)  { hdg = h; }
};

template < class Refs, class Trts >
class HalfedgeDS_face_base< Refs, Tag_false, Trts> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_face_base< Refs, Tag_false, Trts>     Base;
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
    typedef Tag_true                             Supports_face_normal;
    typedef Trts                                 Traits;
    typedef typename Traits::Normal              Normal;
    typedef typename Traits::Plane               Plane;
private:
    Plane   pln;
public:
    HalfedgeDS_face_base() {}
    HalfedgeDS_face_base( const Plane& g) : pln(g) {}
    Plane&                plane()                           { return pln; }
    const Plane&          plane() const                     { return pln; }
    Normal                normal() const { return pln.orthogonal_vector();}
};

template < class Refs, class Trts >
class HalfedgeDS_face_base< Refs, Tag_true, Trts> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_face_base< Refs, Tag_true, Trts>      Base;
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
    typedef Tag_true                             Supports_face_normal;
    typedef Trts                                 Traits;
    typedef typename Traits::Normal              Normal;
    typedef typename Traits::Plane               Plane;
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
    Normal                normal() const { return pln.orthogonal_vector();}
};

#else // CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION //

// Partial specialization doesn't work. We can factor out the
// Plane parameter in a base class with full specialization
// on 'void', but we cannot get rid of the halfedge reference.
// So, we just waste the space and have it always.
//   Furthermore, is is likely to have a non-optimal memory
// price-tag for the base class as well if it is the empty base
// class for traits type 'void', since empty structs probably
// consume at least a byte, probably a word.

template < class Trts>
struct I_HalfedgeDS_face_base_plane {
    typedef Trts                                 Traits;
    typedef typename Traits::Plane               Plane;
    typedef typename Traits::Normal              Normal;
    Plane plane;
    I_HalfedgeDS_face_base_plane() {}
    I_HalfedgeDS_face_base_plane( const Plane& pl) : plane(pl) {}
    typedef Tag_true Supports_plane;
};
template <>
struct I_HalfedgeDS_face_base_plane<void> {
    typedef Tag_false Supports_plane;
    struct Normal_not_supported {};
    typedef Normal_not_supported Normal;
    struct Plane_not_supported {};
    typedef Plane_not_supported Plane;
};

template < class Refs, class T = Tag_true, class Trts = void>
class HalfedgeDS_face_base : public I_HalfedgeDS_face_base_plane<Trts> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef HalfedgeDS_face_base< Refs, T, Trts> Base;
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
    typedef I_HalfedgeDS_face_base_plane<Trts>   Plane_base;
    typedef typename Plane_base::Supports_plane  Supports_face_plane;
    typedef typename Plane_base::Supports_plane  Supports_face_normal;
    typedef Trts                                 Traits;
    typedef typename Plane_base::Normal          Normal;
    typedef typename Plane_base::Plane           Plane;
private:
    Halfedge_handle hdg;
public:
    HalfedgeDS_face_base() {}
    HalfedgeDS_face_base( const Plane& g)
        : I_HalfedgeDS_face_base_plane<Trts>(g) {}
    Halfedge_handle       halfedge()                      { return hdg; }
    Halfedge_const_handle halfedge() const                { return hdg; }
    void                  set_halfedge( Halfedge_handle h){ hdg = h; }
    Plane&                plane()        { return Plane_base::plane; }
    const Plane&          plane()  const { return Plane_base::plane; }
    Normal                normal() const {
                              return Plane_base::plane.orthogonal_vector();
    }
};

#endif // CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION //

CGAL_END_NAMESPACE

#endif // CGAL_HALFEDGEDS_FACE_BASE_H //
// EOF //
