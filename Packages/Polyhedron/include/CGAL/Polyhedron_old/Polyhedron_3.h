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
// file          : Polyhedron_3.h
// chapter       : $CGAL_Chapter: 3D-Polyhedral Surfaces $
// package       : $CGAL_Package: Polyhedron 2.9 (13 Sep 2000) $
// source        : polyhedron.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Polyhedral Surfaces.
// ============================================================================

#ifndef CGAL_POLYHEDRON_OLD_POLYHEDRON_3_H
#define CGAL_POLYHEDRON_OLD_POLYHEDRON_3_H 1

#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Halfedge_data_structure_decorator.h>

#ifdef CGAL_REP_CLASS_DEFINED
#include <CGAL/Polyhedron_default_traits_3.h>
#endif // CGAL_REP_CLASS_DEFINED

#include <CGAL/Polyhedron_iterator_3.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>
#include <CGAL/N_step_adaptor_derived.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/IO/Verbose_ostream.h>

CGAL_BEGIN_NAMESPACE

// Forward declaration of the three element classes.
template <class HDS> class _Polyhedron_vertex;
template <class HDS> class _Polyhedron_halfedge;
template <class HDS> class _Polyhedron_facet;


// A halfedge is an oriented edge. Both orientations exist, i.e.
// an edge is represented by two opposite halfedges. The geometric
// relations are as follows:
//
//              _ _ _   .
//             /        |\.
//                      | \.
//           /             \ next half
//                          \ edge
//         /                 \.
//
//        |                   O  incident vertex
//                facet      ,
//        |                 /| |
//                         / | | opposite
//         \                 | | half edge
//                      half | |
//           \          edge | | /
//                           | |/
//             \_ _ _ _ _ _    '
//

template <class HDS>
class _Polyhedron_halfedge : public HDS::Halfedge {
    public:
        typedef HDS                               Halfedge_data_structure;
    
    protected:
        typedef typename HDS::Vertex              V;
        typedef typename HDS::Halfedge            H;
        typedef typename HDS::Facet               F;
    
    public:
        typedef _Polyhedron_vertex<HDS>           Vertex_;
        typedef _Polyhedron_halfedge<HDS>         Halfedge_;
        typedef _Polyhedron_facet<HDS>            Facet_;
    
        typedef typename  HDS::Size               Size;
        typedef typename  HDS::Difference         Difference;
    
    // The following types denotes the (optionally) associated geometry. If
    // the specific item is not supported the type is `void*'.
    
        typedef typename  HDS::Point              Point;
        typedef typename  HDS::Point              Point_3;
        typedef typename  F::Vector_3             Vector_3;
        typedef typename  F::Plane_3              Plane_3;
    
    // The following types are equal to either `Tag_true' or
    // `Tag_false', dependant whether the named feature is
    // supported or not.
    
        typedef typename  HDS::Supports_vertex_halfedge
                                                      Supports_vertex_halfedge;
        typedef typename  HDS::Supports_halfedge_prev
                                                      Supports_halfedge_prev;
        typedef typename  HDS::Supports_halfedge_vertex
                                                      Supports_halfedge_vertex;
        typedef typename  HDS::Supports_halfedge_facet
                                                      Supports_halfedge_facet;
        typedef typename  HDS::Supports_facet_halfedge
                                                      Supports_facet_halfedge;
    
        typedef typename  HDS::Supports_vertex_point  Supports_vertex_point;
        typedef typename  HDS::Supports_facet_plane   Supports_facet_plane;
        typedef typename  HDS::Supports_facet_normal  Supports_facet_normal;
    
        typedef typename  HDS::Supports_removal       Supports_removal;
    
    // The iterator/circulator categories.
    
        typedef typename  HDS::iterator_category      iterator_category;
        typedef Polyhedron_circulator_traits<Supports_halfedge_prev>
                                                      Circ_traits;
        typedef typename Circ_traits::iterator_category
                                                      circulator_category;
    
    protected:
        // These extra (internal) typedefs are necessary to make
        // SunPro CC 4.2 happy. (And they are used.)
        typedef typename  HDS::Vertex_iterator          TR_VI;
        typedef typename  HDS::Vertex_const_iterator    TR_C_VI;
        typedef typename  HDS::Halfedge_iterator        TR_HI;
        typedef typename  HDS::Halfedge_const_iterator  TR_C_HI;
        typedef typename  HDS::Edge_iterator            TR_EI;
        typedef typename  HDS::Edge_const_iterator      TR_C_EI;
        typedef typename  HDS::Facet_iterator           TR_FI;
        typedef typename  HDS::Facet_const_iterator     TR_C_FI;
    
    public:
        typedef _Polyhedron_iterator<
            TR_VI,
            Vertex_,
            Difference, iterator_category>       Vertex_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_VI, TR_VI,
            Vertex_,
            Difference, iterator_category>       Vertex_const_iterator;
    
        typedef _Polyhedron_iterator<
            TR_HI,
            Halfedge_,
            Difference, iterator_category>       Halfedge_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_HI, TR_HI,
            Halfedge_,
            Difference, iterator_category>       Halfedge_const_iterator;
    
        typedef _Polyhedron_iterator<
            TR_FI,
            Facet_,
            Difference, iterator_category>       Facet_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_FI, TR_FI,
            Facet_,
            Difference, iterator_category>       Facet_const_iterator;
    
    // The circulators around a vertex or around a facet.
    
        typedef _Polyhedron_facet_circ<
            Halfedge_,
            Halfedge_iterator,
            circulator_category>            Halfedge_around_facet_circulator;
    
        typedef _Polyhedron_vertex_circ<
            Halfedge_,
            Halfedge_iterator,
            circulator_category>            Halfedge_around_vertex_circulator;
    
        typedef _Polyhedron_facet_const_circ<
            Halfedge_,
            Halfedge_const_iterator,
            circulator_category>       Halfedge_around_facet_const_circulator;
    
        typedef _Polyhedron_vertex_const_circ<
            Halfedge_,
            Halfedge_const_iterator,
            circulator_category>      Halfedge_around_vertex_const_circulator;
    
    // The handles. They are currently only simple typedef's.
        typedef Vertex_iterator                       Vertex_handle;
        typedef Vertex_const_iterator                 Vertex_const_handle;
        typedef Halfedge_iterator                     Halfedge_handle;
        typedef Halfedge_const_iterator               Halfedge_const_handle;
        typedef Facet_iterator                        Facet_handle;
        typedef Facet_const_iterator                  Facet_const_handle;
    
    // Edge iterator.
    
        typedef N_step_adaptor_derived<Halfedge_iterator, 2>
                                                      Edge_iterator;
        typedef N_step_adaptor_derived<Halfedge_const_iterator, 2>
                                                      Edge_const_iterator;

public:

    typedef _Polyhedron_vertex<HDS>       Vertex;
    typedef _Polyhedron_facet<HDS>        Facet;

    _Polyhedron_halfedge() {}
    _Polyhedron_halfedge( const H& h) : H(h) {}

    Halfedge_handle opposite() { return TR_HI( H::opposite());}
    Halfedge_handle next()     { return TR_HI( H::next());}
    Halfedge_handle prev()     {
        // returns the halfedge previous to this. Uses the `prev()'
        // method if available or performs a search around the facet.
        Halfedge_data_structure_decorator<HDS> decorator;
        return TR_HI(decorator.find_prev((H*)(this)));
    }
    Vertex_handle   vertex()   { return TR_VI( H::vertex());}
    Facet_handle    facet()    { return TR_FI( H::facet());}

    Halfedge_const_handle opposite() const {
                                return TR_C_HI( H::opposite());
    }
    Halfedge_const_handle next()     const {
                                return TR_C_HI( H::next());
    }
    Halfedge_const_handle prev()     const {
        Halfedge_data_structure_decorator<HDS> D;
        return TR_C_HI( D.find_prev(this));
    }
    Vertex_const_handle   vertex()   const {
                                return TR_C_VI( H::vertex());
    }
    Facet_const_handle    facet()    const {
                                return TR_C_FI(H::facet());
    }

// Derived Access Functions (not present in H).

    Halfedge_handle       next_on_vertex()       {
                                return next()->opposite();
    }
    Halfedge_const_handle next_on_vertex() const {
                                return next()->opposite();
    }
        // the next halfedge around the vertex (clockwise). Is equal to
        // `h.next()->opposite()'.

    Halfedge_handle       prev_on_vertex()       {
                                return opposite()->prev();
    }
    Halfedge_const_handle prev_on_vertex() const {
                                return opposite()->prev();
    }
        // the previous halfedge around the vertex (counterclockwise). Is
        // equal to `h.opposite()->prev()'.

    bool is_border_edge() const {
        // is true if `h' or `h.opposite()' is a border halfedge.
        return (opposite()->is_border() || is_border());
    }

    Halfedge_around_vertex_circulator vertex_begin() {
        // a circulator of halfedges around the vertex (clockwise).
        return Halfedge_around_vertex_circulator(TR_HI(this));
    }
    Halfedge_around_vertex_const_circulator vertex_begin() const {
        // a circulator of halfedges around the vertex (clockwise).
        return Halfedge_around_vertex_const_circulator(TR_C_HI(this));
    }

    Halfedge_around_facet_circulator facet_begin() {
        // a circulator of halfedges around the facet (counterclockwise).
        return Halfedge_around_facet_circulator(TR_HI(this));
    }
    Halfedge_around_facet_const_circulator facet_begin() const {
        // a circulator of halfedges around the facet (counterclockwise).
        return Halfedge_around_facet_const_circulator(TR_C_HI(this));
    }

private:
    // Hide some other functions of H.
    void  set_next( Halfedge_* h)  { H::set_next(h);}
    void  set_prev( Halfedge_* h)  { H::set_prev(h);}
    void  set_vertex( Vertex* v)   { H::set_vertex(v);}
    void  set_facet( Facet* f)     { H::set_facet(f);}
};


template <class HDS>
class _Polyhedron_vertex  : public HDS::Vertex  {
    public:
        typedef HDS                               Halfedge_data_structure;
    
    protected:
        typedef typename HDS::Vertex              V;
        typedef typename HDS::Halfedge            H;
        typedef typename HDS::Facet               F;
    
    public:
        typedef _Polyhedron_vertex<HDS>           Vertex_;
        typedef _Polyhedron_halfedge<HDS>         Halfedge_;
        typedef _Polyhedron_facet<HDS>            Facet_;
    
        typedef typename  HDS::Size               Size;
        typedef typename  HDS::Difference         Difference;
    
    // The following types denotes the (optionally) associated geometry. If
    // the specific item is not supported the type is `void*'.
    
        typedef typename  HDS::Point              Point;
        typedef typename  HDS::Point              Point_3;
        typedef typename  F::Vector_3             Vector_3;
        typedef typename  F::Plane_3              Plane_3;
    
    // The following types are equal to either `Tag_true' or
    // `Tag_false', dependant whether the named feature is
    // supported or not.
    
        typedef typename  HDS::Supports_vertex_halfedge
                                                      Supports_vertex_halfedge;
        typedef typename  HDS::Supports_halfedge_prev
                                                      Supports_halfedge_prev;
        typedef typename  HDS::Supports_halfedge_vertex
                                                      Supports_halfedge_vertex;
        typedef typename  HDS::Supports_halfedge_facet
                                                      Supports_halfedge_facet;
        typedef typename  HDS::Supports_facet_halfedge
                                                      Supports_facet_halfedge;
    
        typedef typename  HDS::Supports_vertex_point  Supports_vertex_point;
        typedef typename  HDS::Supports_facet_plane   Supports_facet_plane;
        typedef typename  HDS::Supports_facet_normal  Supports_facet_normal;
    
        typedef typename  HDS::Supports_removal       Supports_removal;
    
    // The iterator/circulator categories.
    
        typedef typename  HDS::iterator_category      iterator_category;
        typedef Polyhedron_circulator_traits<Supports_halfedge_prev>
                                                      Circ_traits;
        typedef typename Circ_traits::iterator_category
                                                      circulator_category;
    
    protected:
        // These extra (internal) typedefs are necessary to make
        // SunPro CC 4.2 happy. (And they are used.)
        typedef typename  HDS::Vertex_iterator          TR_VI;
        typedef typename  HDS::Vertex_const_iterator    TR_C_VI;
        typedef typename  HDS::Halfedge_iterator        TR_HI;
        typedef typename  HDS::Halfedge_const_iterator  TR_C_HI;
        typedef typename  HDS::Edge_iterator            TR_EI;
        typedef typename  HDS::Edge_const_iterator      TR_C_EI;
        typedef typename  HDS::Facet_iterator           TR_FI;
        typedef typename  HDS::Facet_const_iterator     TR_C_FI;
    
    public:
        typedef _Polyhedron_iterator<
            TR_VI,
            Vertex_,
            Difference, iterator_category>       Vertex_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_VI, TR_VI,
            Vertex_,
            Difference, iterator_category>       Vertex_const_iterator;
    
        typedef _Polyhedron_iterator<
            TR_HI,
            Halfedge_,
            Difference, iterator_category>       Halfedge_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_HI, TR_HI,
            Halfedge_,
            Difference, iterator_category>       Halfedge_const_iterator;
    
        typedef _Polyhedron_iterator<
            TR_FI,
            Facet_,
            Difference, iterator_category>       Facet_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_FI, TR_FI,
            Facet_,
            Difference, iterator_category>       Facet_const_iterator;
    
    // The circulators around a vertex or around a facet.
    
        typedef _Polyhedron_facet_circ<
            Halfedge_,
            Halfedge_iterator,
            circulator_category>            Halfedge_around_facet_circulator;
    
        typedef _Polyhedron_vertex_circ<
            Halfedge_,
            Halfedge_iterator,
            circulator_category>            Halfedge_around_vertex_circulator;
    
        typedef _Polyhedron_facet_const_circ<
            Halfedge_,
            Halfedge_const_iterator,
            circulator_category>       Halfedge_around_facet_const_circulator;
    
        typedef _Polyhedron_vertex_const_circ<
            Halfedge_,
            Halfedge_const_iterator,
            circulator_category>      Halfedge_around_vertex_const_circulator;
    
    // The handles. They are currently only simple typedef's.
        typedef Vertex_iterator                       Vertex_handle;
        typedef Vertex_const_iterator                 Vertex_const_handle;
        typedef Halfedge_iterator                     Halfedge_handle;
        typedef Halfedge_const_iterator               Halfedge_const_handle;
        typedef Facet_iterator                        Facet_handle;
        typedef Facet_const_iterator                  Facet_const_handle;
    
    // Edge iterator.
    
        typedef N_step_adaptor_derived<Halfedge_iterator, 2>
                                                      Edge_iterator;
        typedef N_step_adaptor_derived<Halfedge_const_iterator, 2>
                                                      Edge_const_iterator;

public:

    typedef _Polyhedron_halfedge<HDS>     Halfedge;
    typedef _Polyhedron_facet<HDS>        Facet;

    _Polyhedron_vertex() {}
    _Polyhedron_vertex( const V& v) : V(v) {}

    Halfedge_handle       halfedge() {return TR_HI( V::halfedge());}
    Halfedge_const_handle halfedge() const {
                                return TR_C_HI( V::halfedge());
    }

    // Avoids unnecessary matchings with base class. (g++ bug)
    Point&          point()       { return V::point();}
    const Point&    point() const { return V::point();}

// Derived Access Functions (not present in V).

    Halfedge_around_vertex_circulator vertex_begin() {
        // a circulator of halfedges around the vertex (clockwise).
        return Halfedge_around_vertex_circulator( TR_HI(halfedge().ptr()));
    }
    Halfedge_around_vertex_const_circulator vertex_begin() const {
        // a circulator of halfedges around the vertex (clockwise).
        return Halfedge_around_vertex_const_circulator(
                   TR_C_HI(halfedge().ptr()));
    }
private:
    // Hide some other functions of V.
    void      set_halfedge( Halfedge* h) { V::set_halfedge(h);}
};


template <class HDS>
class _Polyhedron_facet : public HDS::Facet {
    public:
        typedef HDS                               Halfedge_data_structure;
    
    protected:
        typedef typename HDS::Vertex              V;
        typedef typename HDS::Halfedge            H;
        typedef typename HDS::Facet               F;
    
    public:
        typedef _Polyhedron_vertex<HDS>           Vertex_;
        typedef _Polyhedron_halfedge<HDS>         Halfedge_;
        typedef _Polyhedron_facet<HDS>            Facet_;
    
        typedef typename  HDS::Size               Size;
        typedef typename  HDS::Difference         Difference;
    
    // The following types denotes the (optionally) associated geometry. If
    // the specific item is not supported the type is `void*'.
    
        typedef typename  HDS::Point              Point;
        typedef typename  HDS::Point              Point_3;
        typedef typename  F::Vector_3             Vector_3;
        typedef typename  F::Plane_3              Plane_3;
    
    // The following types are equal to either `Tag_true' or
    // `Tag_false', dependant whether the named feature is
    // supported or not.
    
        typedef typename  HDS::Supports_vertex_halfedge
                                                      Supports_vertex_halfedge;
        typedef typename  HDS::Supports_halfedge_prev
                                                      Supports_halfedge_prev;
        typedef typename  HDS::Supports_halfedge_vertex
                                                      Supports_halfedge_vertex;
        typedef typename  HDS::Supports_halfedge_facet
                                                      Supports_halfedge_facet;
        typedef typename  HDS::Supports_facet_halfedge
                                                      Supports_facet_halfedge;
    
        typedef typename  HDS::Supports_vertex_point  Supports_vertex_point;
        typedef typename  HDS::Supports_facet_plane   Supports_facet_plane;
        typedef typename  HDS::Supports_facet_normal  Supports_facet_normal;
    
        typedef typename  HDS::Supports_removal       Supports_removal;
    
    // The iterator/circulator categories.
    
        typedef typename  HDS::iterator_category      iterator_category;
        typedef Polyhedron_circulator_traits<Supports_halfedge_prev>
                                                      Circ_traits;
        typedef typename Circ_traits::iterator_category
                                                      circulator_category;
    
    protected:
        // These extra (internal) typedefs are necessary to make
        // SunPro CC 4.2 happy. (And they are used.)
        typedef typename  HDS::Vertex_iterator          TR_VI;
        typedef typename  HDS::Vertex_const_iterator    TR_C_VI;
        typedef typename  HDS::Halfedge_iterator        TR_HI;
        typedef typename  HDS::Halfedge_const_iterator  TR_C_HI;
        typedef typename  HDS::Edge_iterator            TR_EI;
        typedef typename  HDS::Edge_const_iterator      TR_C_EI;
        typedef typename  HDS::Facet_iterator           TR_FI;
        typedef typename  HDS::Facet_const_iterator     TR_C_FI;
    
    public:
        typedef _Polyhedron_iterator<
            TR_VI,
            Vertex_,
            Difference, iterator_category>       Vertex_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_VI, TR_VI,
            Vertex_,
            Difference, iterator_category>       Vertex_const_iterator;
    
        typedef _Polyhedron_iterator<
            TR_HI,
            Halfedge_,
            Difference, iterator_category>       Halfedge_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_HI, TR_HI,
            Halfedge_,
            Difference, iterator_category>       Halfedge_const_iterator;
    
        typedef _Polyhedron_iterator<
            TR_FI,
            Facet_,
            Difference, iterator_category>       Facet_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_FI, TR_FI,
            Facet_,
            Difference, iterator_category>       Facet_const_iterator;
    
    // The circulators around a vertex or around a facet.
    
        typedef _Polyhedron_facet_circ<
            Halfedge_,
            Halfedge_iterator,
            circulator_category>            Halfedge_around_facet_circulator;
    
        typedef _Polyhedron_vertex_circ<
            Halfedge_,
            Halfedge_iterator,
            circulator_category>            Halfedge_around_vertex_circulator;
    
        typedef _Polyhedron_facet_const_circ<
            Halfedge_,
            Halfedge_const_iterator,
            circulator_category>       Halfedge_around_facet_const_circulator;
    
        typedef _Polyhedron_vertex_const_circ<
            Halfedge_,
            Halfedge_const_iterator,
            circulator_category>      Halfedge_around_vertex_const_circulator;
    
    // The handles. They are currently only simple typedef's.
        typedef Vertex_iterator                       Vertex_handle;
        typedef Vertex_const_iterator                 Vertex_const_handle;
        typedef Halfedge_iterator                     Halfedge_handle;
        typedef Halfedge_const_iterator               Halfedge_const_handle;
        typedef Facet_iterator                        Facet_handle;
        typedef Facet_const_iterator                  Facet_const_handle;
    
    // Edge iterator.
    
        typedef N_step_adaptor_derived<Halfedge_iterator, 2>
                                                      Edge_iterator;
        typedef N_step_adaptor_derived<Halfedge_const_iterator, 2>
                                                      Edge_const_iterator;

public:

    typedef _Polyhedron_vertex<HDS>       Vertex;
    typedef _Polyhedron_halfedge<HDS>     Halfedge;

    _Polyhedron_facet() {}
    _Polyhedron_facet( const F& f) : F(f) {}

    Halfedge_handle       halfedge() {return TR_HI( F::halfedge());}
    Halfedge_const_handle halfedge()  const {
                                return TR_C_HI( F::halfedge());
    }

    // Avoids unnecessary matchings with base class. (g++ bug)
    Vector_3        normal() const { return F::normal();}
    Plane_3&        plane()        { return F::plane();}
    const Plane_3&  plane() const  { return F::plane();}

// Derived Access Functions (not present in F).

    Halfedge_around_facet_circulator facet_begin() {
        // a circulator of halfedges around the facet (counterclockwise).
        return Halfedge_around_facet_circulator( TR_HI(halfedge().ptr()));
    }
    Halfedge_around_facet_const_circulator facet_begin() const {
        // a circulator of halfedges around the facet (counterclockwise).
        return Halfedge_around_facet_const_circulator(
                    TR_C_HI(halfedge().ptr()));
    }

private:
    // Hide some other functions of F.
    void      set_halfedge( Halfedge* h) { F::set_halfedge(h);}
};



template < class TR, class HDS 
    = CGAL::Halfedge_data_structure_polyhedron_default_3< TR> >
class Polyhedron_3 {
    //
    // DEFINITION
    //
    // The boundary representation of a 3d-polyhedron P of the type
    // Polyhedron<HDS> consists of vertices, edges and facets. The
    // vertices are points in space. The edges are straight line
    // segments. The facets are planar polygons. We restrict here
    // the facets to be simple planar polygons without holes and the
    // boundary of the polyhedron to be an oriented 2-manifold. Thus
    // facets are consistently oriented and an edge is incident to
    // exactly two facets. We restrict the representation further
    // that an edge has two distinct incident endpoints and
    // following duality that an edge has two distinct incident
    // facets. The class Polyhedron<HDS> is able to guarantee
    // the combinatorial properties, but not all geometric
    // properties. Support functions are provided for testing
    // geometric properties, e.g. test for self intersections which
    // is  too expensive to be guaranteed as a class invariant.

    public:
        typedef HDS                               Halfedge_data_structure;
        typedef HDS                               HalfedgeDS;
    
    protected:
        typedef typename HDS::Vertex              V;
        typedef typename HDS::Halfedge            H;
        typedef typename HDS::Facet               F;
    
    public:
        typedef _Polyhedron_vertex<HDS>           Vertex_;
        typedef _Polyhedron_halfedge<HDS>         Halfedge_;
        typedef _Polyhedron_facet<HDS>            Facet_;
    
        typedef typename  HDS::Size               Size;
        typedef typename  HDS::Difference         Difference;
    
    // The following types denotes the (optionally) associated geometry. If
    // the specific item is not supported the type is `void*'.
    
        typedef typename  HDS::Point              Point;
        typedef typename  HDS::Point              Point_3;
        typedef typename  F::Vector_3             Vector_3;
        typedef typename  F::Plane_3              Plane_3;
    
    // The following types are equal to either `Tag_true' or
    // `Tag_false', dependant whether the named feature is
    // supported or not.
    
        typedef typename  HDS::Supports_vertex_halfedge
                                                      Supports_vertex_halfedge;
        typedef typename  HDS::Supports_halfedge_prev
                                                      Supports_halfedge_prev;
        typedef typename  HDS::Supports_halfedge_vertex
                                                      Supports_halfedge_vertex;
        typedef typename  HDS::Supports_halfedge_facet
                                                      Supports_halfedge_facet;
        typedef typename  HDS::Supports_facet_halfedge
                                                      Supports_facet_halfedge;
    
        typedef typename  HDS::Supports_vertex_point  Supports_vertex_point;
        typedef typename  HDS::Supports_facet_plane   Supports_facet_plane;
        typedef typename  HDS::Supports_facet_normal  Supports_facet_normal;
    
        typedef typename  HDS::Supports_removal       Supports_removal;
    
    // The iterator/circulator categories.
    
        typedef typename  HDS::iterator_category      iterator_category;
        typedef Polyhedron_circulator_traits<Supports_halfedge_prev>
                                                      Circ_traits;
        typedef typename Circ_traits::iterator_category
                                                      circulator_category;
    
    protected:
        // These extra (internal) typedefs are necessary to make
        // SunPro CC 4.2 happy. (And they are used.)
        typedef typename  HDS::Vertex_iterator          TR_VI;
        typedef typename  HDS::Vertex_const_iterator    TR_C_VI;
        typedef typename  HDS::Halfedge_iterator        TR_HI;
        typedef typename  HDS::Halfedge_const_iterator  TR_C_HI;
        typedef typename  HDS::Edge_iterator            TR_EI;
        typedef typename  HDS::Edge_const_iterator      TR_C_EI;
        typedef typename  HDS::Facet_iterator           TR_FI;
        typedef typename  HDS::Facet_const_iterator     TR_C_FI;
    
    public:
        typedef _Polyhedron_iterator<
            TR_VI,
            Vertex_,
            Difference, iterator_category>       Vertex_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_VI, TR_VI,
            Vertex_,
            Difference, iterator_category>       Vertex_const_iterator;
    
        typedef _Polyhedron_iterator<
            TR_HI,
            Halfedge_,
            Difference, iterator_category>       Halfedge_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_HI, TR_HI,
            Halfedge_,
            Difference, iterator_category>       Halfedge_const_iterator;
    
        typedef _Polyhedron_iterator<
            TR_FI,
            Facet_,
            Difference, iterator_category>       Facet_iterator;
    
        typedef _Polyhedron_const_iterator<
            TR_C_FI, TR_FI,
            Facet_,
            Difference, iterator_category>       Facet_const_iterator;
    
    // The circulators around a vertex or around a facet.
    
        typedef _Polyhedron_facet_circ<
            Halfedge_,
            Halfedge_iterator,
            circulator_category>            Halfedge_around_facet_circulator;
    
        typedef _Polyhedron_vertex_circ<
            Halfedge_,
            Halfedge_iterator,
            circulator_category>            Halfedge_around_vertex_circulator;
    
        typedef _Polyhedron_facet_const_circ<
            Halfedge_,
            Halfedge_const_iterator,
            circulator_category>       Halfedge_around_facet_const_circulator;
    
        typedef _Polyhedron_vertex_const_circ<
            Halfedge_,
            Halfedge_const_iterator,
            circulator_category>      Halfedge_around_vertex_const_circulator;
    
    // The handles. They are currently only simple typedef's.
        typedef Vertex_iterator                       Vertex_handle;
        typedef Vertex_const_iterator                 Vertex_const_handle;
        typedef Halfedge_iterator                     Halfedge_handle;
        typedef Halfedge_const_iterator               Halfedge_const_handle;
        typedef Facet_iterator                        Facet_handle;
        typedef Facet_const_iterator                  Facet_const_handle;
    
    // Edge iterator.
    
        typedef N_step_adaptor_derived<Halfedge_iterator, 2>
                                                      Edge_iterator;
        typedef N_step_adaptor_derived<Halfedge_const_iterator, 2>
                                                      Edge_const_iterator;

    typedef TR Traits;

public:
    typedef _Polyhedron_vertex<HDS>       Vertex;
    typedef _Polyhedron_halfedge<HDS>     Halfedge;
    typedef _Polyhedron_facet<HDS>        Facet;

protected:
    HDS hds;  // the boundary representation.
    TR  m_traits;

public:
    typedef Polyhedron_3<TR,HDS> Self;

// CREATION

    Polyhedron_3( const Traits& traits = Traits())
    : m_traits(traits) {
        // the empty polyhedron `P'.
        typedef typename Traits::Point_3  Traits_point;
        assert_equal_types( Point(), Traits_point());
    }
    Polyhedron_3( Size v, Size h, Size f, const Traits&
                      traits = Traits())
    : hds(v,h,f), m_traits(traits) {
        // a polyhedron `P' with storage reserved for v vertices, h
        // halfedges, and f facets. The reservation sizes are a hint for
        // optimizing storage allocation.
        typedef typename Traits::Point_3  Traits_point;
        assert_equal_types( Point(), Traits_point());
    }

    //Polyhedron_3( const Self& p);
        // copy constructor.

    //Self& operator= ( const Self& p);
        // assignment operator.

    // Explicit implementation to help EGCS in compiling ...
    Self& operator= ( const Self& p) {
        hds = p.hds;
        m_traits = p.m_traits;
        return *this;
    }

    void reserve( Size v, Size h, Size f) { hds.reserve(v,h,f);}
        // reserve storage for v vertices, h halfedges, and f facets. The
        // reservation sizes are a hint for optimizing storage allocation.
        // If the `capacity' is already greater than the requested size
        // nothing happens. If the `capacity' changes all iterators and
        // circulators invalidates.

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
protected:
    Halfedge_handle make_tetrahedron( V* v1, V* v2, V* v3, V* v4);

    Halfedge_handle make_triangle( V* v1, V* v2, V* v3);

public:
    Halfedge_handle make_tetrahedron();
        // the combinatorial structure of a tetrahedron is added to the
        // actual polyhedral surface. Returns an arbitrary halfedge of
        // this structure.

    Halfedge_handle make_tetrahedron(const Point& p1,
                                     const Point& p2,
                                     const Point& p3,
                                     const Point& p4);

    Halfedge_handle make_triangle();
        // the combinatorial structure of a single triangle with border
        // edges is added to the actual polyhedral surface. Returns an
        // arbitrary halfedge of this structure.

    Halfedge_handle make_triangle(const Point& p1,
                                  const Point& p2,
                                  const Point& p3);
        // the single triangle p_1, p_2, p_3 with border edges is added to
        // the actual polyhedral surface. Returns an arbitrary halfedge of
        // this structure.
#endif // CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS //

// Access Member Functions

    Size size_of_vertices() const { return hds.size_of_vertices();}
        // number of vertices.

    Size size_of_halfedges() const { return hds.size_of_halfedges();}
        // number of all halfedges (including border halfedges).

    Size size_of_facets() const { return hds.size_of_facets();}
        // number of facets.

    Size capacity_of_vertices() const {return hds.capacity_of_vertices();}
        // space reserved for vertices.

    Size capacity_of_halfedges() const {
        // space reserved for halfedges.
        return hds.capacity_of_halfedges();
    }

    Size capacity_of_facets() const { return hds.capacity_of_facets();}
        // space reserved for facets.

    std::size_t bytes() const {
        // bytes used for the polyhedron.
        return sizeof(Self) - sizeof(HDS) + hds.bytes();
    }

    std::size_t bytes_reserved() const {
        // bytes reserved for the polyhedron.
        return sizeof(Self) - sizeof(HDS) + hds.bytes_reserved();
    }

    Vertex_iterator vertices_begin() { return hds.vertices_begin();}
        // iterator over all vertices.

    Vertex_iterator vertices_end() { return hds.vertices_end();}

    Halfedge_iterator halfedges_begin() { return hds.halfedges_begin();}
        // iterator over all halfedges

    Halfedge_iterator halfedges_end() { return hds.halfedges_end();}

    Facet_iterator facets_begin() { return hds.facets_begin();}
        // iterator over all facets

    Facet_iterator facets_end() { return hds.facets_end();}

    Edge_iterator edges_begin() { return halfedges_begin();}
        // iterator over all edges. The iterator refers to halfedges, but
        // enumerates only one of the two corresponding opposite
        // halfedges.

    Edge_iterator edges_end() { return halfedges_end();}
        // end of the range over all edges.

    // The constant iterators and circulators.

    Vertex_const_iterator vertices_begin() const {
        return hds.vertices_begin();
    }
    Vertex_const_iterator vertices_end() const {
        return hds.vertices_end();
    }

    Halfedge_const_iterator halfedges_begin() const {
      return hds.halfedges_begin();
    }
    Halfedge_const_iterator halfedges_end() const {
        return hds.halfedges_end();
    }
    Facet_const_iterator facets_begin() const {return hds.facets_begin();}
    Facet_const_iterator facets_end()   const {return hds.facets_end();}
    Edge_const_iterator edges_begin()   const {return halfedges_begin();}
    Edge_const_iterator edges_end()     const {return halfedges_end();}

    const Traits& traits() const { return m_traits; }

// Geometric Predicates

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    bool is_triangle( Halfedge_const_handle h) const;
        // returns whether the connected component containing h is a
        // single triangle.

    bool is_tetrahedron( Halfedge_const_handle h) const;
        // returns whether the connected component containing h is a
        // single tetrahedron.

// Euler Operators (Combinatorial Modifications)
//
// The following Euler operations modify consistently the combinatorial
// structure of the polyhedral surface. The geometry remains unchanged.

    Halfedge_handle split_facet( Halfedge_handle h,
                                 Halfedge_handle g);
        // split the facet incident to `h' and `g' into two facets with
        // new diagonal between the two vertices denoted by `h' and `g'
        // respectively. The second (new) facet is a copy of the first
        // facet. It returns the new diagonal. The time is proportional to
        // the distance from `h' to `g' around the facet. Precondition:
        // `h' and `g' are incident to the same facet. `h != g' (no
        // loops). `h->next() != g' and `g->next() != h' (no multi-edges).

    Halfedge_handle join_facet( Halfedge_handle h);
        // join the two facets incident to h. The facet incident to
        // `h->opposite()' gets removed. Both facets might be holes.
        // Returns the predecessor of h. The invariant `join_facet(
        // split_facet( h, g))' returns h and keeps the polyhedron
        // unchanged. The time is proportional to the size of the facet
        // removed and the time to compute `h.prev()'. Precondition:
        // `HDS' supports removal of facets. The degree of both
        // vertices incident to h is at least three (no antennas).

    Halfedge_handle split_vertex( Halfedge_handle h,
                                  Halfedge_handle g);
        // split the vertex incident to `h' and `g' into two vertices and
        // connects them with a new edge. The second (new) vertex is a
        // copy of the first vertex. It returns the new edge. The time is
        // proportional to the distance from `h' to `g' around the vertex.
        // Precondition: `h' and `g' are incident to the same vertex. `h
        // != g' (no antennas). `h->next() != g' and `g->next() != h'.

    Halfedge_handle join_vertex( Halfedge_handle h);

// {join the two vertices incident to h. The vertex denoted by
// `h->opposite()' gets removed. Returns the predecessor of h. The
// invariant `join_vertex( split_vertex( h, g))' returns h and keeps the
// polyhedron unchanged. The time is proportional to the degree of the
// vertex removed and the time to compute `h.prev()'. Precondition:
// `HDS' supports removal of vertices. The size of both facets incident
// to h is at least four (no multi-edges)}
//
// Euler Operators Modifying Genus

    Halfedge_handle split_loop( Halfedge_handle h,
                                Halfedge_handle i,
                                Halfedge_handle j);
        // cut the polyhedron into two parts along the cycle (h,i,j).
        // Three copies of the vertices and two new triangles will be
        // created. h,i,j will be incident to the first new triangle. The
        // returnvalue will be an halfedge iterator denoting the new
        // halfegdes of the second new triangle which was h beforehand.
        // Precondition: h,i,j are distinct, consecutive vertices of the
        // polyhedron and form a cycle: i.e. `h->vertex() == i->opposite()
        // ->vertex()', ..., `j->vertex() == h->opposite()->vertex()'. The
        // six facets incident to h,i,j are all distinct.

    Halfedge_handle join_loop( Halfedge_handle h,
                               Halfedge_handle g);
        // glues the boundary of two facets together. Both facets and the
        // vertices of g gets removed. Returns an halfedge iterator for h.
        // The invariant `join_loop( h, split_loop( h, i, j))' returns h
        // and keeps the polyhedron unchanged. Precondition: `HDS'
        // supports removal of vertices and facets. The facets denoted by
        // h and g have equal size.

#endif // CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS //

// Modifying Facets and Holes

    Halfedge_handle make_hole( Halfedge_handle h) {
        // removes incident facet and makes all halfedges incident to the
        // facet to border edges. Returns h. Precondition: `HDS'
        // supports removal of facets. `! h.is_border()'.
        Halfedge_data_structure_decorator<HDS> D;
        return TR_HI( D.make_hole( hds, h.ptr()));
    }

    Halfedge_handle fill_hole( Halfedge_handle h) {
        // fill a hole with a new created facet. Makes all border
        // halfedges of the hole denoted by h incident to the new facet.
        // Returns h. Precondition: `h.is_border()'.
        Halfedge_data_structure_decorator<HDS> D;
        return TR_HI( D.fill_hole( hds, h.ptr()));
    }

    Halfedge_handle add_vertex_and_facet_to_border(   Halfedge_handle h,
                                                      Halfedge_handle g) {
        // creates a new facet within the hole incident to h and g by
        // connecting the tip of g with the tip of h with two new
        // halfedges and a new vertex and filling this separated part of
        // the hole with a new facet. Returns the new halfedge incident to
        // the new facet and the new vertex. Precondition: `h->is_border(
        // )', `g->is_border()', `h != g', and g can be reached along the
        // same hole starting with h.
        CGAL_precondition( h != g);
        Halfedge_data_structure_decorator<HDS> D;
        H* hh = D.add_facet_to_border( hds, h.ptr(), g.ptr());
        CGAL_assertion( hh == &*(g.ptr()->next()));
        D.split_vertex( hds, g.ptr(), hh->opposite());
        return g->next();
    }

    Halfedge_handle add_facet_to_border( Halfedge_handle h,
                                         Halfedge_handle g) {
        // creates a new facet within the hole incident to h and g by
        // connecting the tip of g with the tip of h with a new halfedge
        // and filling this separated part of the hole with a new facet.
        // Returns the new halfedge incident to the new facet.
        // Precondition: `h->is_border()', `g->is_border()', `h != g',
        // `h->next() != g', and g can be reached along the same hole
        // starting with h.
        CGAL_precondition( h != g);
        CGAL_precondition( h->next() != g);
        Halfedge_data_structure_decorator<HDS> D;
        return TR_HI( D.add_facet_to_border( hds, h.ptr(), g.ptr()));
    }

// Erasing

    void erase_facet( Halfedge_handle h) {
        // removes the incident facet of h and changes all halfedges
        // incident to the facet into border edges or removes them from
        // the polyhedral surface if they were already border edges. See
        // `make_hole(h)' for a more specialized variant. Precondition:
        // `Traits' supports removal.
        Halfedge_data_structure_decorator<HDS> D;
        D.erase_facet( hds, h.ptr());
    }

    void erase_connected_component( Halfedge_handle h) {
        // removes the vertices, halfedges, and facets that belong to the
        // connected component of h. Precondition: `Traits' supports
        // removal.
        Halfedge_data_structure_decorator<HDS> D;
        D.erase_connected_component( hds, h.ptr());
    }

    void clear() {
        // removes all vertices, halfedges, and facets.
        hds.delete_all();
    }
    void erase_all() {
        // equivalent to `clear()'. Depricated.
        clear();
    }

// Special Operations on Polyhedral Surfaces

    void delegate( Modifier_base<HDS>& modifier) {
        // calls the `operator()' of the `modifier'. Precondition: The
        // `modifier' returns a consistent representation.
        modifier( hds);
        CGAL_postcondition( is_valid());
    }

// Operations with Border Halfedges

    Size size_of_border_halfedges() const {
        // number of border halfedges. An edge with no incident facet
        // counts as two border halfedges. Precondition: `normalize_border
        // ()' has been called and no halfedge insertion or removal and no
        // change in border status of the halfedges have occured since
        // then.
        return hds.size_of_border_halfedges();
    }

    Size size_of_border_edges() const {
        // number of border edges. If `size_of_border_edges() ==
        // size_of_border_halfedges()' all border edges are incident to a
        // facet on one side and to a hole on the other side.
        // Precondition: `normalize_border()' has been called and no
        // halfedge insertion or removal and no change in border status of
        // the halfedges have occured since then.
        return hds.size_of_border_edges();
    }

    Halfedge_iterator border_halfedges_begin() {
        // halfedge iterator starting with the border edges. The range [
        // `halfedges_begin(), border_halfedges_begin()') denotes all
        // non-border edges. The range [`border_halfedges_begin(),
        // halfedges_end()') denotes all border edges. Precondition:
        // `normalize_border()' has been called and no halfedge insertion
        // or removal and no change in border status of the halfedges have
        // occured since then.
        return hds.border_halfedges_begin();
    }

    Edge_iterator border_edges_begin() { return border_halfedges_begin();}
        // ... trial to make Edge_iterator obsolete.

    Halfedge_const_iterator border_halfedges_begin() const {
        return hds.border_halfedges_begin();
    }

    Edge_const_iterator border_edges_begin() const {
        return border_halfedges_begin();
    }

    bool normalized_border_is_valid( bool verbose = false) const {
        // checks whether all non-border edges precedes the border edges.
        Halfedge_data_structure_decorator<HDS> decorator;
        bool valid = decorator.normalized_border_is_valid( hds, verbose);
        for ( Halfedge_const_iterator i = border_halfedges_begin();
              valid && (i != halfedges_end()); (++i, ++i)) {
            if ( i->is_border()) {
                Verbose_ostream verr(verbose);
                verr << "    both halfedges of an edge are border "
                        "halfedges." << std::endl;
                        valid = false;
            }
        }
        return valid;
    }

    void normalize_border() {
        // sorts halfedges such that the non-border edges precedes the
        // border edges.
             hds.normalize_border();
             CGAL_postcondition( normalized_border_is_valid());
    }

protected:
                      // Supports: normals,      planes
    void inside_out_geometry( Tag_false, Tag_false) {}
    void inside_out_geometry( Tag_true,  Tag_false);
    void inside_out_geometry( Tag_true,  Tag_true);
    void inside_out_geometry( Tag_false, Tag_true)  {
        inside_out_geometry( Tag_true(), Tag_true());
    }

public:

    void inside_out() {
        // reverse facet orientation.
        Halfedge_data_structure_decorator<HDS> decorator;
        decorator.inside_out( hds);
        inside_out_geometry( Supports_facet_normal(),
                             Supports_facet_plane());
    }

    bool is_valid( bool verb = false, int level = 0) const;
        // checks the combinatorial consistency.
#ifdef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    make_triangle( V* v1, V* v2, V* v3) {
        Halfedge_data_structure_decorator<HDS> decorator;
        H* h  = hds.new_edge();
        h->set_next( hds.new_edge());
        h->next()->set_next( hds.new_edge());
        h->next()->next()->set_next( h);
        decorator.set_prev( h, h->next()->next());
        decorator.set_prev( h->next(), h);
        decorator.set_prev( h->next()->next(), h->next());
        h->opposite()->set_next( h->next()->next()->opposite());
        h->next()->opposite()->set_next( h->opposite());
        h->next()->next()->opposite()->set_next( h->next()->opposite());
        decorator.set_prev( h->opposite(), h->next()->opposite());
        decorator.set_prev( h->next()->opposite(),
                           h->next()->next()->opposite());
        decorator.set_prev( h->next()->next()->opposite(), h->opposite());
        // the vertices
        decorator.set_vertex( h, v1);
        decorator.set_vertex( h->next(), v2);
        decorator.set_vertex( h->next()->next(), v3);
        decorator.set_vertex( h->opposite(), v3);
        decorator.set_vertex( h->next()->opposite(), v1);
        decorator.set_vertex( h->next()->next()->opposite(), v2);
        decorator.set_vertex_halfedge( h);
        decorator.set_vertex_halfedge( h->next());
        decorator.set_vertex_halfedge( h->next()->next());
        // the facet
        F* f = decorator.new_facet( hds);
        decorator.set_facet( h, f);
        decorator.set_facet( h->next(), f);
        decorator.set_facet( h->next()->next(), f);
        decorator.set_facet_halfedge( h);
        return TR_HI(h);
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    make_triangle() {
        Halfedge_data_structure_decorator<HDS> decorator;
        return make_triangle( decorator.new_vertex( hds),
                              decorator.new_vertex( hds),
                              decorator.new_vertex( hds));
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    make_triangle(const Point& p1, const Point& p2, const Point& p3) {
        Halfedge_data_structure_decorator<HDS> decorator;
        return make_triangle( decorator.new_vertex( hds, p1),
                              decorator.new_vertex( hds, p2),
                              decorator.new_vertex( hds, p3));
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    make_tetrahedron( V* v1, V* v2, V* v3, V* v4) {
        Halfedge_data_structure_decorator<HDS> decorator;
        H* h  = make_triangle(v1,v2,v3).ptr();
        // The remaining tip.
        H* g  = hds.new_edge();
        decorator.insert_tip( g->opposite(), h->opposite());
        decorator.close_tip( g);
        decorator.set_vertex( g, v4);
        H* e  = hds.new_edge();
        H* d  = hds.new_edge();
        decorator.insert_tip( e->opposite(), h->next()->opposite());
        decorator.insert_tip( e, g);
        decorator.insert_tip( d->opposite(), h->next()->next()->opposite());
        decorator.insert_tip( d, e);
        decorator.set_vertex_halfedge( g);
        // facets
        F* f = decorator.new_facet( hds);
        decorator.set_facet( h->opposite(), f);
        decorator.set_facet( g, f);
        decorator.set_facet( e->opposite(), f);
        decorator.set_facet_halfedge( g);
        f = decorator.new_facet( hds);
        decorator.set_facet( h->next()->opposite(), f);
        decorator.set_facet( e, f);
        decorator.set_facet( d->opposite(), f);
        decorator.set_facet_halfedge( e);
        f = decorator.new_facet( hds);
        decorator.set_facet( h->next()->next()->opposite(), f);
        decorator.set_facet( d, f);
        decorator.set_facet( g->opposite(), f);
        decorator.set_facet_halfedge( d);
        return TR_HI(h);
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    make_tetrahedron() {
        Halfedge_data_structure_decorator<HDS> decorator;
        return make_tetrahedron( decorator.new_vertex( hds),
                                 decorator.new_vertex( hds),
                                 decorator.new_vertex( hds),
                                 decorator.new_vertex( hds));
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    make_tetrahedron(const Point& p1,
                     const Point& p2,
                     const Point& p3,
                     const Point& p4){
        Halfedge_data_structure_decorator<HDS> decorator;
        return make_tetrahedron( decorator.new_vertex( hds, p1),
                                 decorator.new_vertex( hds, p2),
                                 decorator.new_vertex( hds, p3),
                                 decorator.new_vertex( hds, p4));
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    bool
    Polyhedron_3<TR,HDS>::
    #else
    bool
    #endif
    is_triangle( Halfedge_const_handle h1) const {
        Halfedge_const_handle h2 = h1->next();
        Halfedge_const_handle h3 = h1->next()->next();
        // check halfedge combinatorics.
        // exact two edges at vertices 1, 2, 3.
        if ( h1->opposite()->next() != h3->opposite())    return false;
        if ( h2->opposite()->next() != h1->opposite())    return false;
        if ( h3->opposite()->next() != h2->opposite())    return false;
        // The facet is a triangle.
        if ( h1->next()->next()->next() != h1) return false;
    
        if ( check_tag( Supports_halfedge_facet())
             &&  ! h1->is_border_edge())
            return false;  // implies h2 and h3
        CGAL_assertion( ! h1->is_border() || ! h1->opposite()->is_border());
    
        // Assert consistency.
        CGAL_assertion( h1 != h2);
        CGAL_assertion( h1 != h3);
        CGAL_assertion( h3 != h2);
    
        // check prev pointer.
        CGAL_assertion_code( Halfedge_data_structure_decorator<HDS>
                        decorator;)
        CGAL_assertion( decorator.get_prev( h1.ptr()) == NULL ||
                   decorator.get_prev( h1.ptr()) == h3.ptr());
        CGAL_assertion( decorator.get_prev( h2.ptr()) == NULL ||
                   decorator.get_prev( h2.ptr()) == h1.ptr());
        CGAL_assertion( decorator.get_prev( h3.ptr()) == NULL ||
                   decorator.get_prev( h3.ptr()) == h2.ptr());
    
        // check vertices.
        CGAL_assertion( decorator.get_vertex( h1.ptr()) ==
                   decorator.get_vertex( h2->opposite().ptr()));
        CGAL_assertion( decorator.get_vertex( h2.ptr()) ==
                   decorator.get_vertex( h3->opposite().ptr()));
        CGAL_assertion( decorator.get_vertex( h3.ptr()) ==
                   decorator.get_vertex( h1->opposite().ptr()));
    
        CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
                   decorator.get_vertex( h1.ptr()) !=
                   decorator.get_vertex( h2.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
                   decorator.get_vertex( h1.ptr()) !=
                   decorator.get_vertex( h3.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
                   decorator.get_vertex( h2.ptr()) !=
                   decorator.get_vertex( h3.ptr()));
    
        // check facets.
        CGAL_assertion( decorator.get_facet( h1.ptr()) ==
                   decorator.get_facet( h2.ptr()));
        CGAL_assertion( decorator.get_facet( h1.ptr()) ==
                   decorator.get_facet( h3.ptr()));
    
        return true;
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    bool
    Polyhedron_3<TR,HDS>::
    #else
    bool
    #endif
    is_tetrahedron( Halfedge_const_handle h1) const {
        Halfedge_const_handle h2 = h1->next();
        Halfedge_const_handle h3 = h1->next()->next();
        Halfedge_const_handle h4 = h1->opposite()->next();
        Halfedge_const_handle h5 = h2->opposite()->next();
        Halfedge_const_handle h6 = h3->opposite()->next();
        // check halfedge combinatorics.
        // at least three edges at vertices 1, 2, 3.
        if ( h4 == h3->opposite())    return false;
        if ( h5 == h1->opposite())    return false;
        if ( h6 == h2->opposite())    return false;
        // exact three edges at vertices 1, 2, 3.
        if ( h4->opposite()->next() != h3->opposite())    return false;
        if ( h5->opposite()->next() != h1->opposite())    return false;
        if ( h6->opposite()->next() != h2->opposite())    return false;
        // three edges at v4.
        if ( h4->next()->opposite() != h5) return false;
        if ( h5->next()->opposite() != h6) return false;
        if ( h6->next()->opposite() != h4) return false;
        // All facets are triangles.
        if ( h1->next()->next()->next() != h1) return false;
        if ( h4->next()->next()->next() != h4) return false;
        if ( h5->next()->next()->next() != h5) return false;
        if ( h6->next()->next()->next() != h6) return false;
        // all edges are non-border edges.
        if ( h1->is_border()) return false;  // implies h2 and h3
        CGAL_assertion( ! h2->is_border());
        CGAL_assertion( ! h3->is_border());
        if ( h4->is_border()) return false;
        if ( h5->is_border()) return false;
        if ( h6->is_border()) return false;
    
        // Assert consistency.
        CGAL_assertion( h1 != h2);
        CGAL_assertion( h1 != h3);
        CGAL_assertion( h3 != h2);
        CGAL_assertion( h4 != h5);
        CGAL_assertion( h5 != h6);
        CGAL_assertion( h6 != h4);
    
        // check prev pointer.
        CGAL_assertion_code( Halfedge_data_structure_decorator<HDS>
                        decorator;)
        CGAL_assertion( decorator.get_prev( h1.ptr()) == NULL ||
                   decorator.get_prev( h1.ptr()) == h3.ptr());
        CGAL_assertion( decorator.get_prev( h2.ptr()) == NULL ||
                   decorator.get_prev( h2.ptr()) == h1.ptr());
        CGAL_assertion( decorator.get_prev( h3.ptr()) == NULL ||
                   decorator.get_prev( h3.ptr()) == h2.ptr());
        CGAL_assertion( decorator.get_prev( h4.ptr()) == NULL ||
                   decorator.get_prev( h4.ptr()) == h1->opposite().ptr());
        CGAL_assertion( decorator.get_prev( h5.ptr()) == NULL ||
                   decorator.get_prev( h5.ptr()) == h2->opposite().ptr());
        CGAL_assertion( decorator.get_prev( h6.ptr()) == NULL ||
                   decorator.get_prev( h6.ptr()) == h3->opposite().ptr());
    
        // check vertices.
        CGAL_assertion( decorator.get_vertex( h1.ptr()) ==
                   decorator.get_vertex( h2->opposite().ptr()));
        CGAL_assertion( decorator.get_vertex( h1.ptr()) ==
                   decorator.get_vertex( h5->opposite().ptr()));
        CGAL_assertion( decorator.get_vertex( h2.ptr()) ==
                   decorator.get_vertex( h3->opposite().ptr()));
        CGAL_assertion( decorator.get_vertex( h2.ptr()) ==
                   decorator.get_vertex( h6->opposite().ptr()));
        CGAL_assertion( decorator.get_vertex( h3.ptr()) ==
                   decorator.get_vertex( h1->opposite().ptr()));
        CGAL_assertion( decorator.get_vertex( h3.ptr()) ==
                   decorator.get_vertex( h4->opposite().ptr()));
        CGAL_assertion( decorator.get_vertex( h4.ptr()) ==
                   decorator.get_vertex( h5.ptr()));
        CGAL_assertion( decorator.get_vertex( h4.ptr()) ==
                   decorator.get_vertex( h6.ptr()));
    
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   decorator.get_vertex( h1.ptr()) !=
                   decorator.get_vertex( h2.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   decorator.get_vertex( h1.ptr()) !=
                   decorator.get_vertex( h3.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   decorator.get_vertex( h1.ptr()) !=
                   decorator.get_vertex( h4.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   decorator.get_vertex( h2.ptr()) !=
                   decorator.get_vertex( h3.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   decorator.get_vertex( h2.ptr()) !=
                   decorator.get_vertex( h4.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   decorator.get_vertex( h3.ptr()) !=
                   decorator.get_vertex( h4.ptr()));
    
        // check facets.
        CGAL_assertion( decorator.get_facet(h1.ptr()) ==
                   decorator.get_facet(h2.ptr()));
        CGAL_assertion( decorator.get_facet(h1.ptr()) ==
                   decorator.get_facet(h3.ptr()));
        CGAL_assertion( decorator.get_facet(h4.ptr()) ==
                   decorator.get_facet(h4->next().ptr()));
        CGAL_assertion( decorator.get_facet(h4.ptr()) ==
                   decorator.get_facet(h1->opposite().ptr()));
        CGAL_assertion( decorator.get_facet(h5.ptr()) ==
                   decorator.get_facet(h5->next().ptr()));
        CGAL_assertion( decorator.get_facet(h5.ptr()) ==
                   decorator.get_facet(h2->opposite().ptr()));
        CGAL_assertion( decorator.get_facet(h6.ptr()) ==
                   decorator.get_facet(h6->next().ptr()));
        CGAL_assertion( decorator.get_facet(h6.ptr()) ==
                   decorator.get_facet(h3->opposite().ptr()));
    
        CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
                   decorator.get_facet( h1.ptr()) !=
                   decorator.get_facet( h4.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
                   decorator.get_facet( h1.ptr()) !=
                   decorator.get_facet( h5.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
                   decorator.get_facet( h1.ptr()) !=
                   decorator.get_facet( h6.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
                   decorator.get_facet( h4.ptr()) !=
                   decorator.get_facet( h5.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
                   decorator.get_facet( h4.ptr()) !=
                   decorator.get_facet( h6.ptr()));
        CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
                   decorator.get_facet( h5.ptr()) !=
                   decorator.get_facet( h6.ptr()));
    
        return true;
    }
    
    // Euler Operations
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    split_facet( Halfedge_handle h, Halfedge_handle g) {
        Halfedge_data_structure_decorator<HDS> D;
        CGAL_precondition( D.get_facet(h.ptr()) == D.get_facet(g.ptr()));
        CGAL_precondition( h != g);
        CGAL_precondition( h != g->next());
        CGAL_precondition( h->next() != g);
        return TR_HI( D.split_facet( hds, h.ptr(), g.ptr()));
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    join_facet( Halfedge_handle h) {
        Halfedge_data_structure_decorator<HDS> D;
        CGAL_precondition( circulator_size(h->vertex_begin()) >= Size(3));
        CGAL_precondition( circulator_size(h->opposite()->vertex_begin())
                    >= Size(3));
        return TR_HI( D.join_facet( hds, h.ptr()));
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    split_vertex( Halfedge_handle h, Halfedge_handle g) {
        Halfedge_data_structure_decorator<HDS> D;
        CGAL_precondition( D.get_vertex(h.ptr()) == D.get_vertex(g.ptr()));
        CGAL_precondition( h != g);
        return TR_HI( D.split_vertex( hds, h.ptr(), g.ptr()));
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    join_vertex( Halfedge_handle h) {
        Halfedge_data_structure_decorator<HDS> D;
        CGAL_precondition( circulator_size( h->facet_begin()) >= Size(4));
        CGAL_precondition( circulator_size( h->opposite()->facet_begin())
                    >= Size(4));
        return TR_HI( D.join_vertex( hds, h.ptr()));
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    split_loop( Halfedge_handle h,
                Halfedge_handle i,
                Halfedge_handle j) {
        Halfedge_data_structure_decorator<HDS> D;
        CGAL_precondition( h != i);
        CGAL_precondition( h != j);
        CGAL_precondition( i != j);
        CGAL_precondition( D.get_vertex(h.ptr())
                    == D.get_vertex(i->opposite().ptr()));
        CGAL_precondition( D.get_vertex(i.ptr())
                    == D.get_vertex(j->opposite().ptr()));
        CGAL_precondition( D.get_vertex(j.ptr())
                    == D.get_vertex(h->opposite().ptr()));
        CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                    D.get_facet(h.ptr()) != D.get_facet(i.ptr()));
        CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                    D.get_facet(h.ptr()) != D.get_facet(j.ptr()));
        CGAL_precondition( D.get_facet(i.ptr()) == 0 ||
                    D.get_facet(i.ptr()) != D.get_facet(j.ptr()));
        CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                    D.get_facet(h.ptr())
                    != D.get_facet(h->opposite().ptr()));
        CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                    D.get_facet(h.ptr())
                    != D.get_facet(i->opposite().ptr()));
        CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                    D.get_facet(h.ptr())
                    != D.get_facet(j->opposite().ptr()));
        CGAL_precondition( D.get_facet(i.ptr()) == 0 ||
                    D.get_facet(i.ptr())
                    != D.get_facet(h->opposite().ptr()));
        CGAL_precondition( D.get_facet(i.ptr()) == 0 ||
                    D.get_facet(i.ptr())
                    != D.get_facet(i->opposite().ptr()));
        CGAL_precondition( D.get_facet(i.ptr()) == 0 ||
                    D.get_facet(i.ptr())
                    != D.get_facet(j->opposite().ptr()));
        CGAL_precondition( D.get_facet(j.ptr()) == 0 ||
                    D.get_facet(j.ptr())
                    != D.get_facet(h->opposite().ptr()));
        CGAL_precondition( D.get_facet(j.ptr()) == 0 ||
                    D.get_facet(j.ptr())
                    != D.get_facet(i->opposite().ptr()));
        CGAL_precondition( D.get_facet(j.ptr()) == 0 ||
                    D.get_facet(j.ptr())
                    != D.get_facet(j->opposite().ptr()));
        CGAL_precondition( D.get_facet(h->opposite().ptr()) == 0 ||
                    D.get_facet(h->opposite().ptr())
                    != D.get_facet(i->opposite().ptr()));
        CGAL_precondition( D.get_facet(h->opposite().ptr()) == 0 ||
                    D.get_facet(h->opposite().ptr())
                    != D.get_facet(j->opposite().ptr()));
        CGAL_precondition( D.get_facet(i->opposite().ptr()) == 0 ||
                    D.get_facet(i->opposite().ptr())
                    != D.get_facet(j->opposite().ptr()));
        return TR_HI( D.split_loop( hds, h.ptr(), i.ptr(), j.ptr()));
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class TR, class HDS >  CGAL_LARGE_INLINE
    typename Polyhedron_3<TR,HDS>::Halfedge_handle
    Polyhedron_3<TR,HDS>::
    #else
    Halfedge_handle
    #endif
    join_loop( Halfedge_handle h,
               Halfedge_handle g) {
        Halfedge_data_structure_decorator<HDS> D;
        CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                    D.get_facet(h.ptr()) != D.get_facet(g.ptr()));
        CGAL_precondition( circulator_size( h->facet_begin()) >= Size(3));
        CGAL_precondition( circulator_size( h->facet_begin())
                    == circulator_size( g->facet_begin()));
        return TR_HI( D.join_loop( hds, h.ptr(), g.ptr()));
    }
#endif // CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS //
};

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
make_triangle( V* v1, V* v2, V* v3) {
    Halfedge_data_structure_decorator<HDS> decorator;
    H* h  = hds.new_edge();
    h->set_next( hds.new_edge());
    h->next()->set_next( hds.new_edge());
    h->next()->next()->set_next( h);
    decorator.set_prev( h, h->next()->next());
    decorator.set_prev( h->next(), h);
    decorator.set_prev( h->next()->next(), h->next());
    h->opposite()->set_next( h->next()->next()->opposite());
    h->next()->opposite()->set_next( h->opposite());
    h->next()->next()->opposite()->set_next( h->next()->opposite());
    decorator.set_prev( h->opposite(), h->next()->opposite());
    decorator.set_prev( h->next()->opposite(),
                       h->next()->next()->opposite());
    decorator.set_prev( h->next()->next()->opposite(), h->opposite());
    // the vertices
    decorator.set_vertex( h, v1);
    decorator.set_vertex( h->next(), v2);
    decorator.set_vertex( h->next()->next(), v3);
    decorator.set_vertex( h->opposite(), v3);
    decorator.set_vertex( h->next()->opposite(), v1);
    decorator.set_vertex( h->next()->next()->opposite(), v2);
    decorator.set_vertex_halfedge( h);
    decorator.set_vertex_halfedge( h->next());
    decorator.set_vertex_halfedge( h->next()->next());
    // the facet
    F* f = decorator.new_facet( hds);
    decorator.set_facet( h, f);
    decorator.set_facet( h->next(), f);
    decorator.set_facet( h->next()->next(), f);
    decorator.set_facet_halfedge( h);
    return TR_HI(h);
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
make_triangle() {
    Halfedge_data_structure_decorator<HDS> decorator;
    return make_triangle( decorator.new_vertex( hds),
                          decorator.new_vertex( hds),
                          decorator.new_vertex( hds));
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
make_triangle(const Point& p1, const Point& p2, const Point& p3) {
    Halfedge_data_structure_decorator<HDS> decorator;
    return make_triangle( decorator.new_vertex( hds, p1),
                          decorator.new_vertex( hds, p2),
                          decorator.new_vertex( hds, p3));
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
make_tetrahedron( V* v1, V* v2, V* v3, V* v4) {
    Halfedge_data_structure_decorator<HDS> decorator;
    H* h  = make_triangle(v1,v2,v3).ptr();
    // The remaining tip.
    H* g  = hds.new_edge();
    decorator.insert_tip( g->opposite(), h->opposite());
    decorator.close_tip( g);
    decorator.set_vertex( g, v4);
    H* e  = hds.new_edge();
    H* d  = hds.new_edge();
    decorator.insert_tip( e->opposite(), h->next()->opposite());
    decorator.insert_tip( e, g);
    decorator.insert_tip( d->opposite(), h->next()->next()->opposite());
    decorator.insert_tip( d, e);
    decorator.set_vertex_halfedge( g);
    // facets
    F* f = decorator.new_facet( hds);
    decorator.set_facet( h->opposite(), f);
    decorator.set_facet( g, f);
    decorator.set_facet( e->opposite(), f);
    decorator.set_facet_halfedge( g);
    f = decorator.new_facet( hds);
    decorator.set_facet( h->next()->opposite(), f);
    decorator.set_facet( e, f);
    decorator.set_facet( d->opposite(), f);
    decorator.set_facet_halfedge( e);
    f = decorator.new_facet( hds);
    decorator.set_facet( h->next()->next()->opposite(), f);
    decorator.set_facet( d, f);
    decorator.set_facet( g->opposite(), f);
    decorator.set_facet_halfedge( d);
    return TR_HI(h);
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
make_tetrahedron() {
    Halfedge_data_structure_decorator<HDS> decorator;
    return make_tetrahedron( decorator.new_vertex( hds),
                             decorator.new_vertex( hds),
                             decorator.new_vertex( hds),
                             decorator.new_vertex( hds));
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
make_tetrahedron(const Point& p1,
                 const Point& p2,
                 const Point& p3,
                 const Point& p4){
    Halfedge_data_structure_decorator<HDS> decorator;
    return make_tetrahedron( decorator.new_vertex( hds, p1),
                             decorator.new_vertex( hds, p2),
                             decorator.new_vertex( hds, p3),
                             decorator.new_vertex( hds, p4));
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
bool
Polyhedron_3<TR,HDS>::
#else
bool
#endif
is_triangle( Halfedge_const_handle h1) const {
    Halfedge_const_handle h2 = h1->next();
    Halfedge_const_handle h3 = h1->next()->next();
    // check halfedge combinatorics.
    // exact two edges at vertices 1, 2, 3.
    if ( h1->opposite()->next() != h3->opposite())    return false;
    if ( h2->opposite()->next() != h1->opposite())    return false;
    if ( h3->opposite()->next() != h2->opposite())    return false;
    // The facet is a triangle.
    if ( h1->next()->next()->next() != h1) return false;

    if ( check_tag( Supports_halfedge_facet())
         &&  ! h1->is_border_edge())
        return false;  // implies h2 and h3
    CGAL_assertion( ! h1->is_border() || ! h1->opposite()->is_border());

    // Assert consistency.
    CGAL_assertion( h1 != h2);
    CGAL_assertion( h1 != h3);
    CGAL_assertion( h3 != h2);

    // check prev pointer.
    CGAL_assertion_code( Halfedge_data_structure_decorator<HDS>
                    decorator;)
    CGAL_assertion( decorator.get_prev( h1.ptr()) == NULL ||
               decorator.get_prev( h1.ptr()) == h3.ptr());
    CGAL_assertion( decorator.get_prev( h2.ptr()) == NULL ||
               decorator.get_prev( h2.ptr()) == h1.ptr());
    CGAL_assertion( decorator.get_prev( h3.ptr()) == NULL ||
               decorator.get_prev( h3.ptr()) == h2.ptr());

    // check vertices.
    CGAL_assertion( decorator.get_vertex( h1.ptr()) ==
               decorator.get_vertex( h2->opposite().ptr()));
    CGAL_assertion( decorator.get_vertex( h2.ptr()) ==
               decorator.get_vertex( h3->opposite().ptr()));
    CGAL_assertion( decorator.get_vertex( h3.ptr()) ==
               decorator.get_vertex( h1->opposite().ptr()));

    CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
               decorator.get_vertex( h1.ptr()) !=
               decorator.get_vertex( h2.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
               decorator.get_vertex( h1.ptr()) !=
               decorator.get_vertex( h3.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
               decorator.get_vertex( h2.ptr()) !=
               decorator.get_vertex( h3.ptr()));

    // check facets.
    CGAL_assertion( decorator.get_facet( h1.ptr()) ==
               decorator.get_facet( h2.ptr()));
    CGAL_assertion( decorator.get_facet( h1.ptr()) ==
               decorator.get_facet( h3.ptr()));

    return true;
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
bool
Polyhedron_3<TR,HDS>::
#else
bool
#endif
is_tetrahedron( Halfedge_const_handle h1) const {
    Halfedge_const_handle h2 = h1->next();
    Halfedge_const_handle h3 = h1->next()->next();
    Halfedge_const_handle h4 = h1->opposite()->next();
    Halfedge_const_handle h5 = h2->opposite()->next();
    Halfedge_const_handle h6 = h3->opposite()->next();
    // check halfedge combinatorics.
    // at least three edges at vertices 1, 2, 3.
    if ( h4 == h3->opposite())    return false;
    if ( h5 == h1->opposite())    return false;
    if ( h6 == h2->opposite())    return false;
    // exact three edges at vertices 1, 2, 3.
    if ( h4->opposite()->next() != h3->opposite())    return false;
    if ( h5->opposite()->next() != h1->opposite())    return false;
    if ( h6->opposite()->next() != h2->opposite())    return false;
    // three edges at v4.
    if ( h4->next()->opposite() != h5) return false;
    if ( h5->next()->opposite() != h6) return false;
    if ( h6->next()->opposite() != h4) return false;
    // All facets are triangles.
    if ( h1->next()->next()->next() != h1) return false;
    if ( h4->next()->next()->next() != h4) return false;
    if ( h5->next()->next()->next() != h5) return false;
    if ( h6->next()->next()->next() != h6) return false;
    // all edges are non-border edges.
    if ( h1->is_border()) return false;  // implies h2 and h3
    CGAL_assertion( ! h2->is_border());
    CGAL_assertion( ! h3->is_border());
    if ( h4->is_border()) return false;
    if ( h5->is_border()) return false;
    if ( h6->is_border()) return false;

    // Assert consistency.
    CGAL_assertion( h1 != h2);
    CGAL_assertion( h1 != h3);
    CGAL_assertion( h3 != h2);
    CGAL_assertion( h4 != h5);
    CGAL_assertion( h5 != h6);
    CGAL_assertion( h6 != h4);

    // check prev pointer.
    CGAL_assertion_code( Halfedge_data_structure_decorator<HDS>
                    decorator;)
    CGAL_assertion( decorator.get_prev( h1.ptr()) == NULL ||
               decorator.get_prev( h1.ptr()) == h3.ptr());
    CGAL_assertion( decorator.get_prev( h2.ptr()) == NULL ||
               decorator.get_prev( h2.ptr()) == h1.ptr());
    CGAL_assertion( decorator.get_prev( h3.ptr()) == NULL ||
               decorator.get_prev( h3.ptr()) == h2.ptr());
    CGAL_assertion( decorator.get_prev( h4.ptr()) == NULL ||
               decorator.get_prev( h4.ptr()) == h1->opposite().ptr());
    CGAL_assertion( decorator.get_prev( h5.ptr()) == NULL ||
               decorator.get_prev( h5.ptr()) == h2->opposite().ptr());
    CGAL_assertion( decorator.get_prev( h6.ptr()) == NULL ||
               decorator.get_prev( h6.ptr()) == h3->opposite().ptr());

    // check vertices.
    CGAL_assertion( decorator.get_vertex( h1.ptr()) ==
               decorator.get_vertex( h2->opposite().ptr()));
    CGAL_assertion( decorator.get_vertex( h1.ptr()) ==
               decorator.get_vertex( h5->opposite().ptr()));
    CGAL_assertion( decorator.get_vertex( h2.ptr()) ==
               decorator.get_vertex( h3->opposite().ptr()));
    CGAL_assertion( decorator.get_vertex( h2.ptr()) ==
               decorator.get_vertex( h6->opposite().ptr()));
    CGAL_assertion( decorator.get_vertex( h3.ptr()) ==
               decorator.get_vertex( h1->opposite().ptr()));
    CGAL_assertion( decorator.get_vertex( h3.ptr()) ==
               decorator.get_vertex( h4->opposite().ptr()));
    CGAL_assertion( decorator.get_vertex( h4.ptr()) ==
               decorator.get_vertex( h5.ptr()));
    CGAL_assertion( decorator.get_vertex( h4.ptr()) ==
               decorator.get_vertex( h6.ptr()));

    CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
               decorator.get_vertex( h1.ptr()) !=
               decorator.get_vertex( h2.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
               decorator.get_vertex( h1.ptr()) !=
               decorator.get_vertex( h3.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
               decorator.get_vertex( h1.ptr()) !=
               decorator.get_vertex( h4.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
               decorator.get_vertex( h2.ptr()) !=
               decorator.get_vertex( h3.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
               decorator.get_vertex( h2.ptr()) !=
               decorator.get_vertex( h4.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
               decorator.get_vertex( h3.ptr()) !=
               decorator.get_vertex( h4.ptr()));

    // check facets.
    CGAL_assertion( decorator.get_facet(h1.ptr()) ==
               decorator.get_facet(h2.ptr()));
    CGAL_assertion( decorator.get_facet(h1.ptr()) ==
               decorator.get_facet(h3.ptr()));
    CGAL_assertion( decorator.get_facet(h4.ptr()) ==
               decorator.get_facet(h4->next().ptr()));
    CGAL_assertion( decorator.get_facet(h4.ptr()) ==
               decorator.get_facet(h1->opposite().ptr()));
    CGAL_assertion( decorator.get_facet(h5.ptr()) ==
               decorator.get_facet(h5->next().ptr()));
    CGAL_assertion( decorator.get_facet(h5.ptr()) ==
               decorator.get_facet(h2->opposite().ptr()));
    CGAL_assertion( decorator.get_facet(h6.ptr()) ==
               decorator.get_facet(h6->next().ptr()));
    CGAL_assertion( decorator.get_facet(h6.ptr()) ==
               decorator.get_facet(h3->opposite().ptr()));

    CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
               decorator.get_facet( h1.ptr()) !=
               decorator.get_facet( h4.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
               decorator.get_facet( h1.ptr()) !=
               decorator.get_facet( h5.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
               decorator.get_facet( h1.ptr()) !=
               decorator.get_facet( h6.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
               decorator.get_facet( h4.ptr()) !=
               decorator.get_facet( h5.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
               decorator.get_facet( h4.ptr()) !=
               decorator.get_facet( h6.ptr()));
    CGAL_assertion( ! check_tag( Supports_halfedge_facet()) ||
               decorator.get_facet( h5.ptr()) !=
               decorator.get_facet( h6.ptr()));

    return true;
}

// Euler Operations

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
split_facet( Halfedge_handle h, Halfedge_handle g) {
    Halfedge_data_structure_decorator<HDS> D;
    CGAL_precondition( D.get_facet(h.ptr()) == D.get_facet(g.ptr()));
    CGAL_precondition( h != g);
    CGAL_precondition( h != g->next());
    CGAL_precondition( h->next() != g);
    return TR_HI( D.split_facet( hds, h.ptr(), g.ptr()));
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
join_facet( Halfedge_handle h) {
    Halfedge_data_structure_decorator<HDS> D;
    CGAL_precondition( circulator_size(h->vertex_begin()) >= Size(3));
    CGAL_precondition( circulator_size(h->opposite()->vertex_begin())
                >= Size(3));
    return TR_HI( D.join_facet( hds, h.ptr()));
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
split_vertex( Halfedge_handle h, Halfedge_handle g) {
    Halfedge_data_structure_decorator<HDS> D;
    CGAL_precondition( D.get_vertex(h.ptr()) == D.get_vertex(g.ptr()));
    CGAL_precondition( h != g);
    return TR_HI( D.split_vertex( hds, h.ptr(), g.ptr()));
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
join_vertex( Halfedge_handle h) {
    Halfedge_data_structure_decorator<HDS> D;
    CGAL_precondition( circulator_size( h->facet_begin()) >= Size(4));
    CGAL_precondition( circulator_size( h->opposite()->facet_begin())
                >= Size(4));
    return TR_HI( D.join_vertex( hds, h.ptr()));
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
split_loop( Halfedge_handle h,
            Halfedge_handle i,
            Halfedge_handle j) {
    Halfedge_data_structure_decorator<HDS> D;
    CGAL_precondition( h != i);
    CGAL_precondition( h != j);
    CGAL_precondition( i != j);
    CGAL_precondition( D.get_vertex(h.ptr())
                == D.get_vertex(i->opposite().ptr()));
    CGAL_precondition( D.get_vertex(i.ptr())
                == D.get_vertex(j->opposite().ptr()));
    CGAL_precondition( D.get_vertex(j.ptr())
                == D.get_vertex(h->opposite().ptr()));
    CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                D.get_facet(h.ptr()) != D.get_facet(i.ptr()));
    CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                D.get_facet(h.ptr()) != D.get_facet(j.ptr()));
    CGAL_precondition( D.get_facet(i.ptr()) == 0 ||
                D.get_facet(i.ptr()) != D.get_facet(j.ptr()));
    CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                D.get_facet(h.ptr())
                != D.get_facet(h->opposite().ptr()));
    CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                D.get_facet(h.ptr())
                != D.get_facet(i->opposite().ptr()));
    CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                D.get_facet(h.ptr())
                != D.get_facet(j->opposite().ptr()));
    CGAL_precondition( D.get_facet(i.ptr()) == 0 ||
                D.get_facet(i.ptr())
                != D.get_facet(h->opposite().ptr()));
    CGAL_precondition( D.get_facet(i.ptr()) == 0 ||
                D.get_facet(i.ptr())
                != D.get_facet(i->opposite().ptr()));
    CGAL_precondition( D.get_facet(i.ptr()) == 0 ||
                D.get_facet(i.ptr())
                != D.get_facet(j->opposite().ptr()));
    CGAL_precondition( D.get_facet(j.ptr()) == 0 ||
                D.get_facet(j.ptr())
                != D.get_facet(h->opposite().ptr()));
    CGAL_precondition( D.get_facet(j.ptr()) == 0 ||
                D.get_facet(j.ptr())
                != D.get_facet(i->opposite().ptr()));
    CGAL_precondition( D.get_facet(j.ptr()) == 0 ||
                D.get_facet(j.ptr())
                != D.get_facet(j->opposite().ptr()));
    CGAL_precondition( D.get_facet(h->opposite().ptr()) == 0 ||
                D.get_facet(h->opposite().ptr())
                != D.get_facet(i->opposite().ptr()));
    CGAL_precondition( D.get_facet(h->opposite().ptr()) == 0 ||
                D.get_facet(h->opposite().ptr())
                != D.get_facet(j->opposite().ptr()));
    CGAL_precondition( D.get_facet(i->opposite().ptr()) == 0 ||
                D.get_facet(i->opposite().ptr())
                != D.get_facet(j->opposite().ptr()));
    return TR_HI( D.split_loop( hds, h.ptr(), i.ptr(), j.ptr()));
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class TR, class HDS >  CGAL_LARGE_INLINE
typename Polyhedron_3<TR,HDS>::Halfedge_handle
Polyhedron_3<TR,HDS>::
#else
Halfedge_handle
#endif
join_loop( Halfedge_handle h,
           Halfedge_handle g) {
    Halfedge_data_structure_decorator<HDS> D;
    CGAL_precondition( D.get_facet(h.ptr()) == 0 ||
                D.get_facet(h.ptr()) != D.get_facet(g.ptr()));
    CGAL_precondition( circulator_size( h->facet_begin()) >= Size(3));
    CGAL_precondition( circulator_size( h->facet_begin())
                == circulator_size( g->facet_begin()));
    return TR_HI( D.join_loop( hds, h.ptr(), g.ptr()));
}
#endif // CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS //
// Special Operations on Polyhedral Surfaces

template < class TR, class HDS >  CGAL_LARGE_INLINE
void                               // Supports: normals,      planes
Polyhedron_3<TR,HDS>::inside_out_geometry(Tag_true,Tag_false) {
    typename Traits::Construct_opposite_vector_3 opposite
        = traits().construct_opposite_vector_3_object();
    for ( Facet_iterator i = facets_begin(); i != facets_end(); ++i)
        i->normal() = opposite( i->normal());
}

template < class TR, class HDS >  CGAL_LARGE_INLINE
void                               // Supports: normals,      planes
Polyhedron_3<TR,HDS>::inside_out_geometry(Tag_true,Tag_true) {
    typename Traits::Construct_opposite_plane_3 opposite
        = traits().construct_opposite_plane_3_object();
    for ( Facet_iterator i = facets_begin(); i != facets_end(); ++i)
        i->plane() = opposite( i->plane());
}

template < class TR, class HDS >
bool
Polyhedron_3<TR,HDS>:: is_valid( bool verb, int level) const {
    Verbose_ostream verr(verb);
    verr << "begin Polyhedron_3<TR,HDS>::is_valid( verb=true, "
                      "level = " << level << "):" << std::endl;
    Halfedge_data_structure_decorator<HDS> decorator;
    bool valid = decorator.is_valid( hds, verb, level + 3);
    // All halfedges.
    Halfedge_const_iterator begin = halfedges_begin();
    Halfedge_const_iterator end   = halfedges_end();
    Size  n = 0;
    for( ; valid && (begin != end); begin++) {
        verr << "halfedge " << n << std::endl;
        // At least triangular facets and distinct geometry.
        valid = valid && ( begin->next() != begin);
        valid = valid && ( begin->next()->next() != begin);
        valid = valid && ( ! check_tag( Supports_halfedge_vertex()) ||
                           decorator.get_vertex( begin.ptr()) !=
                           decorator.get_vertex( begin->opposite().ptr()));
        valid = valid && ( ! check_tag( Supports_halfedge_vertex()) ||
                           decorator.get_vertex( begin.ptr()) !=
                           decorator.get_vertex( begin->next().ptr()));
        valid = valid && ( ! check_tag( Supports_halfedge_vertex()) ||
                           decorator.get_vertex( begin.ptr()) !=
                           decorator.get_vertex( begin->next()->
                                                 next().ptr()));
        if ( ! valid) {
            verr << "    incident facet is not at least a triangle."
                 << std::endl;
            break;
        }
        // Distinct facets on each side of an halfegde.
        valid = valid && ( ! check_tag( Supports_halfedge_facet()) ||
                           decorator.get_facet( begin.ptr()) !=
                           decorator.get_facet( begin->opposite().ptr()));
        if ( ! valid) {
            verr << "    both incident facets are equal." << std::endl;
            break;
        }
        ++n;
    }
    valid = valid && (n == size_of_halfedges());
    if ( n != size_of_halfedges())
        verr << "counting halfedges failed." << std::endl;

    verr << "end of Polyhedron_3<TR,HDS>::is_valid(): structure is "
         << ( valid ? "valid." : "NOT VALID.") << std::endl;
    return valid;
}



CGAL_END_NAMESPACE

#endif // CGAL_POLYHEDRON_OLD_POLYHEDRON_3_H //
// EOF //
