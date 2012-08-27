// Copyright (c) 2000, 2006  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>)

// Cuts a connected piece off a Polyhedral Surfaces.

#ifndef CGAL_POLYHEDRON_CUT_PLANE_3_H
#define CGAL_POLYHEDRON_CUT_PLANE_3_H 1

#include <CGAL/basic.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/halfedgeDS_cut_component.h>
#include <CGAL/intersections.h> 


namespace CGAL {

// Auxiliary classes and functions to make polyhedron_cut_plane_3 work.
// See below for the implementations of polyhedron_cut_plane_3.

template < class HDS, class Predicate>
class Polyhedron_cut_component_3 : public Modifier_base<HDS> {
public:
    typedef typename HDS::Halfedge_handle Halfedge_handle;
private:
    Halfedge_handle h;
    Predicate       pred;
public:
    Polyhedron_cut_component_3( Halfedge_handle hh, const Predicate& p)
        : h(hh), pred(p) {}
    void operator()( HDS& target) {
        h = halfedgeDS_cut_component( target, h, pred);
    }
    Halfedge_handle halfedge() const { return h; }
};

template < class Poly, class Predicate >
typename Poly::Halfedge_handle 
I_polyhedron_cut_component( Poly&                          poly,
                            typename Poly::Halfedge_handle h,
                            Predicate                      pred) {
    typedef typename Poly::HalfedgeDS HDS;
    typedef Polyhedron_cut_component_3<HDS,Predicate> Modifier;
    Modifier modifier( h, pred);
    poly.delegate( modifier);
    return modifier.halfedge();
}
    
template < class Poly, class Plane, class Traits>
class I_Polyhedron_cut_plane_predicate {
    const Plane&  plane;
    const Traits& traits;
public:
    typedef typename Poly::Vertex_const_handle Vertex_const_handle;
    I_Polyhedron_cut_plane_predicate( const Plane& pl, const Traits& tr)
        : plane(pl), traits(tr) {}
    bool operator()( Vertex_const_handle v) const {
        return traits.has_on_negative_side_3_object()( plane, v->point());
    }
};

template <class Kernel>
struct I_Construct_point_3 {
    typedef typename Kernel::Point_3 Point_3;
    typedef typename Kernel::Plane_3 Plane_3;
    typedef typename Kernel::Line_3  Line_3;
    Point_3
    operator()( const Plane_3& p, const Plane_3& q, const Plane_3& r) const { 
        Object obj = intersection( p, q);
        Line_3 line;
        if ( assign( line, obj)) {
            obj = intersection( r, line);
            Point_3 pt;
            if ( assign( pt, obj)) {
                return pt;
            }
        }
        std::cerr << "ERROR: coplanar planes used for computing "
                     "intersecting point." << std::endl;
        std::cerr << "       Return ORIGIN. Don't trust result." << std::endl;
        return ORIGIN;
    }
};

template < class Poly, class Plane, class Traits>
typename Poly::Halfedge_handle
polyhedron_cut_plane_3( Poly& poly,
                        typename Poly::Halfedge_handle h,
                        const Plane& plane,
                        const Traits& traits)
// Cuts the polyhedron `poly' at plane `plane' starting at halfedge `h'.
// Traces the intersection curve of `plane' with `poly' starting at `h',
// cuts `poly' along that intersection curve and deletes the (connected)
// component on the negative side of the plane. The hole created along
// the intersection curve is filled with a new face containing the plane
// and the points in the vertices are computed.
{
    typedef typename Poly::Halfedge_handle  Halfedge_handle;
    typedef I_Polyhedron_cut_plane_predicate< Poly, Plane, Traits> Predicate;
    typedef I_Construct_point_3< Traits> Construct_point;
    Predicate pred( plane, traits);
    Construct_point construct_point; // replace this with kernel one day.
    CGAL_precondition( poly.is_valid());
    CGAL_precondition( pred( h->vertex()));
    CGAL_precondition( ! pred( h->opposite()->vertex()));
    h = I_polyhedron_cut_component( poly, h, pred);
    // Assign geometry
    h->facet()->plane() = plane;
    Halfedge_handle start = h;
    do {
        h->vertex()->point() = 
            construct_point( plane,
                             h->next()->opposite()->facet()->plane(),
                             h->opposite()->facet()->plane());
        h = h->next();
    } while ( h != start);
    CGAL_postcondition( poly.is_valid());
    return h;
}


template < class Poly, class Plane>
typename Poly::Halfedge_handle
polyhedron_cut_plane_3( Poly& poly, 
                        typename Poly::Halfedge_handle h,
                        const Plane& plane)
// Same function as above using the kernel that comes with the plane.
{
    typedef CGAL::Kernel_traits<Plane>  KTraits;
    typedef typename KTraits::Kernel    Kernel;
    return polyhedron_cut_plane_3( poly, h, plane, Kernel());
}

} //namespace CGAL
#endif // CGAL_POLYHEDRON_CUT_PLANE_3_H //
// EOF //
