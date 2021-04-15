// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
#ifndef CGAL_NORMAL_VECTOR_NEWELL_3_H
#define CGAL_NORMAL_VECTOR_NEWELL_3_H 1

#include <CGAL/license/Nef_3.h>


#include <CGAL/iterator.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 79
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

namespace internal_nef
{
template <class Handle, class Vector>
CGAL_MEDIUM_INLINE
void newell_single_step_3( const Handle& p, const Handle& q, Vector& n )
{
  n = Vector(
    n.hx()
        * p.hw() * p.hw()
        * q.hw() * q.hw()
      +   n.hw()
        * ( p.hy() * q.hw() - q.hy() * p.hw())
        * ( p.hz() * q.hw() + q.hz() * p.hw()),
    n.hy()
        * p.hw() * p.hw()
        * q.hw() * q.hw()
      +   n.hw()
        * ( p.hz() * q.hw() - q.hz() * p.hw())
        * ( p.hx() * q.hw() + q.hx() * p.hw()),
    n.hz()
        * p.hw() * p.hw()
        * q.hw() * q.hw()
      +   n.hw()
        * ( p.hx() * q.hw() - q.hx() * p.hw())
        * ( p.hy() * q.hw() + q.hy() * p.hw()),
    n.hw()
        * p.hw() * p.hw()
        * q.hw() * q.hw()
    );
}
}

template <class IC, class Vector>
void normal_vector_newell_3( IC first, IC last, Vector& n )
    // compute a facet normal n via the Newell-method for the facet
    // surrounded by the points in the range [`h', `last'), where `IC'
    // denotes either an iterator or a circulator with a point type as
    // value type. The Newell-method computes for non-planar facets the
    // best least-square-fit normal of the plane equation. Precondition:
    // The numbertype for the points and the vector n must be compatible
    // with `double'. The points and the vector n must have dimension
    // three.
{
    CGAL_assertion( !CGAL::is_empty_range( first, last));
    // Compute facet normals via the Newell-method as described in
    // Filippo Tampieri: Newell's Method for Computing the Plane
    // Equation of a Polygon. Graphics Gems III, David Kirk,
    // AP Professional, 1992.
    // The code works with cartesian and with semi-rational points.
    n = Vector( 0, 0, 0);        // init current normal to 0
    IC start_point = first;
    IC prev = first;
    ++first;
    while( first != last) {
        internal_nef::newell_single_step_3( *prev, *first, n);
        prev = first;
        ++first;
    }
    internal_nef::newell_single_step_3( *prev, *start_point, n);
    CGAL_NEF_TRACEN("newell normal vector "<<n);
}

template <class IC, class Vector, class VertexPointMap>
void normal_vector_newell_3( IC first, IC last, VertexPointMap vpm, Vector& n )
{
    CGAL_assertion( !CGAL::is_empty_range( first, last));
    // Compute facet normals via the Newell-method as described in
    // Filippo Tampieri: Newell's Method for Computing the Plane
    // Equation of a Polygon. Graphics Gems III, David Kirk,
    // AP Professional, 1992.
    // The code works with cartesian and with semi-rational points.
    n = Vector( 0, 0, 0);        // init current normal to 0
    IC start_point = first;
    IC prev = first;
    ++first;
    while( first != last) {
        internal_nef::newell_single_step_3( get(vpm,*prev), get(vpm,*first), n);
        prev = first;
        ++first;
    }
    internal_nef::newell_single_step_3( get(vpm,*prev), get(vpm,*start_point), n);
    CGAL_NEF_TRACEN("newell normal vector "<<n);
}

} //namespace CGAL

#endif // CGAL_NORMAL_VECTOR_NEWELL_3_H
