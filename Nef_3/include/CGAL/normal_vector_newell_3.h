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
#include <CGAL/representation_tags.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL {

namespace internal_nef
{
template <class Handle, class Vector>
CGAL_MEDIUM_INLINE
void newell_single_step_3( const Handle& p, const Handle& q, Vector& n, const Homogeneous_tag&)
{
  typedef typename Kernel_traits<Vector>::Kernel::RT RT;
  const RT& phw = p.hw();
  const RT& qhw = q.hw();
  const RT& nhw = n.hw();
  const RT& sq = square(phw) * square(qhw);
  const RT& phyqhw = p.hy() * qhw;
  const RT& qhyphw = q.hy() * phw;
  const RT& phxqhw = p.hx() * qhw;
  const RT& qhxphw = q.hx() * phw;
  const RT& phzqhw = p.hz() * qhw;
  const RT& qhzphw = q.hz() * phw;

  n = Vector(
    n.hx()
        * sq
      + nhw
        * ( phyqhw - qhyphw )
        * ( phzqhw + qhzphw ),
    n.hy()
        * sq
      + nhw
        * ( phzqhw - qhzphw )
        * ( phxqhw + qhxphw ),
    n.hz()
        * sq
      + nhw
        * ( phxqhw - qhxphw )
        * ( phyqhw + qhyphw ),
    n.hw()
        * sq
    );
}

template <class Handle, class Vector>
CGAL_MEDIUM_INLINE
void newell_single_step_3( const Handle& p, const Handle& q, Vector& n, const Cartesian_tag&)
{
  typedef typename Kernel_traits<Vector>::Kernel::FT FT;
  const FT& py = p.y();
  const FT& qy = q.y();
  const FT& px = p.x();
  const FT& qx = q.x();
  const FT& pz = p.z();
  const FT& qz = q.z();

  n = Vector(
    n.x()
        + ( py - qy )
        * ( pz + qz ),
    n.y()
        + ( pz - qz )
        * ( px + qx ),
    n.z()
        + ( px - qx )
        * ( py + qy )
    );
}

template <class IC>
bool is_triangle_3( const IC& first )
{
    return std::next(first,3) == first;
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

    if(internal_nef::is_triangle_3(first)) {
      n = orthogonal_vector(*first,*(std::next(first)),*(std::next(first,2)));
      return;
    }


    typedef typename Kernel_traits<Vector>::Kernel R;
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
        internal_nef::newell_single_step_3( *prev, *first, n, typename R::Kernel_tag());
        prev = first;
        ++first;
    }
    internal_nef::newell_single_step_3( *prev, *start_point, n, typename R::Kernel_tag());
    CGAL_NEF_TRACEN("newell normal vector "<<n);
}

template <class IC, class Vector, class VertexPointMap>
void normal_vector_newell_3( IC first, IC last, VertexPointMap vpm, Vector& n )
{
    CGAL_assertion( !CGAL::is_empty_range( first, last));

    if(internal_nef::is_triangle_3(first)) {

      n = orthogonal_vector(get(vpm,*first),get(vpm,*(std::next(first))),get(vpm,*(std::next(first,2))));
      return;
    }


    typedef typename Kernel_traits<Vector>::Kernel R;
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
        internal_nef::newell_single_step_3( get(vpm,*prev), get(vpm,*first), n, typename R::Kernel_tag());
        prev = first;
        ++first;
    }
    internal_nef::newell_single_step_3( get(vpm,*prev), get(vpm,*start_point), n, typename R::Kernel_tag());
    CGAL_NEF_TRACEN("newell normal vector "<<n);
}

} //namespace CGAL

#endif // CGAL_NORMAL_VECTOR_NEWELL_3_H
