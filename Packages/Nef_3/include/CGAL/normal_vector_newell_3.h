#ifndef CGAL_NORMAL_VECTOR_NEWELL_3_H
#define CGAL_NORMAL_VECTOR_NEWELL_3_H 1

#include <CGAL/iterator.h>

CGAL_BEGIN_NAMESPACE

template <class Handle, class Vector>
CGAL_MEDIUM_INLINE
void newell_single_step_3( Handle p, Handle q, Vector& n )
{
  n = Vector(
    n.hx()
        * p->hw() * p->hw()
        * q->hw() * q->hw()
      +   n.hw()
        * ( p->hy() * q->hw() - q->hy() * p->hw())
        * ( p->hz() * q->hw() + q->hz() * p->hw()),
    n.hy()
        * p->hw() * p->hw()
        * q->hw() * q->hw()
      +   n.hw()
        * ( p->hz() * q->hw() - q->hz() * p->hw())
        * ( p->hx() * q->hw() + q->hx() * p->hw()),
    n.hz()
        * p->hw() * p->hw()
        * q->hw() * q->hw()
      +   n.hw()
        * ( p->hx() * q->hw() - q->hx() * p->hw())
        * ( p->hy() * q->hw() + q->hy() * p->hw()),
    n.hw()
        * p->hw() * p->hw()
        * q->hw() * q->hw()
    );
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
        newell_single_step_3( prev, first, n);
        prev = first;
        ++first;
    }
    newell_single_step_3( prev, start_point, n);
}

CGAL_END_NAMESPACE

#endif // CGAL_NORMAL_VECTOR_NEWELL_3_H
