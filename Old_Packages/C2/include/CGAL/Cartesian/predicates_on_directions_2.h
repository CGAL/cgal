// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_PREDICATES_ON_DIRECTIONS_2_H
#define CGAL_CARTESIAN_PREDICATES_ON_DIRECTIONS_2_H

CGAL_BEGIN_NAMESPACE

template < class R >
inline
bool
equal_direction(const DirectionC2<R CGAL_CTAG>& d1,
                const DirectionC2<R CGAL_CTAG>& d2)
{
  return equal_directionC2(d1.dx(),d1.dy(),d2.dx(),d2.dy());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
compare_angle_with_x_axis(const DirectionC2<R CGAL_CTAG>& d1,
                          const DirectionC2<R CGAL_CTAG>& d2)
{
  return compare_angle_with_x_axisC2(d1.dx(),d1.dy(),d2.dx(),d2.dy());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PREDICATES_ON_DIRECTIONS_2_H
