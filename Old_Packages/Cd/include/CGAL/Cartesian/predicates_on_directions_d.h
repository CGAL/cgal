#ifndef CGAL_PREDICATES_ON_DIRECTIONS_D_H
#define CGAL_PREDICATES_ON_DIRECTIONS_D_H

CGAL_BEGIN_NAMESPACE

template < class R >
inline
bool
equal_direction(const DirectionC2<R CGAL_CTAG>& d1,
                const DirectionC2<R CGAL_CTAG>& d2)
{
  return equal_directionC2(d1.begin(),d1.end(),d2.begin(),d2.end());
}

CGAL_END_NAMESPACE


#endif // CGAL_PREDICATES_ON_DIRECTIONS_D_H
