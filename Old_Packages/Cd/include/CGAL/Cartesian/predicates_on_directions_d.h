#ifndef CGAL_PREDICATES_ON_DIRECTIONS_D_H
#define CGAL_PREDICATES_ON_DIRECTIONS_D_H

#include <CGAL/predicates/kernel_ftCd.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
bool
equal_direction(const DirectionCd<R CGAL_CTAG>& d1,
                const DirectionCd<R CGAL_CTAG>& d2)
{
  return equal_directionCd(d1.begin(),d1.end(),d2.begin(),d2.end());
}

CGAL_END_NAMESPACE


#endif // CGAL_PREDICATES_ON_DIRECTIONS_D_H
