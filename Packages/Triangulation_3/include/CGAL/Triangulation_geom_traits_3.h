// ============================================================================
//
// $Id$
//
// geometric traits for a <=3 D triangulation
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_GEOM_TRAITS_H
#define CGAL_TRIANGULATION_GEOM_TRAITS_H

#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names.h>

template < class R >
class CGAL_Triangulation_geom_traits 
{
public:

  typedef CGAL_Point_3<R>  Point;

  CGAL_Orientation orientation(const Point &p,
			       const Point &q,
			       const Point &r,
			       const Point &s) const
  {
    return CGAL_orientation(p, q, r, s);
  }

};


#endif // CGAL_TRIANGULATION_GEOM_TRAITS_H
