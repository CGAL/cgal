// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_POINTS_3_H
#define CGAL_CARTESIAN_CONSTRUCTIONS_ON_POINTS_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
PointC3<R CGAL_CTAG>
midpoint(PointC3<R CGAL_CTAG> const& p,
         PointC3<R CGAL_CTAG> const& q )
{
  typename R::FT x,y,z;
  midpointC3(p.x(),p.y(),p.z(),q.x(),q.y(),q.z(),x,y,z);
  return PointC3<R CGAL_CTAG>(x,y,z);
}

template < class R >
PointC3<R CGAL_CTAG>
circumcenter( PointC3<R CGAL_CTAG> const& p,
              PointC3<R CGAL_CTAG> const& q,
              PointC3<R CGAL_CTAG> const& r,
              PointC3<R CGAL_CTAG> const& s)
{
  typename R::FT x,y,z;
  circumcenterC3(p.x(),p.y(),p.z(),
                 q.x(),q.y(),q.z(),
                 r.x(),r.y(),r.z(),
                 s.x(),s.y(),s.z(),
                 x,y,z);
  return PointC3<R CGAL_CTAG>(x,y,z);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_POINTS_3_H
