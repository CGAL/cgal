// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_PLANES_D_H
#define CGAL_CARTESIAN_CONSTRUCTIONS_ON_PLANES_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/Cartesian/Point_d.h>
#include <CGAL/Cartesian/Plane_d.h>
#include <CGAL/constructions/kernel_ftCd.h>

CGAL_BEGIN_NAMESPACE

template < class R, class InputIterator >
CGAL_KERNEL_LARGE_INLINE
PlaneCd<R CGAL_CTAG>
plane_from_points(int dim,
                  const InputIterator &first, const InputIterator &last)
{
  CGAL_kernel_precondition(last-first == dim);
  typename R::FT *h = new typename R::FT[dim+1];
  const typename R::FT **p;
  const typename R::FT **q = p = new (const typename R::FT*)[dim];
  PointCd<R CGAL_CTAG> *i;
  for (i=first; i!=last; ++i,++q) {
    CGAL_kernel_precondition( i->dimension() == dim );
    *q = i->begin();
  }
  plane_from_pointsCd(dim, p, q, h);
  return PlaneCd<R CGAL_CTAG>(dim, h, h+dim+1);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PlaneCd<R CGAL_CTAG>
plane_from_point_direction(const PointCd<R CGAL_CTAG>& p,
                           const DirectionCd<R CGAL_CTAG>& d)
{
  CGAL_kernel_precondition( p.dimension() == d.dimension() );
  typename R::FT e = new typename R::FT[p.dimension()+1];
  plane_from_point_directionCd(p.dimension(), p.begin(), p.end(),
                               d.begin(), d.end(), e);
  return PlaneCd<R CGAL_CTAG>(dim, e, e+p.dimension()+1);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointCd<R CGAL_CTAG>
point_on_plane(const PlaneCd<R CGAL_CTAG>& h, int i)
{
  typename R::FT e = new typename R::FT[h.dimension()];
  point_on_planeCd(h.dimension(), h.begin(), h.end(), i, e);
  return PointCd<R CGAL_CTAG>(h.dimension(), e, e+h.dimension());
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointCd<R CGAL_CTAG>
projection_plane(const PointCd<R CGAL_CTAG>& p,
                 const PlaneCd<R CGAL_CTAG>& h)
{
  CGAL_kernel_precondition( p.dimension() == h.dimension() );
  typename R::FT e = new typename R::FT[dim];
  projection_planeCd(h.dimension(), h.begin(), h.end(),
                     p.begin(), p.end(), e);
  return PointCd<R CGAL_CTAG>(h.dimension(), e, e+h.dimension());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_PLANES_D_H
