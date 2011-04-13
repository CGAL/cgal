// revision      : 2.8
// revision_date : 28 Oct 1999 
// author(s)     : Hervé Brönnimann

#ifndef CGAL_DISTANCE_PREDICATES_D_H
#define CGAL_DISTANCE_PREDICATES_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/distance_predicatesHd.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Cartesian/distance_predicates_d.h>
#endif // CGAL_CARTESIAN_H

#include <CGAL/Point_d.h>
#include <CGAL/Plane_d.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Comparison_result
cmp_dist_to_point( const Point_d<R> &p,
                   const Point_d<R> &q,
                   const Point_d<R> &r)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return cmp_dist_to_point((const RPoint_d& )p,
                           (const RPoint_d& )q,
                           (const RPoint_d& )r);
}

template < class R >
inline
bool
has_larger_dist_to_point( const Point_d<R> &p,
                          const Point_d<R> &q,
                          const Point_d<R> &r)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return has_larger_dist_to_point((const RPoint_d& )p,
                                  (const RPoint_d& )q,
                                  (const RPoint_d& )r);
}

template < class R >
inline
bool
has_smaller_dist_to_point( const Point_d<R> &p,
                           const Point_d<R> &q,
                           const Point_d<R> &r)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return has_smaller_dist_to_point((const RPoint_d& )p,
                                   (const RPoint_d& )q,
                                   (const RPoint_d& )r);
}

template < class R >
inline
Comparison_result
cmp_signed_dist_to_plane( const Point_d<R> &p,
                          const Point_d<R> &q,
                          const Point_d<R> &r)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return cmp_signed_dist_to_plane((const RPoint_d& )p,
                                  (const RPoint_d& )q,
                                  (const RPoint_d& )r);
}

template < class R >
inline
bool
has_larger_signed_dist_to_plane( const Point_d<R> &p,
                                 const Point_d<R> &q,
                                 const Point_d<R> &r)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return has_larger_signed_dist_to_plane((const RPoint_d& )p,
                                         (const RPoint_d& )q,
                                         (const RPoint_d& )r);
}

template < class R >
inline
bool
has_smaller_signed_dist_to_plane( const Point_d<R> &p,
                                  const Point_d<R> &q,
                                  const Point_d<R> &r)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return has_smaller_signed_dist_to_plane((const RPoint_d& )p,
                                          (const RPoint_d& )q,
                                          (const RPoint_d& )r);
}

CGAL_END_NAMESPACE

#endif //CGAL_DISTANCE_PREDICATES_D_H
