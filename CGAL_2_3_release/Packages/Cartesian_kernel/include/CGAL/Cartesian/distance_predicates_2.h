// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Cartesian/distance_predicates_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_DISTANCE_PREDICATES_2_H
#define CGAL_CARTESIAN_DISTANCE_PREDICATES_2_H

#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template <class R >
inline
Comparison_result
compare_distance_to_point(const PointC2<R CGAL_CTAG>& p,
                          const PointC2<R CGAL_CTAG>& q,
                          const PointC2<R CGAL_CTAG>& r)
{
  return cmp_dist_to_pointC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template <class R>
inline
bool
has_larger_distance_to_point(const PointC2<R CGAL_CTAG>& p,
                             const PointC2<R CGAL_CTAG>& q,
                             const PointC2<R CGAL_CTAG>& r)
{
  return has_larger_dist_to_pointC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template <class R>
inline
bool
has_smaller_distance_to_point(const PointC2<R CGAL_CTAG>& p,
                              const PointC2<R CGAL_CTAG>& q,
                              const PointC2<R CGAL_CTAG>& r)
{
  return has_smaller_dist_to_pointC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template <class R>
inline
Comparison_result
compare_signed_distance_to_line(const LineC2<R CGAL_CTAG>&  l,
                                const PointC2<R CGAL_CTAG>& p,
                                const PointC2<R CGAL_CTAG>& q)
{
  return cmp_signed_dist_to_directionC2(l.a(), l.b(), p.x(), p.y(),
                                        q.x(), q.y());
}

template <class R>
inline
bool
has_larger_signed_distance_to_line(const LineC2<R CGAL_CTAG>&  l,
                                   const PointC2<R CGAL_CTAG>& p,
                                   const PointC2<R CGAL_CTAG>& q)
{
  return has_larger_signed_dist_to_directionC2(l.a(), l.b(), p.x(), p.y(),
                                               q.x(), q.y());
}

template <class R>
inline
bool
has_smaller_signed_distance_to_line(const LineC2<R CGAL_CTAG>&  l,
                                    const PointC2<R CGAL_CTAG>& p,
                                    const PointC2<R CGAL_CTAG>& q)
{
  return has_smaller_signed_dist_to_directionC2(l.a(), l.b(), p.x(), p.y(),
                                                q.x(), q.y());
}

template <class R>
inline
Comparison_result
compare_signed_distance_to_line(const PointC2<R CGAL_CTAG>& p,
                                const PointC2<R CGAL_CTAG>& q,
                                const PointC2<R CGAL_CTAG>& r,
                                const PointC2<R CGAL_CTAG>& s)
{
  return cmp_signed_dist_to_lineC2(p.x(), p.y(), q.x(), q.y(),
                                   r.x(), r.y(), s.x(), s.y());
}

template <class R>
inline
bool
has_smaller_signed_distance_to_line(const PointC2<R CGAL_CTAG>& p,
                                    const PointC2<R CGAL_CTAG>& q,
                                    const PointC2<R CGAL_CTAG>& r,
                                    const PointC2<R CGAL_CTAG>& s)
{
  return has_smaller_signed_dist_to_lineC2(p.x(), p.y(), q.x(), q.y(),
                                           r.x(), r.y(), s.x(), s.y());
}

template <class R>
inline
bool
has_larger_signed_distance_to_line(const PointC2<R CGAL_CTAG>& p,
                                   const PointC2<R CGAL_CTAG>& q,
                                   const PointC2<R CGAL_CTAG>& r,
                                   const PointC2<R CGAL_CTAG>& s)
{
  return has_larger_signed_dist_to_lineC2(p.x(), p.y(), q.x(), q.y(),
                                          r.x(), r.y(), s.x(), s.y());
}

#ifndef CGAL_NO_DEPRECATED_CODE
template <class R >
inline
Comparison_result
cmp_dist_to_point(const PointC2<R CGAL_CTAG>& p,
                  const PointC2<R CGAL_CTAG>& q,
                  const PointC2<R CGAL_CTAG>& r)
{
  return cmp_dist_to_pointC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template <class R>
inline
bool
has_larger_dist_to_point(const PointC2<R CGAL_CTAG>& p,
                         const PointC2<R CGAL_CTAG>& q,
                         const PointC2<R CGAL_CTAG>& r)
{
  return has_larger_dist_to_pointC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template <class R>
inline
bool
has_smaller_dist_to_point(const PointC2<R CGAL_CTAG>& p,
                          const PointC2<R CGAL_CTAG>& q,
                          const PointC2<R CGAL_CTAG>& r)
{
  return has_smaller_dist_to_pointC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template <class R>
inline
Comparison_result
cmp_signed_dist_to_line(const LineC2<R CGAL_CTAG>&  l,
                        const PointC2<R CGAL_CTAG>& p,
                        const PointC2<R CGAL_CTAG>& q)
{
  return cmp_signed_dist_to_directionC2(l.a(), l.b(), p.x(), p.y(),
                                        q.x(), q.y());
}

template <class R>
inline
bool
has_larger_signed_dist_to_line(const LineC2<R CGAL_CTAG>&  l,
                               const PointC2<R CGAL_CTAG>& p,
                               const PointC2<R CGAL_CTAG>& q)
{
  return has_larger_signed_dist_to_directionC2(l.a(), l.b(), p.x(), p.y(),
                                               q.x(), q.y());
}

template <class R>
inline
bool
has_smaller_signed_dist_to_line(const LineC2<R CGAL_CTAG>&  l,
                                const PointC2<R CGAL_CTAG>& p,
                                const PointC2<R CGAL_CTAG>& q)
{
  return has_smaller_signed_dist_to_directionC2(l.a(), l.b(), p.x(), p.y(),
                                                q.x(), q.y());
}

template <class R>
inline
Comparison_result
cmp_signed_dist_to_line(const PointC2<R CGAL_CTAG>& p,
                        const PointC2<R CGAL_CTAG>& q,
                        const PointC2<R CGAL_CTAG>& r,
                        const PointC2<R CGAL_CTAG>& s)
{
  return cmp_signed_dist_to_lineC2(p.x(), p.y(), q.x(), q.y(),
                                   r.x(), r.y(), s.x(), s.y());
}

template <class R>
inline
bool
has_smaller_signed_dist_to_line(const PointC2<R CGAL_CTAG>& p,
                                const PointC2<R CGAL_CTAG>& q,
                                const PointC2<R CGAL_CTAG>& r,
                                const PointC2<R CGAL_CTAG>& s)
{
  return has_smaller_signed_dist_to_lineC2(p.x(), p.y(), q.x(), q.y(),
                                           r.x(), r.y(), s.x(), s.y());
}

template <class R>
inline
bool
has_larger_signed_dist_to_line(const PointC2<R CGAL_CTAG>& p,
                               const PointC2<R CGAL_CTAG>& q,
                               const PointC2<R CGAL_CTAG>& r,
                               const PointC2<R CGAL_CTAG>& s)
{
  return has_larger_signed_dist_to_lineC2(p.x(), p.y(), q.x(), q.y(),
                                          r.x(), r.y(), s.x(), s.y());
}
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DISTANCE_PREDICATES_2_H
