// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 16
//
// source        : webS2/S2.lw
// file          : include/CGAL/SimpleCartesian/distance_predicatesS2.h
// package       : S2 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.6
// revision_date : 27 Jun 2000
// author(s)     : Stefan Schirra
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================


#ifndef CGAL_DISTANCE_PREDICATESS2_H
#define CGAL_DISTANCE_PREDICATESS2_H

#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template <class FT >
inline
Comparison_result
cmp_dist_to_point(const PointS2<FT>& p,
                  const PointS2<FT>& q,
                  const PointS2<FT>& r)
{
  return cmp_dist_to_pointC2(p.x(),p.y(),q.x(),q.y(),r.x(),r.y());
}

template <class FT>
inline
bool
has_larger_dist_to_point(const PointS2<FT>& p,
                         const PointS2<FT>& q,
                         const PointS2<FT>& r)
{
  return has_larger_dist_to_pointC2(p.x(),p.y(),q.x(),q.y(),r.x(),r.y());
}

template <class FT>
inline
bool
has_smaller_dist_to_point(const PointS2<FT>& p,
                          const PointS2<FT>& q,
                          const PointS2<FT>& r)
{
  return has_smaller_dist_to_pointC2(p.x(),p.y(),q.x(),q.y(),r.x(),r.y());
}
template <class FT>
inline
Comparison_result
cmp_signed_dist_to_line(const LineS2<FT>&  l,
                        const PointS2<FT>& p,
                        const PointS2<FT>& q)
{
  return cmp_signed_dist_to_directionC2
           (l.a(),l.b(),p.x(),p.y(),q.x(),q.y());
}

template <class FT>
inline
bool
has_larger_signed_dist_to_line(const LineS2<FT>&  l,
                               const PointS2<FT>& p,
                               const PointS2<FT>& q)
{
  return has_larger_signed_dist_to_directionC2
           (l.a(),l.b(),p.x(),p.y(),q.x(),q.y());
}

template <class FT>
inline
bool
has_smaller_signed_dist_to_line(const LineS2<FT>&  l,
                                const PointS2<FT>& p,
                                const PointS2<FT>& q)
{
  return has_smaller_signed_dist_to_directionC2
           (l.a(),l.b(),p.x(),p.y(),q.x(),q.y());
}
template <class FT>
inline
Comparison_result
cmp_signed_dist_to_line(const PointS2<FT>& p,
                        const PointS2<FT>& q,
                        const PointS2<FT>& r,
                        const PointS2<FT>& s)
{
  return cmp_signed_dist_to_lineC2
           (p.x(),p.y(),q.x(),q.y(),r.x(),r.y(),s.x(),s.y());
}

template <class FT>
inline
bool
has_smaller_signed_dist_to_line(const PointS2<FT>& p,
                                const PointS2<FT>& q,
                                const PointS2<FT>& r,
                                const PointS2<FT>& s)
{
  return has_smaller_signed_dist_to_lineC2
           (p.x(),p.y(),q.x(),q.y(),r.x(),r.y(),s.x(),s.y());
}

template <class FT>
inline
bool
has_larger_signed_dist_to_line(const PointS2<FT>& p,
                               const PointS2<FT>& q,
                               const PointS2<FT>& r,
                               const PointS2<FT>& s)
{
  return has_larger_signed_dist_to_lineC2
           (p.x(),p.y(),q.x(),q.y(),r.x(),r.y(),s.x(),s.y());
}


CGAL_END_NAMESPACE

#endif //CGAL_DISTANCE_PREDICATESS2_H
