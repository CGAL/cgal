// ============================================================================
//
// Copyright (c) 1997-1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Regular_triangulation_euclidean_traits_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Weighted_point_2.h>

#ifdef CGAL_CARTESIAN_H
#include <CGAL/predicates/Regular_triangulation_ftC2.h>
#endif

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/predicates/Regular_triangulation_rtH2.h>
#endif

CGAL_BEGIN_NAMESPACE 

template < class R, class W >
class Regular_triangulation_euclidean_traits_2
  : public Triangulation_euclidean_traits_2<R>
{
public:
  typedef R                                     Rep;
  typedef W                                     Weight;
  typedef Triangulation_euclidean_traits_2 <R>  Traits;
  typedef Traits::Point                         Bare_point;
  typedef Traits::Segment                       Segment;
  typedef Traits::Triangle                      Triangle;
  typedef Traits::Line                          Line;
  typedef Traits::Direction                     Direction;
  typedef Traits::Ray                           Ray;
  typedef Weighted_point_2 <Bare_point, W>      Weighted_point;
  typedef Weighted_point                        Point;

  // power test for 2 dimension triangulation
  Oriented_side power_test(const Weighted_point &p,
			   const Weighted_point &q,
			   const Weighted_point &r,
			   const Weighted_point &t) const
  {
    CGAL_triangulation_precondition( ! collinear(p, q, r) );
    return CGAL::power_test(p, q, r, t);
  }

  // power test for 1 dimension triangulation
  Oriented_side power_test(const Weighted_point &p,
			   const Weighted_point &q,
			   const Weighted_point &t) const
  {
    CGAL_triangulation_precondition( collinear(p, q, t) );
    CGAL_triangulation_precondition( p.point() != q.point() );
    return CGAL::power_test(p, q, t);
  }
};

#ifdef CGAL_CARTESIAN_H
template < class FT, class Weight >
inline
Oriented_side
power_test(const Weighted_point_2<Point_2<Cartesian<FT> >, Weight> &p,
           const Weighted_point_2<Point_2<Cartesian<FT> >, Weight> &q,
           const Weighted_point_2<Point_2<Cartesian<FT> >, Weight> &r,
           const Weighted_point_2<Point_2<Cartesian<FT> >, Weight> &t)
{
    return power_testC2(p.x(), p.y(), FT(p.weight()),
                        q.x(), q.y(), FT(q.weight()),
                        r.x(), r.y(), FT(r.weight()),
                        t.x(), t.y(), FT(t.weight()));
}

template < class FT, class Weight >
inline
Oriented_side
power_test(const Weighted_point_2<Point_2<Cartesian<FT> >, Weight> &p,
           const Weighted_point_2<Point_2<Cartesian<FT> >, Weight> &q,
           const Weighted_point_2<Point_2<Cartesian<FT> >, Weight> &t)
{
    return power_testC2(p.x(), p.y(), FT(p.weight()),
                        q.x(), q.y(), FT(q.weight()),
                        t.x(), t.y(), FT(t.weight()));
}
#endif // CGAL_CARTESIAN_H

#ifdef CGAL_HOMOGENEOUS_H
template < class RT, class Weight >
inline
Oriented_side
power_test(const Weighted_point_2<Point_2<Homogeneous<RT> >, Weight> &p,
           const Weighted_point_2<Point_2<Homogeneous<RT> >, Weight> &q,
           const Weighted_point_2<Point_2<Homogeneous<RT> >, Weight> &r,
           const Weighted_point_2<Point_2<Homogeneous<RT> >, Weight> &t)
{
    return power_testH2(p.x(), p.y(), p.w(), RT(p.weight()),
                        q.x(), q.y(), q.w(), RT(q.weight()),
                        r.x(), r.y(), r.w(), RT(r.weight()),
                        t.x(), t.y(), t.w(), RT(t.weight()));
}

template < class RT, class Weight >
inline
Oriented_side
power_test(const Weighted_point_2<Point_2<Homogeneous<RT> >, Weight> &p,
           const Weighted_point_2<Point_2<Homogeneous<RT> >, Weight> &q,
           const Weighted_point_2<Point_2<Homogeneous<RT> >, Weight> &t)
{
    return power_testH2(p.x(), p.y(), p.w(), RT(p.weight()),
                        q.x(), q.y(), q.w(), RT(q.weight()),
                        t.x(), t.y(), t.w(), RT(t.weight()));
}
#endif // CGAL_HOMOGENEOUS_H

CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
