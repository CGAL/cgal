// ============================================================================
//
// Copyright (c) 1999  The CGAL Consortium
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
// file          : include/CGAL/Regular_triangulation_euclidean_traits_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H

// This file is based on Regular_triangulation_euclidean_traits_2.h

#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Weighted_point.h>

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_CARTESIAN_H
#include <CGAL/predicates/Regular_triangulation_ftC3.h>
#endif

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/predicates/Regular_triangulation_rtH3.h>
#endif

CGAL_BEGIN_NAMESPACE 

template < class R, class W = typename R::RT>
class Regular_triangulation_euclidean_traits_3
  : public Triangulation_geom_traits_3<R>
{
public:
  typedef R                                     Rep;
  typedef W                                     Weight;
  typedef Triangulation_geom_traits_3 <R>       Traits;
  typedef Traits::Point                         Bare_point;
  typedef Weighted_point <Bare_point, W>        Weighted_point;
  typedef Weighted_point                        Point;

  // power test for 3 dimension triangulation
  Oriented_side power_test(const Weighted_point &p,
			   const Weighted_point &q,
			   const Weighted_point &r,
			   const Weighted_point &s,
			   const Weighted_point &t) const
  {
    CGAL_triangulation_precondition( ! coplanar(p, q, r, s) );
    return CGAL::power_test(p, q, r, s, t);
  }

  // power test for 2 dimension triangulation
  Oriented_side power_test(const Weighted_point &p,
			   const Weighted_point &q,
			   const Weighted_point &r,
			   const Weighted_point &t) const
  {
    CGAL_triangulation_precondition( ! collinear(p, q, r) );
    CGAL_triangulation_precondition( coplanar(p,q,r,t) );
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
power_test(const Weighted_point<Point_3<Cartesian<FT> >, Weight> &p,
           const Weighted_point<Point_3<Cartesian<FT> >, Weight> &q,
           const Weighted_point<Point_3<Cartesian<FT> >, Weight> &r,
           const Weighted_point<Point_3<Cartesian<FT> >, Weight> &s,
           const Weighted_point<Point_3<Cartesian<FT> >, Weight> &t)
{
    return power_testC3(p.x(), p.y(), p.z(), FT(p.weight()),
                        q.x(), q.y(), q.z(), FT(q.weight()),
                        r.x(), r.y(), r.z(), FT(r.weight()),
                        s.x(), s.y(), s.z(), FT(s.weight()),
                        t.x(), t.y(), t.z(), FT(t.weight()));
}

template < class FT, class Weight >
inline
Oriented_side
power_test(const Weighted_point<Point_3<Cartesian<FT> >, Weight> &p,
           const Weighted_point<Point_3<Cartesian<FT> >, Weight> &q,
           const Weighted_point<Point_3<Cartesian<FT> >, Weight> &r,
           const Weighted_point<Point_3<Cartesian<FT> >, Weight> &t)
{
    return power_testC3(p.x(), p.y(), p.z(), FT(p.weight()),
                        q.x(), q.y(), q.z(), FT(q.weight()),
                        r.x(), r.y(), r.z(), FT(r.weight()),
                        t.x(), t.y(), t.z(), FT(t.weight()));
}

template < class FT, class Weight >
inline
Oriented_side
power_test(const Weighted_point<Point_3<Cartesian<FT> >, Weight> &p,
           const Weighted_point<Point_3<Cartesian<FT> >, Weight> &q,
           const Weighted_point<Point_3<Cartesian<FT> >, Weight> &t)
{
    return power_testC3(p.x(), p.y(), p.z(), FT(p.weight()),
                        q.x(), q.y(), q.z(), FT(q.weight()),
                        t.x(), t.y(), t.z(), FT(t.weight()));
}
#endif // CGAL_CARTESIAN_H

#ifdef CGAL_HOMOGENEOUS_H
template < class RT, class Weight >
inline
Oriented_side
power_test(const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &p,
           const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &q,
           const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &r,
           const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &s,
           const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &t)
{
    return power_testH3(p.hx(), p.hy(), p.hz(), p.hw(), RT(p.weight()),
                        q.hx(), q.hy(), q.hz(), q.hw(), RT(q.weight()),
                        r.hx(), r.hy(), r.hz(), r.hw(), RT(r.weight()),
                        s.hx(), s.hy(), s.hz(), s.hw(), RT(s.weight()),
                        t.hx(), t.hy(), t.hz(), t.hw(), RT(t.weight()));
}

template < class RT, class Weight >
inline
Oriented_side
power_test(const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &p,
           const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &q,
           const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &r,
           const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &t)
{
    return power_testH3(p.hx(), p.hy(), p.hz(), p.hw(), RT(p.weight()),
                        q.hx(), q.hy(), q.hz(), q.hw(), RT(q.weight()),
                        r.hx(), r.hy(), r.hz(), r.hw(), RT(r.weight()),
                        t.hx(), t.hy(), t.hz(), t.hw(), RT(t.weight()));
}

template < class RT, class Weight >
inline
Oriented_side
power_test(const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &p,
           const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &q,
           const Weighted_point<Point_3<Homogeneous<RT> >, Weight> &t)
{
    return power_testH3(p.hx(), p.hy(), p.hz(), p.hw(), RT(p.weight()),
                        q.hx(), q.hy(), q.hz(), q.hw(), RT(q.weight()),
                        t.hx(), t.hy(), t.hz(), t.hw(), RT(t.weight()));
}
#endif // CGAL_HOMOGENEOUS_H

CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H
