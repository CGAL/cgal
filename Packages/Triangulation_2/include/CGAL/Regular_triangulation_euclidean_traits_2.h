// ============================================================================
//
// Copyright (c) 1997  The CGAL Consortium
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
// author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Weighted_point.h>

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_CARTESIAN_H
#include <CGAL/predicates/Regular_triangulation_ftC2.h>
#endif

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/predicates/Regular_triangulation_rtH2.h>
#endif

CGAL_BEGIN_NAMESPACE 

template < class R, class W = typename R::FT>
class Regular_triangulation_euclidean_traits_2
  : public Triangulation_euclidean_traits_2<R>
{
public:
  typedef R                                     Rep;
  typedef W                                     Weight;
  typedef Triangulation_euclidean_traits_2 <R>  Traits;
  typedef typename Traits::Point                Bare_point;
  typedef Weighted_point <Bare_point, W>        Weighted_point;
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

// #ifdef CGAL_CARTESIAN_H
// template < class FT, class Weight >
// inline
// Oriented_side
// power_test(const Weighted_point<Point_2<Cartesian<FT> >, Weight> &p,
//            const Weighted_point<Point_2<Cartesian<FT> >, Weight> &q,
//            const Weighted_point<Point_2<Cartesian<FT> >, Weight> &r,
//            const Weighted_point<Point_2<Cartesian<FT> >, Weight> &t)
// {
//     return power_testC2(p.x(), p.y(), FT(p.weight()),
//                         q.x(), q.y(), FT(q.weight()),
//                         r.x(), r.y(), FT(r.weight()),
//                         t.x(), t.y(), FT(t.weight()));
// }

// template < class FT, class Weight >
// inline
// Oriented_side
// power_test(const Weighted_point<Point_2<Cartesian<FT> >, Weight> &p,
//            const Weighted_point<Point_2<Cartesian<FT> >, Weight> &q,
//            const Weighted_point<Point_2<Cartesian<FT> >, Weight> &t)
// {
//     return power_testC2(p.x(), p.y(), FT(p.weight()),
//                         q.x(), q.y(), FT(q.weight()),
//                         t.x(), t.y(), FT(t.weight()));
// }
// #endif // CGAL_CARTESIAN_H

// #ifdef CGAL_HOMOGENEOUS_H
// template < class RT, class Weight >
// inline
// Oriented_side
// power_test(const Weighted_point<Point_2<Homogeneous<RT> >, Weight> &p,
//            const Weighted_point<Point_2<Homogeneous<RT> >, Weight> &q,
//            const Weighted_point<Point_2<Homogeneous<RT> >, Weight> &r,
//            const Weighted_point<Point_2<Homogeneous<RT> >, Weight> &t)
// {
//     return power_testH2(p.hx(), p.hy(), p.hw(), RT(p.weight()),
//                         q.hx(), q.hy(), q.hw(), RT(q.weight()),
//                         r.hx(), r.hy(), r.hw(), RT(r.weight()),
//                         t.hx(), t.hy(), t.hw(), RT(t.weight()));
// }

// template < class RT, class Weight >
// inline
// Oriented_side
// power_test(const Weighted_point<Point_2<Homogeneous<RT> >, Weight> &p,
//            const Weighted_point<Point_2<Homogeneous<RT> >, Weight> &q,
//            const Weighted_point<Point_2<Homogeneous<RT> >, Weight> &t)
// {
//     return power_testH2(p.hx(), p.hy(), p.hw(), RT(p.weight()),
//                         q.hx(), q.hy(), q.hw(), RT(q.weight()),
//                         t.hx(), t.hy(), t.hw(), RT(t.weight()));
// }
// #endif // CGAL_HOMOGENEOUS_H

template < class Point, class Weight >
inline
Oriented_side
power_test(const Weighted_point<Point, Weight> &p,
           const Weighted_point<Point, Weight> &q,
           const Weighted_point<Point, Weight> &r,
           const Weighted_point<Point, Weight> &t,
	   Cartesian_tag )
{
  typedef typename Point::FT  FT;
  return power_testC2(p.x(), p.y(), FT(p.weight()),
		      q.x(), q.y(), FT(q.weight()),
		      r.x(), r.y(), FT(r.weight()),
		      t.x(), t.y(), FT(t.weight()));
}

template < class Point, class Weight >
inline
Oriented_side
power_test(const Weighted_point<Point, Weight> &p,
           const Weighted_point<Point, Weight> &q,
           const Weighted_point<Point, Weight> &r,
           const Weighted_point<Point, Weight> &t,
	   Homogeneous_tag )
{
  typedef typename Point::RT  RT;
  return power_testH2(p.hx(), p.hy(), p.hw(), RT(p.weight()),
		      q.hx(), q.hy(), q.hw(), RT(q.weight()),
		      r.hx(), r.hy(), r.hw(), RT(r.weight()),
		      t.hx(), t.hy(), t.hw(), RT(t.weight()));
}

template < class Point, class Weight >
inline
Oriented_side
power_test(const Weighted_point<Point, Weight> &p,
           const Weighted_point<Point, Weight> &q,
           const Weighted_point<Point, Weight> &r,
           const Weighted_point<Point, Weight> &t)
{
  typedef typename Point::R::Rep_tag Tag;
  return power_test(p, q, r, t, Tag());
}
  
template < class Point, class Weight >
inline
Oriented_side
power_test(const Weighted_point<Point, Weight> &p,
           const Weighted_point<Point, Weight> &q,
           const Weighted_point<Point, Weight> &t,
	   Cartesian_tag )
{
    typedef typename Point::FT  FT;
    return power_testC2(p.x(), p.y(), FT(p.weight()),
                        q.x(), q.y(), FT(q.weight()),
                        t.x(), t.y(), FT(t.weight()));
}


template < class Point, class Weight >
inline
Oriented_side
power_test(const Weighted_point<Point, Weight> &p,
           const Weighted_point<Point, Weight> &q,
           const Weighted_point<Point, Weight> &t,
	   Homogeneous_tag )
{
   typedef typename Point::RT  RT;
    return power_testH2(p.hx(), p.hy(), p.hw(), RT(p.weight()),
                        q.hx(), q.hy(), q.hw(), RT(q.weight()),
                        t.hx(), t.hy(), t.hw(), RT(t.weight()));
}

template < class Point, class Weight >
inline
Oriented_side
power_test(const Weighted_point<Point, Weight> &p,
           const Weighted_point<Point, Weight> &q,
           const Weighted_point<Point, Weight> &t)
{
  typedef typename Point::R::Rep_tag Tag;
  return power_test(p, q, t, Tag());
}
CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
