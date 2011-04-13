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
//                 Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Weighted_point.h>

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/predicates/Regular_triangulation_ftC3.h>
#endif

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/predicates/Regular_triangulation_rtH3.h>
#endif

CGAL_BEGIN_NAMESPACE 

template <class Point, class Weight>
class Power_test_3
{
public:
  typedef Oriented_side  result_type;
  typedef Weighted_point <Point, Weight>        Weighted_point;

  Oriented_side operator() ( const Weighted_point & p,
			     const Weighted_point & q,
			     const Weighted_point & r,
			     const Weighted_point & s,
			     const Weighted_point & t) const
    {
      return CGAL::power_test(p,q,r,s,t);
    }
  Oriented_side operator() ( const Weighted_point & p,
			     const Weighted_point & q,
			     const Weighted_point & r,
			     const Weighted_point & s) const
    {
      return CGAL::power_test(p,q,r,s);
    }

  Oriented_side operator() ( const Weighted_point & p,
			     const Weighted_point & q,
			     const Weighted_point & r) const
    {
      return CGAL::power_test(p,q,r);
    }
};



template < class Repres, class Weight = CGAL_TYPENAME_MSVC_NULL Repres::RT >
class Regular_triangulation_euclidean_traits_3
  : public Repres
{
public:
  typedef typename Repres::Point_3 Bare_point;
  typedef Weighted_point <Bare_point, Weight>   Weighted_point;
  typedef Weighted_point                        Point_3;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_3 Point;

  typedef Power_test_3<Bare_point, Weight> Power_test_3;
  
  Power_test_3 
  power_test_3_object() const
    {  return Power_test_3();}

};

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt, Weight> &p,
           const Weighted_point<pt, Weight> &q,
           const Weighted_point<pt, Weight> &r,
           const Weighted_point<pt, Weight> &s,
           const Weighted_point<pt, Weight> &t,
	   Cartesian_tag)
{
  typedef typename pt::FT FT;
    return power_testC3(p.x(), p.y(), p.z(), FT(p.weight()),
                        q.x(), q.y(), q.z(), FT(q.weight()),
                        r.x(), r.y(), r.z(), FT(r.weight()),
                        s.x(), s.y(), s.z(), FT(s.weight()),
                        t.x(), t.y(), t.z(), FT(t.weight()));
}

template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt, Weight> &p,
           const Weighted_point<pt, Weight> &q,
           const Weighted_point<pt, Weight> &r,
           const Weighted_point<pt, Weight> &t,
	   Cartesian_tag)
{
  typedef typename pt::FT FT;
    return power_testC3(p.x(), p.y(), p.z(), FT(p.weight()),
                        q.x(), q.y(), q.z(), FT(q.weight()),
                        r.x(), r.y(), r.z(), FT(r.weight()),
                        t.x(), t.y(), t.z(), FT(t.weight()));
}

template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt, Weight> &p,
           const Weighted_point<pt, Weight> &q,
           const Weighted_point<pt, Weight> &t,
	   Cartesian_tag)
{
  typedef typename pt::FT FT;
    return power_testC3(p.x(), p.y(), p.z(), FT(p.weight()),
                        q.x(), q.y(), q.z(), FT(q.weight()),
                        t.x(), t.y(), t.z(), FT(t.weight()));
}
#endif // CGAL_CARTESIAN_H

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt, Weight> &p,
           const Weighted_point<pt, Weight> &q,
           const Weighted_point<pt, Weight> &r,
           const Weighted_point<pt, Weight> &s,
           const Weighted_point<pt, Weight> &t,
	   Homogeneous_tag)
{
  typedef typename pt::RT RT;
    return power_testH3(p.hx(), p.hy(), p.hz(), p.hw(), RT(p.weight()),
                        q.hx(), q.hy(), q.hz(), q.hw(), RT(q.weight()),
                        r.hx(), r.hy(), r.hz(), r.hw(), RT(r.weight()),
                        s.hx(), s.hy(), s.hz(), s.hw(), RT(s.weight()),
                        t.hx(), t.hy(), t.hz(), t.hw(), RT(t.weight()));
}

// The 2 following call the cartesian version over FT, because an homogeneous
// special version would be boring to write.

template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt, Weight> &p,
           const Weighted_point<pt, Weight> &q,
           const Weighted_point<pt, Weight> &r,
           const Weighted_point<pt, Weight> &t,
	   Homogeneous_tag)
{
    typedef typename pt::FT FT;
    return power_testC3(p.x(), p.y(), p.z(), FT(p.weight()),
                        q.x(), q.y(), q.z(), FT(q.weight()),
                        r.x(), r.y(), r.z(), FT(r.weight()),
                        t.x(), t.y(), t.z(), FT(t.weight()));
}

template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt, Weight> &p,
           const Weighted_point<pt, Weight> &q,
           const Weighted_point<pt, Weight> &t,
	   Homogeneous_tag)
{
    typedef typename pt::FT FT;
    return power_testC3(p.x(), p.y(), p.z(), FT(p.weight()),
                        q.x(), q.y(), q.z(), FT(q.weight()),
                        t.x(), t.y(), t.z(), FT(t.weight()));
}
#endif // CGAL_HOMOGENEOUS_H

// Kludges for M$.

template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt,Weight> &p,
	   const Weighted_point<pt,Weight> &q,
	   const Weighted_point<pt,Weight> &r,
	   const Weighted_point<pt,Weight> &s,
	   const Weighted_point<pt,Weight> &t)
{
  typedef typename pt::R::Rep_tag Tag;
  return( power_test(p,q,r,s,t, Tag()) );
}

template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt,Weight> &p,
	   const Weighted_point<pt,Weight> &q,
	   const Weighted_point<pt,Weight> &r,
	   const Weighted_point<pt,Weight> &t)
{
  typedef typename pt::R::Rep_tag Tag;
  return( power_test(p,q,r,t, Tag()) );
}

template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt,Weight> &p,
	   const Weighted_point<pt,Weight> &q,
	   const Weighted_point<pt,Weight> &t)
{
  typedef typename pt::R::Rep_tag Tag;
  return( power_test(p,q,t, Tag()) );
}

CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H
