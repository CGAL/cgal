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
#include <CGAL/representation_tags.h>

#include <CGAL/predicates/Regular_triangulation_ftC3.h>
#include <CGAL/predicates/Regular_triangulation_rtH3.h>

CGAL_BEGIN_NAMESPACE 

template <class Point, class Weight>
class Power_test_3
{
public:
  typedef Oriented_side  result_type;
  typedef CGAL::Weighted_point <Point, Weight>        Weighted_point;

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

  Oriented_side operator() ( const Weighted_point & p,
			     const Weighted_point & q) const
    {
      return CGAL::power_test(p,q);
    }
};



template < class K, class Weight = CGAL_TYPENAME_MSVC_NULL K::RT >
class Regular_triangulation_euclidean_traits_3
  : public K
{
public:
  typedef typename K::Point_3                        Bare_point;
  typedef CGAL::Weighted_point<Bare_point, Weight>   Weighted_point;
  typedef Weighted_point                             Point_3;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_3 Point;

  typedef CGAL::Power_test_3<Bare_point, Weight> Power_test_3;
  
  Power_test_3 
  power_test_3_object() const
    {  return Power_test_3();}

};


// Cartesian versions.
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
  typedef typename pt::R::FT FT;
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
  typedef typename pt::R::FT FT;
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
  typedef typename pt::R::FT FT;
    return power_testC3(p.x(), p.y(), p.z(), FT(p.weight()),
                        q.x(), q.y(), q.z(), FT(q.weight()),
                        t.x(), t.y(), t.z(), FT(t.weight()));
}

template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt, Weight> &p,
           const Weighted_point<pt, Weight> &q,
	   Cartesian_tag)
{
    typedef typename pt::R::FT FT;
    return power_testC3(FT(p.weight()),
                        FT(q.weight()));
}


// Homogeneous versions.
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
  typedef typename pt::R::RT RT;
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
    typedef typename pt::R::FT FT;
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
    typedef typename pt::R::FT FT;
    return power_testC3(p.x(), p.y(), p.z(), FT(p.weight()),
                        q.x(), q.y(), q.z(), FT(q.weight()),
                        t.x(), t.y(), t.z(), FT(t.weight()));
}

template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt, Weight> &p,
           const Weighted_point<pt, Weight> &q,
	   Homogeneous_tag)
{
    typedef typename pt::R::FT FT;
    return power_testC3(FT(p.weight()),
                        FT(q.weight()));
}


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

template < class pt, class Weight >
inline
Oriented_side
power_test(const Weighted_point<pt,Weight> &p,
	   const Weighted_point<pt,Weight> &q)
{
  typedef typename pt::R::Rep_tag Tag;
  return( power_test(p,q, Tag()) );
}

CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H
