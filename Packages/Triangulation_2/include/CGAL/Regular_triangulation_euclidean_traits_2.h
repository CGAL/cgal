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

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/predicates/Regular_triangulation_ftC2.h>
#include <CGAL/constructions_on_weighted_points_cartesian_2.h>
#endif

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/predicates/Regular_triangulation_rtH2.h>
#include <CGAL/constructions_on_weighted_points_homogeneous_2.h>
#endif

CGAL_BEGIN_NAMESPACE 

// constructions for DUALITY: weighted_circumcenter and radical 
//  axis
template < class Point, class We >
inline
Point
weighted_circumcenter(const Weighted_point< Point,We >& p,
		      const Weighted_point< Point,We >& q,
		      const Weighted_point< Point,We >& r,
		      Cartesian_tag )
{
  
  typename Point::RT x,y;
  weighted_circumcenterC2(p.x(),p.y(),p.weight(),
			  q.x(),q.y(),q.weight(),
			  r.x(),r.y(),r.weight(),x,y);
  return Point(x,y);
}

template < class Point, class We >
inline
Point
weighted_circumcenter(const Weighted_point< Point,We >& p,
		      const Weighted_point< Point,We >& q,
		      const Weighted_point< Point,We >& r,
		      Homogeneous_tag )
{
  
  typename Point::RT x,y,w;
  weighted_circumcenterH2(p.hx(),p.hy(),p.hw(),p.weight(),
			  q.hx(),q.hy(),q.hw(),q.weight(),
			  r.hx(),r.hy(),r.hw(),r.weight(),
			  x,y,w);
  return Point(x,y,w);
}


template < class Point, class We >
inline
Point
weighted_circumcenter(const Weighted_point< Point,We >& p,
		      const Weighted_point< Point,We >& q,
		      const Weighted_point< Point,We >& r)
{
  typedef typename Point::R::Rep_tag Tag;
  return weighted_circumcenter(p, q, r, Tag()); 
}


template <class Point, class Weight>
class Construct_weighted_circumcenter_2
{
public:
  typedef Weighted_point <Point, Weight>        Weighted_point;
  Point operator() ( Weighted_point p,
		     Weighted_point q,
		     Weighted_point r) 
    {
      CGAL_triangulation_precondition( ! collinear(p, q, r) );
      return CGAL::weighted_circumcenter(p,q,r);
    }
};
 


template < class Point, class We >
inline
Line_2<typename Point::R>
radical_axis(const Weighted_point< Point,We >& p,
	     const Weighted_point< Point,We >& q,
	     Cartesian_tag )
{
  typedef typename Point::RT RT;
  typedef typename Point::R  Rep;
  RT a,b,c;
  radical_axisC2(p.x(),p.y(),p.weight(),q.x(),q.y(),q.weight(),a,b,c);
  return Line_2<Rep>(a,b,c);
}

template < class Point, class We >
inline
Line_2<typename Point::R>
radical_axis(const Weighted_point< Point,We >& p,
	     const Weighted_point< Point,We >& q,
	      Homogeneous_tag)
{
  typedef typename Point::RT RT;
  typedef typename Point::R  Rep;
  RT a,b,c;
  radical_axisH2(p.hx(),p.hy(), p.hw(), p.weight(),
		 q.hx(),q.hy(), q.hw(), q.weight(),a,b,c);
  return Line_2<Rep>(a,b,c);
}

template < class Point, class We >
inline
Line_2<typename Point::R>
radical_axis(const Weighted_point< Point,We >& p,
	     const Weighted_point< Point,We >& q)
{
  typedef typename Point::R::Rep_tag Tag;
  return radical_axis(p, q, Tag()); 
}


template <class Point, class Weight>
class Construct_radical_axis_2
{
public:
  typedef Weighted_point <Point, Weight>   Weighted_point;
  typedef typename Point::R  R;

  Line_2<R> operator() ( Weighted_point p, Weighted_point q) 
  {
    return CGAL::radical_axis(p,q);
  }
};

///-----------------------------------------------------------


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

template <class Point, class Weight>
class Power_test_2
{
public:
  typedef Weighted_point <Point, Weight>        Weighted_point;
  Oriented_side operator() ( Weighted_point p,
			     Weighted_point q,
			     Weighted_point r,
			     Weighted_point s) 
    {
      //CGAL_triangulation_precondition( ! collinear(p, q, r) );
      return CGAL::power_test(p,q,r,s);
    }
};

template <class Point, class Weight>
class Power_test_degenerated_2
{
public:
  typedef Weighted_point <Point, Weight>        Weighted_point;
  Oriented_side operator() ( Weighted_point p,
			     Weighted_point q,
			     Weighted_point r)
    {
      //CGAL_triangulation_precondition( collinear(p, q, r) );
      //CGAL_triangulation_precondition( p.point() != q.point() );
      return CGAL::power_test(p,q,r);
    }
};

template < class R, class W = typename R::FT>
class Regular_triangulation_euclidean_traits_2
  : public Triangulation_euclidean_traits_2<R>
{
public:
  typedef R                                     Rep;
  typedef W                                     Weight;
  typedef Triangulation_euclidean_traits_2 <R>  Traits;
  typedef typename Traits::Point_2              Bare_point;
  typedef Weighted_point <Bare_point, W>        Weighted_point;
  // This is required for the point() function of vertex base class
  // to be correctly return a weighted_point;
  // patch 27/11/00
  typedef Weighted_point                        Point_2;

  typedef CGAL::Power_test_2<Bare_point, W>     Power_test_2;
  typedef CGAL::Power_test_degenerated_2<Bare_point, W>  
                                                Power_test_degenerated_2;
 //concstruction objects
  typedef CGAL::Construct_weighted_circumcenter_2<Bare_point, W> 
                                            Construct_weighted_circumcenter_2;
  typedef CGAL::Construct_radical_axis_2<Bare_point, W> 
                                            Construct_radical_axis_2;
  
  Power_test_2 
  power_test_2_object() const
    {  return Power_test_2();}

  Power_test_degenerated_2
  power_test_degenerated_2_object() const
    {return Power_test_degenerated_2();}

  //constructions for dual:
  Construct_weighted_circumcenter_2
  construct_weighted_circumcenter_2_object() const
    {return Construct_weighted_circumcenter_2();}
  
  Construct_radical_axis_2
  construct_radical_axis_2_object() const
    {return Construct_radical_axis_2();}

};
 
CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
