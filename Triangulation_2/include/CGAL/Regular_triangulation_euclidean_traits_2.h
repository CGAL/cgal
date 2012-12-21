// Copyright (c) 1997   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                 Sylvain Pion

#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

#include <CGAL/Weighted_point.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/predicates/Regular_triangulation_ftC2.h>
#include <CGAL/constructions_on_weighted_points_cartesian_2.h>

#include <CGAL/predicates/Regular_triangulation_rtH2.h>
#include <CGAL/constructions_on_weighted_points_homogeneous_2.h>

namespace CGAL { 

// constructions for DUALITY: weighted_circumcenter and radical 
//  axis
template < class Bare_point, class We >
inline
Bare_point
weighted_circumcenter(const Weighted_point< Bare_point,We >& p,
		      const Weighted_point< Bare_point,We >& q,
		      const Weighted_point< Bare_point,We >& r,
		      Cartesian_tag )
{
  typename Kernel_traits<Bare_point>::Kernel::RT x,y;
  weighted_circumcenterC2(p.x(),p.y(),p.weight(),
			  q.x(),q.y(),q.weight(),
			  r.x(),r.y(),r.weight(),x,y);
  return Bare_point(x,y);
}

template < class Bare_point, class We >
inline
Bare_point
weighted_circumcenter(const Weighted_point< Bare_point,We >& p,
		      const Weighted_point< Bare_point,We >& q,
		      const Weighted_point< Bare_point,We >& r,
		      Homogeneous_tag )
{
  typename Kernel_traits<Bare_point>::Kernel::RT x,y,w;
  weighted_circumcenterH2(p.hx(),p.hy(),p.hw(),p.weight(),
			  q.hx(),q.hy(),q.hw(),q.weight(),
			  r.hx(),r.hy(),r.hw(),r.weight(),
			  x,y,w);
  return Bare_point(x,y,w);
}


template < class Bare_point, class We >
inline
Bare_point
weighted_circumcenter(const Weighted_point< Bare_point,We >& p,
		      const Weighted_point< Bare_point,We >& q,
		      const Weighted_point< Bare_point,We >& r)
{
  typedef typename Kernel_traits<Bare_point>::Kernel::Rep_tag Tag;
  return weighted_circumcenter(p, q, r, Tag()); 
}


template < typename K >
class Construct_weighted_circumcenter_2
{
public:
  typedef typename K::Weighted_point_2         Weighted_point_2;
  typedef typename K::Bare_point               Bare_point;

  typedef Bare_point       result_type;

  Bare_point operator() ( const Weighted_point_2 & p,
		          const Weighted_point_2 & q,
		          const Weighted_point_2 & r) const
    {
      CGAL_triangulation_precondition( ! collinear(p, q, r) );
      return CGAL::weighted_circumcenter(p,q,r);
    }
};
 


template < class Bare_point, class We >
inline
Line_2<typename Kernel_traits<Bare_point>::Kernel>
radical_axis(const Weighted_point< Bare_point,We >& p,
	     const Weighted_point< Bare_point,We >& q,
	     Cartesian_tag )
{
  typedef typename Kernel_traits<Bare_point>::Kernel::RT RT;
  typedef typename Kernel_traits<Bare_point>::Kernel     Rep;
  RT a,b,c;
  radical_axisC2(p.x(),p.y(),p.weight(),q.x(),q.y(),q.weight(),a,b,c);
  return Line_2<Rep>(a,b,c);
}

template < class Bare_point, class We >
inline
Line_2<typename Kernel_traits<Bare_point>::Kernel>
radical_axis(const Weighted_point< Bare_point,We >& p,
	     const Weighted_point< Bare_point,We >& q,
	      Homogeneous_tag)
{
  typedef typename Kernel_traits<Bare_point>::Kernel::RT RT;
  typedef typename Kernel_traits<Bare_point>::Kernel     Rep;
  RT a,b,c;
  radical_axisH2(p.hx(),p.hy(), p.hw(), p.weight(),
		 q.hx(),q.hy(), q.hw(), q.weight(),a,b,c);
  return Line_2<Rep>(a,b,c);
}

template < class Bare_point, class We >
inline
Line_2<typename Kernel_traits<Bare_point>::Kernel>
radical_axis(const Weighted_point< Bare_point,We >& p,
	     const Weighted_point< Bare_point,We >& q)
{
  typedef typename Kernel_traits<Bare_point>::Kernel::Rep_tag Tag;
  return radical_axis(p, q, Tag()); 
}


template < typename K >
class Construct_radical_axis_2
{
public:
  typedef typename K::Weighted_point_2                Weighted_point_2;
  typedef typename K::Line_2                          Line_2;

  typedef Line_2           result_type;

  Line_2
  operator() ( const Weighted_point_2 & p, const Weighted_point_2 & q) const
  {
    return CGAL::radical_axis(p,q);
  }
};

///-----------------------------------------------------------

template < class Bare_point, class Weight >
inline
typename Kernel_traits<Bare_point>::Kernel::Comparison_result
compare_power_distance(const Weighted_point<Bare_point, Weight>& p,
		       const Weighted_point<Bare_point, Weight>& q,
		       const Bare_point& r, Cartesian_tag)
{
  typedef typename Kernel_traits<Bare_point>::Kernel::FT  FT;
  return compare_power_distanceC2(p.x(), p.y(), FT(p.weight()),
				  q.x(), q.y(), FT(q.weight()),
				  r.x(), r.y());
}

template < class Bare_point, class Weight >
inline
typename Kernel_traits<Bare_point>::Kernel::Comparison_result
compare_power_distance(const Weighted_point<Bare_point, Weight>& p,
		       const Weighted_point<Bare_point, Weight>& q,
		       const Bare_point& r, Homogeneous_tag)
{
  typedef typename Kernel_traits<Bare_point>::Kernel::RT  RT;
  return compare_power_distanceH2(p.hx(), p.hy(), p.hw(), RT(p.weight()),
				  q.hx(), q.hy(), q.hw(), RT(q.weight()),
				  r.hx(), r.hy(), r.hw());
}

template < class Bare_point, class Weight >
inline
typename Kernel_traits<Bare_point>::Kernel::Comparison_result
compare_power_distance(const Weighted_point<Bare_point, Weight>& p,
		       const Weighted_point<Bare_point, Weight>& q,
		       const Bare_point& r)
{
  typedef typename Kernel_traits<Bare_point>::Kernel::Rep_tag Tag;
  return compare_power_distance(p, q, r, Tag());
}

template < typename K >
class Compare_power_distance_2
{
public:
  typedef typename K::Weighted_point_2         Weighted_point_2;
  typedef typename K::Point_2                  Point_2;
  typedef typename K::Comparison_result        Comparison_result;

  typedef Comparison_result   result_type;

  Comparison_result operator()(const Point_2& p,
			       const Weighted_point_2& q,
			       const Weighted_point_2& r) const
  {
    return CGAL::compare_power_distance(q, r, p);
  }
};

///-----------------------------------------------------------

template < class Bare_point, class Weight >
inline
typename Kernel_traits<Bare_point>::Kernel::Oriented_side
power_test_2(const Weighted_point<Bare_point, Weight> &p,
           const Weighted_point<Bare_point, Weight> &q,
           const Weighted_point<Bare_point, Weight> &r,
           const Weighted_point<Bare_point, Weight> &t,
	   Cartesian_tag )
{
  typedef typename Kernel_traits<Bare_point>::Kernel::FT  FT;
  return power_testC2(p.x(), p.y(), FT(p.weight()),
		      q.x(), q.y(), FT(q.weight()),
		      r.x(), r.y(), FT(r.weight()),
		      t.x(), t.y(), FT(t.weight()));
}

template < class Bare_point, class Weight >
inline
typename Kernel_traits<Bare_point>::Kernel::Oriented_side
power_test_2(const Weighted_point<Bare_point, Weight> &p,
           const Weighted_point<Bare_point, Weight> &q,
           const Weighted_point<Bare_point, Weight> &r,
           const Weighted_point<Bare_point, Weight> &t,
	   Homogeneous_tag )
{
  typedef typename Kernel_traits<Bare_point>::Kernel::RT  RT;
  return power_testH2(p.hx(), p.hy(), p.hw(), RT(p.weight()),
		      q.hx(), q.hy(), q.hw(), RT(q.weight()),
		      r.hx(), r.hy(), r.hw(), RT(r.weight()),
		      t.hx(), t.hy(), t.hw(), RT(t.weight()));
}

template < class Bare_point, class Weight >
inline
typename Kernel_traits<Bare_point>::Kernel::Oriented_side
power_test_2(const Weighted_point<Bare_point, Weight> &p,
           const Weighted_point<Bare_point, Weight> &q,
           const Weighted_point<Bare_point, Weight> &r,
           const Weighted_point<Bare_point, Weight> &t)
{
  typedef typename Kernel_traits<Bare_point>::Kernel::Rep_tag Tag;
  return power_test_2(p, q, r, t, Tag());
}
  
template < class Bare_point, class Weight >
inline
typename Kernel_traits<Bare_point>::Kernel::Oriented_side
power_test_2(const Weighted_point<Bare_point, Weight> &p,
	     const Weighted_point<Bare_point, Weight> &q,
	     const Weighted_point<Bare_point, Weight> &t,
	     Cartesian_tag )
{
    typedef typename Kernel_traits<Bare_point>::Kernel::FT  FT;
    return power_testC2(p.x(), p.y(), FT(p.weight()),
                        q.x(), q.y(), FT(q.weight()),
                        t.x(), t.y(), FT(t.weight()));
}


template < class Bare_point, class Weight >
inline
typename Kernel_traits<Bare_point>::Kernel::Oriented_side
power_test_2(const Weighted_point<Bare_point, Weight> &p,
	     const Weighted_point<Bare_point, Weight> &q,
	     const Weighted_point<Bare_point, Weight> &t,
	     Homogeneous_tag )
{
   typedef typename Kernel_traits<Bare_point>::Kernel::RT  RT;
    return power_testH2(p.hx(), p.hy(), p.hw(), RT(p.weight()),
                        q.hx(), q.hy(), q.hw(), RT(q.weight()),
                        t.hx(), t.hy(), t.hw(), RT(t.weight()));
}

template < class Bare_point, class Weight >
inline
typename Kernel_traits<Bare_point>::Kernel::Oriented_side
power_test_2(const Weighted_point<Bare_point, Weight> &p,
	     const Weighted_point<Bare_point, Weight> &q,
	     const Weighted_point<Bare_point, Weight> &t)
{
  typedef typename Kernel_traits<Bare_point>::Kernel::Rep_tag Tag;
  return power_test_2(p, q, t, Tag());
}

template < class Bare_point, class Weight >
inline
typename Kernel_traits<Bare_point>::Kernel::Oriented_side
power_test_2(const Weighted_point<Bare_point, Weight> &p,
	     const Weighted_point<Bare_point, Weight> &t)
{
  Comparison_result r = compare(p.weight(), t.weight());
  if(r == LARGER)    return ON_NEGATIVE_SIDE;
  else if (r == SMALLER) return ON_POSITIVE_SIDE;
  return ON_ORIENTED_BOUNDARY;
}


template < typename K >
class Power_test_2
{
public:
  typedef typename K::Weighted_point_2         Weighted_point_2;
  typedef typename K::Oriented_side            Oriented_side;

  typedef Oriented_side    result_type;

  Oriented_side operator() ( const Weighted_point_2 & p,
			     const Weighted_point_2 & q,
			     const Weighted_point_2 & r,
			     const Weighted_point_2 & s) const
    {
      //CGAL_triangulation_precondition( ! collinear(p, q, r) );
      return CGAL::power_test_2(p,q,r,s);
    }

  Oriented_side operator() ( const Weighted_point_2 & p,
			     const Weighted_point_2 & q,
			     const Weighted_point_2 & r) const
    {
      //CGAL_triangulation_precondition( collinear(p, q, r) );
      //CGAL_triangulation_precondition( p.point() != q.point() );
      return CGAL::power_test_2(p,q,r);
    }  

  Oriented_side operator() ( const Weighted_point_2 & p,
			     const Weighted_point_2 & r) const
    {
      //CGAL_triangulation_precondition( p.point() == r.point() );
      return CGAL::power_test_2(p,r);
    }
};

template < class R, class W = typename R::RT>
class Regular_triangulation_euclidean_traits_base_2
  : public R
{
public:
  typedef R                                     Kernel;
  typedef R                                     Rep;
  typedef W                                     Weight;
  typedef R                                     Traits;
  typedef typename Traits::Point_2              Bare_point;
  typedef typename Traits::Point_2              Point_2;
  typedef CGAL::Weighted_point<Bare_point, W>   Weighted_point_2;
  // This is required for the point() function of vertex base class
  // to be correctly return a weighted_point;
  // patch 27/11/00
  // 18/03/03 I put now the same typedef in Regulat_triangulation_2
  // for the need of hierarchy
  // don't know if this is definitive
  //typedef Weighted_point                        Point_2;

  typedef Regular_triangulation_euclidean_traits_base_2<R, W>   Self;

  typedef CGAL::Power_test_2<Self>              Power_test_2;
  typedef CGAL::Compare_power_distance_2<Self>  Compare_power_distance_2;

  // construction objects
  typedef CGAL::Construct_weighted_circumcenter_2<Self> 
                                            Construct_weighted_circumcenter_2;
  typedef CGAL::Construct_radical_axis_2<Self>  Construct_radical_axis_2;
  
  Power_test_2 
  power_test_2_object() const
    {  return Power_test_2();}

  Compare_power_distance_2
  compare_power_distance_2_object() const {
    return Compare_power_distance_2();
  }

  //constructions for dual:
  Construct_weighted_circumcenter_2
  construct_weighted_circumcenter_2_object() const
    {return Construct_weighted_circumcenter_2();}
  
  Construct_radical_axis_2
  construct_radical_axis_2_object() const
    {return Construct_radical_axis_2();}

};
 
// We use a base class here to have the specialization below to work.
// Otherwise there is a cycle in the derivation.
template < class R, class W = typename R::RT, bool Has_filtered_predicates = R::Has_filtered_predicates>
class Regular_triangulation_euclidean_traits_2
  : public Regular_triangulation_euclidean_traits_base_2<R, W>
{};

} //namespace CGAL

// Now specialize for Kernel with filtered predicates, to get
// the filtered traits automatically.
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Filtered_kernel.h>

namespace CGAL {

// This declaration is needed to break the cyclic dependency.
template < typename K >
class Regular_triangulation_filtered_traits_2;


template < typename K, class W >
class Regular_triangulation_euclidean_traits_2 < K, W, true >
  : public Regular_triangulation_filtered_traits_2 < K >
{
public:
  typedef K   Kernel;
};

} //namespace CGAL

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
