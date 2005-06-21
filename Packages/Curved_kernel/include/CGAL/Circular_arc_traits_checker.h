// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/Circular_arc_traits_checker.h

#ifndef CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_CHECKER_H
#define CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_CHECKER_H

// TODO :
// - Similar to Pm_checker_traits, but it doesn't require conversion
//   between the 2 types (especially endpoints).
// - Maybe this should move to TAU's CGAL packages (?)
// - It's not specific to circular arcs...

#include <CGAL/basic.h>
#include <cassert>
#include <utility>

namespace CGAL {

/// Traits checker : compares the output of all predicates by running two
/// traits in parallel.

template < typename P1, typename P2 >
struct Predicate_checker {

  P1 p1;
  P2 p2;

  Predicate_checker(const P1& _p1, const P2& _p2) : p1(_p1), p2(_p2) {}

  typedef typename P1::result_type  result_type;

  template < typename A1 >
  result_type
  operator()(const A1 & a1) const
    {
      std::cout << " call " << __PRETTY_FUNCTION__ << std::endl
		<< " [_1 ] " << a1.first << std::endl
		<< " [_2 ] " << a1.second << std::endl;

      result_type r1 = p1(a1.first);
      result_type r2 = p2(a1.second);
      
      if ( r1 != r2) {
	std::cerr << " Error in " << __PRETTY_FUNCTION__ << std::endl;
	std::cerr << " P1 returns : " << r1 << "  but P2 returns : "
		  << r2 << std::endl;
	std::cerr << " Arguments for P1 are : "
		  << a1.first << std::endl;
	std::cerr << " Arguments for P2 are : "
		  << a1.second << std::endl;
	abort(); // message... a1, a2
      }

      return r1;
    }

  template < typename A1, typename A2 >
  result_type
  operator()(const A1 & a1, const A2 & a2 ) const
    {
      std::cout << " call " << __PRETTY_FUNCTION__ << std::endl
		<< " [_1 ] " << a1.first << " et " << a2.first
		<< std::endl
		<< " [_2 ] " << a1.second << " et " << a2.second 
		<< std::endl;
      
      result_type r1 = p1(a1.first, a2.first);
      result_type r2 = p2(a1.second, a2.second);
      
      if ( r1 != r2) {
	std::cerr << " Error in " << __PRETTY_FUNCTION__ << std::endl;
	std::cerr << " P1 returns : " << r1 << "  but P2 returns : "
		  << r2 << std::endl;
	std::cerr << " Arguments for P1 are : "
		  << a1.first << std::endl
		  << a2.first << std::endl;
	std::cerr << " Arguments for P2 are : "
		  << a1.second << std::endl
		  << a2.second << std::endl;
	abort(); // message... a1, a2
      }

      return r1;
    }

  template < typename A1, typename A2, typename A3 >
  result_type
  operator()(const A1 & a1, const A2 & a2, const A3 & a3 ) const
    {
      std::cout << " call " << __PRETTY_FUNCTION__ << std::endl
		<< " [_1 ] " << a1.first << " et " << a2.first
		<< " et " << a3.first
		<< std::endl
		<< " [_2 ] " << a1.second << " et " << a2.second 
		<< " et " << a3.second 
		<< std::endl;
     
      result_type r1 = p1(a1.first, a2.first, a3.first);
      result_type r2 = p2(a1.second, a2.second, a3.second);
      
      if ( r1 != r2) {
	std::cerr << " Error in " << __PRETTY_FUNCTION__ << std::endl;
	std::cerr << " P1 returns : " << r1 << "  but P2 returns : "
		  << r2 << std::endl;
	std::cerr << " Arguments for P1 are : "
		  << a1.first << std::endl
		  << a2.first << std::endl
		  << a3.first << std::endl;
	std::cerr << " Arguments for P2 are : "
		  << a1.second << std::endl
		  << a2.second << std::endl
		  << a3.second << std::endl;
	abort(); // message... a1, a2
      }

      return r1;
    }
};


template < typename Traits >
struct cmp_less
{
  typedef typename Traits::Point_2             Point_2;
  typedef typename Traits::X_monotone_curve_2  X_monotone_curve_2;

  bool operator()(const Point_2 &p, const Point_2 &q) const
    {
      return Traits().compare_xy_2_object()(p, q) == CGAL::SMALLER;
    }

  bool operator()(const std::pair<Point_2, unsigned> &p,
		  const std::pair<Point_2, unsigned> &q) const
    {
      return this->operator()(p.first, q.first) ||
	(!this->operator()(q.first, p.first) && p.second < q.second);
    }

  bool operator()(const X_monotone_curve_2 &c,
		  const X_monotone_curve_2 &d) const
    {
      const Point_2& cp = Traits().construct_min_vertex_2_object()(c);
      const Point_2& dp = Traits().construct_min_vertex_2_object()(d);

      return this->operator()(cp, dp) ||
	(!this->operator()(dp, cp) &&
	 Traits().compare_y_at_x_right_2_object()(c, d, cp) == CGAL::SMALLER);
    }
};


template < typename Traits_1, typename Traits_2 >
class Circular_arc_traits_checker {

  Traits_1 t1;
  Traits_2 t2;

public:

  typedef Circular_arc_traits_checker<Traits_1, Traits_2 > Self;

  typedef Traits_1                                     First_traits;
  typedef Traits_2                                     Second_traits;

  typedef std::pair<typename Traits_1::X_monotone_curve_2,
                    typename Traits_2::X_monotone_curve_2> X_monotone_curve_2;
  typedef std::pair<typename Traits_1::Curve_2,
                    typename Traits_2::Curve_2>        Curve_2;
  typedef std::pair<typename Traits_1::Point_2,
                    typename Traits_2::Point_2>        Point_2;

  //typedef Curve                                        Curve_2;
  //typedef Curve                                        X_curve;
  //typedef Curve                                        X_monotone_curve_2;
  //typedef Point_2                                      Point;

  typedef typename Traits_1::Has_left_category         Has_left_category;
  typedef typename Traits_1::Has_merge_category        Has_merge_category;

  Circular_arc_traits_checker()
  {
    // "Check" that the categories are the same...
    //typename Traits_1::Has_left_category cat1;
    //typename Traits_2::Has_left_category cat2 = cat1;
  }

  typedef Predicate_checker<typename Traits_1::Compare_x_2,
                            typename Traits_2::Compare_x_2>  Compare_x_2;

  typedef Predicate_checker<typename Traits_1::Compare_xy_2,
                            typename Traits_2::Compare_xy_2>  Compare_xy_2;

  typedef Predicate_checker<typename Traits_1::Compare_y_at_x_2,
                            typename Traits_2::Compare_y_at_x_2>  Compare_y_at_x_2;

  typedef Predicate_checker<typename Traits_1::Compare_y_at_x_right_2,
                            typename Traits_2::Compare_y_at_x_right_2>  Compare_y_at_x_right_2;

  typedef Predicate_checker<typename Traits_1::Equal_2,
                            typename Traits_2::Equal_2>  Equal_2;

  typedef Predicate_checker<typename Traits_1::Is_vertical_2,
                            typename Traits_2::Is_vertical_2>  Is_vertical_2;

  Compare_x_2
  compare_x_2_object() const
    { return Compare_x_2(t1.compare_x_2_object(), t2.compare_x_2_object()); }

  Compare_xy_2
  compare_xy_2_object() const
    { return Compare_xy_2(t1.compare_xy_2_object(), t2.compare_xy_2_object()); }

  Compare_y_at_x_2
  compare_y_at_x_2_object() const
    { return Compare_y_at_x_2(t1.compare_y_at_x_2_object(),
			      t2.compare_y_at_x_2_object()); }

  Compare_y_at_x_right_2
  compare_y_at_x_right_2_object() const
    { return
	Compare_y_at_x_right_2(t1.compare_y_at_x_right_2_object(),
			       t2.compare_y_at_x_right_2_object()); }

  Equal_2
  equal_2_object() const
    { return Equal_2(t1.equal_2_object(), t2.equal_2_object()); }

  Is_vertical_2
  is_vertical_2_object() const
    { return Is_vertical_2(t1.is_vertical_2_object(),
			   t2.is_vertical_2_object()); }

  struct Construct_min_vertex_2
  {
    const Self *t;

    Construct_min_vertex_2(const Self*s) : t(s) {}

    typedef Point_2 result_type;

    Point_2 operator()(const X_monotone_curve_2 & cv) const
    {
      std::cout << " [Construct_min_vertex_2] "
		<< std::endl << "[Point _1]  "
		<< t->t1.construct_min_vertex_2_object()(cv.first)
		<< std::endl << "[Point  _2] "
		<< t->t2.construct_min_vertex_2_object()(cv.second)
		<< std::endl;
      return Point_2(t->t1.construct_min_vertex_2_object()(cv.first),
		     t->t2.construct_min_vertex_2_object()(cv.second));
    }
  };

  Construct_min_vertex_2
  construct_min_vertex_2_object() const
    { return Construct_min_vertex_2(this); }

  struct Construct_max_vertex_2
  {
    const Self *t;

    Construct_max_vertex_2(const Self*s) : t(s) {}

    typedef Point_2 result_type;

    Point_2 operator()(const X_monotone_curve_2 & cv) const
    {
      std::cout << " [Construct_max_vertex_2] "
		<< std::endl << "[Point _1]  "
		<< t->t1.construct_max_vertex_2_object()(cv.first)
		<< std::endl << "[Point  _2] "
		<< t->t2.construct_max_vertex_2_object()(cv.second)
		<< std::endl;
	
      return Point_2(t->t1.construct_max_vertex_2_object()(cv.first),
		     t->t2.construct_max_vertex_2_object()(cv.second));
    }
  };

  Construct_max_vertex_2
  construct_max_vertex_2_object() const
    { return Construct_max_vertex_2(this); }

  struct Split_2 {

    const Self *t;

    Split_2(const Self*s) : t(s) {}

    typedef void result_type;

    result_type
    operator()(const X_monotone_curve_2 &a, const Point_2 &p,
               X_monotone_curve_2 &ca1, X_monotone_curve_2 &ca2) const
    {
      using std::swap;
      t->t1.split_2_object()(a.first, p.first, ca1.first, ca2.first);
      t->t2.split_2_object()(a.second, p.second, ca1.second, ca2.second);
      if (cmp_less<Traits_1>()(ca2.first, ca1.first))
	swap(ca1.first, ca2.first);
      if (cmp_less<Traits_2>()(ca2.second, ca1.second))
	swap(ca1.second, ca2.second);

      std::cout << " [Split_2]"
		<< std::endl << "[Curve #1 _1] " << ca1.first;
      std::cout << CGALi::print(std::cout, ca1.first) << std::endl;
      std::cout << "[Curve #2 _1] " << ca2.first;
      std::cout << CGALi::print(std::cout, ca2.first) << std::endl;
      std::cout << "[Curve #1 _2] " << ca1.second << std::endl;
      std::cout << "[Curve #2 _2] " << ca2.second << std::endl;
    }
  };

  Split_2 split_2_object() const
  { return Split_2(this); }

  struct Make_x_monotone_2 {

    const Self *t;

    Make_x_monotone_2(const Self*s) : t(s) {}

    template < class OutputIterator >
    OutputIterator
    operator()(const Curve_2 &a, OutputIterator res) const
      {
	std::vector<typename Traits_1::X_monotone_curve_2> v1;
	std::vector<typename Traits_2::X_monotone_curve_2> v2;
	
	t->t1.make_x_monotone_2_object()(a.first, std::back_inserter(v1));
	t->t2.make_x_monotone_2_object()(a.second, std::back_inserter(v2));

	assert(v1.size() == v2.size());

	std::sort(v1.begin(), v1.end(), cmp_less<Traits_1>());
	std::sort(v2.begin(), v2.end(), cmp_less<Traits_2>());

	for(unsigned int i=0; i<v1.size(); ++i)
	  *res++ = std::make_pair(v1[i], v2[i]);

	for(unsigned int i=0; i<v1.size(); ++i) {
	  std::cout << " [Monotone_2] "
		    << std::endl << "[Curve  _1] : "
		    << v1[i] << std::endl;
	  std::cout << CGALi::print(std::cout, v1[i]) << std::endl;
	  std::cout << "[Curve  _2] : "
		    << v2[i] << std::endl;
	}

	return res;
      }
  };

  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(this); }

  struct Intersect_2 {

    const Self *t;

    Intersect_2(const Self*s) : t(s) {}

    template < class OutputIterator >
    OutputIterator
    operator()(const X_monotone_curve_2 &a,
	       const X_monotone_curve_2 &b,
	       OutputIterator res) const
      {
	std::vector<CGAL::Object> v1, v2;
	std::vector<typename Traits_1::X_monotone_curve_2>  vc1;
	std::vector<typename Traits_2::X_monotone_curve_2>  vc2;
	std::vector<std::pair<typename Traits_1::Point_2, unsigned> >  vp1;
	std::vector<std::pair<typename Traits_2::Point_2, unsigned> >  vp2;

	t->t1.intersect_2_object()(a.first, b.first, std::back_inserter(v1));
	t->t2.intersect_2_object()(a.second, b.second, std::back_inserter(v2));

	for(unsigned i=0; i<v1.size(); ++i) {
	  typename Traits_1::X_monotone_curve_2 c1;
	  std::pair<typename Traits_1::Point_2, unsigned int> pp1;
	  
	  if (assign(c1, v1[i])) {
	    vc1.push_back(c1);
	  }
	  else if (assign(pp1, v1[i])) {
	    vp1.push_back(pp1);
	  }
	  else
	    abort();
	}

	for(unsigned i=0; i<v2.size(); ++i) {
	  typename Traits_2::X_monotone_curve_2 c2;
	  std::pair<typename Traits_2::Point_2, unsigned int> pp2;
	  
	  if (assign(c2, v2[i])) {
	    vc2.push_back(c2);
	  }
	  else if (assign(pp2, v2[i])) {
	    vp2.push_back(pp2);
	  }
	  else
	    abort();
	}


	std::sort(vp1.begin(), vp1.end(), cmp_less<Traits_1>());
	std::sort(vp2.begin(), vp2.end(), cmp_less<Traits_2>());
	std::sort(vc1.begin(), vc1.end(), cmp_less<Traits_1>());
	std::sort(vc2.begin(), vc2.end(), cmp_less<Traits_2>());

	for(unsigned int i=0; i<vp1.size(); ++i) {
	  std::cout << " [Intersect_2] "
		    << std::endl << "[pair<Pt, int> _1] : "
		    << vp1[i].first << "  , multiplicite = "
		    << vp1[i].second << std::endl;
	  std::cout << "[pair<Pt, int> _2] : "
		    << vp2[i].first << "  , multiplicite = "
		    << vp2[i].second << std::endl;
	  *res++ = CGAL::make_object(std::pair<Point_2, unsigned>
					 (Point_2(vp1[i].first, vp2[i].first),
					  vp1[i].second));
	}

	for(unsigned int i=0; i<vc1.size(); ++i) {
	  std::cout << " [Intersect_2]"
		    << std::endl << "[Curve  _1] : "
		    << vc1[i] << std::endl;
	  std::cout << "[Curve  _2] : "
		    << vc2[i] << std::endl;
	  *res++ = CGAL::make_object(X_monotone_curve_2(vc1[i], vc2[i]));
	}

	assert(vp1.size() == vp2.size());
	assert(vc1.size() == vc2.size());
	assert(v1.size()  == v2.size());

	return res;
      }
  };

  Intersect_2 intersect_2_object() const
  { return Intersect_2(this); }


};

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_CHECKER_H
