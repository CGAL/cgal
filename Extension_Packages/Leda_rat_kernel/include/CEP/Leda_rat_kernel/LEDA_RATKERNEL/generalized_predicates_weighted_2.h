#ifndef CEP_LEDA_RAT_GENERALIZED_PREDICATES_WEIGHTED_2_H
#define CEP_LEDA_RAT_GENERALIZED_PREDICATES_WEIGHTED_2_H

// LEDA rational kernel generalized predicates for weighted objects ...

#include <CGAL/Origin.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Quotient.h>
#include <CGAL/functional_base.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>

#include <LEDA/rat_point.h>

#if !defined(LEDA_NAMESPACE_NAME)
#define LEDA_NAMESPACE_NAME
#endif

CGAL_BEGIN_NAMESPACE

class Predicate_leda_rat_power_test_2
{
public:
  typedef Arity_tag<4> Arity;
  typedef CGAL::Oriented_side    result_type;

  typedef CGAL::Weighted_point<leda_rat_point, leda_rational>        Weighted_point;
  
  CGAL::Oriented_side operator() (Weighted_point p,
			     Weighted_point q,
			     Weighted_point r,
			     Weighted_point t) const
  {
      // perform CGAL power test ...
      return CGAL::power_testC2(p.xcoord(), p.ycoord(), p.weight(),
		                q.xcoord(), q.ycoord(), q.weight(),
		                r.xcoord(), r.ycoord(), r.weight(),
		                t.xcoord(), t.ycoord(), t.weight());
  }
};


class Predicate_leda_rat_power_test_degenerated_2
{
public:
  typedef Arity_tag<3> Arity;
  typedef CGAL::Oriented_side    result_type;

  typedef CGAL::Weighted_point<leda_rat_point, leda_rational>        Weighted_point;
  
  CGAL::Oriented_side operator()(Weighted_point p,
			     Weighted_point q,
			     Weighted_point t) const
    {
      // perform CGAL power test ...
      return power_testC2(p.xcoord(), p.ycoord(), p.weight(),
                          q.xcoord(), q.ycoord(), q.weight(),
                          t.xcoord(), t.ycoord(), t.weight());
    }
};


CGAL_END_NAMESPACE


#endif
