// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.5-I-99 $
// release_date  : $CGAL_Date: 2003/05/23 $
//
// file          : include/CGAL/Weighted_Minkowski_distance.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

// Note: Use p=0 to denote the weighted Linf-distance 
// For 0<p<1 Lp is not a metric

#ifndef CGAL_WEIGHTED_MINKOWSKI_DISTANCE_H
#define CGAL_WEIGHTED_MINKOWSKI_DISTANCE_H
#include <math.h>
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class GeomTraits>
  class Weighted_Minkowski_distance {

    public:

    typedef typename GeomTraits::Point Point;
    typedef typename GeomTraits::NT NT;
    typedef std::vector<NT> Weight_vector;

    private:

    typedef typename GeomTraits::Cartesian_const_iterator Coord_iterator;
    NT power; 

    Weight_vector the_weights;

    public:


    // default constructor
    Weighted_Minkowski_distance()
      : power(2) 
    {}

    Weighted_Minkowski_distance(const int d) 
      : power(2)
    {
      the_weights.reserve(d);
      for (unsigned int i = 0; i < d; ++i) the_weights[i]=NT(1);
    }

    //default copy constructor and destructor
    

    Weighted_Minkowski_distance (NT pow, int dim,
				 const Weight_vector& weights) 
      : power(pow)
    {
      assert(power >= NT(0));
      assert(dim==weights.size());
      for (unsigned int i = 0; i < weights.size(); ++i)
	assert(weights[i]>=NT(0));
      the_weights.resize(weights.size());
      the_weights = weights;
    }

    template <class InputIterator>
    Weighted_Minkowski_distance (NT pow, int dim,
				 InputIterator begin, InputIterator end) 
      : power(pow)
    {
      assert(power >= NT(0));
      the_weights.resize(dim);
      std::copy(begin, end, the_weights.begin());
      for (unsigned int i = 0; i < dim; ++i){
	the_weights[i] = *begin;
	++begin;
	assert(the_weights[i]>=NT(0));
      }
      assert(begin == end);
    }


    inline 
    NT 
    distance(const Point& q, const Point& p) 
    {
      NT distance = NT(0);
      typename GeomTraits::Construct_cartesian_const_iterator construct_it;
      Coord_iterator qit = construct_it(q),
	             qe = construct_it(q,1), 
	             pit = construct_it(p);
      if (power == NT(0)) {
	for (unsigned int i = 0; qit != qe; ++qit, ++i)
	  if (the_weights[i] * fabs((*qit) - (*pit)) > distance)
	    distance = the_weights[i] * fabs((*qit)-(*pit));
      }
      else
	for (unsigned int i = 0; qit != qe; ++qit, ++i)
	  distance += 
	    the_weights[i] * pow(fabs((*qit)-(*pit)),power);
      return distance;
    }
    

    inline 
    NT 
    min_distance_to_queryitem(const Point& q,
			      const Kd_tree_rectangle<GeomTraits>& r) const 
    {
      NT distance = NT(0);
      typename GeomTraits::Construct_cartesian_const_iterator construct_it;
      Coord_iterator qit = construct_it(q), qe = construct_it(q,1);
      if (power == NT(0))
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if (the_weights[i]*(r.min_coord(i) - 
				(*qit)) > distance)
	      distance = the_weights[i] * (r.min_coord(i)-
					   (*qit));
	    if (the_weights[i] * ((*qit) - r.max_coord(i)) > 
		distance)
	      distance = the_weights[i] * 
		((*qit)-r.max_coord(i));
	  }
	}
      else
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) < r.min_coord(i))
	      distance += the_weights[i] * 
		pow(r.min_coord(i)-(*qit),power);
	    if ((*qit) > r.max_coord(i))
	      distance += the_weights[i] * 
		pow((*qit)-r.max_coord(i),power);
	  }
	};
      return distance;
    }

    inline 
    NT
    max_distance_to_queryitem(const Point& q,
			      const Kd_tree_rectangle<GeomTraits>& r) const {
      NT distance=NT(0);
      typename GeomTraits::Construct_cartesian_const_iterator construct_it;
      Coord_iterator qit = construct_it(q), qe = construct_it(q,1);
      if (power == NT(0))
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) >= (r.min_coord(i) + 
			 r.max_coord(i))/NT(2.0))
	      if (the_weights[i] * ((*qit) - 
				    r.min_coord(i)) > distance)
		distance = the_weights[i] * 
		  ((*qit)-r.min_coord(i));
	      else
		if (the_weights[i] * 
		    (r.max_coord(i) - (*qit)) > distance)
		  distance = the_weights[i] * 
		    ( r.max_coord(i)-(*qit));
	  }
	}
      else
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) <= (r.min_coord(i)+r.max_coord(i))/NT(2.0))
	      distance += the_weights[i] * pow(r.max_coord(i)-(*qit),power);
	    else
	      distance += the_weights[i] * pow((*qit)-r.min_coord(i),power);
	  }
	};
      return distance;
    }
    
    inline 
    NT 
    new_distance(NT dist, NT old_off, NT new_off,
		 int cutting_dimension)  const 
    {
      NT new_dist;
      if (power == NT(0))
	{
	  if (the_weights[cutting_dimension]*fabs(new_off) 
	      > dist) 
	    new_dist= 
	      the_weights[cutting_dimension]*fabs(new_off);
	  else new_dist=dist;
	}
      else
	{
	  new_dist = dist + the_weights[cutting_dimension] * 
	    (pow(fabs(new_off),power)-pow(fabs(old_off),power));
	}
      return new_dist;
    }
    
    inline 
    NT 
    transformed_distance(NT d) const 
    {
      if (power <= NT(0)) return d;
      else return pow(d,power);
      
    }
    
    inline 
    NT 
    inverse_of_transformed_distance(NT d) const 
    {
      if (power <= NT(0)) return d;
      else return pow(d,1/power);
      
    }

  }; // class Weighted_Minkowski_distance

} // namespace CGAL
#endif // WEIGHTED_MINKOWSKI_DISTANCE_H
