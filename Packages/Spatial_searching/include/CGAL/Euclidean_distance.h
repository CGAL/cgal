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
// file          : include/CGAL/Euclidean_distance.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 2.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================



#ifndef CGAL_EUCLIDEAN_DISTANCE_H
#define CGAL_EUCLIDEAN_DISTANCE_H
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class SearchTraits>
  class Euclidean_distance {

    public:

    typedef typename SearchTraits::Point_d Point_d;
    typedef typename SearchTraits::FT    FT;
    typedef Point_d Query_item;
    
    public:

    	// default constructor
    	Euclidean_distance() {}


	inline FT transformed_distance(const Point_d& q, const Point_d& p) const {
	        FT distance = FT(0);
		typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
                typename SearchTraits::Cartesian_const_iterator_d qit = construct_it(q),
		  qe = construct_it(q,1), pit = construct_it(p);
		for(; qit != qe; qit++, pit++){
		  distance += ((*qit)-(*pit))*((*qit)-(*pit));
		}
        	return distance;
	}


	inline FT min_distance_to_rectangle(const Point_d& q,
					    const Kd_tree_rectangle<SearchTraits>& r) const {
		FT distance = FT(0);
		typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
                typename SearchTraits::Cartesian_const_iterator_d qit = construct_it(q),
		  qe = construct_it(q,1);
		for(unsigned int i = 0;qit != qe; i++, qit++){
		  if((*qit) < r.min_coord(i))
				distance += 
				(r.min_coord(i)-(*qit))*(r.min_coord(i)-(*qit));
			if ((*qit) > r.max_coord(i))
				distance +=  
				((*qit)-r.max_coord(i))*((*qit)-r.max_coord(i));
			
		}
		return distance;
	}

	inline FT max_distance_to_rectangle(const Point_d& q,
					     const Kd_tree_rectangle<SearchTraits>& r) const {
		FT distance=FT(0);
		typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
                typename SearchTraits::Cartesian_const_iterator_d qit = construct_it(q),
		  qe = construct_it(q,1);
		for(unsigned int i = 0;qit != qe; i++, qit++){
				if ((*qit) <= (r.min_coord(i)+r.max_coord(i))/FT(2.0))
					distance += (r.max_coord(i)-(*qit))*(r.max_coord(i)-(*qit));
				else
					distance += ((*qit)-r.min_coord(i))*((*qit)-r.min_coord(i));
		};
		return distance;
	}

	inline FT new_distance(FT dist, FT old_off, FT new_off,
			int cutting_dimension)  const {
		
		FT new_dist = dist + new_off*new_off - old_off*old_off;
                return new_dist;
	}

  inline FT transformed_distance(FT d) const {
		return d*d;
	}

  inline FT inverse_of_transformed_distance(FT d) const {
		return CGAL::sqrt(d);
	}

  }; // class Euclidean_distance

} // namespace CGAL
#endif // EUCLIDEAN_DISTANCE_H
