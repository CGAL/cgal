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
// release       :
// release_date  :
//
// file          : include/CGAL/L1_distance_rectangle_point.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

// Note: this implementation is based on a weight vector,
// another traits class for weighted Minkowski distances based on a 
// weight matrix will be added to ASPAS in the future

// Note: Use p=0 to denote the weighted Linf-distance 
// For 0<p<1 Lp is not a metric

#ifndef CGAL_L1_DISTANCE_RECTANGLE_POINT_H
#define CGAL_L1_DISTANCE_RECTANGLE_POINT_H

#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class Query_item, class Item>
  class L1_distance_rectangle_point {

    public:


    typedef typename Item::R::FT NT;

    private:

    unsigned int The_dimension;

    public:

    
    // default constructor
    L1_distance_rectangle_point() : The_dimension(2) {}
    
    L1_distance_rectangle_point(const int d) : The_dimension(d) {}

    ~L1_distance_rectangle_point() {
    };

    inline NT distance(const Query_item& q, const Item& p) {
		NT distance = 0.0;
		for (unsigned int i = 0; i < The_dimension; ++i) {
			if (p[i]>q.upper(i)) distance += (p[i]-q.upper(i)); 
			if (p[i]<q.lower(i)) distance += (q.lower(i)-p[i]);	
		}
        	return distance;
    }


    inline NT min_distance_to_queryitem(const Query_item& q,
					      const Kd_tree_rectangle<NT>& r) {
		NT distance = 0.0;
		for (unsigned int i = 0; i < The_dimension; ++i)  {
			if (r.lower(i)>q.upper(i)) distance += (r.lower(i)-q.upper(i)); 
			if (r.upper(i)<q.lower(i)) distance += (q.lower(i)-r.upper(i));
	        }
		return distance;
	}

    inline NT max_distance_to_queryitem(const Query_item& q,
					      const Kd_tree_rectangle<NT>& r) {
		NT distance=0.0;
		for (unsigned int i = 0; i < The_dimension; ++i)
			if ( r.upper(i)-q.lower(i) > q.upper(i)-r.lower(i) )  
				distance += (r.upper(i)-q.lower(i));
			else 
				distance += (q.upper(i)-r.lower(i));
		return distance;
	}

	
  inline NT transformed_distance(NT d) {

		return d;

	}

  inline NT inverse_of_transformed_distance(NT d) {

		return d;

	}

  }; // class L1_distance_rectangle_point

} // namespace CGAL
#endif // L1_DISTANCE_RECTANGLE_POINT_H
