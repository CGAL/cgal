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
// revision      : 2.4 
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
		NT distance = NT(0);
		for (unsigned int i = 0; i < The_dimension; ++i) {
			if (p[i]>q.max_coord(i)) distance += 
			(p[i]-q.max_coord(i)); 
			if (p[i]<q.min_coord(i)) distance += 
			(q.min_coord(i)-p[i]);	
		}
        	return distance;
    }


    inline NT min_distance_to_queryitem(const Query_item& q,
					      const Kd_tree_rectangle<NT>& r) {
		NT distance = NT(0);
		for (unsigned int i = 0; i < The_dimension; ++i)  {
			if (r.min_coord(i)>q.max_coord(i)) distance += 
			(r.min_coord(i)-q.max_coord(i)); 
			if (r.max_coord(i)<q.min_coord(i)) distance += 
			(q.min_coord(i)-r.max_coord(i));
	        }
		return distance;
	}

    inline NT max_distance_to_queryitem(const Query_item& q,
					      const Kd_tree_rectangle<NT>& r) {
		NT distance=NT(0);
		for (unsigned int i = 0; i < The_dimension; ++i)
			if ( r.max_coord(i)-q.min_coord(i) > 
			     q.max_coord(i)-r.min_coord(i) )  
				distance += (r.max_coord(i)-q.min_coord(i));
			else 
				distance += (q.max_coord(i)-r.min_coord(i));
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
