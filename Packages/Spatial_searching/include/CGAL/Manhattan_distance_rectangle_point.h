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
// file          : include/CGAL/Manhattan_distance_rectangle_point.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


#ifndef CGAL_MANHATTAN_DISTANCE_RECTANGLE_POINT_H
#define CGAL_MANHATTAN_DISTANCE_RECTANGLE_POINT_H

#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class QueryItem, class Point>
  class Manhattan_distance_rectangle_point {

    public:

    typedef typename Kernel_traits<Point>::Kernel::FT NT;

    private:

    unsigned int the_dimension;

    public:

    
    // default constructor
    Manhattan_distance_rectangle_point() {
		Point p;
		the_dimension=p.dimension();
		assert(the_dimension>0);
    }
    
    Manhattan_distance_rectangle_point(const int d) : the_dimension(d) {}

    //copy constructor
    Manhattan_distance_rectangle_point(const Manhattan_distance_rectangle_point& d) : 
    the_dimension(d.the_dimension) {}

    ~Manhattan_distance_rectangle_point() {
    };

    inline NT distance(const QueryItem& q, const Point& p) {
		NT distance = NT(0);
		for (unsigned int i = 0; i < the_dimension; ++i) {
			if (p[i]>q.max()[i]) distance += 
			(p[i]-q.max()[i]); 
			if (p[i]<q.min()[i]) distance += 
			(q.min()[i]-p[i]);	
		}
        	return distance;
    }


    inline NT min_distance_to_queryitem(const QueryItem& q,
					      const Kd_tree_rectangle<NT>& r) {
		NT distance = NT(0);
		for (unsigned int i = 0; i < the_dimension; ++i)  {
			if (r.min_coord(i)>q.max()[i]) distance += 
			(r.min_coord(i)-q.max()[i]); 
			if (r.max_coord(i)<q.min()[i]) distance += 
			(q.min()[i]-r.max_coord(i));
	        }
		return distance;
	}

    inline NT max_distance_to_queryitem(const QueryItem& q,
					      const Kd_tree_rectangle<NT>& r) {
		NT distance=NT(0);
		for (unsigned int i = 0; i < the_dimension; ++i)
			if ( r.max_coord(i)-q.min()[i] > 
			     q.max()[i]-r.min_coord(i) )  
				distance += (r.max_coord(i)-q.min()[i]);
			else 
				distance += (q.max()[i]-r.min_coord(i));
		return distance;
	}

	
  inline NT transformed_distance(NT d) {

		return d;

	}

  inline NT inverse_of_transformed_distance(NT d) {

		return d;

	}

  }; // class Manhattan_distance_rectangle_point

} // namespace CGAL
#endif // MANHATTAN_DISTANCE_RECTANGLE_POINT_H
