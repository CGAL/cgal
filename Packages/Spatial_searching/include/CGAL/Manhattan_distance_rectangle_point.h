// Copyright (c) 2002  Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Authors       : Hans Tangelder (<hanst@cs.uu.nl>)


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
