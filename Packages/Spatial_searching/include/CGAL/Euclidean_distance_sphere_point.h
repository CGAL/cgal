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



#ifndef CGAL_EUCLIDEAN_DISTANCE_SPHERE_POINT_H
#define CGAL_EUCLIDEAN_DISTANCE_SPHERE_POINT_H
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class QueryItem, class Point>
  class Euclidean_distance_sphere_point {

    public:

    typedef typename Kernel_traits<Point>::Kernel::FT NT;
    
    private:

    unsigned int the_dimension;

    public:


    	// default constructor
    	Euclidean_distance_sphere_point() {
		Point p;
		the_dimension=p.dimension();
		assert(the_dimension>0);
    	}


        Euclidean_distance_sphere_point(const int d) : the_dimension(d) {}

	~Euclidean_distance_sphere_point() {}

	inline NT distance(const QueryItem& q, const Point& p) const {
                Point c=q.center();
	        NT distance = NT(0);
		for (unsigned int i = 0; i < the_dimension; ++i)
			distance += (c[i]-p[i])*(c[i]-p[i]);
                distance += -q.squared_radius();
                if (distance<0) distance=NT(0);
        	return distance;
	}


	inline NT min_distance_to_queryitem(const QueryItem& q,
					    const Kd_tree_rectangle<NT>& r) const {
                Point c=q.center();
		NT distance = NT(0);
		for (unsigned int i = 0; i < the_dimension; ++i) {
			if (c[i] < r.min_coord(i))
				distance += 
				(r.min_coord(i)-c[i])*(r.min_coord(i)-c[i]);
			if (c[i] > r.max_coord(i))
				distance +=  
				(c[i]-r.max_coord(i))*(c[i]-r.max_coord(i));
			
		};
                distance += -q.squared_radius();
                if (distance<0) distance=NT(0);
		return distance;
	}

	inline NT max_distance_to_queryitem(const QueryItem& q,
					      const Kd_tree_rectangle<NT>& r) const {
                Point c=q.center();
		NT distance=NT(0);
		for (unsigned int i = 0; i < the_dimension; ++i) {
				if (c[i] <= (r.min_coord(i)+r.max_coord(i))/NT(2.0))
					distance += (r.max_coord(i)-c[i])*(r.max_coord(i)-c[i]);
				else
					distance += (c[i]-r.min_coord(i))*(c[i]-r.min_coord(i));
		};
		distance += -q.squared_radius();
                if (distance<0) distance=NT(0);
		return distance;
	}



  inline NT transformed_distance(NT d) const {
		return d*d;
	}

  inline NT inverse_of_transformed_distance(NT d) const {
		return CGAL::sqrt(d);
	}

  }; // class Euclidean_distance_sphere_point

} // namespace CGAL
#endif // EUCLIDEAN_DISTANCE_SPHERE_POINT_H
