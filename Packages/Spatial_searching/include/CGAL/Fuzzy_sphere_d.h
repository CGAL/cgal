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



#ifndef CGAL_FUZZY_SPHERE_D_H
#define CGAL_FUZZY_SPHERE_D_H
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class Point>
  class Fuzzy_sphere_d{

    public:

    typedef typename Kernel_traits<Point>::Kernel::FT NT;
    
    private:

    Point c;
    NT r;
    NT eps;
    unsigned int dim;

    public:

    	// default constructor
    	Fuzzy_sphere_d() {}
		

	// constructor
	Fuzzy_sphere_d(Point center, NT radius, NT epsilon=NT(0)) : 
	c(center), r(radius), eps(epsilon), dim(c.dimension()) 
	{ 	// avoid problems if eps > r
		if (eps>r) eps=r; 
	} 
        	
        bool contains(const Point& p) const {
		// test whether the squared distance 
		// between P and c 
		// is at most the squared_radius
		NT squared_radius = r*r;
		NT distance=NT(0);		  
		for (unsigned int i = 0; 
		(i < dim) && (distance <= squared_radius); ++i) {
			distance += 
			(c[i]-p[i])*(c[i]-p[i]);
		}
		return (distance < squared_radius); 
        }

        
	bool inner_range_intersects(const Kd_tree_rectangle<NT>* rectangle) const {                          
                // test whether the interior of a sphere
		// with radius (r-eps) intersects r, i.e.
                // if the minimal distance of r to c is less than r-eps
		NT distance = NT(0);
		NT squared_radius = (r-eps)*(r-eps);
		for (unsigned int i = 0; (i < dim) && (distance < squared_radius); ++i) {
			if (c[i] < rectangle->min_coord(i))
				distance += 
				(rectangle->min_coord(i)-c[i])*(rectangle->min_coord(i)-c[i]);
			if (c[i] > rectangle->max_coord(i))
				distance += 
				(c[i]-rectangle->max_coord(i))*(c[i]-rectangle->max_coord(i));
		}
		return (distance < squared_radius);
	}


	bool outer_range_is_contained_by(const Kd_tree_rectangle<NT>* rectangle) const { 
        // test whether the interior of a sphere
	// with radius (r+eps) is contained by r, i.e.
        // if the minimal distance of the boundary of r 
        // to c is less than r+eps                         
	NT distance=NT(0);
	NT squared_radius = (r+eps)*(r+eps);	
	for (unsigned int i = 0; (i < dim) && (distance < squared_radius) ; ++i) {
		if (c[i] <= (rectangle->min_coord(i)+rectangle->max_coord(i))/NT(2))
			distance += 
			(rectangle->max_coord(i)-c[i])*(rectangle->max_coord(i)-c[i]);
		else
			distance += (c[i]-rectangle->min_coord(i))*(c[i]-rectangle->min_coord(i));
		}
		return (distance < squared_radius);
	}

	~Fuzzy_sphere_d() {}

	

  }; // class Fuzzy_sphere_d

} // namespace CGAL
#endif // FUZZY_SPHERE_D_H
