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
// file          : include/CGAL/Fuzzy_sphere.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0 
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================



#ifndef CGAL_FUZZY_SPHERE_H
#define CGAL_FUZZY_SPHERE_H
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class SearchTraits>
  class Fuzzy_sphere{

    public:

    typedef typename SearchTraits::FT FT;
    
    typedef typename SearchTraits::Point_d Point_d;
    private:

    Point_d c;
    FT r;
    FT eps;
    unsigned int dim;

    public:

    	// default constructor
    	Fuzzy_sphere() {}
		

	// constructor
	Fuzzy_sphere(const Point_d& center, FT radius, FT epsilon=FT(0)) : 
	c(center), r(radius), eps(epsilon), dim(c.dimension()) 
	{ 	// avoid problems if eps > r
		if (eps>r) eps=r; 
	} 
        	
        bool contains(const Point_d& p) const {
		// test whether the squared distance 
		// between P and c 
		// is at most the squared_radius
		FT squared_radius = r*r;
		FT distance=FT(0);		  
		for (unsigned int i = 0; 
		(i < dim) && (distance <= squared_radius); ++i) {
			distance += 
			(c[i]-p[i])*(c[i]-p[i]);
		}
		return (distance < squared_radius); 
        }

        
	bool inner_range_intersects(const Kd_tree_rectangle<SearchTraits>& rectangle) const {                          
                // test whether the interior of a sphere
		// with radius (r-eps) intersects r, i.e.
                // if the minimal distance of r to c is less than r-eps
		FT distance = FT(0);
		FT squared_radius = (r-eps)*(r-eps);
		for (unsigned int i = 0; (i < dim) && (distance < squared_radius); ++i) {
			if (c[i] < rectangle.min_coord(i))
				distance += 
				(rectangle.min_coord(i)-c[i])*(rectangle.min_coord(i)-c[i]);
			if (c[i] > rectangle.max_coord(i))
				distance += 
				(c[i]-rectangle.max_coord(i))*(c[i]-rectangle.max_coord(i));
		}
		return (distance < squared_radius);
	}


	bool outer_range_is_contained_by(const Kd_tree_rectangle<SearchTraits>& rectangle) const { 
        // test whether the interior of a sphere
	// with radius (r+eps) is contained by r, i.e.
        // if the minimal distance of the boundary of r 
        // to c is less than r+eps                         
	FT distance=FT(0);
	FT squared_radius = (r+eps)*(r+eps);	
	for (unsigned int i = 0; (i < dim) && (distance < squared_radius) ; ++i) {
		if (c[i] <= (rectangle.min_coord(i)+rectangle.max_coord(i))/FT(2))
			distance += 
			(rectangle.max_coord(i)-c[i])*(rectangle.max_coord(i)-c[i]);
		else
			distance += (c[i]-rectangle.min_coord(i))*(c[i]-rectangle.min_coord(i));
		}
		return (distance < squared_radius);
	}

	~Fuzzy_sphere() {}

	

  }; // class Fuzzy_sphere

} // namespace CGAL
#endif // FUZZY_SPHERE_H
