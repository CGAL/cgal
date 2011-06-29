// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)


#ifndef CGAL_FUZZY_SPHERE_H
#define CGAL_FUZZY_SPHERE_H
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class SearchTraits>
  class Fuzzy_sphere{
    SearchTraits traits;
    public:

    typedef typename SearchTraits::FT FT;
    typedef typename SearchTraits::Point_d Point_d;
    private:

    Point_d c;
    FT r;
    FT eps;

    public:


    	// default constructor
    	Fuzzy_sphere(const SearchTraits& traits_=SearchTraits()):traits(traits_) {}
		

	// constructor
	Fuzzy_sphere(const Point_d& center, FT radius, FT epsilon=FT(0),const SearchTraits& traits_=SearchTraits()) : 
	traits(traits_),c(center), r(radius), eps(epsilon) 
	{ 	// avoid problems if eps > r
		if (eps>r) eps=r; 
	}
        	
        bool contains(const Point_d& p) const {
		// test whether the squared distance 
		// between P and c 
		// is at most the squared_radius
		FT squared_radius = r*r;
		FT distance=FT(0);
		typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
                  traits.construct_cartesian_const_iterator_d_object();
                typename SearchTraits::Cartesian_const_iterator_d cit = construct_it(c),
		                                                  pit = construct_it(p),
                                                                  end = construct_it(c, 0);
		for (; cit != end
                       && (distance <= squared_radius); ++cit, ++pit) {
		  distance += 
			((*cit)-(*pit))*((*cit)-(*pit));
		}
		
		return (distance < squared_radius); 
        }

        
	bool inner_range_intersects(const Kd_tree_rectangle<FT>& rectangle) const {                          
                // test whether the interior of a sphere
		// with radius (r-eps) intersects r, i.e.
                // if the minimal distance of r to c is less than r-eps
		FT distance = FT(0);
		FT squared_radius = (r-eps)*(r-eps);
		typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
                  traits.construct_cartesian_const_iterator_d_object();
                typename SearchTraits::Cartesian_const_iterator_d cit = construct_it(c),
                                                                  end = construct_it(c, 0);
		for (int i = 0; cit != end && (distance < squared_radius); ++cit, ++i) {
			if ((*cit) < rectangle.min_coord(i))
				distance += 
				(rectangle.min_coord(i)-(*cit))*(rectangle.min_coord(i)-(*cit));
			if ((*cit) > rectangle.max_coord(i))
				distance += 
				((*cit)-rectangle.max_coord(i))*((*cit)-rectangle.max_coord(i));
		}
		
		return (distance < squared_radius);
	}


	bool outer_range_contains(const Kd_tree_rectangle<FT>& rectangle) const { 
        // test whether the interior of a sphere
	// with radius (r+eps) is contained by r, i.e.
        // if the minimal distance of the boundary of r 
        // to c is less than r+eps                         
	FT distance=FT(0);
	FT squared_radius = (r+eps)*(r+eps);	
	typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
          traits.construct_cartesian_const_iterator_d_object();
	typename SearchTraits::Cartesian_const_iterator_d cit = construct_it(c),
                                                          end = construct_it(c, 0);
        for (int i = 0; cit != end && (distance < squared_radius) ; ++cit,++i) {
		if ((*cit) <= (rectangle.min_coord(i)+rectangle.max_coord(i))/FT(2))
			distance += 
			(rectangle.max_coord(i)-(*cit))*(rectangle.max_coord(i)-(*cit));
		else
			distance += ((*cit)-rectangle.min_coord(i))*((*cit)-rectangle.min_coord(i));
		}
		
		return (distance < squared_radius);
	}
	
  

	~Fuzzy_sphere() {}

	

  }; // class Fuzzy_sphere

} // namespace CGAL
#endif // FUZZY_SPHERE_H
