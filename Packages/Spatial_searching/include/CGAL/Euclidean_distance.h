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
#include <math.h>
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

  template <class Item>
  class Euclidean_distance {

    public:

    typedef typename Item::R::FT NT;
    
    private:

    unsigned int The_dimension;

    public:

    	// default constructor
    	Euclidean_distance() {
		Item P;
		The_dimension=P.dimension();
    	}

	~Euclidean_distance() {}

	inline NT distance(const Item& p1, const Item& p2) const {
	        NT distance = NT(0);
		for (unsigned int i = 0; i < The_dimension; ++i)
			distance += (p1[i]-p2[i])*(p1[i]-p2[i]);
        	return distance;
	}


	inline NT min_distance_to_queryitem(const Item& Point,
					    const Kd_tree_rectangle<NT>& r) const {
		NT distance = NT(0);
		for (unsigned int i = 0; i < The_dimension; ++i) {
			if (Point[i] < r.min_coord(i))
				distance += 
				(r.min_coord(i)-Point[i])*(r.min_coord(i)-Point[i]);
			if (Point[i] > r.max_coord(i))
				distance +=  
				(Point[i]-r.max_coord(i))*(Point[i]-r.max_coord(i));
			
		};
		return distance;
	}

	inline NT max_distance_to_queryitem(const Item& Point,
					      const Kd_tree_rectangle<NT>& r) const {
		NT distance=NT(0);
		for (unsigned int i = 0; i < The_dimension; ++i) {
				if (Point[i] <= (r.min_coord(i)+r.max_coord(i))/NT(2.0))
					distance += (r.max_coord(i)-Point[i])*(r.max_coord(i)-Point[i]);
				else
					distance += (Point[i]-r.min_coord(i))*(Point[i]-r.min_coord(i));
		};
		return distance;
	}

	inline NT new_distance(NT dist, NT old_off, NT new_off,
			int cutting_dimension)  const {
		
		NT new_dist = dist + new_off*new_off - old_off*old_off;
                return new_dist;
	}

  inline NT transformed_distance(NT d) const {
		return d*d;
	}

  inline NT inverse_of_transformed_distance(NT d) const {
		return sqrt(d);
	}

  }; // class Euclidean_distance

} // namespace CGAL
#endif // EUCLIDEAN_DISTANCE_H
