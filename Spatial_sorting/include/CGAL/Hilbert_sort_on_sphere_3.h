// Copyright (c) 2013  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Olivier Devillers
//                 Pedro Machado Manhaes de Castro

#ifndef CGAL_HILBERT_SORT_ON_SPHERE_3_H
#define CGAL_HILBERT_SORT_ON_SPHERE_3_H

#include <CGAL/Hilbert_sort_2.h>
#include <CGAL/internal/Transform_coordinates_traits_3.h>
#include <algorithm>
#include <vector>

namespace CGAL {

template <class K,  class Hilbert_policy >
class Hilbert_sort_on_sphere_3 {
	
	static const double _sqrt_of_one_over_three = 0.57735026919;
	
	// Face 1, x > sqrt(1/3)
	// Face 2, y > sqrt(1/3)
	// Face 3, x < -sqrt(1/3)
	// Face 4, z > sqrt(1/3)
	// Face 5, y < -sqrt(1/3)
	// Face 6, z < -sqrt(1/3)
	
	typedef internal::Transform_coordinates_traits_3<K,0,1,1,0> Face_1_traits_3;  // +y +z
	typedef internal::Transform_coordinates_traits_3<K,-1,0,1,0> Face_2_traits_3; // -x +z
	typedef internal::Transform_coordinates_traits_3<K,0,1,-1,1> Face_3_traits_3; // +z -y
	typedef internal::Transform_coordinates_traits_3<K,-1,1,0,1> Face_4_traits_3;  // -y +x
	typedef internal::Transform_coordinates_traits_3<K,-1,0,1,1> Face_5_traits_3;  // -z +x
	typedef internal::Transform_coordinates_traits_3<K,1,1,0,0> Face_6_traits_3;  // +x +y
	
	Hilbert_sort_2<Face_1_traits_3, Hilbert_policy > _hs_1_object;
	Hilbert_sort_2<Face_2_traits_3, Hilbert_policy > _hs_2_object;
	Hilbert_sort_2<Face_3_traits_3, Hilbert_policy > _hs_3_object;
	Hilbert_sort_2<Face_4_traits_3, Hilbert_policy > _hs_4_object;
	Hilbert_sort_2<Face_5_traits_3, Hilbert_policy > _hs_5_object;
	Hilbert_sort_2<Face_6_traits_3, Hilbert_policy > _hs_6_object;
	K _k;
	
public:
	Hilbert_sort_on_sphere_3 (const K &k=K() , std::ptrdiff_t limit=1 ) 
	: _k(k), 
	  _hs_1_object(Face_1_traits_3(),limit), 
	  _hs_2_object(Face_2_traits_3(),limit),
	  _hs_3_object(Face_3_traits_3(),limit),
	  _hs_4_object(Face_4_traits_3(),limit),
	  _hs_5_object(Face_5_traits_3(),limit),
	  _hs_6_object(Face_6_traits_3(),limit) {}
	
	
	template <class RandomAccessIterator>
	void operator()(RandomAccessIterator begin, RandomAccessIterator end) const {
		std::vector<typename K::Point_3> vec[6];
		static const typename K::FT sqrt_of_one_over_three = typename K::FT(_sqrt_of_one_over_three);
		for(RandomAccessIterator i = begin; i != end; ++i) {
			const typename K::Point_3 &p = *i;
			const typename K::FT x = _k.compute_x_3_object()(p);
			const typename K::FT y = _k.compute_y_3_object()(p);
			const typename K::FT z = _k.compute_z_3_object()(p);
			if(x > sqrt_of_one_over_three) vec[0].push_back(p);             // Face 1, x > sqrt(1/3)
			else if(y > sqrt_of_one_over_three) vec[1].push_back(p);        // Face 2, y > sqrt(1/3)
			else if(x < -sqrt_of_one_over_three) vec[2].push_back(p);       // Face 3, x < -sqrt(1/3)
			else if(z > sqrt_of_one_over_three) vec[3].push_back(p);        // Face 4, z > sqrt(1/3)
			else if(y < -sqrt_of_one_over_three) vec[4].push_back(p);       // Face 5, y < -sqrt(1/3)
			else vec[5].push_back(p);                              			// Face 6, z < -sqrt(1/3)
		}
		_hs_1_object(vec[0].begin(), vec[0].end());
		_hs_2_object(vec[1].begin(), vec[1].end());
		_hs_3_object(vec[2].begin(), vec[2].end());
		_hs_4_object(vec[3].begin(), vec[3].end());
		_hs_5_object(vec[4].begin(), vec[4].end());
		_hs_6_object(vec[5].begin(), vec[5].end());
		
		// this is the order that set of points in a face should appear 
		// after sorting points wrt each face
		for(int i=0; i<6; i++)
			for(int j=0; j<vec[i].size(); j++) 
				*begin++ = vec[i][j];
	}
};


} // namespace CGAL

#endif//CGAL_HILBERT_SORT_ON_SPHERE_3_H
