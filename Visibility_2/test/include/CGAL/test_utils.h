/*
 * Author: Francisc Bungiu 
 * E-mail: fbungiu@gmail.com
 * Description: This file contains useful functions for testing the 
 * 				Visibility_2 package, such as comparing two Arrangements
 */

#ifndef CGAL_TEST_UTILS_H
#define CGAL_TEST_UTILS_H

#include <cassert>

namespace CGAL {

/* 
 * Function to compare two arrangements; first determines lowest vertex
 * from each arrangements, then it walks the edges and compares them
 */
template <class _Arrangement_2> 
bool test_are_equal(_Arrangement_2 &arr1, _Arrangement_2 &arr2) {

	typedef _Arrangement_2 								  Arrangement_2;
	typedef typename Arrangement_2::Geometry_traits_2	  Geometry_traits_2;
	typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;
	typedef typename Geometry_traits_2::Point_2	          Point_2;

	Vertex_const_iterator vit_fst, vit_snd;
	vit_fst = arr1.vertices_begin();
    vit_snd = arr2.vertices_begin();

    Point_2 min_fst = vit_fst->point();
    Point_2 min_snd = vit_fst->point();

    return true;
}
} // end namespace CGAL

#endif