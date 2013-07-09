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
    typedef typename Arrangement_2::Edge_const_iterator   Edge_const_iterator;
    typedef typename Geometry_traits_2::Segment_2         Segment_2;
	typedef typename Geometry_traits_2::Point_2	          Point_2;

    // First make sure they have the same size
    if (arr1.number_of_vertices() != arr2.number_of_vertices()) {
        return false;
    }
    if (arr1.number_of_edges() != arr2.number_of_edges()) {
        return false;
    }

	Vertex_const_iterator vit_fst, vit_snd;
	vit_fst = arr1.vertices_begin();
    vit_snd = arr2.vertices_begin();

    Point_2 min_fst = vit_fst->point();
    Point_2 min_snd = vit_snd->point();
    vit_fst++;
    vit_snd++;

    for (vit_fst ; vit_fst != arr1.vertices_end(); ++vit_fst) {
        if (vit_fst->point().y() < min_fst.y()) {
            min_fst = vit_fst->point();
        }
        else if (vit_fst->point().y() == min_fst.y() &&
                vit_fst->point().x() < min_fst.x()) {
            min_fst = vit_fst->point();
        }
    }

    for (vit_snd ; vit_snd != arr2.vertices_end() ; ++vit_snd) {
        if (vit_snd->point().y() < min_snd.y()) {
            min_snd = vit_snd->point();
        }
        else if (vit_snd->point().y() == min_snd.y() &&
                vit_snd->point().x() < min_snd.y()) {
            min_snd = vit_snd->point();
        }
    }

    // Determine from which edge to start from 
    // (the one that contains the minimum vertex)
    Edge_const_iterator eit_fst, eit_snd, eit_fst_mid, eit_snd_mid;

    for (eit_fst = arr1.edges_begin(); eit_fst != arr1.edges_end(); ++eit_fst) {
        Segment_2 curr = eit_fst->curve();
        if (curr.source() == min_fst || curr.target() == min_fst) {
            break;
        }
    }

    for (eit_snd = arr2.edges_begin() ; eit_snd != arr2.edges_end(); ++eit_snd){
        Segment_2 curr = eit_snd->curve();
        if (curr.source() == min_snd || curr.target() == min_snd) {
            break;
        }
    }
    eit_fst_mid = eit_fst;
    eit_snd_mid = eit_snd;
    // Now we know where to start the edge comparisons from
    for (eit_fst, eit_snd ; eit_fst != arr1.edges_end(), 
         eit_snd != arr2.edges_end() ; ++eit_fst, ++eit_snd) {

        Segment_2 seg_fst = eit_fst->curve();
        Segment_2 seg_snd = eit_snd->curve();
        if (seg_fst != seg_snd) {
            return false;
        }
    }
    // Check first edges as well
    for (eit_fst = arr1.edges_begin(), eit_snd = arr2.edges_begin() ;
         eit_fst != eit_fst_mid, eit_snd != eit_snd_mid ; ++eit_fst, ++eit_snd){

        Segment_2 seg_fst = eit_fst->curve();
        Segment_2 seg_snd = eit_snd->curve();
        if (seg_fst != seg_snd) {
            return false;
        }
    }
    return true;
}
} // end namespace CGAL

#endif 