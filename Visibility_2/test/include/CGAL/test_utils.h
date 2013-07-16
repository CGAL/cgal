/*
 * Author: Francisc Bungiu, Kan Huang 
 * E-mail: fbungiu@gmail.com, huangkandiy@gmail.com
 * Description: This file contains useful functions for testing the 
 * 				Visibility_2 package, such as comparing two Arrangements
 */

#ifndef CGAL_TEST_UTILS_H
#define CGAL_TEST_UTILS_H

#include <cassert>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <CGAL/Gmpq.h>


namespace CGAL {

/* 
 * Function to compare two arrangements; first determines lowest vertex
 * from each arrangements, then it walks the edges and compares them
 */
template <class _Arrangement_2> 
bool test_are_equal(const _Arrangement_2 &arr1, const _Arrangement_2 &arr2) {

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
        std::cout << seg_fst << std::endl;
        std::cout << seg_snd << std::endl;
        if ((seg_fst.source() != seg_snd.source() || seg_fst.target() != seg_snd.target()) &&
            (seg_fst.source() != seg_snd.target() || seg_fst.target() != seg_snd.source())) {
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



template<class Number_type>
Number_type string2num(const std::string& s) {
    int i;
    if (s.find("/") != std::string::npos) {
        i = s.find("/");
        std::string p = s.substr(0, i);
        std::string q = s.substr(i+1);
        std::stringstream convert(p);
        int n, d;
        convert >> n;
        std::stringstream convert2(q);
        convert2 >> d;        
        return Number_type(n)/Number_type(d);

    }
    else {
        std::stringstream convert(s);
        double n;
        convert >> n;
        return Number_type(n);
    }
}

template <class _Arrangement_2>
void create_arrangement_from_file(_Arrangement_2 &arr, std::ifstream& input) {
    typedef _Arrangement_2 								  Arrangement_2;
    typedef typename Arrangement_2::Geometry_traits_2	  Geometry_traits_2;
    typedef typename Geometry_traits_2::Segment_2         Segment_2;
    typedef typename Geometry_traits_2::Point_2	          Point_2;
    typedef typename Geometry_traits_2::FT                Number_type;
    if (input) {
        std::string line;
        while (std::getline(input, line)) {
            if (line[0] != '#' && line[0] != '/')
                break;
        }
        std::vector<Point_2> points;
        std::vector<Segment_2> segments;
        std::stringstream convert(line);
        int number_of_points;
        convert >> number_of_points;

        for (int i = 0; i != number_of_points; i++) {
            std::getline(input, line);
            std::string n1, n2;
            std::istringstream iss(line);
            iss>> n1 >> n2;
            points.push_back(Point_2(string2num<Number_type>(n1), string2num<Number_type>(n2)));
        }
        int number_of_edges;
        input >> number_of_edges;
        for (int i = 0; i != number_of_edges; i++) {
            unsigned i1,i2;
            input >> i1 >> i2;
            segments.push_back(Segment_2(points[i1], points[i2]));
        }
        CGAL::insert(arr, segments.begin(), segments.end());
    }
    else {
        std::cout<<"Can't open the file. Check the file name.";
    }
}

template <class _Arrangement_2>
void create_polygons_from_file(_Arrangement_2 &arr, std::ifstream& input) {
    typedef _Arrangement_2 								  Arrangement_2;
    typedef typename Arrangement_2::Geometry_traits_2	  Geometry_traits_2;
    typedef typename Geometry_traits_2::Segment_2         Segment_2;
    typedef typename Geometry_traits_2::Point_2	          Point_2;
    typedef typename Geometry_traits_2::FT                Number_type;
    if (input) {
        std::string line;
        while (std::getline(input, line)) {
            if (line[0] != '#' && line[0] != '/')
                break;
        }
        std::stringstream convert(line);
        int number_of_polygons;
        convert >> number_of_polygons;
        for (int i = 0; i != number_of_polygons; i++) {
            std::vector<Point_2> points;
            std::vector<Segment_2> segments;
            int number_of_vertex;
            input >> number_of_vertex;
            for (int j = 0; j != number_of_vertex-1; j++) {
                std::getline(input, line);
                std::string n1, n2;
                std::istringstream iss(line);
                iss >> n1 >> n2;
                points.push_back(Point_2(string2num<Number_type>(n1), string2num<Number_type>(n2)));
            }
            for (int j = 0; j != number_of_vertex-1; j++) {

                segments.push_back(Segment_2(points[j], points[j+1]));
            }
            segments.push_back(Segment_2(points.front(), points.back()));
            CGAL::insert(arr, segments.begin(), segments.end());
        }

    }
    else {
        std::cout<<"Can't open the file. Check the file name.";
    }

}

} // end namespace CGAL

#endif 
