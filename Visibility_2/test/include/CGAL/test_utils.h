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
#include <set>

namespace CGAL {

template <class Arrangement_2> 
typename Arrangement_2::Halfedge_handle get_initial_halfedge(const Arrangement_2 &arr) {
  
  typedef typename Arrangement_2::Vertex Vertex;
  typedef typename Arrangement_2::Vertex_handle Vertex_handle;
  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Arrangement_2::Halfedge Halfedge; 
  typedef typename Arrangement_2::Halfedge_handle Halfedge_handle; 
  
  // find the min vertex 
  Vertex v = *arr.vertices_begin(); 
  for(Vertex_const_iterator vit = arr.vertices_begin(); vit !=  arr.vertices_end(); vit++){
    if(arr.traits()->compare_xy_2_object()((*vit).point(),v.point()) == CGAL::SMALLER){
      v = *vit;
    }
  }
  
  // take the edge with the smallest source 
  Halfedge_handle he_final = v.incident_halfedges(); 
  Halfedge_handle he1 = v.incident_halfedges(); 
  he1=he1->next()->twin();
  
  while(he1 != v.incident_halfedges()){
    if(arr.traits()->compare_xy_2_object()(
           he1->source()->point(),
           he_final->source()->point()) == CGAL::SMALLER){
      he_final = he1;      
    }
    he1=he1->next()->twin();
  }
  
  // as this may be on a needle, continue until faces on both sides differ 
  while(he_final->face() == he_final->twin()->face())
    he_final = he_final->next(); 
  
  return he_final;  
  
}

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
  typedef typename Arrangement_2::Halfedge              Halfedge;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Geometry_traits_2::Segment_2         Segment_2;
  typedef typename Geometry_traits_2::Point_2	          Point_2;


  // First make sure they have the same size
  if (arr1.number_of_vertices() != arr2.number_of_vertices()) {
    return false;
  }
  if (arr1.number_of_edges() != arr2.number_of_edges()) {
    return false;
  }
  if (arr1.number_of_isolated_vertices() != arr2.number_of_isolated_vertices()) {
    return false;
  }
  if (arr1.number_of_faces() != arr2.number_of_faces()) {
    return false;
  }
  
  // currently checking for closed for visibility region 
  assert(arr1.number_of_faces() == 2);
  assert(arr2.number_of_faces() == 2);


  // get unique halfedge 
  Halfedge_handle he_start_1 = get_initial_halfedge(arr1);
  Halfedge_handle he_start_2 = get_initial_halfedge(arr2);
  
  // run on first loop and compare sources 
  assert(arr1.traits()->compare_xy_2_object()(
             he_start_1->source()->point(),
             he_start_2->source()->point()) == CGAL::EQUAL);

  Halfedge_handle he_run_1 = he_start_1->next(); 
  Halfedge_handle he_run_2 = he_start_2->next(); 

  while(he_run_1 != he_start_1){
    
    assert(arr1.traits()->compare_xy_2_object()(
               he_run_1->source()->point(),
               he_run_2->source()->point()) == CGAL::EQUAL);
    
    he_run_1 = he_run_1->next(); 
    he_run_2 = he_run_2->next(); 
  }
  assert(he_run_2 == he_start_2); 


  // run on second loop and compare sources. 
  he_start_1 =  he_start_1->twin(); 
  he_start_2 =  he_start_2->twin(); 
  he_run_1 = he_start_1->next(); 
  he_run_2 = he_start_2->next(); 
  while(he_run_1 != he_start_1){
    
    assert(arr1.traits()->compare_xy_2_object()(
               he_run_1->source()->point(),
               he_run_2->source()->point()) == CGAL::EQUAL);
    
    he_run_1 = he_run_1->next(); 
    he_run_2 = he_run_2->next(); 
  }
  assert(he_run_2 == he_start_2); 

  
  return true; 
  
//   Edge_const_iterator eit_fst, eit_snd;
//   std::vector<Halfedge> halfedges_fst, halfedges_snd;

//   for (eit_fst = arr1.edges_begin(), eit_snd = arr2.edges_begin() ; 
//        eit_fst != arr1.edges_end(), eit_snd != arr2.edges_end() ; 
//        ++eit_fst, ++eit_snd) {

//     halfedges_fst.push_back(*eit_fst);
//     halfedges_snd.push_back(*eit_snd);
//   }

//   // Compare the two vectors
//   for (unsigned int i = 0 ; i < halfedges_fst.size() ; i++) {
//     Halfedge he_curr = halfedges_fst[i];
//     bool found = false;
//     for (unsigned int j = 0 ; j < halfedges_snd.size() ; j++) {
//       if (he_curr.source()->point() == halfedges_snd[j].source()->point() &&
//           he_curr.target()->point() == halfedges_snd[j].target()->point()) {
//         found = true;
//         break;
//       }
//     }
//     if (found == false) {
//       return false;
//     }
//   }

//   for (unsigned int i = 0 ; i < halfedges_snd.size() ; i++) {
//     Halfedge he_curr = halfedges_snd[i];
//     bool found = false;
//     for (unsigned int j = 0 ; j < halfedges_fst.size() ; j++) {
//       if (he_curr.source()->point() == halfedges_fst[j].source()->point() &&
//           he_curr.target()->point() == halfedges_fst[j].target()->point()) {
//         found = true;
//         break;
//       }
//     }
//     if (found == false) {
//       return false;
//     }
//   }
//   return true;
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

template<class Number_type>
std::string num2string(Number_type& n) {
  std::stringstream ss;
  ss<<n;
  return ss.str();
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

template<typename Arrangement_2>
bool compare_arr_by_edges(const Arrangement_2& arr1, const Arrangement_2& arr2) {
  std::set<std::string> s1;
  typedef typename Arrangement_2::Edge_const_iterator Edge_const_iterator;
  for (Edge_const_iterator eit = arr1.edges_begin(); eit != arr1.edges_end(); ++eit) {
    s1.insert(edge2string(eit->target()->point(), eit->source()->point()));
  }
  std::set<std::string> s2;
  for (Edge_const_iterator eit = arr2.edges_begin(); eit != arr2.edges_end(); ++eit) {
    s2.insert(edge2string(eit->target()->point(), eit->source()->point()));
  }
  return s1==s2;
}

template<typename Point_2>
bool is_ahead(Point_2& p1, Point_2& p2) {
  if (p1.x() > p2.x()) {
    return true;
  }
  if (p1.x() == p2.x() && p1.y() > p2.y()) {
    return true;
  }
  return false;
}


template<typename Point_2>
std::string edge2string(const Point_2& p1, const Point_2& p2) {
  Point_2 q1, q2;
  if (is_ahead(p1, p2)) {
    q1 = p1;
    q2 = p2;
  }
  else {
    q1 = p2;
    q2 = p1;
  }
  return num2string(q1.x()) + num2string(q1.y()) + num2string(q2.x()) + num2string(q2.y());
}



} // end namespace CGAL

#endif 
