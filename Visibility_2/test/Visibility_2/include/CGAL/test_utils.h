// Copyright (c) 2013 Technical University Braunschweig (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s):  Kan Huang <huangkandiy@gmail.com>
//             Francisc Bungiu <fbungiu@gmail.com>
//             Michael Hemmer <michael.hemmer@cgal.org>

#ifndef CGAL_TEST_UTILS_H
#define CGAL_TEST_UTILS_H

//#define RESET   "\033[0m"
//#define RED     "\033[31m"
//#define GREEN   "\033[32m"

#define RESET   ""
#define RED     ""
#define GREEN   ""

#include <cassert>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <CGAL/Gmpq.h>
#include <CGAL/Timer.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Visibility_2/visibility_utils.h>

namespace CGAL {

enum Query_choice {VERTEX, EDGE, FACE};  

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

template<class Arrangement_2>
CGAL::Object get_location(
                  const Arrangement_2 &arr, 
                  const typename Arrangement_2::Geometry_traits_2::Point_2 &q) {

  typedef CGAL::Arr_naive_point_location<Arrangement_2> Naive_pl;
  Naive_pl naive_pl(arr);

  // Perform the point-location query.
  CGAL::Object obj = naive_pl.locate(q);
  return obj;
}   

template<class Arrangement_2> 
bool is_inside_face(
                  const Arrangement_2 &arr, 
                  const typename Arrangement_2::Face_const_handle face,
                  const typename Arrangement_2::Geometry_traits_2::Point_2 &q) {

  CGAL::Object obj = get_location<Arrangement_2>(arr, q);
  typename Arrangement_2::Face_const_handle f;
  if (CGAL::assign (f, obj)) {
    // q is located inside a face:
    if (!f->is_unbounded()) {
      if (f == face) {
        return true;
      }
    }
  }
  return false;
}

template <class _Arrangement_2> 
bool create_arrangement_from_dat_file(std::ifstream &input,
                                       _Arrangement_2 &arr) {
  arr.clear();

  typedef _Arrangement_2                                      Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2           Geometry_traits_2;
  typedef typename Arrangement_2::Face_handle                 Face_handle;
  typedef typename Geometry_traits_2::Segment_2               Segment_2;
  typedef typename Geometry_traits_2::Point_2                 Point_2;
  typedef typename Geometry_traits_2::FT                      Number_type;

  if (input) {
    std::string curr_line;
    std::vector<Point_2> isolated_vertices;
    std::getline(input, curr_line);
    std::stringstream convert(curr_line);
    int number_of_isolated_vertices;
    convert >> number_of_isolated_vertices;
    Face_handle uface = arr.unbounded_face();
    for (int i = 0 ; i < number_of_isolated_vertices ; i++) {
      std::getline(input, curr_line);
      std::istringstream iss(curr_line);
      std::string x, y;
      iss >> x >> y;
      arr.insert_in_face_interior(Point_2(string2num<Number_type>(x), 
                                          string2num<Number_type>(y)),
                                  uface);
    }
    std::vector<Segment_2> edges;
    int number_of_edges;
    std::getline(input, curr_line);
    std::stringstream convert2(curr_line);
    convert2 >> number_of_edges;
    for (int i = 0 ; i < number_of_edges ; i++) {
      std::getline(input, curr_line);
      std::string x1, y1, x2, y2;
      std::istringstream iss(curr_line);
      iss >> x1 >> y1 >> x2 >> y2;
      edges.push_back(Segment_2(Point_2(string2num<Number_type>(x1), 
                                        string2num<Number_type>(y1)), 
                                Point_2(string2num<Number_type>(x2), 
                                        string2num<Number_type>(y2))));
    }
    CGAL::insert(arr, edges.begin(), edges.end());
    return true;
  }
  else {
    return false;
  }
}

template <class Arrangement_2>
void regularize(Arrangement_2& arr){
  //std::cout << "\n regularize" << std::endl; 
  // remove all edges with the same face on both sides 
  typedef typename Arrangement_2::Edge_iterator EIT; 
  for(EIT eit = arr.edges_begin(); eit != arr.edges_end();){
    if(eit->face() == eit->twin()->face()){
      // arr.remove_edge(eit++,false,false); did not compile 
      EIT eh = eit; 
      ++eit;
      arr.remove_edge(eh,false,false);
    }else{
      ++eit; 
    }
  }
  // remove all isolated vertices (also those left from prvious step)
  typedef typename Arrangement_2::Vertex_iterator VIT; 
  for(VIT vit = arr.vertices_begin(); vit != arr.vertices_end();){
    if(vit->degree()== 0){
      VIT vh = vit;
      vit++; 
      arr.remove_isolated_vertex(vh);
    }else{
      vit++; 
    }
  }
  //std::cout << "regularize done" << std::endl; 
}

template <class Visibility_2>
bool run_test_case_from_file(Visibility_2 visibility, std::ifstream &input) {
  typedef typename Visibility_2::Input_arrangement_2          Input_arrangement_2;
  typedef typename Visibility_2::Output_arrangement_2         Output_arrangement_2;
  typedef typename Input_arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Geometry_traits_2::Point_2                 Point_2; 
  typedef typename Geometry_traits_2::FT                      Number_type;
  typedef typename Input_arrangement_2::Halfedge_around_vertex_const_circulator
                                        Halfedge_around_vertex_const_circulator;
  typedef typename Input_arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Input_arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Input_arrangement_2::Vertex_const_handle   Vertex_const_handle;

  Input_arrangement_2 arr_in;
  Output_arrangement_2 arr_correct_out;
  Output_arrangement_2 arr_out;
    
  std::string curr_line;
  while (std::getline(input, curr_line)) {
    if (curr_line[0] != '#' && curr_line[0] != '/')
      break;
  }
  std::stringstream convert(curr_line);
  std::string x, y;
  convert >> x >> y;
  Point_2 query_pt(string2num<Number_type>(x),
                   string2num<Number_type>(y));
  std::getline(input, curr_line);
  std::string x1, y1;
  std::istringstream iss(curr_line);
  iss >> x1 >> y1;
  Point_2 reference_pt(string2num<Number_type>(x1),
                       string2num<Number_type>(y1));

  std::getline(input, curr_line);
  if (!create_arrangement_from_dat_file<Input_arrangement_2>(input, arr_in)) {
    return false;
  }
  visibility.attach(arr_in);
  std::getline(input, curr_line);
  if (!create_arrangement_from_dat_file<Output_arrangement_2>
                                                   (input, arr_correct_out)) {
    return false;
  }

  if(Visibility_2::Regularization_tag::value){
    regularize(arr_correct_out);   
  }


  CGAL::Object obj = CGAL::get_location<Input_arrangement_2>(arr_in, query_pt);
  Face_const_handle f;
  Halfedge_const_handle e;
  Vertex_const_handle v;

  if (CGAL::assign (f, obj)) {
    if (!f->is_unbounded()) {
      visibility.compute_visibility(query_pt, f, arr_out);
    }
  }
  else if (CGAL::assign(e, obj)) {
    if (e->source()->point() == reference_pt) {
      visibility.compute_visibility(query_pt, e, arr_out);
    }
    else {
      visibility.compute_visibility(query_pt, e->twin(), arr_out);
    }
  }
  else if (CGAL::assign(v, obj)) {
    // Check which halfedge we want in the query
    Halfedge_around_vertex_const_circulator he_circ = v->incident_halfedges();
    Halfedge_around_vertex_const_circulator he_curr = he_circ;
    do {
      if (he_curr->source()->point() == reference_pt) {
        visibility.compute_visibility(query_pt, he_curr, arr_out);
      }
    } while (++he_curr != he_circ);
  }


  if (!test_are_equal<Output_arrangement_2>(arr_out, arr_correct_out)) {
    std::cout<<"the result is:\n";
    CGAL::Visibility_2::print_arrangement(arr_out);
    std::cout<<"however, the expected answer is:\n";
    CGAL::Visibility_2::print_arrangement(arr_correct_out);
    return false;
  }
  visibility.detach();
  return true;
}

template <class Visibility_2> 
void run_tests(int case_number_simple, int case_number_non_simple) {
  
  Visibility_2 visibility;
  bool one_failed = false; 
  if (Visibility_2::Supports_simple_polygon_tag::value 
      && case_number_simple >= 1) {
    int cnt = 0;
    int cnt_passed = 0;
    std::cout << "    Running simple polygon test cases...\n";
    for (int i = 1 ; i <= case_number_simple ; i++) {
      std::string input_arr_file("data/test_simple_polygon_");
      input_arr_file += num2string<int>(i);
      input_arr_file += ".dat";
      std::cout << "        Running test " 
                << GREEN << input_arr_file << RESET 
                << " - ";
      std::ifstream input(input_arr_file.c_str());
      if (run_test_case_from_file<Visibility_2>(visibility, input)) {
        cnt_passed++;
        std::cout << GREEN << "Done!" << RESET << std::endl;
      }
      else {
        one_failed = true;
        std::cout << RED << "Failed!" << RESET << std::endl; 
      }
      cnt++;
    }
    std::cout << "    Visibility_2 object passed " << cnt_passed 
              << "/" << cnt << " tests"      
              << " (";
    double result = (double)cnt_passed/cnt*100;
    if (result > 99.9) {
      std::cout << GREEN << result << "%" << RESET;
    }
    else {
      std::cout << RED << result << "%" << RESET;
    }
    std::cout << ")" << std::endl;
  }
  if (Visibility_2::Supports_general_polygon_tag::value
      && case_number_non_simple >= 1) {

    int cnt = 0;
    int cnt_passed = 0;
    std::cout << "    Running non-simple polygon test cases...\n";
    for (int i = 1 ; i <= case_number_non_simple ; i++) {
      std::string input_arr_file("data/test_non_simple_polygon_");
      input_arr_file += num2string<int>(i);
      input_arr_file += ".dat";
      std::cout << "        Running test "
                << GREEN << input_arr_file << RESET 
                << " - ";      
      std::ifstream input(input_arr_file.c_str());
      if (run_test_case_from_file<Visibility_2>(visibility, input)) {
        cnt_passed++;
        std::cout << GREEN << "Done!" << RESET << std::endl;
      }
      else {
        one_failed = true;
        std::cout << RED << "Failed!" << RESET << std::endl;
      }
      cnt++;
    }
    std::cout << "    Visibility_2 object passed " << cnt_passed
              << "/" << cnt  << " tests"     
              << " (";
    double result = (double)cnt_passed/cnt*100;
    if (result > 99.9) {
      std::cout << GREEN << result << "%" << RESET;
    }
    else {
      std::cout << RED << result << "%" << RESET;
    }
    std::cout << ")" << std::endl;
  }
  if (one_failed) {
    assert(false);
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
      points.push_back(Point_2(string2num<Number_type>(n1), 
                               string2num<Number_type>(n2)));
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
void create_arrangement_from_env_file(_Arrangement_2 &arr, std::ifstream& input) {
  typedef _Arrangement_2 								                Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2	    Geometry_traits_2;
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
    for (int i = 0 ; i < number_of_polygons ; i++) {
      std::vector<Point_2> points;
      std::vector<Segment_2> segments;

      std::getline(input, line);
      std::stringstream convert(line);
      int number_of_vertices;
      convert >> number_of_vertices;
      for (int j = 0; j < number_of_vertices; j++) {
        std::getline(input, line);
        std::string n1, n2;
        std::istringstream iss(line);
        iss >> n1 >> n2;
        points.push_back(Point_2(string2num<Number_type>(n1),   
                                 string2num<Number_type>(n2)));
      }
      for (int j = 0; j < number_of_vertices-1 ; j++) {
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
bool compare_arr_by_edges(const Arrangement_2& arr1, const Arrangement_2& arr2){
  std::set<std::string> s1;
  typedef typename Arrangement_2::Edge_const_iterator Edge_const_iterator;
  for (Edge_const_iterator eit = arr1.edges_begin(); 
                           eit != arr1.edges_end(); ++eit) {
    s1.insert(edge2string(eit->target()->point(), eit->source()->point()));
  }
  std::set<std::string> s2;
  for (Edge_const_iterator eit = arr2.edges_begin();  
                           eit != arr2.edges_end(); ++eit) {
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
  return num2string(q1.x()) + num2string(q1.y()) 
  + num2string(q2.x()) + num2string(q2.y());
}

template<class Point_2, class Number_type> 
Point_2 random_linear_interpolation(const Point_2 &p, const Point_2 &q) {

//  srand(time(NULL));
  Number_type min_x, max_x;
  Number_type y0, y1, x0, x1;

  if (p.x() < q.x()) {
    min_x = p.x();
    max_x = q.x();
    x0 = p.x();
    x1 = q.x();
    y0 = p.y();
    y1 = q.y();
  }
  else {
    min_x = q.x();
    max_x = p.x();
    x0 = q.x();
    x1 = p.x();
    y0 = q.y();
    y1 = p.y();
  }

  Number_type x_normalized = rand()/static_cast<Number_type>(RAND_MAX);
  Number_type x = min_x + static_cast<Number_type>(x_normalized*(max_x - min_x));
  if (x == max_x && x == min_x) {
    Number_type min_y, max_y;
    if (p.y() < q.y()) {
      min_y = p.y();
      max_y = q.y();
    }
    else {
      min_y = q.y();
      max_y = p.y();
    }
    Number_type y_normalized = rand()/static_cast<Number_type>(RAND_MAX);
    Number_type y = min_y + static_cast<Number_type>(y_normalized*(max_y - min_y));
    return Point_2(x, y);
  }
  else {
    Number_type y = y0 + (y1 - y0)*(x - x0)/(x1 - x0);
    return Point_2(x, y);
  }
}

//make sure q is in fh or on the bound.
template<class Visibility_2>
bool is_star_shape(const typename Visibility_2::Point_2& q,
                   const typename Visibility_2::Face_handle fh) {
  typedef typename Visibility_2::Output_arrangement_2
                                                            Output_arrangement_2;
  typedef typename Output_arrangement_2::Face_const_handle   Face_const_handle;
  typedef typename Output_arrangement_2::Halfedge_handle     Halfedge_handle;
  typedef typename Output_arrangement_2::Geometry_traits_2   Geometry_traits_2;
  typedef typename Output_arrangement_2::Edge_iterator       Edge_iterator;
  typedef typename Geometry_traits_2::Point_2               Point_2;
  typedef typename Geometry_traits_2::Segment_2             Segment_2;
  if (fh->is_unbounded())
    return false;
  if (fh->has_outer_ccb()) {
    typename Output_arrangement_2::Ccb_halfedge_const_circulator curr, circ;
    curr = circ = fh->outer_ccb();
    do {
      if (CGAL::right_turn(q, curr->source()->point(), curr->target()->point())) {
        return false;
      }
//      Point_2 p = curr->target()->point();
//      typename Output_arrangement_2::Ccb_halfedge_const_circulator curr1, circ1;
//      curr1 = circ1 = fh->outer_ccb();
//      do {
//        Segment_2 intersect_s;
//        Point_2 intersect_p;
//        Point_2 source = curr1->source()->point();
//        Point_2 target = curr1->target()->point();
//        int i = intersect_seg<Segment_2, Point_2>(Segment_2(p, q), Segment_2(source, target), intersect_s, intersect_p);
//        if (i == 1 && intersect_p != source && intersect_p != target && intersect_p != q)
//          return false;
//      } while (++curr1 != circ1);
    } while (++curr != circ);
  }
  return true;
}

template <class Visibility_2_fst, class Visibility_2_snd>
void simple_benchmark_one_unit(
          typename Visibility_2_fst::Input_arrangement_2 &arr,
          const Query_choice &choice,
          typename Visibility_2_fst::Input_arrangement_2::Face_const_handle &fit,
          Visibility_2_fst visibility_fst,
          Visibility_2_snd visibility_snd,
          double& qtime1,
          double& qtime2,
          double& ptime1,
          double& ptime2,
          int& case_cnt) {

  typedef typename Visibility_2_fst::Input_arrangement_2    Input_arrangement_2;
  typedef typename Input_arrangement_2::Face_const_handle   Face_const_handle;
  typedef typename Visibility_2_fst::Output_arrangement_2   Output_arrangement_2;
  typedef typename Input_arrangement_2::Halfedge_const_handle
                                                            Halfedge_const_handle;
  typedef typename Input_arrangement_2::Geometry_traits_2   Geometry_traits_2;
  typedef typename Input_arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;

  typedef typename Output_arrangement_2::Face_handle        Face_handle;
  typedef typename Geometry_traits_2::Point_2               Point_2;
  typedef typename Geometry_traits_2::FT                    Number_type;
  typedef Timer Benchmark_timer;

  Benchmark_timer timer;
//  std::cout << "INPUT:\n";
//  CGAL::Visibility_2::print_arrangement<Input_arrangement_2>(arr);
//  std::cout << "END INPUT\n";

  timer.start();
  visibility_fst.attach(arr);
  timer.stop();
  ptime1 += timer.time();
//  std::cout << "    Time to attach to first object: "
//            << GREEN << timer.time() << " sec" << RESET << std::endl;

  timer.reset();
  timer.start();
  visibility_snd.attach(arr);
  timer.stop();
  ptime2 += timer.time();
//  std::cout << "    Time to attach to second object: "
//            << GREEN << timer.time() << " sec" << RESET << std::endl;

  Ccb_halfedge_const_circulator circ = fit->outer_ccb();
  Ccb_halfedge_const_circulator curr = circ;

  do {
    Halfedge_const_handle he = curr;
    Point_2 curr_query_pt;
    bool selected_query_pt = true;
    switch (choice) {
      case VERTEX:
        curr_query_pt = he->target()->point();
        break;
      case EDGE:
//        std::cout << "Generating qpoint on: " << he->curve() << std::endl;
        curr_query_pt = random_linear_interpolation<Point_2, Number_type>
                      (he->source()->point(), he->target()->point());
        break;
      case FACE:
        Ccb_halfedge_const_circulator curr_next = curr;
        curr_next++;
        Halfedge_const_handle he_next = curr_next;
        Point_2 p1 = he->source()->point();
        Point_2 p2 = he->target()->point();
        Point_2 p3 = he_next->target()->point();
        Point_2 avg((p1.x() + p2.x() + p3.x())/3, (p1.y() + p2.y() + p3.y())/3);
        if (is_inside_face<Input_arrangement_2>(arr, fit, avg)) {
          curr_query_pt = avg;
        }
        else {
          selected_query_pt = false;
        }
        break;
    }
    if (!selected_query_pt) {
      curr++;
      if (curr == circ) {
        break;
      }
      continue;
    }
//    std::cout << "    Running with qpoint: "
//              << RED << curr_query_pt << RESET <<  std::endl;
    Input_arrangement_2 out_arr_fst, out_arr_snd;
    timer.reset();
    timer.start();
    Face_handle f_fst;
    if (choice == FACE) {
      f_fst = visibility_fst.compute_visibility(curr_query_pt, fit, out_arr_fst);
    }
    else {
      f_fst = visibility_fst.compute_visibility(curr_query_pt, he, out_arr_fst);
    }
    timer.stop();
    qtime1 += timer.time();
//    std::cout << "        Time to compute visibility region using first object: "
//              << GREEN << timer.time() << " sec" << RESET << std::endl;
    timer.reset();
    if ( !is_star_shape<Visibility_2_fst>(curr_query_pt, f_fst) ) {
      std::cout << RED << "         Warning: the first output is not star-shape." << RESET << std::endl;
    }
    timer.start();
    Face_handle f_snd;
    if (choice == FACE) {
      f_snd = visibility_snd.compute_visibility(curr_query_pt, fit, out_arr_snd);
    }
    else {
      f_snd = visibility_snd.compute_visibility(curr_query_pt, he, out_arr_snd);
    }
    timer.stop();
    if ( !is_star_shape<Visibility_2_snd>(curr_query_pt, f_snd) ) {
      std::cout << RED << "         Warning: the second output is not star-shape." << RESET << std::endl;
    }
    qtime2 += timer.time();
//    std::cout << "        Time to compute visibility region using second object: "
//              << GREEN << timer.time() << " sec" << RESET << std::endl;

 //   CGAL::Visibility_2::print_arrangement<Output_arrangement_2>(out_arr_snd);
    if (! CGAL::test_are_equal<Output_arrangement_2> (out_arr_fst, out_arr_snd)) {
      if (choice == FACE) {
        std::cout << "query in Face:\n";
        CGAL::Visibility_2::print_simple_face<Face_const_handle, Ccb_halfedge_const_circulator>(fit);
      }
      std::cout << "the query point is" << curr_query_pt << std::endl;
      std::cout << RED << "two outputs are different:\n"
                << "first output is:\n" << RESET;
      CGAL::Visibility_2::print_arrangement(out_arr_fst);
      std::cout << RED << "second output is:\n"<< RESET;
      CGAL::Visibility_2::print_arrangement(out_arr_snd);
      assert(false);
    }
    case_cnt++;
  } while (++curr != circ);
}

template<class Visibility_2_fst, class Visibility_2_snd>
void simple_benchmark(Visibility_2_fst &visibility_fst,
               Visibility_2_snd &visibility_snd,
               const Query_choice &choice,
               std::ifstream &input) {

  typedef typename Visibility_2_fst::Input_arrangement_2
                                                            Input_arrangement_2;
  typedef typename Input_arrangement_2::Halfedge_const_handle
                                                          Halfedge_const_handle;
  typedef typename Input_arrangement_2::Geometry_traits_2   Geometry_traits_2;
  typedef typename Input_arrangement_2::Face_const_iterator Face_const_iterator;
  typedef typename Input_arrangement_2::Face_const_handle   Face_const_handle;
  typedef typename Input_arrangement_2::Hole_const_iterator Hole_const_iterator;
  typedef typename Input_arrangement_2::Halfedge_const_handle
                                                          Halfedge_const_handle;
  typedef typename Input_arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;
  typedef typename Geometry_traits_2::Point_2               Point_2;
  typedef typename Geometry_traits_2::Segment_2             Segment_2;

  assert(Visibility_2_fst::Regularization_tag::value == Visibility_2_snd::Regularization_tag::value);

  Input_arrangement_2 arr;
  create_arrangement_from_env_file<Input_arrangement_2>(arr, input);
//  std::cout << "Input arrangement has: "
//            << GREEN << arr.number_of_faces()-1 << RESET
//            << " faces." << std::endl;
  int query_cnt(0);
  double qtime1(0), qtime2(0), ptime1(0), ptime2(0);
  if (Visibility_2_fst::Supports_general_polygon_tag::value
    && Visibility_2_snd::Supports_general_polygon_tag::value) {
    int cnt(1);
    Face_const_iterator fit;

    for (fit = arr.faces_begin() ; fit != arr.faces_end() ; fit++) {
      if (!fit->is_unbounded()) {
//        std::cout << "Benchmarking with face "
//                  << GREEN << cnt << RESET << " ..." << std::endl;
        simple_benchmark_one_unit<Visibility_2_fst, Visibility_2_snd>(arr,
                                                               choice,
                                                               fit,
                                                               visibility_fst,
                                                               visibility_snd,
                                                                      qtime1,
                                                                      qtime2,
                                                                      ptime1,
                                                                      ptime2,
                                                                      query_cnt);
      }
      cnt++;
    }
    std::cout << "Preprocessing: "  << std::endl
              << "Model 1 uses " << ptime1/cnt << "  sec" << std::endl
              << "Model 2 uses " << ptime2/cnt << "  sec" << std::endl;
    std::cout << query_cnt << " queries are done.\n"
              << "Model 1 uses " << qtime1 << "  sec" << std::endl
              << "Model 2 uses " << qtime2 << "  sec" << std::endl;
    std::cout << "total times are:" << std::endl
              << "Model 1 uses " << ptime1/cnt + qtime1 << "  sec" << std::endl
              << "Model 2 uses " << qtime1/cnt + qtime2 << "  sec" << std::endl;
  }
  else {  // Only run the benchmark on the outer loop of the arrangement
    Face_const_iterator fit;
    // See which face has holes
    int cnt(1);
    for (fit = arr.faces_begin() ; fit != arr.faces_end() ; fit++) {
      if (!fit->is_unbounded()) {
//        std::cout << "Benchmarking with face "
//                  << GREEN << cnt << RESET << " ..." << std::endl;
        Hole_const_iterator hit;
        bool has_holes = false;
        for (hit = fit->holes_begin() ; hit != fit->holes_end() ; hit++) {
          has_holes = true;
          break;
        }
        if (has_holes && fit->has_outer_ccb()) {
          Input_arrangement_2 arr_trimmed;
          std::vector<Segment_2> segments;
          Ccb_halfedge_const_circulator circ = fit->outer_ccb();
          Ccb_halfedge_const_circulator curr = circ;
          do {
            Halfedge_const_handle he = curr;
            segments.push_back(Segment_2(he->source()->point(), he->target()->point()));
          } while (++curr != circ);
          CGAL::insert(arr_trimmed, segments.begin(), segments.end());
          Face_const_handle fch;
          if (arr_trimmed.faces_begin()->is_unbounded()) {
            fch = ++arr_trimmed.faces_begin();
          }
          else {
            fch = arr_trimmed.faces_begin();
          }
          simple_benchmark_one_unit<Visibility_2_fst, Visibility_2_snd>(arr_trimmed,
                                                                 choice,
                                                                 fch,
                                                                 visibility_fst,
                                                                 visibility_snd,
                                                                        qtime1,
                                                                        qtime2,
                                                                        ptime1,
                                                                        ptime2,
                                                                        query_cnt
                                                                 );

          //CGAL::Visibility_2::print_arrangement(arr_trimmed);
        }
        else if (!has_holes) {
          simple_benchmark_one_unit<Visibility_2_fst, Visibility_2_snd>(arr,
                                                                 choice,
                                                                 fit,
                                                                 visibility_fst,
                                                                 visibility_snd,
                                                                        qtime1,
                                                                        qtime2,
                                                                        ptime1,
                                                                        ptime2,
                                                                        query_cnt
                                                                 );

        }
        cnt++;
      }
    }
    std::cout << "attach " << cnt << " times." << std::endl
              << "Model 1 uses " << ptime1 << "  sec" << std::endl
              << "Model 2 uses " << ptime2 << "  sec" << std::endl;
    std::cout << query_cnt << " queries are done.\n"
              << "Model 1 uses " << qtime1 << "  sec" << std::endl
              << "Model 2 uses " << qtime2 << "  sec" << std::endl;
    std::cout << "total times are:" << std::endl
              << "Model 1 uses " << ptime1 + qtime1 << "  sec" << std::endl
              << "Model 2 uses " << qtime1 + qtime2 << "  sec" << std::endl;
  }
}
template <class Visibility_2_fst, class Visibility_2_snd>
void benchmark_one_unit(
          typename Visibility_2_fst::Input_arrangement_2 &arr,
          const Query_choice &choice,
          typename Visibility_2_fst::Input_arrangement_2::Face_const_handle &fit,
          Visibility_2_fst visibility_fst, 
          Visibility_2_snd visibility_snd) {

  typedef typename Visibility_2_fst::Input_arrangement_2    Input_arrangement_2;
  typedef typename Input_arrangement_2::Face_const_handle   Face_const_handle;
  typedef typename Visibility_2_fst::Output_arrangement_2   Output_arrangement_2;
  typedef typename Input_arrangement_2::Halfedge_const_handle     
                                                            Halfedge_const_handle;
  typedef typename Input_arrangement_2::Geometry_traits_2   Geometry_traits_2;
  typedef typename Input_arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;

  typedef typename Output_arrangement_2::Face_handle        Face_handle;
  typedef typename Geometry_traits_2::Point_2               Point_2;
  typedef typename Geometry_traits_2::FT                    Number_type;
  typedef Timer Benchmark_timer;

  Benchmark_timer timer;
//  std::cout << "INPUT:\n";
//  CGAL::Visibility_2::print_arrangement<Input_arrangement_2>(arr);
//  std::cout << "END INPUT\n";

  timer.start();
  visibility_fst.attach(arr); 
  timer.stop();
  std::cout << "    Time to attach to first object: " 
            << GREEN << timer.time() << " sec" << RESET << std::endl;

  timer.reset();
  timer.start();
  visibility_snd.attach(arr);
  timer.stop();
  std::cout << "    Time to attach to second object: " 
            << GREEN << timer.time() << " sec" << RESET << std::endl;

  Ccb_halfedge_const_circulator circ = fit->outer_ccb();
  Ccb_halfedge_const_circulator curr = circ;

  do {
    Halfedge_const_handle he = curr;
    Point_2 curr_query_pt;
    bool selected_query_pt = true;
    switch (choice) {
      case VERTEX:
        curr_query_pt = he->target()->point();        
        break;
      case EDGE:
        std::cout << "Generating qpoint on: " << he->curve() << std::endl;
        curr_query_pt = random_linear_interpolation<Point_2, Number_type>
                      (he->source()->point(), he->target()->point());
        break;
      case FACE:
        Ccb_halfedge_const_circulator curr_next = curr;
        curr_next++;
        Halfedge_const_handle he_next = curr_next;
        Point_2 p1 = he->source()->point();
        Point_2 p2 = he->target()->point();
        Point_2 p3 = he_next->target()->point();
        Point_2 avg((p1.x() + p2.x() + p3.x())/3, (p1.y() + p2.y() + p3.y())/3);
        if (is_inside_face<Input_arrangement_2>(arr, fit, avg)) {
          curr_query_pt = avg;
        }
        else {
          selected_query_pt = false;
        }
        break;
    }
    if (!selected_query_pt) {
      curr++;
      if (curr == circ) {
        break;
      }
      continue;
    }
    std::cout << "    Running with qpoint: " 
              << RED << curr_query_pt << RESET <<  std::endl;
    Input_arrangement_2 out_arr_fst, out_arr_snd;
    timer.reset();
    timer.start();
    Face_handle f_fst;
    if (choice == FACE) {
      f_fst = visibility_fst.compute_visibility(curr_query_pt, fit, out_arr_fst);
    }
    else {
      f_fst = visibility_fst.compute_visibility(curr_query_pt, he, out_arr_fst);
    }
    timer.stop();

    std::cout << "        Time to compute visibility region using first object: " 
              << GREEN << timer.time() << " sec" << RESET << std::endl;
    timer.reset();
    if ( !is_star_shape<Visibility_2_fst>(curr_query_pt, f_fst) ) {
      std::cout << RED << "         Warning: the first output is not star-shape." << RESET << std::endl;
    }
    timer.start();
    Face_handle f_snd;
    if (choice == FACE) {
      f_snd = visibility_snd.compute_visibility(curr_query_pt, fit, out_arr_snd);
    }
    else {
      f_snd = visibility_snd.compute_visibility(curr_query_pt, he, out_arr_snd);
    }
    timer.stop();
    if ( !is_star_shape<Visibility_2_snd>(curr_query_pt, f_snd) ) {
      std::cout << RED << "         Warning: the second output is not star-shape." << RESET << std::endl;
    }

    std::cout << "        Time to compute visibility region using second object: " 
              << GREEN << timer.time() << " sec" << RESET << std::endl;

 //   CGAL::Visibility_2::print_arrangement<Output_arrangement_2>(out_arr_snd);
    if (! CGAL::test_are_equal<Output_arrangement_2> (out_arr_fst, out_arr_snd)) {
      if (choice == FACE) {
        std::cout << "query in Face:\n";
        CGAL::Visibility_2::print_simple_face<Face_const_handle, Ccb_halfedge_const_circulator>(fit);
      }
      std::cout << RED << "two outputs are different:\n"
                << "first output is:\n" << RESET;
      CGAL::Visibility_2::print_arrangement(out_arr_fst);
      std::cout << RED << "second output is:\n"<< RESET;
      CGAL::Visibility_2::print_arrangement(out_arr_snd);
      assert(false);
    }
  } while (++curr != circ);
}



template<class Visibility_2_fst, class Visibility_2_snd>
void benchmark(Visibility_2_fst &visibility_fst, 
               Visibility_2_snd &visibility_snd,
               const Query_choice &choice,
               std::ifstream &input) {

  typedef typename Visibility_2_fst::Input_arrangement_2
                                                            Input_arrangement_2;
  typedef typename Input_arrangement_2::Halfedge_const_handle       
                                                          Halfedge_const_handle;
  typedef typename Input_arrangement_2::Geometry_traits_2   Geometry_traits_2;
  typedef typename Input_arrangement_2::Face_const_iterator Face_const_iterator;
  typedef typename Input_arrangement_2::Face_const_handle   Face_const_handle;
  typedef typename Input_arrangement_2::Hole_const_iterator Hole_const_iterator;
  typedef typename Input_arrangement_2::Halfedge_const_handle     
                                                          Halfedge_const_handle;
  typedef typename Input_arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;
  typedef typename Geometry_traits_2::Point_2               Point_2;
  typedef typename Geometry_traits_2::Segment_2             Segment_2;

  assert(Visibility_2_fst::Regularization_tag::value == Visibility_2_snd::Regularization_tag::value);

  Input_arrangement_2 arr;
  create_arrangement_from_env_file<Input_arrangement_2>(arr, input);
  std::cout << "Input arrangement has: " 
            << GREEN << arr.number_of_faces()-1 << RESET
            << " faces." << std::endl;
  if (Visibility_2_fst::Supports_general_polygon_tag::value 
    && Visibility_2_snd::Supports_general_polygon_tag::value) {
    int cnt(1);
    Face_const_iterator fit;
    for (fit = arr.faces_begin() ; fit != arr.faces_end() ; fit++) {
      if (!fit->is_unbounded()) {
        std::cout << "Benchmarking with face " 
                  << GREEN << cnt << RESET << " ..." << std::endl;
        benchmark_one_unit<Visibility_2_fst, Visibility_2_snd>(arr,
                                                               choice,
                                                               fit,
                                                               visibility_fst,    
                                                               visibility_snd);
      }
      cnt++;
    }
  }
  else {  // Only run the benchmark on the outer loop of the arrangement
    Face_const_iterator fit;
    // See which face has holes
    int cnt(1);
    for (fit = arr.faces_begin() ; fit != arr.faces_end() ; fit++) {
      if (!fit->is_unbounded()) {
        std::cout << "Benchmarking with face " 
                  << GREEN << cnt << RESET << " ..." << std::endl;
        Hole_const_iterator hit;
        bool has_holes = false;
        for (hit = fit->holes_begin() ; hit != fit->holes_end() ; hit++) {
          has_holes = true;
          break;
        }
        if (has_holes && fit->has_outer_ccb()) {
          Input_arrangement_2 arr_trimmed;
          std::vector<Segment_2> segments;
          Ccb_halfedge_const_circulator circ = fit->outer_ccb();
          Ccb_halfedge_const_circulator curr = circ;
          do {
            Halfedge_const_handle he = curr;
            segments.push_back(Segment_2(he->source()->point(), he->target()->point()));
          } while (++curr != circ);
          CGAL::insert(arr_trimmed, segments.begin(), segments.end());
          Face_const_handle fch;
          if (arr_trimmed.faces_begin()->is_unbounded()) {
            fch = ++arr_trimmed.faces_begin();
          } 
          else {
            fch = arr_trimmed.faces_begin();
          }
          benchmark_one_unit<Visibility_2_fst, Visibility_2_snd>(arr_trimmed,
                                                                 choice,
                                                                 fch,
                                                                 visibility_fst,  
                                                                 visibility_snd
                                                                 );
          //CGAL::Visibility_2::print_arrangement(arr_trimmed);
        }
        else if (!has_holes) {
          benchmark_one_unit<Visibility_2_fst, Visibility_2_snd>(arr,
                                                                 choice,
                                                                 fit,
                                                                 visibility_fst,  
                                                                 visibility_snd
                                                                 );
        }
        cnt++;
      }
    }
  }
}

template<class Segment_2, class Point_2>
int intersect_seg(const Segment_2& seg1, const Segment_2& seg2, Segment_2& seg_out, Point_2& p_out)
{
    CGAL::Object result = CGAL::intersection(seg1, seg2);
    if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        p_out = *ipoint;
        return 1;
    } else
        if (const Segment_2 *iseg = CGAL::object_cast<Segment_2>(&result)) {
            seg_out = *iseg;
            return 2;
        } else {
            return 0;
        }
}


template<class Visibility_2>
void test_star_shape_one_face(  typename Visibility_2::Input_arrangement_2 &arr,
                                const Query_choice &choice,
                                typename Visibility_2::Input_arrangement_2::Face_const_handle &fit,
                                Visibility_2 visibility)
{

  typedef typename Visibility_2::Input_arrangement_2        Input_arrangement_2;
  typedef typename Visibility_2::Output_arrangement_2       Output_arrangement_2;
  typedef typename Input_arrangement_2::Halfedge_const_handle
                                                            Halfedge_const_handle;
  typedef typename Input_arrangement_2::Geometry_traits_2   Geometry_traits_2;
  typedef typename Input_arrangement_2::Face_const_handle   Face_const_handle;


  typedef typename Input_arrangement_2::Ccb_halfedge_const_circulator
                                                            Ccb_halfedge_const_circulator;

  typedef typename Output_arrangement_2::Face_handle        Face_handle;
  typedef typename Geometry_traits_2::Point_2               Point_2;
  typedef typename Geometry_traits_2::FT                    Number_type;
  typedef Timer Benchmark_timer;

  visibility.attach(arr);

  Ccb_halfedge_const_circulator circ = fit->outer_ccb();
  Ccb_halfedge_const_circulator curr = circ;
  do {
    Halfedge_const_handle he = curr;
    Point_2 curr_query_pt;
    bool selected_query_pt = true;
    switch (choice) {
      case VERTEX:
        curr_query_pt = he->target()->point();
        break;
      case EDGE:
        curr_query_pt = random_linear_interpolation<Point_2, Number_type>
                      (he->source()->point(), he->target()->point());
        break;
      case FACE:
        Ccb_halfedge_const_circulator curr_next = curr;
        curr_next++;
        Halfedge_const_handle he_next = curr_next;
        Point_2 p1 = he->source()->point();
        Point_2 p2 = he->target()->point();
        Point_2 p3 = he_next->target()->point();
        Point_2 avg((p1.x() + p2.x() + p3.x())/3, (p1.y() + p2.y() + p3.y())/3);
        if (is_inside_face<Input_arrangement_2>(arr, fit, avg)) {
          curr_query_pt = avg;
        }
        else {
          selected_query_pt = false;
        }
        break;
    }
    if (!selected_query_pt) {
      curr++;
      continue;
    }
    std::cout << "    Running with qpoint: "
              << RED << curr_query_pt << RESET <<  std::endl;
    Output_arrangement_2 out_arr;
    Face_handle fh;
    if (choice == FACE) {
      fh = visibility.compute_visibility(curr_query_pt, fit, out_arr);
    }
    else {
      fh = visibility.compute_visibility(curr_query_pt, he, out_arr);
    }
    if ( !is_star_shape<Visibility_2>(curr_query_pt, fh)) {
      std::cout << RED << "     The face is not a star shape to qpoint." << RESET <<  std::endl;
    }
  } while (++curr != circ);
}

template<class Visibility_2>
void test_star_shape(Visibility_2 &visibility,
               const Query_choice &choice,
               std::ifstream &input) {

  typedef typename Visibility_2::Output_arrangement_2
                                                            Output_arrangement_2;
  typedef typename Visibility_2::Input_arrangement_2
                                                            Input_arrangement_2;
  typedef typename Output_arrangement_2::Halfedge_const_handle
                                                          Halfedge_const_handle;
  typedef typename Output_arrangement_2::Geometry_traits_2   Geometry_traits_2;
  typedef typename Output_arrangement_2::Face_const_iterator Face_const_iterator;
  typedef typename Output_arrangement_2::Face_const_handle   Face_const_handle;
  typedef typename Output_arrangement_2::Hole_const_iterator Hole_const_iterator;
  typedef typename Output_arrangement_2::Halfedge_const_handle
                                                          Halfedge_const_handle;
  typedef typename Output_arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;
  typedef typename Geometry_traits_2::Point_2               Point_2;
  typedef typename Geometry_traits_2::Segment_2             Segment_2;

  Input_arrangement_2 arr;
  create_arrangement_from_env_file<Input_arrangement_2>(arr, input);
  std::cout << "Input arrangement has: "
            << GREEN << arr.number_of_faces()-1 << RESET
            << " faces." << std::endl;
  if (Visibility_2::Supports_general_polygon_tag::value) {
    int cnt(1);
    Face_const_iterator fit;
    for (fit = arr.faces_begin() ; fit != arr.faces_end() ; fit++) {
      if (!fit->is_unbounded()) {
        std::cout << "Test star-shape with face "
                  << GREEN << cnt << RESET << " ..." << std::endl;
        test_star_shape_one_face<Visibility_2>(  arr,
                                                 choice,
                                                 fit,
                                                 visibility);
      }
    }
  }
  else {  // Only run the test_star_shape_one_face() on the outer loop of the arrangement
    Face_const_iterator fit;
    // See which face has holes
    int cnt(1);
    for (fit = arr.faces_begin() ; fit != arr.faces_end() ; fit++) {
      if (!fit->is_unbounded()) {
        std::cout << "Test star-shape with face "
                  << GREEN << cnt << RESET << " ..." << std::endl;
        Hole_const_iterator hit;
        bool has_holes = false;
        for (hit = fit->holes_begin() ; hit != fit->holes_end() ; hit++) {
          has_holes = true;
          break;
        }
        if (has_holes && fit->has_outer_ccb()) {
          Output_arrangement_2 arr_trimmed;
          std::vector<Segment_2> segments;
          Ccb_halfedge_const_circulator circ = fit->outer_ccb();
          Ccb_halfedge_const_circulator curr = circ;
          do {
            Halfedge_const_handle he = curr;
            segments.push_back(Segment_2(he->source()->point(), he->target()->point()));
          } while (++curr != circ);
          CGAL::insert(arr_trimmed, segments.begin(), segments.end());
          Face_const_handle fch;
          if (arr_trimmed.faces_begin()->is_unbounded()) {
            fch = ++arr_trimmed.faces_begin();
          }
          else {
            fch = arr_trimmed.faces_begin();
          }
          test_star_shape_one_face<Visibility_2>(arr_trimmed,
                                                 choice,
                                                 fch,
                                                 visibility
                                                 );
          //CGAL::Visibility_2::print_arrangement(arr_trimmed);
        }
        else if (!has_holes) {
          test_star_shape_one_face<Visibility_2>(arr,
                                                 choice,
                                                 fit,
                                                 visibility);
        }
        cnt++;
      }
    }
  }
}

} // end namespace CGAL

#endif 
