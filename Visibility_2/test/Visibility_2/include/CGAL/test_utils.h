// Copyright (c) 2013 Technical University Braunschweig (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
#include <CGAL/Timer.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Visibility_2/visibility_utils.h>

namespace CGAL {

enum Query_choice {VERTEX, EDGE, FACE};

template <class Arrangement_2>
typename Arrangement_2::Halfedge_handle get_initial_halfedge(const Arrangement_2 &arr) {

  typedef typename Arrangement_2::Vertex Vertex;
  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;
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
template <class ARR1, class ARR2>
bool test_are_equal(const ARR1 &arr1, const ARR2 &arr2) {
  typedef typename ARR1::Halfedge_handle       HE1;
  typedef typename ARR2::Halfedge_handle       HE2;


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
  HE1 he_start_1 = get_initial_halfedge(arr1);
  HE2 he_start_2 = get_initial_halfedge(arr2);

  // run on first loop and compare sources
  assert(arr1.traits()->compare_xy_2_object()(
             he_start_1->source()->point(),
             he_start_2->source()->point()) == CGAL::EQUAL);

  HE1 he_run_1 = he_start_1->next();
  HE2 he_run_2 = he_start_2->next();

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

  if (input) {
    std::string curr_line;
    std::getline(input, curr_line);
    std::stringstream convert(curr_line);
    int number_of_isolated_vertices;
    convert >> number_of_isolated_vertices;
    Face_handle uface = arr.unbounded_face();
    for (int i = 0 ; i < number_of_isolated_vertices ; i++) {
      std::getline(input, curr_line);
      std::istringstream iss(curr_line);
      Point_2 p;
      iss >> p;
      arr.insert_in_face_interior(p, uface);
    }
    std::vector<Segment_2> edges;
    int number_of_edges;
    std::getline(input, curr_line);
    std::stringstream convert2(curr_line);
    convert2 >> number_of_edges;
    for (int i = 0 ; i < number_of_edges ; i++) {
      std::getline(input, curr_line);
      std::istringstream iss(curr_line);
      Segment_2 seg;
      iss >> seg;
      edges.push_back(seg);
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


template <class Visibility_2, class Visibility_arrangement_2>
bool run_test_case_from_file(Visibility_2& visibility, std::ifstream &input) {
  typedef typename Visibility_2::Arrangement_2          Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Geometry_traits_2::Point_2                 Point_2;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
                                        Halfedge_around_vertex_const_circulator;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;

  Arrangement_2 arr_in;
  Visibility_arrangement_2 arr_correct_out;
  Visibility_arrangement_2 arr_out;

  std::string curr_line;
  while (std::getline(input, curr_line)) {
    if (curr_line[0] != '#' && curr_line[0] != '/')
      break;
  }
  std::stringstream ss(curr_line);
  Point_2 query_pt;
  ss >> query_pt;

  std::getline(input, curr_line);
  std::istringstream iss(curr_line);
  Point_2 reference_pt;
  iss >> reference_pt;

  std::getline(input, curr_line);
  if (!create_arrangement_from_dat_file<Arrangement_2>(input, arr_in)) {
    return false;
  }
  visibility.attach(arr_in);
  std::getline(input, curr_line);
  if (!create_arrangement_from_dat_file<Visibility_arrangement_2>
                                                   (input, arr_correct_out)) {
    return false;
  }

  if(Visibility_2::Regularization_category::value){
    regularize(arr_correct_out);
  }


  CGAL::Object obj = CGAL::get_location<Arrangement_2>(arr_in, query_pt);
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


  if (!test_are_equal<Visibility_arrangement_2, Visibility_arrangement_2>(arr_out, arr_correct_out)) {
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
void test_interface() {
    typedef typename Visibility_2::Arrangement_2              Arrangement_2;
    typedef typename Arrangement_2::Geometry_traits_2         Geometry_traits_2;
    typedef typename Arrangement_2::Face_const_handle         Face_const_handle;
    typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Geometry_traits_2::Segment_2             Segment_2;
    typedef typename Geometry_traits_2::Point_2               Point_2;

    Point_2 query(0.1, 0.1);

    std::vector<Point_2> vertices;
    vertices.push_back(Point_2(0, 0));
    vertices.push_back(Point_2(1, 0));
    vertices.push_back(Point_2(1, 1));
    vertices.push_back(Point_2(0, 1));


    Arrangement_2 arr_out;
    Arrangement_2 square;

    for(unsigned int i = 0; i < vertices.size(); ++i) {
        CGAL::insert(square, Segment_2(vertices[i],
                                       vertices[(i+1) % vertices.size()]));
    }

    Face_const_handle location;
    CGAL::assign(location, get_location(square, query));



    const Arrangement_2& arr = square;

    // Constructor and attach method must accept a const arrangement.
    Visibility_2 visibility(arr);

    visibility.detach();
    visibility.attach(arr);


    const Visibility_2& vis = visibility;

    // compute_visibility must be const
    vis.compute_visibility(query, location, arr_out);

    Halfedge_const_handle he = arr.edges_begin();

    if(he->face()->is_unbounded())
        he = he->twin();

    vis.compute_visibility(he->target()->point(), he, arr_out);

    // must have const arrangement_2();
    const Arrangement_2& a = vis.arrangement_2();

    visibility.attach(a);

    // must have const is_attached();
    vis.is_attached();




}

template <class Visibility_2>
void run_tests_with_changes_to_arr() {

  std::cout << "    Testing changes to attached arrangement:";

  typedef typename Visibility_2::Arrangement_2                Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2           Geometry_traits_2;
  typedef typename Arrangement_2::Face_const_handle           Face_const_handle;
  typedef typename Geometry_traits_2::Segment_2               Segment_2;
  typedef typename Geometry_traits_2::Point_2                 Point_2;

  bool all_passed = true;


  Point_2 query(0.1, 0.1);

  std::vector<Point_2> vertices;
  vertices.push_back(Point_2(0, 0));
  vertices.push_back(Point_2(1, 0));
  vertices.push_back(Point_2(1, 1));
  vertices.push_back(Point_2(0, 1));


  Arrangement_2 square;
  for(unsigned int i = 0; i < vertices.size(); ++i) {
      CGAL::insert(square, Segment_2(vertices[i],
                                     vertices[(i+1) % vertices.size()]));
  }

  Arrangement_2 lower_tri;
  CGAL::insert(lower_tri, Segment_2(vertices[0], vertices[1]));
  CGAL::insert(lower_tri, Segment_2(vertices[1], vertices[3]));
  CGAL::insert(lower_tri, Segment_2(vertices[3], vertices[0]));

  Visibility_2 visibility;
  Arrangement_2 arr;
  Arrangement_2 arr_out;

  // Attach empty arr and fill it afterwards
  visibility.attach(arr);


  for(unsigned int i = 0; i < vertices.size(); ++i) {
      CGAL::insert(arr, Segment_2(vertices[i],
                                  vertices[(i+1) % vertices.size()]));
  }

  Face_const_handle location;
  CGAL::assign(location, get_location(arr, query));

  visibility.compute_visibility(query, location, arr_out);

  all_passed &= test_are_equal(arr_out, square);


  // Change attached arrangement and query again

  //arr.clear();

  //CGAL::insert(arr, Segment_2(vertices[0], vertices[1]));
  CGAL::insert(arr, Segment_2(vertices[1], vertices[3]));
  //CGAL::insert(arr, Segment_2(vertices[3], vertices[0]));

  CGAL::assign(location, get_location(arr, query));

  visibility.compute_visibility(query, location, arr_out);

  all_passed &= test_are_equal(arr_out, lower_tri);


  // Detach and attach again

  visibility.detach();
  visibility.attach(arr);

  CGAL::assign(location, get_location(arr, query));

  visibility.compute_visibility(query, location, arr_out);

  all_passed &= test_are_equal(arr_out, lower_tri);


  // Attach another arrangement without detaching the old one first.

  visibility.attach(square);

  CGAL::assign(location, get_location(square, query));

  visibility.compute_visibility(query, location, arr_out);

  all_passed &= test_are_equal(arr_out, square);


  if (!all_passed) {
    std::cout << "\tFailed: Modifying attached arrangement causes wrong output.\n";
    assert(false);
  } else {
    std::cout << "\tPassed.\n" ;
  }
}


template <class Visibility_2, class Visibility_arrangement_2>
void run_tests(int case_number_simple, int case_number_non_simple) {

  // Make sure the code only compiles with a conforming interface
  test_interface<Visibility_2>();


  Visibility_2 visibility;
  bool one_failed = false;
  if (Visibility_2::Supports_simple_polygon_category::value
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
      if (run_test_case_from_file<Visibility_2,Visibility_arrangement_2>(visibility, input)) {
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
  if (Visibility_2::Supports_general_polygon_category::value
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
      if (run_test_case_from_file<Visibility_2,Visibility_arrangement_2>(visibility, input)) {
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

  run_tests_with_changes_to_arr<Visibility_2>();

}


template <class _Arrangement_2>
void create_arrangement_from_file(_Arrangement_2 &arr, std::ifstream& input) {
  typedef _Arrangement_2                                                                   Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2          Geometry_traits_2;
  typedef typename Geometry_traits_2::Segment_2         Segment_2;
  typedef typename Geometry_traits_2::Point_2                  Point_2;

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
      std::istringstream iss(line);
      Point_2 p;
      iss >> p;
      points.push_back(p);
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
  typedef _Arrangement_2                                                                                 Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2            Geometry_traits_2;
  typedef typename Geometry_traits_2::Segment_2         Segment_2;
  typedef typename Geometry_traits_2::Point_2                  Point_2;

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
      typename std::vector<Point_2>::size_type number_of_vertices;
      convert >> number_of_vertices;
      for (typename std::vector<Point_2>::size_type j = 0;
           j < number_of_vertices; j++) {
        std::getline(input, line);
        std::istringstream iss(line);
        Point_2 p;
        iss >> p;
        points.push_back(p);
      }
      for (typename std::vector<Point_2>::size_type j = 0;
           j < number_of_vertices-1 ; j++)
      {
        segments.push_back(Segment_2(points[j], points[j+1]));
      }
      segments.push_back(Segment_2(points.front(), points.back()));
  //    std::cout << "before insertion\n";
      CGAL::insert(arr, segments.begin(), segments.end());
 //     std::cout << "after\n";
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
template<class Arrangement_2>
bool is_star_shape(
    const typename Arrangement_2::Point_2& q,
    const Arrangement_2& arr) {

  typedef typename Arrangement_2::Face_const_handle   Face_const_handle;

  // this test is written for an arr that contains on star shaped polygon
  if (arr.number_of_faces()!=2){
    return false;
  }

  // get the bounded face
  Face_const_handle fh;
  if(arr.faces_begin()->is_unbounded()){
    fh = arr.faces_begin();
  }else{
    fh = ++(arr.faces_begin());
  }
  assert(fh->is_unbounded());

  if (fh->has_outer_ccb()) {
    typename Arrangement_2::Ccb_halfedge_const_circulator curr, circ;
    curr = circ = fh->outer_ccb();
    do {
      if (CGAL::right_turn(q, curr->source()->point(), curr->target()->point())) {
        return false;
      }
//      Point_2 p = curr->target()->point();
//      typename Arrangement_2::Ccb_halfedge_const_circulator curr1, circ1;
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

template <class Arrangement_2>
int count_edges_in_face(typename Arrangement_2::Face_const_handle &fch) {
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;

  Ccb_halfedge_const_circulator circ = fch->outer_ccb();
  Ccb_halfedge_const_circulator curr = circ;

  int edge_cnt(0);
  do {
    edge_cnt++;
  } while (++curr != circ);
  return edge_cnt;
}

template <class Arrangement_2>
typename Arrangement_2::Face_const_handle construct_biggest_arr_with_no_holes(
                                                    Arrangement_2 &arr_in,
                                                    Arrangement_2 &arr_out) {

  typedef typename Arrangement_2::Face_const_iterator   Face_const_iterator;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Geometry_traits_2::Segment_2         Segment_2;

  int curr_max(0);
  Ccb_halfedge_const_circulator curr_max_circ;
  Ccb_halfedge_const_circulator circ;
  Ccb_halfedge_const_circulator curr;
  Face_const_iterator fit;

  for (fit = arr_in.faces_begin() ; fit != arr_in.faces_end() ; fit++) {
    if (fit->has_outer_ccb()) {
      Ccb_halfedge_const_circulator circ = fit->outer_ccb();
      int edge_cnt = count_edges_in_face<Arrangement_2>(fit);
      if (edge_cnt > curr_max) {
        curr_max = edge_cnt;
        curr_max_circ = circ;
      }
    }
  }

  std::vector<Segment_2> segments;
  curr = curr_max_circ;
  Halfedge_const_handle he;
  do {
    he = curr;
    segments.push_back(Segment_2(he->source()->point(), he->target()->point()));
  } while (++curr != curr_max_circ);

  arr_out.clear();
//  std::cout << "before insertion\n";
  CGAL::insert_non_intersecting_curves(arr_out, segments.begin(), segments.end());
//  std::cout << "after\n";

  Face_const_handle fch;
  curr_max = 0;
  for (fit = arr_out.faces_begin() ; fit != arr_out.faces_end() ; fit++) {
    if (fit->has_outer_ccb()) {
      int edge_cnt = count_edges_in_face<Arrangement_2>(fit);
      if (edge_cnt > curr_max) {
        curr_max = edge_cnt;
        fch = fit;
      }
    }
  }

  return fch;
}

template <class Visibility_2_fst, class Visibility_2_snd>
void simple_benchmark_one_unit(
          typename Visibility_2_fst::Arrangement_2 &arr,
          const Query_choice &choice,
          typename Visibility_2_fst::Arrangement_2::Face_const_handle &fit,
          Visibility_2_fst& visibility_fst,
          Visibility_2_snd& visibility_snd,
          double& qtime1,
          double& qtime2,
          int& query_cnt) {

  typedef typename Visibility_2_fst::Arrangement_2      Arrangement_2;
  typedef typename Visibility_2_fst::Arrangement_2      Visibility_arrangement_2;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;

  typedef typename Visibility_arrangement_2::Face_handle        Face_handle;
  typedef typename Geometry_traits_2::Point_2               Point_2;
  typedef typename Geometry_traits_2::FT                    Number_type;
  typedef Timer Benchmark_timer;

  Benchmark_timer timer;

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
        if (is_inside_face<Arrangement_2>(arr, fit, avg)) {
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

    Visibility_arrangement_2 out_arr_fst, out_arr_snd;
    timer.reset();
    timer.start();
    Face_handle f_fst;
    if (choice == CGAL::FACE) {
      f_fst = visibility_fst.compute_visibility(curr_query_pt, fit, out_arr_fst);
      query_cnt++;
    }
    else {
      f_fst = visibility_fst.compute_visibility(curr_query_pt, he, out_arr_fst);
      query_cnt++;
    }
    timer.stop();
    qtime1 += timer.time();

    timer.reset();
    if ( !is_star_shape(curr_query_pt, out_arr_fst) ) {
      std::cout << RED << "         Warning: the first output is not star-shape." << RESET << std::endl;
    }
    timer.start();
    Face_handle f_snd;
    if (choice == CGAL::FACE) {
      f_snd = visibility_snd.compute_visibility(curr_query_pt, fit, out_arr_snd);
    }
    else {
      f_snd = visibility_snd.compute_visibility(curr_query_pt, he, out_arr_snd);
    }
    timer.stop();
    if ( !is_star_shape(curr_query_pt, out_arr_snd) ) {
      std::cout << RED << "         Warning: the second output is not star-shape." << RESET << std::endl;
    }
    qtime2 += timer.time();

    if (! CGAL::test_are_equal<Visibility_arrangement_2> (out_arr_fst, out_arr_snd)) {
      if (choice == CGAL::FACE) {
        std::cout << "query in Face:\n";
        CGAL::Visibility_2::print_simple_face<Face_const_handle, Ccb_halfedge_const_circulator>(fit);
      }
      std::cout << "the query point is " << curr_query_pt << std::endl;
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
void simple_benchmark(Visibility_2_fst &visibility_fst,
                      Visibility_2_snd &visibility_snd,
                      const Query_choice &choice,
                      std::ifstream &input) {

  typedef typename Visibility_2_fst::Arrangement_2    Arrangement_2;

  typedef typename Arrangement_2::Face_const_iterator Face_const_iterator;
  typedef typename Arrangement_2::Face_const_handle   Face_const_handle;

  assert(Visibility_2_fst::Regularization_category::value
      == Visibility_2_snd::Regularization_category::value);

  Arrangement_2 arr;
  create_arrangement_from_env_file<Arrangement_2>(arr, input);

  int query_cnt(0);
  double qtime1(0), qtime2(0), ptime1(0), ptime2(0);
  if (Visibility_2_fst::Supports_general_polygon_category::value
    && Visibility_2_snd::Supports_general_polygon_category::value) {

    Face_const_iterator fit;
    Timer timer;

    timer.start();
    visibility_fst.attach(arr);
    timer.stop();
    ptime1 = timer.time();

    timer.reset();
    timer.start();
    visibility_snd.attach(arr);
    timer.stop();
    ptime2 = timer.time();

    for (fit = arr.faces_begin() ; fit != arr.faces_end() ; fit++) {
      if (!fit->is_unbounded()) {

        simple_benchmark_one_unit<Visibility_2_fst, Visibility_2_snd>(arr,
                                                               choice,
                                                               fit,
                                                               visibility_fst,
                                                               visibility_snd,
                                                               qtime1,
                                                               qtime2,
                                                               query_cnt);
      }
    }
  }
  else {
    Arrangement_2 arr_trimmed;
    Face_const_handle fch = construct_biggest_arr_with_no_holes
                                        <Arrangement_2>(arr, arr_trimmed);
    Timer timer;

    timer.start();
    visibility_fst.attach(arr_trimmed);
    timer.stop();
    ptime1 = timer.time();

    timer.reset();
    timer.start();
    visibility_snd.attach(arr_trimmed);
    timer.stop();
    ptime2 = timer.time();
    simple_benchmark_one_unit<Visibility_2_fst, Visibility_2_snd>(arr_trimmed,
                                                                 choice,
                                                                 fch,
                                                                 visibility_fst,
                                                                 visibility_snd,
                                                                 qtime1,
                                                                 qtime2,
                                                                 query_cnt);
  }
  std::cout << "Preprocessing: "  << std::endl
            << "Model 1 uses " << ptime1 << "  sec" << std::endl
            << "Model 2 uses " << ptime2 << "  sec" << std::endl;
  std::cout << query_cnt << " queries are done.\n"
            << "Model 1 uses " << qtime1 << "  sec" << std::endl
            << "Model 2 uses " << qtime2 << "  sec" << std::endl;
  std::cout << "total times are:" << std::endl
            << "Model 1 uses " << ptime1 + qtime1 << "  sec" << std::endl
            << "Model 2 uses " << ptime2 + qtime2 << "  sec" << std::endl;
}

template<class Visibility_2>
void pure_benchmark_one_unit(
          typename Visibility_2::Arrangement_2 &arr,
          const Query_choice &choice,
          typename Visibility_2::Arrangement_2::Face_const_handle &fit,
          Visibility_2& visibility,
          double& qtime,
          int& query_cnt) {

  typedef typename Visibility_2::Arrangement_2                    Arrangement_2;
  typedef typename Visibility_2::Arrangement_2         Visibility_arrangement_2;
  typedef typename Arrangement_2::Halfedge_const_handle
                                                          Halfedge_const_handle;
  typedef typename Arrangement_2::Geometry_traits_2           Geometry_traits_2;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;

  typedef typename Visibility_arrangement_2::Face_handle            Face_handle;
  typedef typename Geometry_traits_2::Point_2                       Point_2;
  typedef typename Geometry_traits_2::FT                            Number_type;
  typedef Timer Benchmark_timer;

  Benchmark_timer timer;

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
        if (is_inside_face<Arrangement_2>(arr, fit, avg)) {
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

    Visibility_arrangement_2 out_arr;
    timer.reset();
    timer.start();
    Face_handle fh;
    if (choice == CGAL::FACE) {
      fh = visibility.compute_visibility(curr_query_pt, fit, out_arr);
      query_cnt++;
    }
    else {
      fh = visibility.compute_visibility(curr_query_pt, he, out_arr);
      query_cnt++;
    }
    timer.stop();
    qtime += timer.time();

    timer.reset();
  } while (++curr != circ);
}

template<class Visibility_2>
void pure_benchmark(  Visibility_2 &visibility,
                      const Query_choice &choice,
                      std::ifstream &input) {

  typedef typename Visibility_2::Arrangement_2        Arrangement_2;

  typedef typename Arrangement_2::Face_const_iterator Face_const_iterator;
  typedef typename Arrangement_2::Face_const_handle   Face_const_handle;

  Arrangement_2 arr;
  create_arrangement_from_env_file<Arrangement_2>(arr, input);

  int query_cnt(0);
  double qtime(0), ptime(0);
  if (Visibility_2::Supports_general_polygon_category::value) {

    Face_const_iterator fit;
    Timer timer;

    timer.start();
    visibility.attach(arr);
    timer.stop();
    ptime = timer.time();

    for (fit = arr.faces_begin() ; fit != arr.faces_end() ; fit++) {
      if (!fit->is_unbounded()) {

        pure_benchmark_one_unit<Visibility_2>( arr,
                                               choice,
                                               fit,
                                               visibility,
                                               qtime,
                                               query_cnt);
      }
    }
  }
  else {
    Arrangement_2 arr_trimmed;
    Face_const_handle fch = construct_biggest_arr_with_no_holes
                                        <Arrangement_2>(arr, arr_trimmed);
    Timer timer;

    timer.start();
    visibility.attach(arr_trimmed);
    timer.stop();
    ptime = timer.time();

    pure_benchmark_one_unit<Visibility_2>( arr_trimmed,
                                           choice,
                                           fch,
                                           visibility,
                                           qtime,
                                           query_cnt);
  }

  // std::cout << "NAME TAG  PreProTime NQueries TimeQueries TotalTime QAVE TAVE" << std::endl;
  std::cout << " " << visibility.name()
            << " " << Visibility_2::Regularization_category::value
            << "  " << ptime
            << " " << query_cnt
            << " " << qtime
            << " " << ptime+qtime
            << "  " << qtime/query_cnt
            << " " << (ptime+qtime)/query_cnt
            << " " << std::endl;
//   std::cout << "Preprocessing: "  << std::endl
//             << "cost " << ptime << "  sec" << std::endl;
//   std::cout << query_cnt << " queries are done.\n"
//             << "cost " << qtime << "  sec" << std::endl;
//   std::cout << "total time is:" << ptime + qtime << "  sec" << std::endl;
}


template<class Visibility_2>
void test_star_shape_one_face(  typename Visibility_2::Arrangement_2 &arr,
                                const Query_choice &choice,
                                typename Visibility_2::Arrangement_2::Face_const_handle &fit,
                                Visibility_2& visibility)
{

  typedef typename Visibility_2::Arrangement_2          Arrangement_2;
  typedef typename Visibility_2::Arrangement_2          Visibility_arrangement_2;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;


  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                            Ccb_halfedge_const_circulator;

  typedef typename Visibility_arrangement_2::Face_handle        Face_handle;
  typedef typename Geometry_traits_2::Point_2               Point_2;
  typedef typename Geometry_traits_2::FT                    Number_type;

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
        if (is_inside_face<Arrangement_2>(arr, fit, avg)) {
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
    Visibility_arrangement_2 out_arr;
    Face_handle fh;
    if (choice == FACE) {
      fh = visibility.compute_visibility(curr_query_pt, fit, out_arr);
    }
    else {
      fh = visibility.compute_visibility(curr_query_pt, he, out_arr);
    }
    assert(out_arr.number_of_faces()==2);
    assert(fh->is_unbounded());
    if ( !is_star_shape(curr_query_pt, out_arr)) {
      std::cout << RED << "     The face is not a star shape to qpoint." << RESET <<  std::endl;
    }
  } while (++curr != circ);
}

template<class Visibility_2>
void test_star_shape(Visibility_2 &visibility,
               const Query_choice &choice,
               std::ifstream &input) {

  typedef typename Visibility_2::Arrangement_2 Arrangement_2;
  typedef typename Visibility_2::Arrangement_2 Visibility_arrangement_2;
  typedef typename Visibility_arrangement_2::Halfedge_const_handle
                                                          Halfedge_const_handle;
  typedef typename Visibility_arrangement_2::Geometry_traits_2   Geometry_traits_2;
  typedef typename Visibility_arrangement_2::Face_const_iterator Face_const_iterator;
  typedef typename Visibility_arrangement_2::Face_const_handle   Face_const_handle;
  typedef typename Visibility_arrangement_2::Hole_const_iterator Hole_const_iterator;
  typedef typename Visibility_arrangement_2::Halfedge_const_handle
                                                          Halfedge_const_handle;
  typedef typename Visibility_arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;
  typedef typename Geometry_traits_2::Segment_2             Segment_2;

  Arrangement_2 arr;
  create_arrangement_from_env_file<Arrangement_2>(arr, input);
  std::cout << "Input arrangement has: "
            << GREEN << arr.number_of_faces()-1 << RESET
            << " faces." << std::endl;
  if (Visibility_2::Supports_general_polygon_category::value) {
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
          Visibility_arrangement_2 arr_trimmed;
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
