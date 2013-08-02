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

#include <cassert>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <CGAL/Gmpq.h>
#include <CGAL/Timer.h>
#include <CGAL/Arr_naive_point_location.h>

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
void lazy_create_arrangement_from_file(std::ifstream &input,
                                       _Arrangement_2 &arr) {

  typedef _Arrangement_2                                      Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2           Geometry_traits_2;
  typedef typename Arrangement_2::Face_const_handle           Face_const_handle;
  typedef typename Geometry_traits_2::Segment_2               Segment_2;
  typedef typename Geometry_traits_2::Point_2                 Point_2;
  typedef typename Geometry_traits_2::FT                      Number_type;

  if (input) {
    std::string curr_line;
    std::vector<Point_2> isolated_vertices;
    std::stringstream convert(curr_line);
    int number_of_isolated_vertices;
    convert >> number_of_isolated_vertices;
    Face_const_handle uface = arr.unbounded_face();

    for (int i = 0 ; i < number_of_isolated_vertices ; i++) {
      std::getline(input, curr_line);
      std::istringstream iss(curr_line);
      std::string x, y;
      iss >> x >> y;
      arr.insert_in_face_interior(Point_2(string2num<Number_type>(x), 
                                          string2num<Number_type>(y)));
    }

    std::vector<Segment_2> edges;
    int number_of_edges;
    input >> number_of_edges;
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
  }
}

template <class _Visibility_2> 
void run_tests(_Visibility_2 visibility, int case_number) {
  typedef _Visibility_2                                       Visibility_2;
  typedef typename Visibility_2::Input_arrangement_2          Input_arrangement_2;
  typedef typename Visibility_2::Output_arrangement_2         Output_arrangement_2;
  typedef typename Input_arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Geometry_traits_2::Point_2                 Point_2; 
  typedef typename Geometry_traits_2::FT                      Number_type;

  for (int i=1 ; i <= case_number ; i++) {

    std::cout<<"Test "<<i<<" begins"<<std::endl;
    std::string input_arr_file("data/test");
    input_arr_file += num2string<int>(i);
    std::ifstream input(input_arr_file.c_str());
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

    lazy_create_arrangement_from_file<Input_arrangement_2>(input, arr_in);
    lazy_create_arrangement_from_file<Output_arrangement_2>(input, arr_correct_out);
    typename Input_arrangement_2::Face_const_iterator fit;

    for (fit = arr_in.faces_begin(); fit != arr_in.faces_end(); ++fit) {
        if (!fit->is_unbounded()) {
            visibility.visibility_region(query_pt, fit, arr_out);
        }
    }
    assert(true == test_are_equal<Output_arrangement_2>(arr_out, arr_correct_out));
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
        points.push_back(Point_2(string2num<Number_type>(n1),   
                                 string2num<Number_type>(n2)));
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

enum QueryChoice {VERTEX, EDGE, FACE};

template<class Point_2, class Number_type> 
Point_2 random_linear_interpolation(const Point_2 &p, const Point_2 &q) {

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

  Number_type y = y0 + (y1 - y0)*(x - x0)/(x1 - x0);
  return Point_2(x, y);
}

template<class Arrangement_2> 
bool is_inside_face(const Arrangement_2 &arr, 
                    const typename Arrangement_2::Face_const_handle face,
                    const typename Arrangement_2::Geometry_traits_2::Point_2 &q) {

  typedef CGAL::Arr_naive_point_location<Arrangement_2> Naive_pl;
  Naive_pl naive_pl(arr);

  // Perform the point-location query.
  CGAL::Object obj = naive_pl.locate(q);

  // Print the result.
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

template<class Visibility_2_fst, class Visibility_2_snd, class Arrangement_2>
void benchmark(Visibility_2_fst &visibility_fst, 
                             Visibility_2_snd &visibility_snd,
                             const Arrangement_2 &arr, 
                             typename Arrangement_2::Face_const_handle face,
                             QueryChoice choice) {

  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;
  typedef typename Geometry_traits_2::Point_2           Point_2;
  typedef typename Geometry_traits_2::FT                Number_type;
  typedef Timer Benchmark_timer;

  std::cout << "here" << std::endl;
  Benchmark_timer timer;
  timer.start();
  visibility_fst.attach(arr);
  timer.stop();
  std::cout << "Time to attach to first object: " << timer.time() << std::endl;

  timer.reset();
  timer.start();
  visibility_snd.attach(arr);
  timer.stop();
  std::cout << "Time to attach to second object: " << timer.time() << std::endl;

  Ccb_halfedge_const_circulator circ = face->outer_ccb();
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
        Ccb_halfedge_const_circulator curr_next = circ;
        Halfedge_const_handle he_next = curr_next;
        Point_2 p1 = he->source()->point();
        Point_2 p2 = he->target()->point();
        Point_2 p3 = he_next->target()->point();
        Point_2 avg((p1.x() + p2.x() + p3.x())/3, (p1.y() + p2.y() + p3.y())/3);
        if (is_inside_face<Arrangement_2>(arr, face, avg)) {
          std::cout << "Inside face" << std::endl;
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
    std::cout << "Running with qpoint: " << curr_query_pt << std::endl;
    Arrangement_2 out_arr_fst, out_arr_snd;
    timer.reset();
    timer.start();

    if (choice == FACE) {
      visibility_fst.visibility_region(curr_query_pt, face, out_arr_fst);
    }
    else {
 //     visibility_fst.visibility_region(curr_query_pt, he, out_arr_fst);
    }
    timer.stop();

    std::cout << "Time to compute visibility region using first object for " 
              << curr_query_pt << " : " << timer.time() << std::endl;
    timer.reset();

    timer.start();
    if (choice == FACE) {
      visibility_snd.visibility_region(curr_query_pt, face, out_arr_snd);
    }
    else {
//      visibility_snd.visibility_region(curr_query_pt, he, out_arr_snd); 
    }
    timer.stop();
    
    std::cout << "Time to compute visibility region using second object for " 
              << curr_query_pt << " : " << timer.time() << std::endl;
    visibility_fst.print_arrangement(out_arr_fst);
    visibility_fst.print_arrangement(out_arr_snd);
    assert(true == (CGAL::test_are_equal<Arrangement_2>
                          (out_arr_fst, out_arr_snd)));  
  } while (++curr != circ);
}

} // end namespace CGAL

#endif 
