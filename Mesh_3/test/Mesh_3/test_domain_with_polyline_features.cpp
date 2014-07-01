// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include "test_utilities.h"

#include <vector>
#include <list>
#include <utility>

template <typename K>
struct Dummy_domain
{
  typedef K R;
  typedef typename K::Point_3 Point_3;
  typedef typename K::FT FT;
  typedef int Index;
  typedef int Surface_patch_index;
};

typedef Dummy_domain<K_e_i> Smooth_domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Smooth_domain> Mesh_domain;
typedef Mesh_domain::Point_3 Point;
typedef Mesh_domain::FT FT;


class Domain_with_polyline_tester
{
  typedef std::vector<Point>  Polyline;
  typedef std::list<Polyline> Polylines;
  
  typedef Mesh_domain::Corner_index         Ci;
  typedef Mesh_domain::Curve_segment_index  Csi;
  typedef Mesh_domain::Index                Index;
  
  typedef std::vector<std::pair<Ci, Point> >        Corners_vector;
  typedef std::pair<Point, Index>                   P_and_i;
  typedef CGAL::cpp11::tuple<Csi,P_and_i,P_and_i>   Curve_tuple;
  typedef std::vector<Curve_tuple>                 Curves_vector;
  
public:
  Domain_with_polyline_tester()
    : p1_(1,0,0), p2_(1,1,0), p3_(1,2,0.1), p4_(0.9, 0.9, 1)
  { }
  
  void build_curve_segment()
  {
    Polylines polylines (1);
    Polyline& polyline = polylines.front();
    
    polyline.push_back(p1_);
    polyline.push_back(p2_);
    polyline.push_back(p3_);
    
    domain_.add_features(polylines.begin(),polylines.end());    
  }
  
  void build_cycle()
  {
    Polylines polylines (1);
    Polyline& polyline = polylines.front();
    
    polyline.push_back(p1_);
    polyline.push_back(p2_);
    polyline.push_back(p3_);
    polyline.push_back(p4_);
    polyline.push_back(p1_);
    
    domain_.add_features(polylines.begin(),polylines.end());    
  }
  
  void test_curve_segment_corners() const
  {
    Corners_vector corners;
    domain_.get_corners(std::back_inserter(corners));
    
    assert(corners.size() == 2);
    assert(   (corners.front().second == p1_ && corners.back().second == p3_)
           || (corners.front().second == p3_ && corners.back().second == p1_) );
  }
  
  void test_cycle_corners() const
  {
    Corners_vector corners;
    domain_.get_corners(std::back_inserter(corners));
    
    assert(corners.empty());
  }
  
  void test_curve_segments() const
  {
    std::pair<Point,Point> extremities = get_extremities();
    
    assert (   ( extremities.first == p1_ && extremities.second == p3_ )
            || ( extremities.second == p3_ && extremities.first == p1_ ) ); 
  }
  
  void test_cycles() const
  {
    std::pair<Point,Point> extremities = get_extremities();
    
    assert (   extremities.first == extremities.second
            && extremities.first == p1_ );
  }
  
  void test_geodesic_distance() const
  {
    std::pair<Point,Point> extremities = get_extremities();
    const Point& p = extremities.first;
    const Point& q = extremities.second;
    const Csi& curve_index = get_curve_segment_index();
    
    const FT geod (1.3);
    Point r = domain_.construct_point_on_curve_segment(p,curve_index,geod);
    const FT& pq_geo = domain_.geodesic_distance(p,q,curve_index);
    
    assert(CGAL::squared_distance(p,r) < CGAL::square(geod));
    assert(near_equal(domain_.geodesic_distance(p,r,curve_index), geod));
    assert(near_equal(domain_.geodesic_distance(r,q,curve_index), pq_geo - geod));
    assert(near_equal(domain_.geodesic_distance(q,r,curve_index), geod - pq_geo));
  }
  
private:
  std::pair<Point,Point> get_extremities() const
  {
    Curves_vector curve_segments;
    domain_.get_curve_segments(std::back_inserter(curve_segments));
    assert(curve_segments.size() == 1);
    
    const Point& p = get_first_point(curve_segments.front());
    const Point& q = get_second_point(curve_segments.front());
    
    return std::make_pair(p,q);
  }
  
  Point get_first_point(const Curve_tuple& tuple) const
  {
    return CGAL::cpp11::get<1>(tuple).first;
  }
  
  Point get_second_point(const Curve_tuple& tuple) const
  {
    return CGAL::cpp11::get<2>(tuple).first;
  }
  
  Csi get_curve_segment_index() const
  {
    Curves_vector curve_segments;
    domain_.get_curve_segments(std::back_inserter(curve_segments));
    assert(curve_segments.size() == 1);
    
    return get_curve_segment_index(curve_segments.front());
  }
  
  Csi get_curve_segment_index(const Curve_tuple& tuple) const
  {
    return CGAL::cpp11::get<0>(tuple);
  }
  
  bool near_equal(const double d1, const double d2) const
  {
    const double epsilon = 1e-8;
    return ( d1-d2 < epsilon || d2-d1 < epsilon );
  }

private:
  Mesh_domain domain_;
  Point p1_, p2_, p3_, p4_;
};

int main()
{
  std::cout << "Test curve segments" << std::endl;
  Domain_with_polyline_tester domain_tester;
  domain_tester.build_curve_segment();
  
  domain_tester.test_curve_segment_corners();
  domain_tester.test_curve_segments();
  domain_tester.test_geodesic_distance();
  
  std::cout << "Test cycles" << std::endl;
  Domain_with_polyline_tester domain_cycle_tester;
  domain_cycle_tester.build_cycle();
  
  domain_cycle_tester.test_cycle_corners();
  domain_cycle_tester.test_cycles();
  domain_cycle_tester.test_geodesic_distance();
  
  return EXIT_SUCCESS;
}
