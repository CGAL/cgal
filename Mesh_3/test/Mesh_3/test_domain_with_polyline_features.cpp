// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
  typedef unsigned short Subdomain_index;
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
  typedef Mesh_domain::Curve_index          Csi;
  typedef Mesh_domain::Surface_patch_index  Spi;
  typedef Mesh_domain::Index                Index;

  typedef std::vector<std::pair<Ci, Point> >        Corners_vector;
  typedef std::pair<Point, Index>                   P_and_i;
  typedef std::tuple<Csi,P_and_i,P_and_i>   Curve_tuple;
  typedef std::vector<Curve_tuple>                  Curves_vector;

public:
  Domain_with_polyline_tester()
    : p1_(1,0,0), p2_(1,1,0), p3_(1,2,0.1), p4_(0.9, 0.9, 1)
  { }

  void build_corners()
  {
    domain_.add_corner(p1_);

    std::vector<Point> corners;
    corners.push_back(p2_);
    domain_.add_corners(corners.begin(), corners.end());

    Csi dummy_curve_index = 12;
    domain_.register_corner(p3_, dummy_curve_index);

    Spi dummy_surface_patch_index = 21;
    domain_.add_corner_with_context(p4_, dummy_surface_patch_index);
  }

  void build_curve()
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

  void test_corners() const
  {
    Corners_vector corners;
    domain_.get_corners(std::back_inserter(corners));
    assert(corners.size() == 4);

    Corners_vector::const_iterator cit = corners.begin(), end = corners.end();
    for(; cit!=end; ++cit)
    {
      if(cit->second == p3_)
      {
        std::vector<Csi> incident_curves;
        domain_.get_corner_incident_curves(cit->first, std::back_inserter(incident_curves));
        assert(incident_curves.size() == 1 && incident_curves.front() == 12);
      }
      else if(cit->second == p4_)
      {
        std::vector<Spi> incident_surface_patchs;
        domain_.get_corner_incidences(cit->first, std::back_inserter(incident_surface_patchs));
        assert(incident_surface_patchs.size() == 1 && incident_surface_patchs.front() == 21);
      }
      else
      {
        assert(cit->second == p1_ || cit->second == p2_);
      }
    }
  }

  void test_curve_corners() const
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

    assert(corners.size() == 1);
  }

  void test_curves() const
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
    const Csi& curve_index = get_curve_index();

    const FT geod (1.3);
    Point r = domain_.construct_point_on_curve(p,curve_index,geod);
    const FT& pq_geo = domain_.signed_geodesic_distance(p,q,curve_index);

    assert(CGAL::squared_distance(p,r) < CGAL::square(geod));
    assert(near_equal(domain_.signed_geodesic_distance(p,r,curve_index), geod));
    assert(near_equal(domain_.signed_geodesic_distance(r,q,curve_index), pq_geo - geod));
    assert(near_equal(domain_.signed_geodesic_distance(q,r,curve_index), geod - pq_geo));
  }

private:
  std::pair<Point,Point> get_extremities() const
  {
    Curves_vector curves;
    domain_.get_curves(std::back_inserter(curves));
    assert(curves.size() == 1);

    const Point& p = get_first_point(curves.front());
    const Point& q = get_second_point(curves.front());

    return std::make_pair(p,q);
  }

  Point get_first_point(const Curve_tuple& tuple) const
  {
    return std::get<1>(tuple).first;
  }

  Point get_second_point(const Curve_tuple& tuple) const
  {
    return std::get<2>(tuple).first;
  }

  Csi get_curve_index() const
  {
    Curves_vector curves;
    domain_.get_curves(std::back_inserter(curves));
    assert(curves.size() == 1);

    return get_curve_index(curves.front());
  }

  Csi get_curve_index(const Curve_tuple& tuple) const
  {
    return std::get<0>(tuple);
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
  std::cout << "Test corners" << std::endl;
  Domain_with_polyline_tester domain_corner_tester;
  domain_corner_tester.build_corners();
  domain_corner_tester.test_corners();

  std::cout << "Test curve segments" << std::endl;
  Domain_with_polyline_tester domain_tester;
  domain_tester.build_curve();

  domain_tester.test_curve_corners();
  domain_tester.test_curves();
  domain_tester.test_geodesic_distance();

  std::cout << "Test cycles" << std::endl;
  Domain_with_polyline_tester domain_cycle_tester;
  domain_cycle_tester.build_cycle();

  domain_cycle_tester.test_cycle_corners();
  domain_cycle_tester.test_cycles();
  domain_cycle_tester.test_geodesic_distance();

  return EXIT_SUCCESS;
}
