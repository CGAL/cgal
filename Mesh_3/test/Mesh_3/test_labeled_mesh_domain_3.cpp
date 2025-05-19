// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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

#include "test_meshing_utilities.h"
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>


template <typename K>
struct LM3_tester
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::FT FT;

  struct Function
  {
    typedef Point_3 Point;
    typedef FT (*Func)(const Point_3&);

    Function (Func f) : f_(f) {}
    FT operator() (const Point_3& p) const { return f_(p); }

  private:
    Func f_;
  };

  typedef CGAL::Implicit_to_labeling_function_wrapper<Function, K> Function_wrapper;
  typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

  static FT shape_function (const Point_3& p)
  {
    if (p.x() < 0)
      return -1;
    const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
        return x2 + y2 + z2 - 1;
  }

  static FT sphere_function (const Point_3& p)
  {
    const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
        return x2 + y2 + z2 - 1;
  }

  void operator() () const
  {
    typedef typename K::Sphere_3 Sphere_3;
    typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

    test_domain(Sphere_3(CGAL::ORIGIN, 4.));
    test_domain(CGAL::Bbox_3(-2.,-2.,-2., 2.,2.,2.));
    test_domain(Iso_cuboid_3(Point_3(-2.,-2.,-2.), Point_3(2.,2.,2.)));
  }

private:
  template <class BoundingShape>
  void test_domain (const BoundingShape& bounding_shape) const
  {
    FT error_bound(1e-3);

    Function f_sphere(&sphere_function);
    Function_wrapper wrapper_1(f_sphere);
    Mesh_domain domain(wrapper_1, bounding_shape, CGAL::parameters::relative_error_bound(error_bound));
    test_construct_initial_points(domain, error_bound);

    Function f_shape(&shape_function);
    Function_wrapper wrapper_2(f_shape);
    Mesh_domain domain_2(wrapper_2, bounding_shape, CGAL::parameters::relative_error_bound(error_bound));
    test_is_in_domain(domain_2);
    test_do_intersect_surface(domain_2);
    test_construct_intersection(domain_2);
  }

  void test_construct_initial_points (const Mesh_domain& domain, FT error_bound) const
  {
    typedef typename Mesh_domain::Construct_initial_points Construct_initial_points;
    typedef typename Mesh_domain::Index Index;
    typedef typename std::vector<std::pair<Point_3, Index> >::const_iterator Points_const_iterator;
    typedef typename K::Compute_squared_distance_3 Compute_squared_distance_3;

    Compute_squared_distance_3 squared_distance = K().compute_squared_distance_3_object();

    Construct_initial_points construct_initial_points = domain.construct_initial_points_object();
    std::vector<std::pair<Point_3, Index> > points;
    int point_count = 12;
    construct_initial_points(std::back_inserter(points), point_count);

    for (Points_const_iterator iter = points.begin(), end_iter = points.end(); iter != end_iter; ++iter)
    {
      const Point_3& p = iter->first;

      FT sd = squared_distance(CGAL::ORIGIN, p);
      FT diff = sd - 1;
      if (diff < FT(0.))
              diff = -diff;
      assert(diff <= error_bound);
    }
  }

  void test_is_in_domain (const Mesh_domain& domain) const
  {
    typedef typename Mesh_domain::Is_in_domain Is_in_domain;
    typedef typename Mesh_domain::Subdomain Subdomain;
    typedef typename Mesh_domain::Subdomain_index Subdomain_index;

    Is_in_domain is_in_domain = domain.is_in_domain_object();

    {
      Subdomain subdomain = is_in_domain(Point_3(CGAL::ORIGIN));
      assert(subdomain);
      Subdomain_index subdomain_index = *subdomain;
      assert(subdomain_index == 1);
    }

    {
      Subdomain subdomain = is_in_domain(Point_3(1.5, 0., 0.));
      assert(!subdomain);
    }
  }

  void test_do_intersect_surface (const Mesh_domain& domain) const
  {
    typedef typename Mesh_domain::Do_intersect_surface Do_intersect_surface;
    typedef typename Mesh_domain::Surface_patch Surface_patch;
    typedef typename Mesh_domain::Surface_patch_index Surface_patch_index;
    typedef typename Mesh_domain::Segment_3 Segment_3;
    typedef typename Mesh_domain::Ray_3 Ray_3;
    typedef typename Mesh_domain::Line_3 Line_3;
    typedef typename Mesh_domain::Vector_3 Vector_3;

    Do_intersect_surface do_intersect_surface = domain.do_intersect_surface_object();

    {
      Segment_3 s(Point_3(CGAL::ORIGIN), Point_3(1.5, 0., 0.));
      Surface_patch p = do_intersect_surface(s);
      assert(p);
      Surface_patch_index pi = *p;
      assert(pi.first == 0 && pi.second == 1);
    }

    {
      Segment_3 s(Point_3(1.5, 1.5, 0.), Point_3(1.5, 0., 0.));
      Surface_patch p = do_intersect_surface(s);
      assert(!p);
    }

    {
      Ray_3 r(Point_3(CGAL::ORIGIN), Vector_3(1., 0., 0.));
      Surface_patch p = do_intersect_surface(r);
      assert(p);
      Surface_patch_index pi = *p;
      assert(pi.first == 0 && pi.second == 1);
    }

    {
      Ray_3 r(Point_3(1.5, 0., 0.), Vector_3(0., 1., 0.));
      Surface_patch p = do_intersect_surface(r);
      assert(!p);
    }

    {
      Line_3 l(Point_3(CGAL::ORIGIN), Point_3(1.5, 0., 0.));
      Surface_patch p = do_intersect_surface(l);
      assert(p);
      Surface_patch_index pi = *p;
      assert(pi.first == 0 && pi.second == 1);
    }

    {
      Line_3 l(Point_3(1.5, 0., 0.), Point_3(1.5, 0.5, 0.));
      Surface_patch p = do_intersect_surface(l);
      assert(!p);
    }
  }


  void test_construct_intersection (const Mesh_domain& domain) const
  {
    typedef typename Mesh_domain::Construct_intersection Construct_intersection;
    typedef typename Mesh_domain::Intersection Intersection;
    typedef typename Mesh_domain::Subdomain_index Subdomain_index;
    typedef typename Mesh_domain::Surface_patch_index Surface_patch_index;
    typedef typename Mesh_domain::Index Index;
    typedef typename Mesh_domain::Segment_3 Segment_3;
    typedef typename Mesh_domain::Ray_3 Ray_3;
    typedef typename Mesh_domain::Line_3 Line_3;
    typedef typename Mesh_domain::Vector_3 Vector_3;

    Construct_intersection construct_intersection = domain.construct_intersection_object();

    {
      Segment_3 s(Point_3(CGAL::ORIGIN), Point_3(1.5, 0., 0.));
      Intersection i = construct_intersection(s);
      assert(std::get<0>(i) != Point_3(1., 0., 0.));
      Index ii = std::get<1>(i);
      assert(std::get_if<Surface_patch_index>(&ii));
      assert(std::get<2>(i) == 2);
    }

    {
      Segment_3 s(Point_3(1.5, 1.5, 0.), Point_3(1.5, 0., 0.));
      Intersection i = construct_intersection(s);
      Index ii = std::get<1>(i);
      assert(std::get_if<Subdomain_index>(&ii));
      assert(std::get<2>(i) == 0);
    }

    {
      Ray_3 r(Point_3(CGAL::ORIGIN), Vector_3(1., 0., 0.));
      Intersection i = construct_intersection(r);
      assert(std::get<0>(i) != Point_3(1., 0., 0.));
      Index ii = std::get<1>(i);
      assert(std::get_if<Surface_patch_index>(&ii));
      assert(std::get<2>(i) == 2);
    }

    {
      Ray_3 r(Point_3(1.5, 0., 0.), Vector_3(0., 1., 0.));
      Intersection i = construct_intersection(r);
      Index ii = std::get<1>(i);
      assert(std::get_if<Subdomain_index>(&ii));
      assert(std::get<2>(i) == 0);
    }

    {
      Line_3 l(Point_3(CGAL::ORIGIN), Point_3(1.5, 0., 0.));
      Intersection i = construct_intersection(l);
      assert(std::get<0>(i) != Point_3(1., 0., 0.));
      Index ii = std::get<1>(i);
      assert(std::get_if<Surface_patch_index>(&ii));
      assert(std::get<2>(i) == 2);
    }

    {
      Line_3 l(Point_3(1.5, 0., 0.), Point_3(1.5, 0.5, 0.));
      Intersection i = construct_intersection(l);
      Index ii = std::get<1>(i);
      assert(std::get_if<Subdomain_index>(&ii));
      assert(std::get<2>(i) == 0);
    }
  }
};


int main ()
{
  LM3_tester<K_e_i> test_epic;
  test_epic();

  return EXIT_SUCCESS;
}
