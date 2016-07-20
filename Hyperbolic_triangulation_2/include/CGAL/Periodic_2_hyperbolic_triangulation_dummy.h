// Copyright (c) 2010-2016  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: 
// $Id: 
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_PERIODIC_2_HYPERBOLIC_TRIANGULATION_DUMMY_H
#define CGAL_PERIODIC_2_HYPERBOLIC_TRIANGULATION_DUMMY_H

#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Aff_transformation_2.h>

namespace CGAL {

// this class might to be moved to hyperbolic traits
template<class GT>
class Construct_reflection  
{
  typedef typename GT::Point_2 Point_2;
  typedef typename GT::Circle_2 Circle_2;
  typedef typename GT::Vector_2 Vector_2;
  typedef typename GT::FT FT;
  
  typedef typename GT::Compute_squared_distance_2 Compute_squared_distance_2;
  typedef	typename GT::Construct_vector_2 Construct_vector_2;
  typedef	typename GT::Construct_scaled_vector_2 Construct_scaled_vector_2;
  typedef	typename GT::Construct_translated_point_2 Construct_translated_point_2;
  
  typedef CGAL::Regular_triangulation_filtered_traits_2<GT> Rtr_2;
  typedef typename Rtr_2::Construct_weighted_circumcenter_2 Construct_weighted_circumcenter_2;
  typedef typename Rtr_2::Weighted_point_2 Weighted_point_2;
  typedef typename Rtr_2::Bare_point Bare_point;
  
public:
  Construct_reflection(const Circle_2& c = Circle_2(Point_2(0.0, 0.0), 1.)): _unit_circle(c)
  {
  }
  
  // exact arithmetic! very nice!
  Point_2 operator ()(const Point_2& p, const Point_2& q, const Point_2& r) const
  {					
    typedef typename GT::Collinear_2 Collinear_2;
    if(Collinear_2()(p, q, _unit_circle.center())){
      assert(false);
    }
    Weighted_point_2 wp(p);
    Weighted_point_2 wq(q);
    Weighted_point_2 wo(_unit_circle.center(), _unit_circle.squared_radius());
    
    Bare_point center = Construct_weighted_circumcenter_2()(wp, wo, wq);
    FT radius = Compute_squared_distance_2()(p, center);
    FT dist = Compute_squared_distance_2()(r,center);
    Vector_2 vect = Construct_vector_2()(center,r);
    vect = Construct_scaled_vector_2()(vect,radius/dist);
    Point_2 image = Construct_translated_point_2()(center,vect);
    return image;
  }
  
private:
  const Circle_2& _unit_circle;
};

template<class GT>
typename GT::Point_2 apply_rotation(const typename GT::Point_2& p)
{
  CGAL::Aff_transformation_2<GT> rotate(CGAL::ROTATION, std::sqrt(0.5), std::sqrt(0.5));
	return rotate(p);
}
  
template<class GT>
void compute_redundant_dummy_points(std::vector<typename GT::Point_2>& inner_points, std::vector<typename GT::Point_2>& points_on_boundary, std::vector<typename GT::Point_2>& points_on_vertex)
{
  assert(inner_points.size() == 0);
  assert(points_on_boundary.size() == 0);
  assert(points_on_vertex.size() == 0);
  
  typedef typename GT::Kernel K;
  typedef typename GT::Point_2 Point_2;
  
  double phi = CGAL_PI / 8.;
  double psi = CGAL_PI / 3.;
  double rho = std::sqrt(cos(psi)*cos(psi) - sin(phi)*sin(phi));
  
  const Point_2 o(0.0, 0.0);
  const Point_2 a(cos(phi)*cos(phi + psi)/rho, sin(phi)*cos(phi + psi)/rho);
  const Point_2 b(a.x(), -a.y());
  
  inner_points.push_back(a);
  Point_2 c = Construct_reflection<K>()(a, b, o);
  Point_2 d = Construct_reflection<K>()(a, c, b);
  inner_points.push_back(d);
  Point_2 e = Construct_reflection<K>()(d, c, a);
  Point_2 f = Construct_reflection<K>()(d, e, c);
  
  int size = inner_points.size();
  for(int i = 0; i < 7; i++) {
    for(int j = 0; j < size; j++) {
      inner_points.push_back(apply_rotation<K>(inner_points[i*size + j]));
    }
  }
  inner_points.push_back(o);
  
  points_on_boundary.push_back(c);
  points_on_vertex.push_back(f);
  for(int i = 1; i < 8; i++) {
    points_on_boundary.push_back(apply_rotation<K>(points_on_boundary[i-1]));
    points_on_vertex.push_back(apply_rotation<K>(points_on_vertex[i-1]));
  }  
}

/*
template < class GT, class TDS >
void Periodic_2_Delaunay_hyperbolic_triangulation_2<GT, TDS>::insert_dummy_points() {
  clear();
  
  std::vector<Point_2> inner_points, points_on_boundary, points_on_vertex;
  
  compute_redundant_dummy_points<GT>(inner_points, points_on_boundary, points_on_vertex);
  Base::insert(inner_points.begin(), inner_points.end());
  
  size_t on_boundary_nb = points_on_boundary.size();
  std::vector<Vertex_handle> vertices_on_boundary(on_boundary_nb);
  for(size_t i = 0; i < on_boundary_nb; i++) {
    vertices_on_boundary[i] = Base::insert(points_on_boundary[i]);
  }
  
  size_t on_vertex_nb = points_on_vertex.size();
  std::vector<Vertex_handle> vertices_on_vertex(on_vertex_nb);
  for(size_t i = 0; i < on_vertex_nb; i++) {
    vertices_on_vertex[i] = Base::insert(points_on_vertex[i]);
  }
  
}*/

} // namespace CGAL 
  
#endif // CGAL_PERIODIC_2_HYPERBOLIC_TRIANGULATION_DUMMY_H
  
