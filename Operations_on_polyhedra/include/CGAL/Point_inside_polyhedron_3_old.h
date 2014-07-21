// Copyright (c) 2011 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot

/************************************************************************/
/* WARNING: This file is not a part of public and tested API + code     */
/************************************************************************/
#ifndef CGAL_POINT_INSIDE_POLYHEDRON_OLD_H
#define CGAL_POINT_INSIDE_POLYHEDRON_OLD_H

#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/internal/Operations_on_polyhedra/Ray_3_Triangle_3_traversal_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>
#include <boost/math/special_functions/round.hpp>

namespace CGAL {
  
template <class Polyhedron_3,class Kernel>
class Point_inside_polyhedron_3 {

  typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron_3> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
  typedef typename Traits::Bounding_box Bounding_box;
  typedef CGAL::AABB_tree<Traits> Tree;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Ray_3 Ray;
  typedef typename Kernel::Segment_3 Segment;
  
  Kernel m_kernel;
  Polyhedron_3& polyhedron;
  Tree tree;
  int N;
  Point grid_base;
  std::vector<std::pair<Point,bool> > grid;
  double grid_dx, grid_dy, grid_dz;
  static const bool m_on_boundary_point_inside=false;
  static const bool m_use_a_vertical_ray=true;
  double max_diagonal;

public:
  
  Point_inside_polyhedron_3(Polyhedron_3& polyhedron, int N = 0, const Kernel& k=Kernel())
    : m_kernel(k), polyhedron(polyhedron), N(N)
  {
    tree.insert(polyhedron.facets_begin(),polyhedron.facets_end());
    tree.build();
    CGAL::Bbox_3 tree_bbox = bbox();
    max_diagonal = std::sqrt( 
      CGAL::squared_distanceC3( tree_bbox.xmin(), tree_bbox.ymin(), tree_bbox.zmin(),
                                tree_bbox.xmax(), tree_bbox.ymax(), tree_bbox.zmax() ) );
    if(N>0){
      initialize_grid();
    }
  }

	
  void initialize_grid()
  {
    Random rg;
    Bounding_box bbox = tree.bbox();
    grid_dx = bbox.xmax() - bbox.xmin();
    grid_dy = bbox.ymax() - bbox.ymin();
    grid_dz = bbox.zmax() - bbox.zmin();
 	
    bbox = Bounding_box(bbox.xmin()-0.01*grid_dx, bbox.ymin()-0.01*grid_dy, bbox.zmin()-0.01*grid_dz,
                        bbox.xmax()+0.01*grid_dx, bbox.ymax()+0.01*grid_dy, bbox.zmax()+0.01*grid_dz);
    grid_dx = (bbox.xmax() - bbox.xmin())/(N-1);
    grid_dy = (bbox.ymax() - bbox.ymin())/(N-1);
    grid_dz = (bbox.zmax() - bbox.zmin())/(N-1);
    grid_base = Point(bbox.xmin(), bbox.ymin(), bbox.zmin());
 

    grid.reserve(N*N*N);
    int points_inside = 0, points_outside=0;
    for(int i=0; i < N; i++){
      for(int j=0; j < N; j++){
        for(int k=0; k < N; k++){
          Point p(bbox.xmin()+i*grid_dx, bbox.ymin()+j*grid_dy, bbox.zmin()+k*grid_dz);
          if(i==0 || j==0 || k==0 || i==N-1 || j==N-1 || k==N-1){
            double eps = grid_dx/100;
            Vector v(rg.get_double(-eps,eps), rg.get_double(-eps,eps), rg.get_double(-eps,eps));
            p = p+ v;
            grid.push_back(std::make_pair(p,false));
              points_outside++;	
          } else {
            double eps = grid_dx/100;
            Vector v(rg.get_double(-eps,eps), rg.get_double(-eps,eps), rg.get_double(-eps,eps));
            p = p+ v;
            const Segment segment(p, grid.back().first);
            std::size_t M = (grid.back().second)? 0 : 1;
            bool inside = (tree.number_of_intersected_primitives(segment)&1) == M;
            if(inside){
              points_inside++;
            }else{
              points_outside++;
            }
            grid.push_back(std::make_pair(p,inside));
          }
        }
      }
    }
    std::cerr << points_inside << " points inside" << std::endl;
    std::cerr << points_outside << " points outside" << std::endl;
  }
  
  Bounding_box bbox() const {
    return tree.bbox();
  }

private:
  template <class Query,bool ray_is_vertical>
  boost::logic::tribool 
  is_inside_ray_tree_traversal(const Query& query) const {
    std::pair<boost::logic::tribool,std::size_t> status(
      boost::logic::tribool(boost::logic::indeterminate), 0);
    internal::Ray_3_Triangle_3_traversal_traits<Traits,Kernel,Boolean_tag<ray_is_vertical> > traversal_traits(status);
    tree.traversal(query, traversal_traits);
    if ( !boost::logic::indeterminate(status.first) ){
      if (status.first) return (status.second&1) == 1;
      //otherwise the point is on the facet
      return m_on_boundary_point_inside;
    }
    return boost::logic::indeterminate;
  }

public:

  bool operator()(const Point& p) const
  {
    const Bounding_box& bbox = tree.bbox();

    if(   p.x() < bbox.xmin() || p.x() > bbox.xmax()
       || p.y() < bbox.ymin() || p.y() > bbox.ymax()
       || p.z() < bbox.zmin() || p.z() > bbox.zmax() )
    {
      return false;
    }

    if(N>0)
    {
      Vector v = p - grid_base;
      int i = boost::math::iround(v.x() / grid_dx);  
      int j =  boost::math::iround(v.y() / grid_dy);  
      int k =  boost::math::iround(v.z() / grid_dz);

      if(i>N-1)i=N-1;
      if(j>N-1)j=N-1;
      if(k>N-1)k=N-1;
      int index = i*N*N + j*N + k;
       
      const std::pair<Point,bool>& close_point = grid[index];
      typename Kernel::Construct_segment_3 segment = m_kernel.construct_segment_3_object();
      const Segment query = segment(p, close_point.first);
      if(p == close_point.first){
        std::cerr << "error"  << std::endl;
      }
      std::size_t M = (close_point.second)? 0 : 1; 
      bool res = (tree.number_of_intersected_primitives(query)&1) == M;
      return res;
    } 
    else 
    {
      typename Kernel::Construct_ray_3 make_ray = m_kernel.construct_ray_3_object();
      typename Kernel::Construct_vector_3 make_vector = m_kernel.construct_vector_3_object();
     
      //start with a vertical ray
      //~ Ray query = ray(p, vector(CGAL::ORIGIN,*random_point));
      
      Random_points_on_sphere_3<Point> random_point(1.);
      //the direction of the vertical ray depends on the position of the point in the bbox
      //in order to limit the expected number of nodes visited.
      Ray query = 
        m_use_a_vertical_ray ?
        make_ray(p, make_vector(0,0,(2*p.z() <  tree.bbox().zmax()+tree.bbox().zmin()?-1:1))) :
        make_ray(p, make_vector(CGAL::ORIGIN,*random_point));

      //double max_distance( max_diagonal / std::sqrt(query.to_vector().squared_length()));
      //const Vector& scaled_direction = m_kernel.construct_scaled_vector_3_object()(query.to_vector(), max_distance);
      //const Vector& target_vector = m_kernel.construct_sum_of_vectors_3_object()( Vector(Point(ORIGIN), p), scaled_direction);
      //const Point&  target_point = m_kernel.construct_translated_point_3_object()(Point(ORIGIN), target_vector);
      //Segment segment(p, target_point); 

      boost::logic::tribool res=is_inside_ray_tree_traversal<Ray,true>(query);
      while (boost::logic::indeterminate(res)){
        //retry with a random ray
        query = make_ray(p, make_vector(CGAL::ORIGIN,*random_point++));
        res=is_inside_ray_tree_traversal<Ray,false>(query);
      }
      return res;
    }
  }
  
  //the original version
  bool operator()(const Point& p,bool) const
  {
    typename Kernel::Construct_ray_3 ray = m_kernel.construct_ray_3_object();
    typename Kernel::Construct_vector_3 vector = m_kernel.construct_vector_3_object();
   
    Random_points_on_sphere_3<Point> random_point(1.);
    
    //const Ray query = ray(p, vector(CGAL::ORIGIN,*random_point));
    
    //the direction of the vertical ray depends on the position of the point in the bbox
    //in order to limit the expected number of nodes visited.
    const Ray query = ray(p, vector(0,0, (2*p.z() <  tree.bbox().zmax()+tree.bbox().zmin()?-1:1) ));

    return (tree.number_of_intersected_primitives(query)&1) == 1;
  }
  
};

} // namespace CGAL

#endif //CGAL_POINT_INSIDE_POLYHEDRON_H
