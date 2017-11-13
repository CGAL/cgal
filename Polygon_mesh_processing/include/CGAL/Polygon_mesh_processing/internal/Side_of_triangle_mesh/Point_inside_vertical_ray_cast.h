// Copyright (c) 2013 GeometryFactory (France).
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
// Author(s)     : Ilker O. Yaz and Sebastien Loriot


#ifndef CGAL_POINT_INSIDE_POLYHEDRON_POINT_INSIDE_VERTICAL_RAY_CAST_H
#define CGAL_POINT_INSIDE_POLYHEDRON_POINT_INSIDE_VERTICAL_RAY_CAST_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>


#include <CGAL/Polygon_mesh_processing/internal/Side_of_triangle_mesh/Ray_3_Triangle_3_traversal_traits.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Origin.h>

#include <boost/optional.hpp>
#include <boost/none.hpp>

namespace CGAL {
namespace internal {

// internal class for point inside test, using existing AABB tree
template<class Kernel, class AABBTree>
class Point_inside_vertical_ray_cast 
{
  typedef typename Kernel::Point_3       Point;
  typedef typename Kernel::Ray_3         Ray;
  typedef typename AABBTree::AABB_traits Traits;

  static const unsigned int seed = 1340818006;

public:
  Bounded_side operator()(
    const Point& point,
    const AABBTree& tree,
    typename Kernel::Construct_ray_3 ray_functor = Kernel().construct_ray_3_object(),
    typename Kernel::Construct_vector_3 vector_functor = Kernel().construct_vector_3_object() ) const
  {
    const typename Traits::Bounding_box& bbox = tree.bbox();

    if(   point.x() < bbox.xmin() || point.x() > bbox.xmax()
      || point.y() < bbox.ymin() || point.y() > bbox.ymax()
      || point.z() < bbox.zmin() || point.z() > bbox.zmax() )
    {
      return ON_UNBOUNDED_SIDE;
    }

    //the direction of the vertical ray depends on the position of the point in the bbox
    //in order to limit the expected number of nodes visited.
    Ray query = ray_functor(point, vector_functor(0,0,(2*point.z() <  tree.bbox().zmax()+tree.bbox().zmin()?-1:1)));
    boost::optional<Bounded_side> res = is_inside_ray_tree_traversal<true>(query, tree);

    if(res == boost::none)
    {
      CGAL::Random rg(seed); // seed some value for make it easy to debug
      Random_points_on_sphere_3<Point> random_point(1.,rg);

      do { //retry with a random ray
        query = ray_functor(point, vector_functor(CGAL::ORIGIN,*random_point++));
        res = is_inside_ray_tree_traversal<false>(query, tree);
      } while (res == boost::none);
    }
    return *res;
  }

private:
  template<bool ray_is_vertical>
  boost::optional<Bounded_side>
  is_inside_ray_tree_traversal(const Ray& ray, const AABBTree& tree) const
  {
    std::pair<boost::logic::tribool,std::size_t>
      status( boost::logic::tribool(boost::logic::indeterminate), 0);

    Ray_3_Triangle_3_traversal_traits<Traits, Kernel, Boolean_tag<ray_is_vertical> >
      traversal_traits(status, tree.traits());

    tree.traversal(ray, traversal_traits);

    if ( !boost::logic::indeterminate(status.first) )
    {
      if (status.first) {
        return (status.second&1) == 1 ? ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
      }
      //otherwise the point is on the facet
      return ON_BOUNDARY;
    }
    return boost::optional<Bounded_side>(); // indeterminate
  }
};

}// namespace internal
}// namespace CGAL

#endif
