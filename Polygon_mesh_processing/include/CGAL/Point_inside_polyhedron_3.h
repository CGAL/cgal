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
// Author(s)     : Sebastien Loriot and Ilker O. Yaz


#ifndef CGAL_POINT_INSIDE_POLYHEDRON_H
#define CGAL_POINT_INSIDE_POLYHEDRON_H

#include <CGAL/internal/Point_inside_polyhedron_3/Point_inside_vertical_ray_cast.h>
#include <CGAL/internal/Point_inside_polyhedron_3/AABB_triangle_accessor_3_primitive.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Triangle_accessor_3.h>

namespace CGAL {

/** 
 * \ingroup PkgPolygonMeshProcessing
 * This class provides an efficient point location functionality with respect to a domain bounded
 * by one or several disjoint triangulated closed polyhedral (manifold) surfaces.
 * In case several polyhedral surface are provided as input, a point is said to be inside the domain
 * if an odd number of surfaces is crossed when walking from infinity to the point.
 * The implementation depends on the package \ref PkgAABB_treeSummary.
 * @tparam Polyhedron a triangulated polyhedral surface
 * @tparam Kernel a \cgal kernel
 * @tparam TriangleAccessor a model of the concept `TriangleAccessor_3`, with `TriangleAccessor_3::Triangle_3` being `Kernel::Triangle_3`. 
 *         If `Polyhedron` is a \cgal Polyhedron, a default is provided.
 * \todo Code: Use this class as an implementation detail of Mesh_3's Polyhedral_mesh_domain_3
 * \todo Code: current version puts all polyhedra under one AABB, more proper approach might be using separate AABB for each polyhedron 
 *       and filtering query point with bboxes of polyhedra...
 * \todo `Polyhedron` should be a model of `FaceListGraph`
 * \todo Remove `TriangleAccessor_3` as well as the concept in Mesh_3 since making `Polyhedron` a model of `FaceListGraph` will make it useless
 * \todo Add a constructor from AABB-tree (once Polyhedron is a FaceListGraph, the type is hardcoded)
 * \todo check the implementation
 */
template <class Polyhedron, 
          class Kernel,
          class TriangleAccessor_3 = Triangle_accessor_3<Polyhedron, typename Polyhedron::Traits> 
>
class Point_inside_polyhedron_3{
  // typedefs
  typedef CGAL::internal::AABB_triangle_accessor_3_primitive<Kernel, TriangleAccessor_3> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;
  typedef typename Kernel::Point_3 Point;

  //members
  typename Kernel::Construct_ray_3     ray_functor;
  typename Kernel::Construct_vector_3  vector_functor;
  Tree tree;

public:
  /**
   * Default constructor. The domain is considered to be empty.
   */
  Point_inside_polyhedron_3(const Kernel& kernel=Kernel())
  : ray_functor(kernel.construct_ray_3_object()),
  vector_functor(kernel.construct_vector_3_object())
  { }
 
  /** 
   * Constructor with one polyhedral surface. `polyhedron` must be closed and triangulated.
   */
  Point_inside_polyhedron_3(const Polyhedron& polyhedron, const Kernel& kernel=Kernel()) 
  : ray_functor(kernel.construct_ray_3_object()),
  vector_functor(kernel.construct_vector_3_object())
  {
    add_polyhedron(polyhedron);
  }

  /** 
   * Constructor with several polyhedral surfaces. All the polyhedral surfaces must be closed, triangulated and disjoint.
   * \tparam InputIterator is an input iterator with `Polyhedron` or `cpp11::reference_wrapper<Polyhedron>` as value type.
   */
  template <class InputIterator>
  Point_inside_polyhedron_3(InputIterator begin, InputIterator beyond, const Kernel& kernel=Kernel()) 
  : ray_functor(kernel.construct_ray_3_object()),
  vector_functor(kernel.construct_vector_3_object())
  {
    add_polyhedra(begin, beyond);
  }

  /** 
   * Builds internal AABB tree. Optional to call, since the tree is automatically built at the time of first query.
   */ 
  void build() { tree.build(); }

  /**
   * `polyhedron` is added as input
   */
  void add_polyhedron(const Polyhedron& polyhedron) 
  {
    CGAL_assertion(polyhedron.is_pure_triangle());
    CGAL_assertion(polyhedron.is_closed());

    tree.insert(TriangleAccessor_3().triangles_begin(polyhedron),
                TriangleAccessor_3().triangles_end(polyhedron));
  }
 
  /**
   * The polyhedral surfaces in the range `[begin,beyond[` are added as input
   * \tparam InputIterator is an input iterator with `Polyhedron` or `cpp11::reference_wrapper<Polyhedron>` as value type.
   */
  template<class InputIterator>
  void add_polyhedra(InputIterator begin, InputIterator beyond) 
  {
    for(; begin != beyond; ++begin) {
      add_polyhedron(*begin);
    }
  }
 
  /**
   * Query function to determine point location.
   * @return 
   *   - CGAL::ON_BOUNDED_SIDE if the point is inside the polyhedron
   *   - CGAL::ON_BOUNDARY if the point is on polyhedron
   *   - CGAL::ON_UNBOUNDED_SIDE if the point is outside polyhedron
   */
  Bounded_side operator()(const Point& point) const
  {
    return internal::Point_inside_vertical_ray_cast<Kernel, Tree>()(point, tree, ray_functor, vector_functor);
  }
};

} // namespace CGAL

#endif //CGAL_POINT_INSIDE_POLYHEDRON_H
