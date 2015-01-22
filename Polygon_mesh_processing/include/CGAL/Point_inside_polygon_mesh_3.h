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


#ifndef CGAL_POINT_INSIDE_POLYGON_MESH_H
#define CGAL_POINT_INSIDE_POLYGON_MESH_H

#include <CGAL/internal/Point_inside_polygon_mesh/Point_inside_vertical_ray_cast.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

namespace CGAL {

/** 
 * \ingroup PkgPolygonMeshProcessing
 * This class provides an efficient point location functionality with respect to a domain bounded
 * by one or several disjoint triangulated closed polyhedral (manifold) surfaces.
 * In case several polyhedral surface are provided as input, a point is said to be inside the domain
 * if an odd number of surfaces is crossed when walking from infinity to the point.
 * The implementation depends on the package \ref PkgAABB_treeSummary.
 * @tparam TriangleMesh a triangulated polyhedral surface
 * @tparam Kernel a \cgal kernel
 * @tparam TriangleAccessor a model of the concept `TriangleAccessor_3`, with `TriangleAccessor_3::Triangle_3` being `Kernel::Triangle_3`. 
 *         If `TriangleMesh` is a \cgal Polyhedron, a default is provided.
 * \todo Code: Use this class as an implementation detail of Mesh_3's Polyhedral_mesh_domain_3.
       Remove `TriangleAccessor_3` as well as the concept in Mesh_3 since making `TriangleMesh`
       a model of `FaceListGraph` will make it useless
 * \todo `TriangleMesh` should be a model of `FaceListGraph`
 * \todo check the implementation
 */
template <class TriangleMesh,
          class Kernel>
class Point_inside_polygon_mesh
{
  // typedefs
  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> AABB_tree;
  typedef typename Kernel::Point_3 Point;

  //members
  typename Kernel::Construct_ray_3     ray_functor;
  typename Kernel::Construct_vector_3  vector_functor;
  const AABB_tree* tree_ptr;
  bool own_tree;

public:
   /**
   * Constructor with one surface polygon mesh.
   * @pre `mesh` must be closed and triangulated.
   */
  Point_inside_polygon_mesh(const TriangleMesh& mesh,
                            const Kernel& kernel=Kernel())
  : ray_functor(kernel.construct_ray_3_object())
  , vector_functor(kernel.construct_vector_3_object())
  , own_tree(true)
  {
    CGAL_assertion(mesh.is_pure_triangle());
    CGAL_assertion(mesh.is_closed());

    tree_ptr = new AABB_tree(faces(mesh).first,
                             faces(mesh).second,
                             mesh);
  }

  /**
  * Constructor that takes a \cgal `AABB_tree` with 
  * `AABB_face_graph_triangle_primitive` as Primitive type.
  * Note the domain should be closed.
  */
  Point_inside_polygon_mesh(const AABB_tree& tree,
    const Kernel& kernel = Kernel())
  : ray_functor(kernel.construct_ray_3_object())
  , vector_functor(kernel.construct_vector_3_object())
  , tree_ptr(&tree)
  , own_tree(false)
  {
  }

  ~Point_inside_polygon_mesh()
  {
    if (own_tree)
      delete tree_ptr;
  }

public:
  /**
   * Query function to determine point location.
   * @return 
   *   - CGAL::ON_BOUNDED_SIDE if the point is inside the mesh
   *   - CGAL::ON_BOUNDARY if the point is on mesh
   *   - CGAL::ON_UNBOUNDED_SIDE if the point is outside mesh
   */
  Bounded_side operator()(const Point& point) const
  {
    return internal::Point_inside_vertical_ray_cast<Kernel, AABB_tree>()(
      point, *tree_ptr, ray_functor, vector_functor);
  }

};

} // namespace CGAL

#endif //CGAL_POINT_INSIDE_POLYGON_MESH_H
