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

#include <CGAL/Polygon_mesh_processing/internal/Point_inside_polygon_mesh/Point_inside_vertical_ray_cast.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL {

/** 
 * \ingroup PkgPolygonMeshProcessing
 * This class provides an efficient point location functionality with respect to a domain bounded
 * by one or several disjoint closed triangle mesh.
 * In case the triangle mesh as several connected components, a point is said to be inside the domain
 * if an odd number of surfaces is crossed when walking from the point to infinity.
 *
 * This class depends on the package \ref PkgAABB_treeSummary.

 * @tparam TriangleMesh a triangulated polyhedral surface, a model of `FaceListGraph`
 * @tparam GeomTraits a geometric traits class, model of `Kernel`
 * @tparam VertexPointMap a model of `ReadablePropertyMap` with
         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
         `GeomTraits::Point_3` as value type.
 *   The default is `typename boost::property_map<TriangleMesh,vertex_point_t>::%type`.
 
 * \todo Code: Use this class as an implementation detail of Mesh_3's Polyhedral_mesh_domain_3.
       Remove `TriangleAccessor_3` as well as the concept in Mesh_3 since making `TriangleMesh`
       a model of `FaceListGraph` will make it useless
 * \todo check the implementation
 */
template <class TriangleMesh,
          class GeomTraits,
          class VertexPointMap = Default >
class Point_inside_polygon_mesh
{
  // typedefs
  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, VertexPointMap> Primitive;
  typedef CGAL::AABB_traits<GeomTraits, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> AABB_tree;
  typedef typename GeomTraits::Point_3 Point;

  //members
  typename GeomTraits::Construct_ray_3     ray_functor;
  typename GeomTraits::Construct_vector_3  vector_functor;
  const AABB_tree* tree_ptr;
  bool own_tree;

public:
   /**
   * Constructor with one surface triangle mesh.
   * @param mesh the triangle mesh bounding the domain to be tested
   * @param vpmap the property map with the points associated to the vertices of `mesh`
   * @param gt an instance of the geometric traits class
   *
   * @pre `CGAL::is_closed(mesh) && CGAL::is_pure_triangle(mesh)`
   */
  Point_inside_polygon_mesh(const TriangleMesh& mesh,
                            VertexPointMap vpmap,
                            const GeomTraits& gt=GeomTraits())
  : ray_functor(gt.construct_ray_3_object())
  , vector_functor(gt.construct_vector_3_object())
  , own_tree(true)
  {
    CGAL_assertion(CGAL::is_pure_triangle(mesh));
    CGAL_assertion(CGAL::is_closed(mesh));

    tree_ptr = new AABB_tree(faces(mesh).first,
                             faces(mesh).second,
                             mesh, vpmap);
  }

  /**
  * Constructor with one surface triangle mesh, using `get(boost::vertex_point, mesh)` as
  * vertex point property map.
  * @param mesh the triangle mesh bounding the domain to be tested
  * @param gt an instance of the geometric traits class
  *
  * @pre `mesh` must be closed and triangulated.
  */
  Point_inside_polygon_mesh(const TriangleMesh& mesh,
                            const GeomTraits& gt=GeomTraits())
  : ray_functor(gt.construct_ray_3_object())
  , vector_functor(gt.construct_vector_3_object())
  , own_tree(true)
  {
    CGAL_assertion(CGAL::is_pure_triangle(mesh));
    CGAL_assertion(CGAL::is_closed(mesh));

    tree_ptr = new AABB_tree(faces(mesh).first,
                             faces(mesh).second,
                             mesh);
  }

  /**
  * Constructor that takes a pre-built \cgal `AABB_tree`
  * of the triangle mesh primitives.
  * Note the domain described by these primitives should be closed.

  * @param tree a \cgal `AABB_tree` with `AABB_face_graph_triangle_primitive` as `Primitive` type
  * @param gt an instance of the geometric traits class
  */
  Point_inside_polygon_mesh(const AABB_tree& tree,
    const GeomTraits& gt = GeomTraits())
  : ray_functor(gt.construct_ray_3_object())
  , vector_functor(gt.construct_vector_3_object())
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
   * returns the location of a query point
   * @param point the query point to be located with respect to the input
            polyhedral surface
   * @return 
   *   - `CGAL::ON_BOUNDED_SIDE` if the point is inside the triangle mesh
   *   - `CGAL::ON_BOUNDARY` if the point is on triangle mesh
   *   - `CGAL::ON_UNBOUNDED_SIDE` if the point is outside triangle mesh
   */
  Bounded_side operator()(const Point& point) const
  {
    return internal::Point_inside_vertical_ray_cast<GeomTraits, AABB_tree>()(
      point, *tree_ptr, ray_functor, vector_functor);
  }

};

} // namespace CGAL

#endif //CGAL_POINT_INSIDE_POLYGON_MESH_H
