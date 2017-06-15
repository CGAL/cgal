// Copyright (c) 2013, 2014, 2015 GeometryFactory (France).
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


#ifndef CGAL_SIDE_OF_TRIANGLE_MESH_H
#define CGAL_SIDE_OF_TRIANGLE_MESH_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>


#include <CGAL/Polygon_mesh_processing/internal/Side_of_triangle_mesh/Point_inside_vertical_ray_cast.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL {

/** 
 * \ingroup PkgPolygonMeshProcessing
 * This class provides an efficient point location functionality with respect to a domain bounded
 * by one or several disjoint closed triangle meshes.
 *
 * A point is said to be on the bounded side of the domain
 * if an odd number of surfaces is crossed when walking from the point to infinity.
 *
 * The input triangle mesh is expected to contain no self-intersections
 * and to be free from self-inclusions.
 *
 * In case the triangle mesh has several connected components,
 * the same test is performed and returns correct results.
 * In case of self-inclusions,
 * the user should be aware that the predicate called
 * inside every other sub-volume bounded by a nested surface
 * will return in turns `CGAL::ON_BOUNDED_SIDE` and `CGAL::ON_UNBOUNDED_SIDE`,
 * following the aforementioned parity criterion.
 *
 * This class depends on the package \ref PkgAABB_treeSummary.
 *
 * @tparam TriangleMesh a triangulated surface mesh, model of `FaceListGraph`
 * @tparam GeomTraits a geometric traits class, model of `Kernel`
 * @tparam VertexPointMap a model of `ReadablePropertyMap` with
         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
         `GeomTraits::Point_3` as value type.
 *   The default is `typename boost::property_map<TriangleMesh,vertex_point_t>::%type`.
 *
 * \cgalHeading{Implementation Details}
 * The current implementation is based on the number of triangles intersected by a ray
 * having the query point as source. The do-intersect predicate used to detect if a triangle
 * is intersected is able to detect if a triangle is intersected in its
 * interior or on its boundary. In case it is intersected on its boundary, another ray is picked.
 * In order to speed queries, the first ray used is an axis aligned one that depends on the extents of the
 * bbox of the input mesh. In case other rays are needed to conclude, the rays are generated
 * from a random uniform sampling of a sphere.
 */
template <class TriangleMesh,
          class GeomTraits,
          class VertexPointMap = Default >
class Side_of_triangle_mesh
{
  // typedefs
  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, VertexPointMap> Primitive;
  typedef CGAL::AABB_traits<GeomTraits, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> AABB_tree_;
  typedef typename GeomTraits::Point_3 Point;

  //members
  typename GeomTraits::Construct_ray_3     ray_functor;
  typename GeomTraits::Construct_vector_3  vector_functor;
  const AABB_tree_* tree_ptr;
  bool own_tree;

public:

   #ifndef DOXYGEN_RUNNING
   typedef AABB_tree_ AABB_tree;
   #else
   /// AABB-tree accepting faces of `TriangleMesh`
   typedef unspecified_type AABB_tree;
   #endif

   /**
   * Constructor with one triangulated surface mesh.
   * @param tmesh the triangulated surface mesh bounding the domain to be tested
   * @param vpmap the property map with the points associated to the vertices of `tmesh`
   * @param gt an instance of the geometric traits class
   *
   * @pre `CGAL::is_closed(tmesh) && CGAL::is_triangle_mesh(tmesh)`
   */
  Side_of_triangle_mesh(const TriangleMesh& tmesh,
                        VertexPointMap vpmap,
                        const GeomTraits& gt=GeomTraits())
  : ray_functor(gt.construct_ray_3_object())
  , vector_functor(gt.construct_vector_3_object())
  , own_tree(true)
  {
    CGAL_assertion(CGAL::is_triangle_mesh(tmesh));
    CGAL_assertion(CGAL::is_closed(tmesh));

    tree_ptr = new AABB_tree(faces(tmesh).first,
                             faces(tmesh).second,
                             tmesh, vpmap);
  }

  /**
  * Constructor with one surface triangle mesh, using `get(boost::vertex_point, tmesh)` as
  * vertex point property map.
  * @param tmesh the triangulated surface mesh bounding the domain to be tested
  * @param gt an instance of the geometric traits class
  *
  * @pre `CGAL::is_closed(tmesh) && CGAL::is_triangle_mesh(tmesh)`
  */
  Side_of_triangle_mesh(const TriangleMesh& tmesh,
                        const GeomTraits& gt=GeomTraits())
  : ray_functor(gt.construct_ray_3_object())
  , vector_functor(gt.construct_vector_3_object())
  , own_tree(true)
  {
    CGAL_assertion(CGAL::is_triangle_mesh(tmesh));
    CGAL_assertion(CGAL::is_closed(tmesh));

    tree_ptr = new AABB_tree(faces(tmesh).first,
                             faces(tmesh).second,
                             tmesh);
  }

  /**
  * Constructor that takes a pre-built \cgal `AABB_tree`
  * of the triangulated surface mesh primitives.
  *
  * @param tree a \cgal `AABB_tree` with `AABB_face_graph_triangle_primitive` as `Primitive` type
  * @param gt an instance of the geometric traits class
  *
  * @pre `CGAL::is_closed(tmesh) && CGAL::is_triangle_mesh(tmesh)`
  */
  Side_of_triangle_mesh(const AABB_tree& tree,
                        const GeomTraits& gt = GeomTraits())
  : ray_functor(gt.construct_ray_3_object())
  , vector_functor(gt.construct_vector_3_object())
  , tree_ptr(&tree)
  , own_tree(false)
  {
  }

  ~Side_of_triangle_mesh()
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
   *   - `CGAL::ON_BOUNDED_SIDE` if the point is inside the volume bounded by the input triangle mesh
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

#endif //CGAL_SIDE_OF_TRIANGLE_MESH_H
