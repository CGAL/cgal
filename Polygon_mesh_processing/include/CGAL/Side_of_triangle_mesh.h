// Copyright (c) 2013, 2014, 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// 
//
// Author(s)     : Sebastien Loriot and Ilker O. Yaz


#ifndef CGAL_SIDE_OF_TRIANGLE_MESH_H
#define CGAL_SIDE_OF_TRIANGLE_MESH_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_mesh_processing/internal/Side_of_triangle_mesh/Point_inside_vertical_ray_cast.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL {

/** 
 * \ingroup PkgPolygonMeshProcessingRef
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
 * This class depends on the package \ref PkgAABBTree.
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
          class VertexPointMap_ = Default
#ifndef DOXYGEN_RUNNING
          , class AABBTree = Default
#endif
          >
class Side_of_triangle_mesh
{
  // typedefs
  template <typename TriangleMesh_,
            typename GeomTraits_,
            typename VertexPointMap__>
  struct AABB_tree_default {
    typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh_,
                                                     VertexPointMap__> Primitive;
    typedef CGAL::AABB_traits<GeomTraits_, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> type;
  };
  typedef typename Default::Lazy_get<AABBTree,
                                     AABB_tree_default<TriangleMesh,
                                                       GeomTraits,
                                                       VertexPointMap_>
                                     >::type AABB_tree_;
  typedef typename Default::Get<VertexPointMap_,
                                typename boost::property_map<TriangleMesh,
                                                             vertex_point_t>::const_type>::type
                                  VertexPointMap;
  typedef typename GeomTraits::Point_3 Point;

  //members
  typename GeomTraits::Construct_ray_3     ray_functor;
  typename GeomTraits::Construct_vector_3  vector_functor;
  mutable const AABB_tree_* tree_ptr;
  const TriangleMesh* tm_ptr;
  boost::optional<VertexPointMap> opt_vpm;
  bool own_tree;
  CGAL::Bbox_3 box;
#ifdef CGAL_HAS_THREADS
  mutable CGAL_MUTEX tree_mutex;
#endif

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
  , tree_ptr(nullptr)
  , tm_ptr(&tmesh)
  , opt_vpm(vpmap)
  , own_tree(true)
  {
    CGAL_assertion(CGAL::is_triangle_mesh(tmesh));
    CGAL_assertion(CGAL::is_closed(tmesh));
    box = Polygon_mesh_processing::bbox(tmesh, parameters::vertex_point_map(vpmap));
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
  : Side_of_triangle_mesh(tmesh, get(vertex_point, tmesh), gt)
  {}

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
    box = tree.bbox();
  }

  ~Side_of_triangle_mesh()
  {
    if (own_tree && tree_ptr!=nullptr)
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
    if(point.x() < box.xmin()
       || point.x() > box.xmax()
       || point.y() < box.ymin()
       || point.y() > box.ymax()
       || point.z() < box.zmin()
       || point.z() > box.zmax())
    {
      return CGAL::ON_UNBOUNDED_SIDE;
    }
    else
    {
      // Lazily build the tree only when needed
      if (tree_ptr==nullptr)
      {
#ifdef CGAL_HAS_THREADS
        CGAL_SCOPED_LOCK(tree_mutex);
#endif
        CGAL_assertion(tm_ptr != nullptr && opt_vpm!=boost::none);
        if (tree_ptr==nullptr)
        {
          tree_ptr = new AABB_tree(faces(*tm_ptr).first,
                                   faces(*tm_ptr).second,
                                   *tm_ptr, *opt_vpm);
          const_cast<AABB_tree_*>(tree_ptr)->build();
        }
      }

      return internal::Point_inside_vertical_ray_cast<GeomTraits, AABB_tree>()(
            point, *tree_ptr, ray_functor, vector_functor);
    }
  }

};

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_SIDE_OF_TRIANGLE_MESH_H
