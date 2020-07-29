// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno and Sebastien Loriot


#ifndef CGAL_RIGID_TRIANGLE_MESH_COLLISION_DETECTION_H
#define CGAL_RIGID_TRIANGLE_MESH_COLLISION_DETECTION_H

#include <CGAL/license/Polygon_mesh_processing/collision_detection.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polygon_mesh_processing/internal/AABB_traversal_traits_with_transformation.h>
#include <CGAL/Polygon_mesh_processing/internal/Side_of_triangle_mesh/Point_inside_vertical_ray_cast.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/property_map.h>

#include <boost/iterator/counting_iterator.hpp>

#ifndef CGAL_RMCD_CACHE_BOXES
#define CGAL_RMCD_CACHE_BOXES 0
#endif

#if CGAL_RMCD_CACHE_BOXES
#include <boost/dynamic_bitset.hpp>
#endif

namespace CGAL {

/*!
 * \ingroup PkgPolygonMeshProcessing
 * This class provides methods to perform some intersection tests between triangle meshes
 * that undergo affine transformations (rotation, translation, and scaling).
 * Meshes are added to an internal set and are referenced using an id assigned when added to the set.
 * Note that the exact predicate framework applies on the meshes after having applied the transformation
 * to the coordinates of the points of the vertices of each mesh.
 *
 * @tparam TriangleMesh a model of `HalfedgeListGraph` and `FaceListGraph`
 * @tparam VertexPointMap a model of `ReadablePropertyMap` with the vertex descriptor of `TriangleMesh` as key type,
 *                        and a point from a CGAL Kernel as value type. %Default is the internal point property map
 *                        of `TriangleMesh` if it exists.
 * @tparam Kernel a model of CGAL Kernel. %Default is the kernel of the value type of `VertexPointMap` retrieved using
 *                `Kernel_traits`.
 * @tparam AABBTree an `AABB_tree` that can containing faces of `TriangleMesh`. %Default is using `AABB_traits` with
 *                       `AABB_face_graph_triangle_primitive` as primitive type.
 * @tparam Has_rotation tag indicating whether the transformations applied to meshes may contain rotations (`Tag_true`)
 *                      or if only translations and scalings are applied (`Tag_false`). Some optimizations are
 *                      switch on in case there are no rotations.
 */
template <class TriangleMesh,
          class VertexPointMap = Default,
          class Kernel = Default,
          class AABBTree = Default,
          class Has_rotation = CGAL::Tag_true>
class Rigid_triangle_mesh_collision_detection
{
// Vertex point map type
  typedef typename property_map_selector<TriangleMesh, boost::vertex_point_t
    >::const_type                                                   Default_vpm;
  typedef typename Default::Get<VertexPointMap, Default_vpm>::type          Vpm;

// Kernel type
  typedef typename Kernel_traits<
    typename boost::property_traits<Vpm>::value_type>::Kernel    Default_kernel;
  typedef typename Default::Get<Kernel, Default_kernel>::type                 K;

// AABB-tree type
  typedef AABB_face_graph_triangle_primitive<TriangleMesh,
                                             Vpm>             Default_primitive;
  typedef AABB_traits<K, Default_primitive>                 Default_tree_traits;
  typedef CGAL::AABB_tree<Default_tree_traits>                     Default_tree;
  typedef typename Default::Get<AABBTree, Default_tree>::type              Tree;
  typedef typename Tree::AABB_traits                                Tree_traits;

// Transformed Tree traversal traits
  typedef Do_intersect_traversal_traits_with_transformation<Tree_traits,
                                                            K,
                                                            Has_rotation>
                                                               Traversal_traits;

// Data members
  std::vector<bool> m_own_aabb_trees;
  std::vector<Tree*> m_aabb_trees;
  std::vector<bool> m_is_closed;
  std::vector< std::vector<typename K::Point_3> > m_points_per_cc;
  std::vector<Traversal_traits> m_traversal_traits;
  std::size_t m_free_id; // position in m_id_pool of the first free element
  std::vector<std::size_t> m_id_pool; // 0-> m_id_pool-1 are valid mesh ids
#if CGAL_RMCD_CACHE_BOXES
  boost::dynamic_bitset<> m_bboxes_is_invalid;
  std::vector<Bbox_3> m_bboxes;
#endif

// internal functions
  std::size_t get_id_for_new_mesh()
  {
    if (m_free_id==m_id_pool.size())
    {
      m_id_pool.push_back(m_free_id);
      ++m_free_id;
      m_own_aabb_trees.resize(m_free_id);
      m_aabb_trees.resize(m_free_id, nullptr);
      m_is_closed.resize(m_free_id);
      m_points_per_cc.resize(m_free_id);
      m_traversal_traits.resize(m_free_id);
  #if CGAL_RMCD_CACHE_BOXES
      m_bboxes.resize(m_free_id);
      m_bboxes_is_invalid.resize(m_free_id, true);
  #endif
      return m_id_pool.back();
    }
    return m_id_pool[m_free_id++];
  }

  template <class NamedParameters>
  void add_cc_points(const TriangleMesh& tm, std::size_t id, const NamedParameters& np)
  {
    collect_one_point_per_connected_component(tm, m_points_per_cc[id], np);
  }

  // precondition A and B does not intersect
  bool does_A_contains_a_CC_of_B(std::size_t id_A, std::size_t id_B) const
  {
    typename K::Construct_ray_3     ray_functor;
    typename K::Construct_vector_3  vector_functor;
    typedef typename Traversal_traits::Transformed_tree_helper Helper;

    for(const typename K::Point_3& q : m_points_per_cc[id_B])
    {
      if( internal::Point_inside_vertical_ray_cast<K, Tree, Helper>(m_traversal_traits[id_A].get_helper())(
            m_traversal_traits[id_B].transformation()( q ), *m_aabb_trees[id_A],
            ray_functor, vector_functor) == CGAL::ON_BOUNDED_SIDE)
      {
        return true;
      }
    }
    return false;
  }

  // this function expects a protector was initialized
  bool does_A_intersect_B(std::size_t id_A, std::size_t id_B) const
  {
#if CGAL_RMCD_CACHE_BOXES
    if (!do_overlap(m_bboxes[id_B], m_bboxes[id_A])) continue;
#endif

    Do_intersect_traversal_traits_for_two_trees<Tree_traits, K, Has_rotation> traversal_traits(
      m_aabb_trees[id_B]->traits(), m_traversal_traits[id_B].transformation(), m_traversal_traits[id_A]);
    m_aabb_trees[id_B]->traversal(*m_aabb_trees[id_A], traversal_traits);
    return traversal_traits.is_intersection_found();
  }

public:
#ifdef DOXYGEN_RUNNING
  /// The AABB_tree type representing the triangles of each input mesh
  typedef unspecified_type AABB_tree;
  /// The vertex point map type used with `TriangleMesh`
  typedef unspecified_type Vertex_point_map;
#else
  typedef Tree AABB_tree;
  typedef Vpm Vertex_point_map;
#endif
  /// Point type
  typedef typename boost::property_traits<Vertex_point_map>::value_type Point_3;

  Rigid_triangle_mesh_collision_detection()
    : m_free_id(0)
  {}

  ~Rigid_triangle_mesh_collision_detection()
  {
    for (std::size_t k=0; k<m_free_id; ++k)
    {
      std::size_t id = m_id_pool[k];
      if (m_own_aabb_trees[id]) delete m_aabb_trees[id];
    }
  }

 /*!
  * adds mesh `tm` to the set of meshes to be considered for intersection.
  *
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * \return the id of `tm` used to refer to that mesh.
  *
  * @param tm triangulated surface mesh to add
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
  *                     should be available for the vertices of `tm`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{face_index_map}
  *     \cgalParamDescription{a property map associating to each face of `tm` a unique index between `0` and `num_faces(tm) - 1`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and `std::size_t` as value type}
  *     \cgalParamDefault{an automatically indexed internal map}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{apply_per_connected_component}
  *     \cgalParamDescription{If `false`, `tm` is assumed to have only one connected component, avoiding the extraction of connected components.}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`true`}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  */
  template <class NamedParameters>
  std::size_t add_mesh(const TriangleMesh& tm,
                       const NamedParameters& np)
  {
    // handle vpm
    typedef typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Local_vpm;
    CGAL_USE_TYPE(Local_vpm);

    CGAL_assertion_code(
      static const bool same_vpm = (boost::is_same<Local_vpm,Vpm>::value); )
    CGAL_static_assertion(same_vpm);

    Vpm vpm =
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                   get_const_property_map(boost::vertex_point, tm) );
  // now add the mesh
    std::size_t id = get_id_for_new_mesh();
    CGAL_assertion( m_aabb_trees[id] == nullptr );
    m_is_closed[id] = is_closed(tm);
    m_own_aabb_trees[id] = true;
    Tree* t = new Tree(boost::begin(faces(tm)), boost::end(faces(tm)), tm, vpm);
    t->build();
    m_aabb_trees[id] = t;
    m_traversal_traits[id] = Traversal_traits(m_aabb_trees[id]->traits());
    add_cc_points(tm, id, np);

    return id;
  }

 /*!
  * adds an instance of a triangulated surface mesh using an external tree of its faces.
  * \warning The tree is not copied and the lifetime of `tree` must be longer than that of this class.
  *
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * \return the id of `tm` used to refer to that mesh.
  *
  * @param tree an AABB-tree of faces of a mesh
  * @param tm triangulated surface mesh
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
  *                     should be available for the vertices of `tm`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{face_index_map}
  *     \cgalParamDescription{a property map associating to each face of `tm` a unique index between `0` and `num_faces(tm) - 1`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and `std::size_t` as value type}
  *     \cgalParamDefault{an automatically indexed internal map}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{apply_per_connected_component}
  *     \cgalParamDescription{If `false`, `tm` is assumed to have only one connected component, avoiding the extraction of connected components}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`true`}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  */
  template <class NamedParameters>
  std::size_t add_mesh(const AABB_tree& tree, const TriangleMesh& tm, const NamedParameters& np)
  {
    std::size_t id = get_id_for_new_mesh();
    CGAL_assertion( m_aabb_trees[id] == nullptr );
    m_is_closed[id] = is_closed(tm);
    m_own_aabb_trees[id] = false ;
    m_aabb_trees[id] = const_cast<Tree*>(&tree);
    m_traversal_traits[id] = Traversal_traits(m_aabb_trees[id]->traits());
    collect_one_point_per_connected_component(tm, m_points_per_cc[id], np);
    return id;
  }

  /*!
   * sets the transformation associated to a mesh identified by its id in the set.
   */
  void set_transformation(std::size_t mesh_id, const Aff_transformation_3<K>& aff_trans)
  {
    CGAL_assertion(m_aabb_trees[mesh_id] != nullptr);
    m_traversal_traits[mesh_id].set_transformation(aff_trans);
#if CGAL_RMCD_CACHE_BOXES
    m_bboxes_is_invalid.set(mesh_id);
#endif
  }

#if CGAL_RMCD_CACHE_BOXES
  void update_bboxes()
  {
    // protector is supposed to have been set
    for (boost::dynamic_bitset<>::size_type i = m_bboxes_is_invalid.find_first();
                                            i != m_bboxes_is_invalid.npos;
                                            i = m_bboxes_is_invalid.find_next(i))
    {
      m_bboxes[i]=m_traversal_traits[i].get_helper().get_tree_bbox(*m_aabb_trees[i]);
    }
    m_bboxes_is_invalid.reset();
  }
#endif

 /*!
  * returns a vector of the ids of meshes within `ids` that have at least a face
  * intersecting a face of the mesh with id `mesh_id`.
  * If `mesh_id` is in `ids` it is not reported.
  * \tparam MeshIdRange a range of ids convertible to `std::size_t`.
  */
  template <class MeshIdRange>
  std::vector<std::size_t>
  get_all_intersections(std::size_t mesh_id, const MeshIdRange& ids) const
  {
    CGAL_assertion(m_aabb_trees[mesh_id] != nullptr);
    CGAL::Interval_nt_advanced::Protector protector;
#if CGAL_RMCD_CACHE_BOXES
    update_bboxes();
#endif
    std::vector<std::size_t> res;

    // TODO: use a non-naive version
    for(std::size_t k : ids)
    {
      CGAL_assertion(m_aabb_trees[k] != nullptr);
      if(k==mesh_id) continue;

      if (does_A_intersect_B(mesh_id, k))
        res.push_back(k);
    }
    return res;
  }

 /*!
  * returns a vector of the ids of meshes in the set that have at least a face
  * intersecting a face of the mesh with id `mesh_id`
  */
  std::vector<std::size_t>
  get_all_intersections(std::size_t mesh_id) const
  {
    return get_all_intersections(
      mesh_id,
      make_range(m_id_pool.begin(),m_id_pool.begin()+m_free_id) );
  }

 /*!
  * returns a vector of the ids of meshes within `ids` that are intersecting with the mesh with id `mesh_id`,
  * considering volume inclusions for closed meshes.
  * More precisely, if at least one face of a mesh with id `i` intersects a face
  * of the mesh with id `mesh_id`, the pair `(i, false)` is put in the output vector.
  * If there is no face intersection, but at least one of the meshes with ids `i` and `mesh_id` is closed,
  * and at least one connected component is included in the bounded volume defined by a closed mesh then the pair
  * `(i, true)` is put in the output vector (independently of mesh `i` or `mesh_id` being the one including the other).
  * The inclusion test is done using `Side_of_triangle_mesh`, in particular surface orientation is ignored and only the
  * nesting level of connected components defines a bounded volume. If a mesh has some self-intersection the inclusion
  * test may return incorrect results.
  * If `mesh_id` is in `ids` it is not reported.
  *
  * \tparam MeshIdRange a range of ids convertible to `std::size_t`.
  *
  * \note If a mesh is made of several connected components and at least one component is not closed,
  *       then no inclusion test will be made even if some components are closed.
  *
  */
  template <class MeshIdRange>
  std::vector<std::pair<std::size_t, bool> >
  get_all_intersections_and_inclusions(std::size_t mesh_id, const MeshIdRange& ids) const
  {
    CGAL_assertion(m_aabb_trees[mesh_id] != nullptr);
    CGAL::Interval_nt_advanced::Protector protector;
#if CGAL_RMCD_CACHE_BOXES
    update_bboxes();
#endif
    std::vector<std::pair<std::size_t, bool> > res;

    // TODO: use a non-naive version
    for(std::size_t k : ids)
    {
      CGAL_assertion(m_aabb_trees[k] != nullptr);
      if(k==mesh_id) continue;

      if (does_A_intersect_B(mesh_id, k))
        res.push_back(std::make_pair(k, false));
      else{
        if (m_is_closed[mesh_id])
        {
          if ( does_A_contains_a_CC_of_B(mesh_id, k) )
          {
            res.push_back(std::make_pair(k, true));
            continue;
          }
        }
        if (m_is_closed[k])
        {
          if ( does_A_contains_a_CC_of_B(k, mesh_id) )
          {
            res.push_back(std::make_pair(k, true));
            continue;
          }
        }
      }
    }
    return res;
  }

 /*!
  * returns a vector of the ids of meshes in the set that are intersecting with the mesh with id `mesh_id`,
  * considering volume inclusions for closed meshes.
  * See the previous overload for details.
  */
  std::vector<std::pair<std::size_t, bool> >
  get_all_intersections_and_inclusions(std::size_t mesh_id) const
  {
    return get_all_intersections_and_inclusions(mesh_id,
      make_range(m_id_pool.begin(),m_id_pool.begin()+m_free_id) );
  }

/// \name Memory Management
 /*!
  * increases the capacity of data structures used internally, `size` being the number of meshes expected to be added.
  */
  void reserve(std::size_t size)
  {
    m_own_aabb_trees.reserve(size);
    m_aabb_trees.reserve(size);
    m_is_closed.reserve(size);
    m_points_per_cc.reserve(size);
    m_traversal_traits.reserve(size);
#if CGAL_RMCD_CACHE_BOXES
    m_bboxes.reserve(size);
#endif
  }

 /*!
  * removes the mesh with id `mesh_id` from the set, the indices of other meshes are kept unchanged.
  */
  void remove_mesh(std::size_t mesh_id)
  {
    std::vector<std::size_t>::iterator itf =
      std::find(m_id_pool.begin(), m_id_pool.begin()+m_free_id, mesh_id);
    if (itf == m_id_pool.begin()+m_free_id) return;

    if (m_own_aabb_trees[mesh_id]) delete m_aabb_trees[mesh_id];
    m_points_per_cc[mesh_id].clear();
    m_aabb_trees[mesh_id] = nullptr;
    if (m_id_pool[m_free_id-1]!=mesh_id)
      std::swap(m_id_pool[m_free_id-1], *itf);
    --m_free_id;
  }

 /*!
  * returns the number of meshes in the set
  */
  std::size_t size() const
  {
    return m_free_id;
  }


#ifndef DOXYGEN_RUNNING
 /*!
  * returns the number of times `add_mesh()` was called minus the number of times `remove_mesh()` (with a valid mesh id) was called
  */
  std::size_t size_of_garbage()
  {
    return m_id_pool.size() - m_free_id;
  }

 /*!
  * returns `false` if `mesh_id` corresponds to an id not used (removed and/or not affected), and `true` otherwise.
  */
  bool is_valid_index(std::size_t mesh_id)
  {
    if (mesh_id >= m_id_pool.size()) return false;
    return m_aabb_trees[mesh_id] != nullptr;
  }

/// \name Helper Static Function

 /*!
  * fills `points` with one point per connected component of `tm`. This is a helper function
  * intended to be used before calling the `add_mesh()` overload taking an AABB-tree instead of a mesh
  * as input parameter.
  *
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param tm input triangulated surface mesh
  * @param [out] points will contain one point per connected component of `tm`
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
  *                     should be available for the vertices of `tm`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{face_index_map}
  *     \cgalParamDescription{a property map associating to each face of `tm` a unique index between `0` and `num_faces(tm) - 1`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and `std::size_t` as value type}
  *     \cgalParamDefault{an automatically indexed internal map}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{apply_per_connected_component}
  *     \cgalParamDescription{If `false`, `tm` is assumed to have only one connected component, avoiding the extraction of connected components}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`true`}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  */
  template <class NamedParameters>
  static
  void collect_one_point_per_connected_component(
    const TriangleMesh& tm,
          std::vector<Point_3>& points,
    const NamedParameters& np)
  {
    const bool maybe_several_cc =
      parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::apply_per_connected_component), true);

    typedef typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Local_vpm;
    CGAL_USE_TYPE(Local_vpm);

    CGAL_assertion_code(
      static const bool same_vpm = (boost::is_same<Local_vpm,Vpm>::value); )
    CGAL_static_assertion(same_vpm);

    Vpm vpm =
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                   get_const_property_map(boost::vertex_point, tm) );

    if (maybe_several_cc)
    {
      // used for face cc id map
      std::vector<std::size_t> cc_ids(num_faces(tm));

      // face index map
      typedef typename GetInitializedFaceIndexMap<TriangleMesh, NamedParameters>::const_type FaceIndexMap;
      FaceIndexMap fid_map = CGAL::get_initialized_face_index_map(tm, np);

      std::size_t nb_cc =
        Polygon_mesh_processing::connected_components(
          tm, bind_property_maps(fid_map, make_property_map(cc_ids)),
          parameters::face_index_map(fid_map));
      if (nb_cc != 1)
      {
        typedef boost::graph_traits<TriangleMesh> GrTr;
        std::vector<typename GrTr::vertex_descriptor>
          vertex_per_cc(nb_cc, GrTr::null_vertex());

        for(typename GrTr::face_descriptor f : faces(tm))
        {
          std::size_t cc_id = cc_ids[get(fid_map, f)];
          if  (vertex_per_cc[cc_id] == GrTr::null_vertex())
          {
            vertex_per_cc[cc_id] = target( halfedge(f, tm), tm);
            points.push_back( get(vpm, vertex_per_cc[cc_id]) );
          }
        }
        return;
      }
    }
    // only one CC
    points.push_back( get(vpm, *boost::begin(vertices(tm))) );
  }

 /*!
  * adds an instance of a triangulated surface mesh using an external tree of its faces.
  * \warning The tree is not copied and the lifetime of `tree` must be longer than that of this class.
  *
  * \return the id of the instance used to refer to that mesh.
  *
  * @param tree an AABB-tree of faces of a mesh
  * @param is_closed `true` is the mesh in `tree` is closed, and `false` otherwise.
  *                  \link is_closed() `CGAL::is_closed()` \endlink can be used for that purpose.
  * @param points_per_cc a vector containing one point of a vertex for each connected
  *                      component of the triangle surface mesh in `tree`
  */
  std::size_t add_mesh(const AABB_tree& tree,
                       bool is_closed,
                       const std::vector<Point_3>& points_per_cc)
  {
    std::size_t id = get_id_for_new_mesh();
    CGAL_assertion( m_aabb_trees[id] == nullptr );
    m_is_closed[id] = is_closed;
    m_own_aabb_trees[id] = false ;
    m_aabb_trees[id] = const_cast<Tree*>(&tree);
    m_traversal_traits[id] = Traversal_traits(m_aabb_trees[id]->traits());
    m_points_per_cc[id] = points_per_cc;

    return id;
  }

  // versions without NP
  static
  void collect_one_point_per_connected_component(
    const TriangleMesh& tm,
          std::vector<typename K::Point_3>& points)
  {
    collect_one_point_per_connected_component(tm, points, parameters::all_default());
  }

  std::size_t add_mesh(const TriangleMesh& tm)
  {
    return add_mesh(tm, parameters::all_default());
  }

  std::size_t add_mesh(const AABB_tree& tree, const TriangleMesh& tm)
  {
    return add_mesh(tree, tm, parameters::all_default());
  }
#endif
};

} // end of CGAL namespace

#undef CGAL_RMCD_CACHE_BOXES

#endif // CGAL_RIGID_TRIANGLE_MESH_COLLISION_DETECTION_H
