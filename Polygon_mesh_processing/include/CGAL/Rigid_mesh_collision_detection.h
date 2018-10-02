// Copyright (c) 2018 GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Maxime Gimeno and Sebastien Loriot


#ifndef CGAL_RIGID_MESH_COLLISION_DETECTION_H
#define CGAL_RIGID_MESH_COLLISION_DETECTION_H

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

template <class TriangleMesh,
          class VertexPointMap = Default,
          class Kernel_ = Default,
          class HAS_ROTATION = CGAL::Tag_true,
          class AABBTree_type = Default>
class Rigid_mesh_collision_detection
{
// Vertex point map type
  typedef typename property_map_selector<TriangleMesh, boost::vertex_point_t
    >::const_type                                                   Default_vpm;
  typedef typename Default::Get<VertexPointMap, Default_vpm>::type          Vpm;

// Kernel type
  typedef typename Kernel_traits<
    typename boost::property_traits<Vpm>::value_type>::Kernel    Default_kernel;
  typedef typename Default::Get<Kernel_, Default_kernel>::type           Kernel;

// AABB-tree type
  typedef AABB_face_graph_triangle_primitive<TriangleMesh,
                                             Vpm>             Default_primitive;
  typedef AABB_traits<Kernel, Default_primitive>            Default_tree_traits;
  typedef CGAL::AABB_tree<Default_tree_traits>                     Default_tree;
  typedef typename Default::Get<AABBTree_type, Default_tree>::type         Tree;
  typedef typename Tree::AABB_traits                                Tree_traits;

// Transformed Tree traversal traits
  typedef Do_intersect_traversal_traits_with_transformation<Tree_traits,
                                                            Kernel,
                                                            HAS_ROTATION>
                                                               Traversal_traits;

// Data members
  std::vector<bool> m_own_aabb_trees;
  std::vector<Tree*> m_aabb_trees;
  std::vector<bool> m_is_closed;
  std::vector< std::vector<typename Kernel::Point_3> > m_points_per_cc;
  std::vector<Traversal_traits> m_traversal_traits;
#if CGAL_RMCD_CACHE_BOXES
  boost::dynamic_bitset<> m_bboxes_is_invalid;
  std::vector<Bbox_3> m_bboxes;
#endif

  void clear_trees()
  {
    int i=0;
    BOOST_FOREACH(Tree* tree, m_aabb_trees){
      if (m_own_aabb_trees[i]) delete tree;
      ++i;
    }
    m_aabb_trees.clear();
  }

public:
  typedef Tree AABB_tree;
  typedef Vpm Vertex_point_map;

  static
  void get_one_point_per_cc(const TriangleMesh& tm,
                                  Vertex_point_map vpm,
                                  std::vector<typename Kernel::Point_3>& points,
                            const bool assume_one_CC=false)
  {
    if (!assume_one_CC)
    {
      std::vector<std::size_t> CC_ids(num_faces(tm));

      // TODO use dynamic property if no defaut fid is available
      typename boost::property_map<TriangleMesh, boost::face_index_t>::type fid_map
        = get(boost::face_index, tm);

      std::size_t nb_cc =
        Polygon_mesh_processing::connected_components(
          tm, bind_property_maps(fid_map, make_property_map(CC_ids)) );
      if (nb_cc != 1)
      {
        typedef boost::graph_traits<TriangleMesh> GrTr;
        std::vector<typename GrTr::vertex_descriptor> vertex_per_cc(nb_cc, GrTr::null_vertex());

        BOOST_FOREACH(typename GrTr::face_descriptor f, faces(tm))
        {
          if  (vertex_per_cc[get(fid_map, f)]!=GrTr::null_vertex())
          {
            points.push_back(
              get(vpm, target( halfedge(f, tm), tm)) );
          }
        }
        return;
      }
    }
    // only one CC
    points.push_back( get(boost::vertex_point, tm, *boost::begin(vertices(tm))) );
  }

  static
  void get_one_point_per_cc(const TriangleMesh& tm,
                                  std::vector<typename Kernel::Point_3>& points,
                            const bool assume_one_CC=false)
  {
    get_one_point_per_cc(tm, get(boost::vertex_point, tm), points, assume_one_CC);
  }

private:
  void add_cc_points(const TriangleMesh& tm, bool assume_one_CC)
  {
    m_points_per_cc.resize(m_points_per_cc.size()+1);
    get_one_point_per_cc(tm, m_points_per_cc.back(), assume_one_CC);
  }

  // precondition A and B does not intersect
  bool does_A_contains_a_CC_of_B(std::size_t id_A, std::size_t id_B)
  {
    typename Kernel::Construct_ray_3     ray_functor;
    typename Kernel::Construct_vector_3  vector_functor;
    typedef typename Traversal_traits::Transformed_tree_helper Helper;

    BOOST_FOREACH(const typename Kernel::Point_3& q, m_points_per_cc[id_B])
    {
      if( internal::Point_inside_vertical_ray_cast<Kernel, Tree, Helper>(m_traversal_traits[id_A].get_helper())(
            m_traversal_traits[id_B].transformation()( q ), *m_aabb_trees[id_A],
            ray_functor, vector_functor) == CGAL::ON_BOUNDED_SIDE)
      {
        return true;
      }
    }
    return false;
  }

public:

  // TODO: document that the default vertex point map will be used
  template <class MeshRange>
  Rigid_mesh_collision_detection(const MeshRange& triangle_meshes, bool assume_one_CC_per_mesh = false)
  {
    init(triangle_meshes, assume_one_CC_per_mesh);
  }

  ~Rigid_mesh_collision_detection()
  {
    clear_trees();
  }

  // TODO: document that the default vertex point map will be used
  template <class MeshRange>
  void init(const MeshRange& triangle_meshes, bool assume_one_CC)
  {
    std::size_t nb_meshes = triangle_meshes.size();
    m_own_aabb_trees.clear();
    m_own_aabb_trees.reserve(nb_meshes);
    m_points_per_cc.clear();
    m_points_per_cc.reserve(nb_meshes);
    clear_trees();
    m_aabb_trees.reserve(nb_meshes);
    m_is_closed.clear();
    m_is_closed.resize(nb_meshes, false);
    m_traversal_traits.reserve(m_aabb_trees.size());
#if CGAL_RMCD_CACHE_BOXES
    m_bboxes_is_invalid.clear();
    m_bboxes_is_invalid.resize(nb_meshes, true);
    m_bboxes.clear();
    m_bboxes.resize(nb_meshes);
#endif
    BOOST_FOREACH(const TriangleMesh& tm, triangle_meshes)
    {
      if (is_closed(tm))
        m_is_closed[m_aabb_trees.size()]=true;
      m_own_aabb_trees.push_back( true );
      Tree* t = new Tree(faces(tm).begin(), faces(tm).end(), tm);
      m_aabb_trees.push_back(t);
      m_traversal_traits.push_back( Traversal_traits(m_aabb_trees.back()->traits()) );
      add_cc_points(tm, assume_one_CC);
    }
  }

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

  std::size_t add_mesh(const TriangleMesh& tm,
                             Vertex_point_map vpm,
                       const bool assume_one_CC_per_mesh = false)
  {
    std::size_t id = m_aabb_trees.size();
    m_is_closed.push_back(is_closed(tm));
    m_own_aabb_trees.push_back( true );
    Tree* t = new Tree(faces(tm).begin(), faces(tm).end(), tm, vpm);
    m_aabb_trees.push_back(t);
    m_traversal_traits.push_back( Traversal_traits(m_aabb_trees.back()->traits()) );
#if CGAL_RMCD_CACHE_BOXES
    m_bboxes.push_back(Bbox_3());
    m_bboxes_is_invalid.resize(id+1, true);
#endif
    add_cc_points(tm, assume_one_CC_per_mesh);

    return id;
  }

  std::size_t add_mesh(const TriangleMesh& tm,
                       const bool assume_one_CC_per_mesh = false)
  {
    return add_mesh(tm, get(boost::vertex_point, tm), assume_one_CC_per_mesh);
  }

  std::size_t add_mesh(const AABB_tree& tree,
                       bool is_closed,
                       const std::vector<typename Kernel::Point_3>& points_per_cc)
  {
    std::size_t id = m_aabb_trees.size();
    m_is_closed.push_back(is_closed);
    m_own_aabb_trees.push_back( false );
    m_aabb_trees.push_back( const_cast<Tree*>(&tree));
    m_traversal_traits.push_back( Traversal_traits(m_aabb_trees.back()->traits()) );
#if CGAL_RMCD_CACHE_BOXES
    m_bboxes.push_back(Bbox_3());
    m_bboxes_is_invalid.resize(id+1, true);
#endif
    m_points_per_cc.push_back(points_per_cc);

    return id;
  }

  void remove_mesh(std::size_t mesh_id)
  {
    if(mesh_id >= m_aabb_trees.size()) return;
    if (m_own_aabb_trees[mesh_id]) delete m_aabb_trees[mesh_id];
    m_own_aabb_trees.erase( m_own_aabb_trees.begin()+mesh_id );
    m_aabb_trees.erase( m_aabb_trees.begin()+mesh_id);
    m_is_closed.erase(m_is_closed.begin()+mesh_id);
    m_points_per_cc.erase(m_points_per_cc.begin()+mesh_id);
    m_traversal_traits.erase(m_traversal_traits.begin()+mesh_id);
#if CGAL_RMCD_CACHE_BOXES
    // TODO this is a lazy approach that is not optimal
    m_bboxes.pop_back();
    m_bboxes_is_invalid.set();
    m_bboxes_is_invalid.resize(m_aabb_trees.size());
#endif
  }

  void set_transformation(std::size_t mesh_id, const Aff_transformation_3<Kernel>& aff_trans)
  {
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

  template <class MeshRangeIds>
  std::vector<std::size_t>
  get_all_intersections(std::size_t mesh_id, const MeshRangeIds& ids)
  {
    CGAL::Interval_nt_advanced::Protector protector;
#if CGAL_RMCD_CACHE_BOXES
    update_bboxes();
#endif
    std::vector<std::size_t> res;

    // TODO: use a non-naive version
    BOOST_FOREACH(std::size_t k, ids)
    {
      if(k==mesh_id) continue;
#if CGAL_RMCD_CACHE_BOXES
      if (!do_overlap(m_bboxes[k], m_bboxes[mesh_id])) continue;
#endif

      Do_intersect_traversal_traits_for_two_trees<Tree_traits, Kernel, HAS_ROTATION> traversal_traits(
        m_aabb_trees[k]->traits(), m_traversal_traits[k].transformation(), m_traversal_traits[mesh_id]);
      m_aabb_trees[k]->traversal(*m_aabb_trees[mesh_id], traversal_traits);
      if (traversal_traits.is_intersection_found())
        res.push_back(k);
    }
    return res;
  }

  std::vector<std::size_t>
  get_all_intersections(std::size_t mesh_id)
  {
    return get_all_intersections(
      mesh_id,
      make_range(boost::make_counting_iterator<std::size_t>(0),
                 boost::make_counting_iterator<std::size_t>(m_aabb_trees.size())));
  }

  std::vector<std::size_t>
  set_transformation_and_get_all_intersections(std::size_t mesh_id,
                                               const Aff_transformation_3<Kernel>& aff_trans)
  {
    CGAL::Interval_nt_advanced::Protector protector;
    set_transformation(mesh_id, aff_trans);
    return get_all_intersections(mesh_id);
  }

  // TODO: document that if a model is composed of several CC on one of them is not closed,
  // no inclusion test will be made
  // TODO: document that the inclusion can be partial in case there are several CC
  template <class MeshRangeIds>
  std::vector<std::pair<std::size_t, bool> >
  get_all_intersections_and_inclusions(std::size_t mesh_id, const MeshRangeIds& ids)
  {
    CGAL::Interval_nt_advanced::Protector protector;
#if CGAL_RMCD_CACHE_BOXES
    update_bboxes();
#endif
    std::vector<std::pair<std::size_t, bool> > res;

    // TODO: use a non-naive version
    BOOST_FOREACH(std::size_t k, ids)
    {
      if(k==mesh_id) continue;
#if CGAL_RMCD_CACHE_BOXES
      if (!do_overlap(m_bboxes[k], m_bboxes[mesh_id])) continue;
#endif

      Do_intersect_traversal_traits_for_two_trees<Tree_traits, Kernel, HAS_ROTATION> traversal_traits(
        m_aabb_trees[k]->traits(), m_traversal_traits[k].transformation(), m_traversal_traits[mesh_id]);
      m_aabb_trees[k]->traversal(*m_aabb_trees[mesh_id], traversal_traits);
      if (traversal_traits.is_intersection_found())
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

  std::vector<std::pair<std::size_t, bool> >
  get_all_intersections_and_inclusions(std::size_t mesh_id)
  {
    return get_all_intersections_and_inclusions(
      mesh_id,
      make_range(boost::make_counting_iterator<std::size_t>(0),
                 boost::make_counting_iterator<std::size_t>(m_aabb_trees.size())));
  }

  std::vector<std::pair<std::size_t, bool> >
  set_transformation_and_get_all_intersections_and_inclusions(std::size_t mesh_id,
                                           const Aff_transformation_3<Kernel>& aff_trans)
  {
    CGAL::Interval_nt_advanced::Protector protector;
    set_transformation(mesh_id, aff_trans);
    return get_all_intersections_and_inclusions(mesh_id);
  }
};

} // end of CGAL namespace

#undef CGAL_RMCD_CACHE_BOXES

#endif // CGAL_RIGID_MESH_COLLISION_DETECTION_H
