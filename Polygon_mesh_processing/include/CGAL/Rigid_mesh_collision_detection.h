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
#include <CGAL/Polygon_mesh_processing/internal/AABB_do_intersect_transform_traits.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/property_map.h>

#include <boost/iterator/counting_iterator.hpp>

#ifndef CGAL_CACHE_BOXES
#define CGAL_CACHE_BOXES 0
#endif

#if CGAL_CACHE_BOXES
#include <boost/dynamic_bitset.hpp>
#endif

namespace CGAL {

//TODO handle vertex point point in the API
template <class TriangleMesh, class Kernel, class HAS_ROTATION = CGAL::Tag_true>
class Rigid_mesh_collision_detection
{
  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh> Primitive;
  typedef CGAL::AABB_do_intersect_transform_traits<Kernel, Primitive, HAS_ROTATION> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;
  typedef Side_of_triangle_mesh<TriangleMesh, Kernel, Default, Tree> Side_of_tm;


  std::vector<const TriangleMesh*> m_triangle_mesh_ptrs;
  // TODO: we probably want an option with external trees
  std::vector<Tree*> m_aabb_trees;
  std::vector<bool> m_is_closed;
  std::vector< std::vector<typename Kernel::Point_3> > m_points_per_cc;

#if CGAL_CACHE_BOXES
  boost::dynamic_bitset<> m_bboxes_is_invalid;
  std::vector<Bbox_3> m_bboxes;
#endif

  void clear_trees()
  {
    BOOST_FOREACH(Tree* tree, m_aabb_trees){
      delete tree;
    }
    m_aabb_trees.clear();
  }

  void add_cc_points(const TriangleMesh& tm, bool assume_one_CC)
  {
    m_points_per_cc.resize(m_points_per_cc.size()+1);
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
        typedef boost::graph_traits<TriangleMesh> GrT;
        std::vector<typename GrT::vertex_descriptor> vertex_per_cc(nb_cc, GrT::null_vertex());

        BOOST_FOREACH(typename GrT::face_descriptor f, faces(tm))
        {
          if  (vertex_per_cc[get(fid_map, f)]!=GrT::null_vertex())
          {
            m_points_per_cc.back().push_back(
              get(boost::vertex_point, tm, target( halfedge(f, tm), tm)) );
          }
        }
        return;
      }
    }
    // only one CC
    m_points_per_cc.back().push_back( get(boost::vertex_point, tm, *boost::begin(vertices(tm))) );
  }

public:
  template <class MeshRange>
  Rigid_mesh_collision_detection(const MeshRange& triangle_meshes, bool assume_one_CC_per_mesh = false)
  {
    init(triangle_meshes, assume_one_CC_per_mesh);
  }

  ~Rigid_mesh_collision_detection()
  {
    clear_trees();
  }
  template <class MeshRange>
  void init(const MeshRange& triangle_meshes, bool assume_one_CC)
  {
    std::size_t nb_meshes = triangle_meshes.size();
    m_triangle_mesh_ptrs.clear();
    m_triangle_mesh_ptrs.reserve(nb_meshes);
    m_points_per_cc.clear();
    m_points_per_cc.reserve(nb_meshes);
    clear_trees();
    m_aabb_trees.reserve(nb_meshes);
    m_is_closed.clear();
    m_is_closed.resize(nb_meshes, false);
#if CGAL_CACHE_BOXES
    m_bboxes_is_invalid.clear();
    m_bboxes_is_invalid.resize(nb_meshes, true);
    m_bboxes.clear();
    m_bboxes.resize(nb_meshes);
#endif
    BOOST_FOREACH(const TriangleMesh& tm, triangle_meshes)
    {
      if (is_closed(tm))
        m_is_closed[m_triangle_mesh_ptrs.size()]=true;
      m_triangle_mesh_ptrs.push_back( &tm );
      Tree* t = new Tree(faces(tm).begin(), faces(tm).end(), tm);
      m_aabb_trees.push_back(t);
      add_cc_points(tm, assume_one_CC);
    }
  }

  void add_mesh(const TriangleMesh& tm, bool assume_one_CC_per_mesh = false)
  {
    m_is_closed.push_back(is_closed(tm));
    m_triangle_mesh_ptrs.push_back( &tm );
    Tree* t = new Tree(faces(tm).begin(), faces(tm).end(), tm);
    m_aabb_trees.push_back(t);
#if CGAL_CACHE_BOXES
    m_bboxes.push_back(Bbox_3());
    m_bboxes_is_invalid.resize(m_bboxes_is_invalid.size()+1, true);
#endif
    add_cc_points(tm, assume_one_CC_per_mesh);
  }

  void remove_mesh(std::size_t mesh_id)
  {
    if(mesh_id >= m_triangle_mesh_ptrs.size()) return;
    m_triangle_mesh_ptrs.erase( m_triangle_mesh_ptrs.begin()+mesh_id );
    delete m_aabb_trees[mesh_id];
    m_aabb_trees.erase( m_aabb_trees.begin()+mesh_id);
    m_is_closed.erase(m_is_closed.begin()+mesh_id);
    m_points_per_cc.erase(m_points_per_cc.begin()+mesh_id);
#if CGAL_CACHE_BOXES
    // TODO this is a lazy approach that is not optimal
    m_bboxes.pop_back();
    m_bboxes_is_invalid.set();
    m_bboxes_is_invalid.resize(m_triangle_mesh_ptrs.size());
#endif
  }

  void set_transformation(std::size_t mesh_id, const Aff_transformation_3<Kernel>& aff_trans)
  {
    m_aabb_trees[mesh_id]->traits().set_transformation(aff_trans);
#if CGAL_CACHE_BOXES
    m_bboxes_is_invalid.set(mesh_id);
#endif
  }

#if CGAL_CACHE_BOXES
  void update_bboxes()
  {
    // protector is supposed to have been set
    for (boost::dynamic_bitset<>::size_type i = m_bboxes.find_first();
                                            i != m_bboxes_is_invalid.npos;
                                            i = m_bboxes_is_invalid.find_next(i))
    {
      m_bboxes[i]=internal::get_tree_bbox(*m_aabb_trees[i]);
    }
    m_bboxes_is_invalid.reset();
  }
#endif

  template <class MeshRangeIds>
  std::vector<std::size_t>
  get_all_intersections(std::size_t mesh_id, const MeshRangeIds& ids)
  {
    CGAL::Interval_nt_advanced::Protector protector;
#if CGAL_CACHE_BOXES
    update_bboxes();
#endif
    std::vector<std::size_t> res;

    BOOST_FOREACH(std::size_t k, ids)
    {
      if(k==mesh_id) continue;
#if CGAL_CACHE_BOXES
       if (!do_overlap(m_bboxes[k], m_bboxes[mesh_id])) continue;
#endif
      // TODO: think about an alternative that is using a traversal traits
      if ( m_aabb_trees[k]->do_intersect( *m_aabb_trees[mesh_id] ) )
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

  // TODO: document that is a model is composed of several CC on one of them is not closed,
  // no inclusion test will be made
  // TODO: document that the inclusion can be partial in case there are several CC
  template <class MeshRangeIds>
  std::vector<std::pair<std::size_t, bool> >
  get_all_intersections_and_inclusions(std::size_t mesh_id, const MeshRangeIds& ids)
  {
    CGAL::Interval_nt_advanced::Protector protector;
#if CGAL_CACHE_BOXES
    update_bboxes();
#endif
    std::vector<std::pair<std::size_t, bool> > res;

    // TODO: use a non-naive version
    BOOST_FOREACH(std::size_t k, ids)
    {
      if(k==mesh_id) continue;
#if CGAL_CACHE_BOXES
      if (!do_overlap(m_bboxes[k], m_bboxes[mesh_id])) continue;
#endif
      // TODO: think about an alternative that is using a traversal traits
      if ( m_aabb_trees[k]->do_intersect( *m_aabb_trees[mesh_id] ) )
        res.push_back(std::make_pair(k, false));
      else{
        if (m_is_closed[mesh_id])
        {
          Side_of_tm side_of_mid(*m_aabb_trees[mesh_id]);
          typename Kernel::Point_3 q = get(boost::vertex_point, *m_triangle_mesh_ptrs[k], *boost::begin(vertices(*m_triangle_mesh_ptrs[k])));
          if(side_of_mid(m_aabb_trees[k]->traits().transformation()( q )) == CGAL::ON_BOUNDED_SIDE)
          {
            res.push_back(std::make_pair(k, true));
            continue;
          }
        }
        if (m_is_closed[k])
        {
          Side_of_tm side_of_mk(*m_aabb_trees[k]);
          BOOST_FOREACH(const typename Kernel::Point_3 q, m_points_per_cc[mesh_id])
          {
            if(side_of_mk(m_aabb_trees[mesh_id]->traits().transformation()( q )) == CGAL::ON_BOUNDED_SIDE)
            {
              res.push_back(std::make_pair(k, true));
              break;
            }
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


#endif // CGAL_RIGID_MESH_COLLISION_DETECTION_H
