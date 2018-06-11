// Copyright (c) 2016 GeometryFactory (France).
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
// Author(s)     : Jane Tournois
//

#ifndef CGAL_LIPSCHITZ_SIZING_H
#define CGAL_LIPSCHITZ_SIZING_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/Mesh_3/experimental/AABB_filtered_projection_traits.h>
#include <CGAL/Mesh_3/experimental/Get_facet_patch_id.h>

#include <CGAL/Mesh_3/experimental/Lipschitz_sizing_parameters.h>

#include <CGAL/Default.h>
#include <CGAL/array.h>
#include <CGAL/Bbox_3.h>

#include <boost/shared_ptr.hpp>

#include <list>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

namespace CGAL
{
namespace Mesh_3
{

template <class Kernel, class MeshDomain
        , typename AABBTreeTemplate = CGAL::Default
#ifdef CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
        , typename Get_facet_patch_id_ = CGAL::Default
        , typename Patches_ids_ = CGAL::Default
#endif
>
class Lipschitz_sizing
{
public:
  typedef Kernel                      K;
  typedef typename Kernel::FT         FT;
  typedef typename Kernel::Triangle_3 Triangle;
  typedef typename Kernel::Point_3    Point_3;

  typedef typename std::list<Triangle>::iterator        Tr_iterator;
  typedef CGAL::AABB_triangle_primitive<K, Tr_iterator> Primitive;
  typedef CGAL::AABB_traits<K, Primitive>               AABB_tr_traits;
  typedef CGAL::AABB_tree<AABB_tr_traits>               AABB_tree;

  typedef typename CGAL::Default::Get<AABBTreeTemplate, AABB_tree>::type Tree;

  typedef typename MeshDomain::Index                    Index;
  typedef typename MeshDomain::Corner_index             Corner_index;
  typedef typename MeshDomain::Subdomain_index          Subdomain_index;
  typedef typename MeshDomain::Surface_patch_index      Surface_patch_index;

  typedef CGAL::Lipschitz_sizing_parameters<MeshDomain, FT> Parameters;

#ifdef CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
  typedef typename CGAL::Default::Get<Patches_ids_,
    typename MeshDomain::Surface_patch_index_set>::type Patches_ids;
  typedef std::vector<Patches_ids>                      Patches_ids_map;

  typedef typename CGAL::Default::Get<
    Get_facet_patch_id_,
    CGAL::Mesh_3::Get_facet_patch_id<typename Tree::Primitive>
  >::type                                               Get_facet_patch_id;

  typedef CGAL::Mesh_3::Filtered_projection_traits<
    typename Tree::AABB_traits, Get_facet_patch_id>     AABB_filtered_traits;

private:
  typedef CGAL::Search_traits_3<Kernel>                 KdTreeTraits;
  typedef CGAL::Orthogonal_k_neighbor_search<KdTreeTraits> Neighbor_search;
  typedef typename Neighbor_search::Tree                Kd_tree;
#endif

private:
  //only one of these aabb_trees is needed
  const Tree* m_ptree;
  boost::shared_ptr<Tree> m_own_ptree;

  const MeshDomain& m_domain;
  Parameters m_params;

  const CGAL::cpp11::array<double, 3>& m_vxyz;
  const CGAL::Bbox_3& m_bbox;
  const bool m_domain_is_a_box;

#ifdef CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
  //help to accelerate aabb_tree queries in m_ptree
  boost::shared_ptr<Kd_tree> m_kd_tree;

  Get_facet_patch_id m_get_facet_patch_id;
  const Patches_ids_map& patches_ids_map;
#endif

public:
  Lipschitz_sizing(const MeshDomain& domain)
    : m_ptree(NULL)
    , m_own_ptree()
    , m_domain(domain)
    , m_params(domain)
  {
  }

  Lipschitz_sizing(const MeshDomain& domain
    , const Tree* ptree
    , const CGAL::cpp11::array<double, 3>& vxyz
    , const CGAL::Bbox_3& bbox
    , const bool domain_is_a_box
#ifdef CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
    , const Patches_ids_map& patches_ids_map
#endif
    )
    : m_ptree(ptree)
    , m_own_ptree()
    , m_domain(domain)
    , m_params(domain)
    , m_vxyz(vxyz)
    , m_bbox(bbox)
    , m_domain_is_a_box(domain_is_a_box)
#ifdef CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
    , m_get_facet_patch_id()
    , patches_ids_map(patches_ids_map)
#endif
  {
  }

  FT operator()(const Point_3& p, const int dim, const Index& index) const
  {
    CGAL_assertion(!m_params.empty());
#ifdef CGAL_MESH_3_LIPSCHITZ_SIZING_VERBOSE
    std::cout << "D = " << dim << "\t";
#endif
    if (dim == 3)
    {
      return size_in_subdomain(p, m_domain.subdomain_index(index));
    }
    else if (dim == 2)
    {
#ifdef CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
      Surface_patch_index sp_index = m_domain.surface_patch_index(index);

      if(!is_on_cube_boundary(sp_index)
        && !is_on_cube_boundary(p))
      {
#ifdef CGAL_MESH_3_LIPSCHITZ_SIZING_VERBOSE
        std::cout << "          \n";
#endif
        FT size_max;
        m_params.get_parameters(sp_index, size_max);
        return size_max;
      }
      else
      {
#ifdef CGAL_MESH_3_LIPSCHITZ_SIZING_VERBOSE
        std::cout << "(on cube) ";
#endif
        const std::pair<Subdomain_index, Subdomain_index>& index
          = m_params.incident_subdomains(sp_index);

#ifdef CGAL_MESH_3_IMAGE_EXAMPLE
        if (index.first == INT_MIN)
          return size_in_subdomain(p, index.second);
        else
          return size_in_subdomain(p, index.first);
#else //==POLYHEDRAL EXAMPLE
        if (!is_in_domain(index.first))
          return size_in_subdomain(p, index.second);
        else
          return size_in_subdomain(p, index.first);
#endif //CGAL_MESH_3_IMAGE_EXAMPLE
      }

#else  //CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
      CGAL_assertion(false);
      return 0.;
#endif //CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
    }

#ifndef CGAL_MESH_3_IMAGE_EXAMPLE
    else if (dim == 1)
    {
#ifdef CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
      const typename MeshDomain::Curve_index& curve_id =
        m_domain.curve_segment_index(index);
      const Patches_ids& ids = patches_ids_map[curve_id];
      
      if (m_domain_is_a_box && ids.size() == 2)
      {
        //we are on an edge of the box
        //same code as when dim == 2
        Surface_patch_index spi = *(ids.begin());
        const std::pair<Subdomain_index, Subdomain_index>& subdomains
          = m_params.incident_subdomains(spi);
        if (!is_in_domain(subdomains.first))
          return size_in_subdomain(p, subdomains.second);
        else
          return size_in_subdomain(p, subdomains.first);
      }
      return min_size_in_incident_subdomains(ids);
#else
      CGAL_assertion(false);//should not be used for dimension 1
      return 0.;
#endif
    }
    else if (dim == 0)
    {
#ifdef CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
      const Corner_index cid = m_domain.corner_index(index);
      const Patches_ids& ids = m_domain.corners_incidences_map().find(cid)->second;

      if (m_domain_is_a_box && ids.size() == 3)
      {
        //we are on a corner of the box
        //same code as when dim == 2
        Surface_patch_index spi = *(ids.begin());
        const std::pair<Subdomain_index, Subdomain_index>& subdomains
          = m_params.incident_subdomains(spi);
        if (!is_in_domain(subdomains.first))
          return size_in_subdomain(p, subdomains.second);
        else
          return size_in_subdomain(p, subdomains.first);
      }

      return min_size_in_incident_subdomains(ids);;
#else
      CGAL_assertion(false);//should not be used for dimension 0
      return 0;
#endif
    }
#endif

    CGAL_assertion(false);
    return 0.;
  }

private:
#ifdef CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
  std::vector<Subdomain_index> incident_subdomains(const Patches_ids& ids) const
  {
    std::vector<Subdomain_index> vec;
    BOOST_FOREACH(Surface_patch_index spi, ids)
    {
      const std::pair<Subdomain_index, Subdomain_index>& subdomains
        = m_params.incident_subdomains(spi);

      if (is_in_domain(subdomains.first))
        vec.push_back(subdomains.first);

      if (is_in_domain(subdomains.second))
        vec.push_back(subdomains.second);
    }
    return vec;
  }

  FT min_size_in_incident_subdomains(const Patches_ids& ids) const
  {
    FT size = static_cast<FT>((std::numeric_limits<double>::max)());
    BOOST_FOREACH(Surface_patch_index spi, ids)
    {
      const std::pair<Subdomain_index, Subdomain_index>& subdomains
        = m_params.incident_subdomains(spi);

      FT k, size_min, size_max;
      if (is_in_domain(subdomains.first))
      {
        m_params.get_parameters(subdomains.first, k, size_min, size_max);
        size = (std::min)(size, size_min);
      }
      if (is_in_domain(subdomains.second))
      {
        m_params.get_parameters(subdomains.second, k, size_min, size_max);
        size = (std::min)(size, size_min);
      }
    }
    return size;
  }
#endif

public:
  template <typename C3T3>
  void init_aabb_tree_from_c3t3(const C3T3* p_c3t3)
  {
    static std::list<Triangle> triangles;
    for (typename C3T3::Facets_in_complex_iterator
      fit = p_c3t3->facets_in_complex_begin();
      fit != p_c3t3->facets_in_complex_end();
      ++fit)
    {
      if (!is_on_cube_boundary(*fit))
        triangles.push_back(p_c3t3->triangulation().triangle(*fit));
    }

    m_own_ptree.reset(new Tree(triangles.begin(), triangles.end()));
    m_own_ptree->build();
    m_own_ptree->accelerate_distance_queries();
  }

private:
  template <typename Facet>
  bool is_on_cube_boundary(const Facet& f) const
  {
    return is_on_cube_boundary(f.first->surface_patch_index(f.second));
  }

  bool is_on_cube_boundary(const Surface_patch_index& sp_index) const
  {
    const std::pair<Subdomain_index, Subdomain_index>& index
      = m_params.incident_subdomains(sp_index);

#ifdef CGAL_MESH_3_IMAGE_EXAMPLE
    return (index.first == INT_MIN || index.second == INT_MIN);
#else //POLYHEDRAL EXAMPLE

    if (m_domain_is_a_box)
      return !is_in_domain(index.first) || !is_in_domain(index.second);
    else
      return false;

#endif //CGAL_MESH_3_IMAGE_EXAMPLE
  }

  bool is_on_cube_boundary(const Point_3& p) const
  {
#ifndef CGAL_MESH_3_IMAGE_EXAMPLE//POLYHEDRAL EXAMPLE

    if (m_domain_is_a_box)
      //checks that p is in the outer 'shell' of voxels
      return p.x() < m_bbox.xmin() + m_vxyz[0]
        || p.x() > m_bbox.xmax() - m_vxyz[0]
        || p.y() < m_bbox.ymin() + m_vxyz[1]
        || p.y() > m_bbox.ymax() - m_vxyz[1]
        || p.z() < m_bbox.zmin() + m_vxyz[2]
        || p.z() > m_bbox.zmax() - m_vxyz[2];
    else
      return false;

#else //IMAGE EXAMPLE
    CGAL_USE(p);
    return false;
#endif //IMAGE_EXAMPLE
  }

  bool is_in_domain(const Subdomain_index& index) const
  {
#ifdef CGAL_MESH_3_IMAGE_EXAMPLE
    return (index != 0 && index != INT_MIN);
#else //POLYHEDRAL EXAMPLE
    return (index != 0);
#endif
  }

  FT size_in_subdomain(const Point_3& p, const Subdomain_index& index) const
  {
    FT k, size_min, size_max;
    m_params.get_parameters(index, k, size_min, size_max);

    FT sqdist = 0.;
    if(m_ptree == NULL)
    {
      sqdist = m_own_ptree->squared_distance(p);
    }
    else
    {
      Point_3 closest = compute_closest_point(p);
      sqdist = CGAL::squared_distance(p, closest);
    }

    FT size = k * CGAL::sqrt(sqdist) + size_min;
    return (std::min)(size, size_max);
  }

#ifdef CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
  void kd_tree()
  {
    typedef typename MeshDomain::Polyhedron Polyhedron;
    if(m_kd_tree.get() == 0) {
      m_kd_tree.reset(new Kd_tree);
      BOOST_FOREACH(std::size_t poly_id, m_domain.inside_polyhedra()) {
        const Polyhedron& poly = m_domain.polyhedra()[poly_id];
        BOOST_FOREACH(typename Polyhedron::Vertex_handle v, vertices(poly))
        {
          m_kd_tree->insert(v->point());
        }
      }
      BOOST_FOREACH(std::size_t poly_id, m_domain.boundary_polyhedra()) {
        const Polyhedron& poly = m_domain.polyhedra()[poly_id];
        BOOST_FOREACH(typename Polyhedron::Vertex_handle v, vertices(poly))
        {
          if(!is_on_cube_boundary(v->point()))
            m_kd_tree->insert(v->point());
        }
      }
      m_kd_tree->build();
    }
  }
#endif

  Point_3 compute_closest_point(const Point_3& p) const
  {
#ifndef CGAL_MESH_3_IMAGE_EXAMPLE //POLYHEDRAL_EXAMPLE

#ifdef CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
    const std::vector<Surface_patch_index>& boundary_ids =
      m_domain.boundary_patches();

    CGAL_STATIC_THREAD_LOCAL_VARIABLE_4(AABB_filtered_traits,
                                        projection_traits,
                                        boundary_ids.begin(),
                                        boundary_ids.end(),
                                        m_ptree->traits(),
                                        m_get_facet_patch_id);
    kd_tree();//build it if needed
    Neighbor_search search(*m_kd_tree, p, 1);
    projection_traits.reset(search.begin()->first);

    m_ptree->traversal(p, projection_traits);
    CGAL_assertion(projection_traits.found());
    return projection_traits.closest_point();

#else //CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS
    return m_ptree->closest_point(p);
#endif //CGAL_MESH_3_EXPERIMENTAL_USE_PATCHES_IDS

#else
    CGAL_assertion(false);
    return CGAL::ORIGIN;//not used
#endif
  }

public:
  void add_parameters_for_subdomain(const Subdomain_index& id
                                  , const FT& k
                                  , const FT& size_min
                                  , const FT& size_max)
  {
    m_params.add_subdomain(id, k, size_min, size_max);
  }

};

}//namespace Mesh_3
}//namespace CGAL

#endif // CGAL_LIPSCHITZ_SIZING_H
