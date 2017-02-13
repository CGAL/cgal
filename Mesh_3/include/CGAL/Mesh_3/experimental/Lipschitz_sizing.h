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
//
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************


#ifndef _LIPSCHITZ_SIZING_
#define _LIPSCHITZ_SIZING_

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/Mesh_3/experimental/AABB_filtered_projection_traits.h>
#include <CGAL/Mesh_3/experimental/Get_facet_patch_id.h>

#include <CGAL/Mesh_3/Lipschitz_sizing_parameters.h>

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

template <class Kernel, class C3T3, class MeshDomain
        , typename AABBTreeTemplate = CGAL::Default
        , typename Get_facet_patch_id_ = CGAL::Default
        , typename Patches_ids_ = CGAL::Default
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

  typedef typename C3T3::Triangulation                  Tr;
  typedef typename C3T3::Triangulation::Facet           Facet;

  typedef typename MeshDomain::Index                    Index;
  typedef typename MeshDomain::Subdomain_index          Subdomain_index;
  typedef typename MeshDomain::Surface_patch_index      Surface_patch_index;

  typedef typename CGAL::Default::Get<Patches_ids_,
    typename MeshDomain::Surface_patch_index_set>::type Patches_ids;
  typedef std::vector<Patches_ids>                      Patches_ids_map;

  typedef typename CGAL::Default::Get<
    Get_facet_patch_id_,
    CGAL::Mesh_3::Get_facet_patch_id<typename Tree::Primitive>
  >::type                                               Get_facet_patch_id;

  typedef CGAL::Mesh_3::Filtered_projection_traits<
    typename Tree::AABB_traits, Get_facet_patch_id>     AABB_filtered_traits;

  typedef CGAL::Lipschitz_sizing_parameters<MeshDomain, FT> Parameters;

private:
  typedef CGAL::Search_traits_3<Kernel>                 KdTreeTraits;
  typedef CGAL::Orthogonal_k_neighbor_search<KdTreeTraits> Neighbor_search;
  typedef typename Neighbor_search::Tree                Kd_tree;

private:
  //only one of these aabb_trees is needed
  const Tree* m_ptree;
  boost::shared_ptr<Tree> m_own_ptree;

  //help to accelerate aabb_tree queries in m_ptree
  boost::shared_ptr<Kd_tree> m_kd_tree;

  C3T3* m_pc3t3;
  const MeshDomain& m_domain;
  Parameters m_params;

#ifndef CGAL_MESH_3_IMAGE_EXAMPLE//POLYHEDRAL EXAMPLE
  Get_facet_patch_id m_get_facet_patch_id;
  const Patches_ids_map& patches_ids_map;
  const CGAL::cpp11::array<double, 3>& m_vxyz;
  const CGAL::Bbox_3& m_bbox;
  const bool m_domain_is_a_box;
#endif

public:
#ifdef CGAL_MESH_3_IMAGE_EXAMPLE
  Lipschitz_sizing(C3T3* pc3t3
    , const MeshDomain& domain
    )
    : m_ptree(NULL)
    , m_own_ptree()
    , m_pc3t3(pc3t3)
    , m_domain(domain)
    , m_params(domain)
  {
    init_aabb_tree();
  }

#else //POLYHEDRAL EXAMPLE
  Lipschitz_sizing(const MeshDomain& domain
    , const Tree* ptree
    , const Patches_ids_map& patches_ids_map
    , const CGAL::cpp11::array<double, 3>& vxyz
    , const CGAL::Bbox_3& bbox
    , const bool domain_is_a_box
    )
    : m_ptree(ptree)
    , m_own_ptree()
    , m_pc3t3(NULL)
    , m_domain(domain)
    , m_params(domain)
    , m_get_facet_patch_id()
    , patches_ids_map(patches_ids_map)
    , m_vxyz(vxyz)
    , m_bbox(bbox)
    , m_domain_is_a_box(domain_is_a_box)
  {
#  ifdef CGAL_MESH_3_LIPSCHITZ_SIZING_VERBOSE
    const std::vector<Surface_patch_index>&
      ids = m_domain.boundary_patches();

    std::cout << "Boundary patch ids : ";
    if (ids.empty()) std::cout << "empty";
    BOOST_FOREACH(Surface_patch_index spi, ids)
      std::cout << spi << " ";
    std::cout << std::endl;
#  endif // CGAL_MESH_3_LIPSCHITZ_SIZING_VERBOSE
    kd_tree();
  }
#endif // not CGAL_MESH_3_IMAGE_EXAMPLE


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
    }
#ifndef CGAL_MESH_3_IMAGE_EXAMPLE
    else if (dim == 1)
    {
      const typename MeshDomain::Curve_segment_index& curve_id =
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
    }
    else if (dim == 0)
    {
      const Patches_ids& ids =
        (m_domain.corners_incidences_map().find(p)->second);

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
    }
#endif

    CGAL_assertion(false);
    return 0.;
  }

private:
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

      FT k, layer_thickness, size_in_layer, size_max;
      if (is_in_domain(subdomains.first))
      {
        m_params.get_parameters(subdomains.first,
          k, layer_thickness, size_in_layer, size_max);
        size = (std::min)(size, size_in_layer);
      }
      if (is_in_domain(subdomains.second))
      {
        m_params.get_parameters(subdomains.second,
          k, layer_thickness, size_in_layer, size_max);
        size = (std::min)(size, size_in_layer);
      }
    }
    return size;
  }

  void init_aabb_tree()
  {
#ifdef CGAL_MESH_3_IMAGE_EXAMPLE
    static std::list<Triangle> triangles;
    for (typename C3T3::Facets_in_complex_iterator
      fit = m_pc3t3->facets_in_complex_begin();
      fit != m_pc3t3->facets_in_complex_end();
      ++fit)
    {
      if (!is_on_cube_boundary(*fit))
        triangles.push_back(m_pc3t3->triangulation().triangle(*fit));
    }

    m_own_ptree.reset(new Tree(triangles.begin(), triangles.end()));
    m_own_ptree->build();
    m_own_ptree->accelerate_distance_queries();
#endif
  }

#ifdef CGAL_MESH_3_IMAGE_EXAMPLE
  bool is_on_cube_boundary(const Facet& f) const
  {
    return is_on_cube_boundary(m_pc3t3->surface_patch_index(f));
  }
#endif //CGAL_MESH_3_IMAGE_EXAMPLE

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
    FT k, layer_thickness, size_in_layer, size_max;
    m_params.get_parameters(index,
      k, layer_thickness, size_in_layer, size_max);

#ifdef CGAL_MESH_3_IMAGE_EXAMPLE
    FT sqdist = (m_ptree->empty())
      ? m_own_ptree->squared_distance(p)
      : m_ptree->squared_distance(p);

#else //POLYHEDRAL EXAMPLE

    Point_3 closest = compute_closest_point(p);
    FT sqdist = CGAL::squared_distance(p, closest);

#endif

    if (sqdist < layer_thickness * layer_thickness)//inside the constant layer
      return size_in_layer;

    FT size = k * (CGAL::sqrt(sqdist) - layer_thickness) + size_in_layer;
    return (std::min)(size, size_max);
  }

  void kd_tree()
  {
#ifndef CGAL_MESH_3_IMAGE_EXAMPLE
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
#endif // not CGAL_MESH_3_IMAGE_EXAMPLE
  }

#ifndef CGAL_MESH_3_IMAGE_EXAMPLE
  Point_3 compute_closest_point(const Point_3& p) const
  {
    const std::vector<Surface_patch_index>& boundary_ids =
      m_domain.boundary_patches();

    CGAL_STATIC_THREAD_LOCAL_VARIABLE_4(AABB_filtered_traits,
                                        projection_traits,
                                        boundary_ids.begin(),
                                        boundary_ids.end(),
                                        m_ptree->traits(),
                                        m_get_facet_patch_id);
    Neighbor_search search(*m_kd_tree, p, 1);
    projection_traits.reset(search.begin()->first);

    m_ptree->traversal(p, projection_traits);
    CGAL_assertion(projection_traits.found());
    return projection_traits.closest_point();
  }
#endif // not CGAL_MESH_3_IMAGE_EXAMPLE

public:
  void add_parameters_for_subdomain(const Subdomain_index& id
                                  , const FT& k
                                  , const FT& layer_thickness
                                  , const FT& size_in_layer
                                  , const FT& size_max)
  {
    m_params.add_subdomain(id, k, layer_thickness, size_in_layer, size_max);
  }

  bool read_parameters(const char* filename)
  {
    std::ifstream infile(filename);
    if (!infile)
    {
      std::cerr << "ERROR : " << filename << " is not a valid file name" << std::endl;
      return false;
    }
    else
      std::cout << "Reading Lipschitz parameters...";

    std::string line;
    while (std::getline(infile, line))
    {
      if (line.empty() || line.at(0) == '#')
        continue; //skip comments

      std::istringstream iss(line);
      Subdomain_index id;
      FT size_in_layer, size_max, layer_thickness, k;
      if (!(iss >> id >> size_in_layer >> size_max >> layer_thickness >> k))
      {
        std::cerr << "ERROR while reading line : [" << line << "]" << std::endl;
        return false;
      }

      add_parameters_for_subdomain(id, k, layer_thickness, size_in_layer, size_max);
    }
    std::cout << "\b\b\b done (" << m_params.size() << " subdomains)." << std::endl;
    return !m_params.empty();
  }

};

}//namespace CGAL

#endif // _LIPSCHITZ_SIZING_
