// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// Copyright (c) 2010-2013 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Mael Rouxel-Labbé, Maxime Gimeno, Jane Tournois
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_SMDS_3_TETRAHEDRON_SOUP_TO_C3T3_H
#define CGAL_SMDS_3_TETRAHEDRON_SOUP_TO_C3T3_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>
#include <CGAL/SMDS_3/internal/SMDS_3_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/unordered_map.hpp>
#include <boost/container/flat_set.hpp>

#include <array>
#include <queue>
#include <set>
#include <vector>

namespace CGAL {

  /**
  * \ingroup PkgSMDS3Functions
  *
  * \brief returns `true` if the soup of tetrahedra defines a valid triangulation
  * that can be handled by `tetrahedron_soup_to_triangulation_3()`.
  *
  * It checks that each facet has at most two incident cells,
  * no tetrahedron has twice the same vertex,
  * and the tetrahedron soup describes a manifold volume.
  * This function does not require a range of points as an argument
  * since the check is purely topological.
  *
  * @tparam TetrahedronRange a model of the concept `RandomAccessContainer`
  * whose `value_type` is a model of the concept `RandomAccessContainer`
  * whose `value_type` is `std::size_t`.
  *
  * @param tets each element in the range describes a tetrahedron
  * using the indices of the vertices.
  */
  template <typename TetrahedronRange,
            typename NamedParameters = parameters::Default_named_parameters>
  bool is_tetrahedron_soup_a_triangulation(const TetrahedronRange& tets,
                                           const NamedParameters& np = parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::is_default_parameter;
    using parameters::get_parameter;

    typedef typename boost::range_value<TetrahedronRange>::type           Tetrahedron;
    typedef typename boost::range_value<Tetrahedron>::type                V_ID;

    const bool verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);

    if(std::begin(tets) == std::end(tets))
      return true;

    // check that no tetrahedron has twice the same vertex
    // and build mapping of facets to incident cells
    V_ID max_id = 0;

    typedef std::array<V_ID, 3> Facet_vvv;
    typedef std::vector<std::size_t> Incident_tet_indices;
    typedef boost::unordered_map<Facet_vvv, Incident_tet_indices> Incident_cells_map;
    Incident_cells_map incident_cells_map;

    for(std::size_t tet_idx = 0; tet_idx < tets.size(); ++tet_idx)
    {
      const Tetrahedron& tet = tets[tet_idx];

      if(tet.size() != 4) {
        if (verbose) {
          std::cerr << "Tetrahedron #" << tet_idx << " has " << tet.size() << " vertices" << std::endl;
        }
        return false;
      }

      boost::container::flat_set<V_ID> tet_vertices;
      std::vector<V_ID> tet_vertices_vec;

      for(V_ID id : tet)
      {
        if(max_id < id)
          max_id = id;

        if(!tet_vertices.insert(id).second) {
          if (verbose) {
            std::cerr << "Tetrahedron #" << tet_idx << " has duplicate vertex " << id << std::endl;
          }
          return false; // vertex met twice in the same tetrahedron
        }

        tet_vertices_vec.push_back(id);
      }

      for(int i=0; i<4; ++i)
      {
        Facet_vvv facet;
        int idx = 0;
        for(int j = 0; j < 4; ++j)
        {
          if(i != j)
            facet[idx++] = tet_vertices_vec[j];
        }

        // Sort the facet vertices to create a canonical ordering
        if(facet[1] < facet[0]) std::swap(facet[0], facet[1]);
        if(facet[2] < facet[1]) std::swap(facet[1], facet[2]);
        if(facet[1] < facet[0]) std::swap(facet[0], facet[1]);

        incident_cells_map[facet].push_back(tet_idx);

        // check that each facet has at most two incident cells
        if(incident_cells_map[facet].size() > 2) {
          if (verbose) {
            std::cerr << "Facet with > 2 incident cells" << std::endl;
          }
          return false;
        }
      }
    }

    // check that there is only a single connected component of tetrahedra
    // Build adjacency graph: two tetrahedra are adjacent if they share a facet
    std::vector<std::vector<std::size_t> > tet_adjacency(tets.size());

    for(const auto& facet_and_incident : incident_cells_map)
    {
      const Incident_tet_indices& incident = facet_and_incident.second;
      if(incident.size() == 2)
      {
        std::size_t tet1 = incident[0];
        std::size_t tet2 = incident[1];
        tet_adjacency[tet1].push_back(tet2);
        tet_adjacency[tet2].push_back(tet1);
      }
    }

    // Check connectivity via BFS
    if(tets.size() > 0)
    {
      std::vector<bool> visited(tets.size(), false);
      std::queue<std::size_t> q;
      q.push(0);
      visited[0] = true;
      std::size_t component_size = 1;

      while(!q.empty())
      {
        std::size_t current = q.front();
        q.pop();

        for(std::size_t neighbor : tet_adjacency[current])
        {
          if(!visited[neighbor])
          {
            visited[neighbor] = true;
            q.push(neighbor);
            component_size++;
          }
        }
      }

      if(component_size != tets.size())
      {
        if (verbose) {
          std::cerr << component_size << " volume connected components" << std::endl;
        }
        return false;
      }
    }

    // check that there is only a single volume connected component of tetrahedra around each vertex
    // and around each edge
    typedef std::pair<V_ID, V_ID> E_ID;
    typedef std::map<V_ID, std::vector<std::size_t> > Vertex_to_tets_map;
    typedef std::map<E_ID, std::vector<std::size_t> > Edge_to_tets_map;

    Vertex_to_tets_map vertex_to_tets;
    Edge_to_tets_map edge_to_tets;

    // Build vertex-to-tetrahedra and edge-to-tetrahedra maps
    for(std::size_t tet_idx=0; tet_idx<tets.size(); ++tet_idx)
    {
      const Tetrahedron& tet = tets[tet_idx];

      for(V_ID v : tet)
        vertex_to_tets[v].push_back(tet_idx);

      for(int i=0; i<4; ++i)
      {
        for(int j=i+1; j<4; ++j)
        {
          V_ID v1 = tet[i];
          V_ID v2 = tet[j];
          if(v1 > v2)
            std::swap(v1, v2);
          edge_to_tets[std::make_pair(v1, v2)].push_back(tet_idx);
        }
      }
    }

    // Check that around each vertex, tetrahedra form a single connected component
    for(const auto& v_and_tets : vertex_to_tets)
    {
      const std::vector<std::size_t>& inc_tets = v_and_tets.second;
      CGAL_assertion(inc_tets.size() != 0);

      std::set<std::size_t> remaining(inc_tets.begin(), inc_tets.end());
      std::queue<std::size_t> local_q;

      std::size_t seed = *remaining.begin();
      local_q.push(seed);
      remaining.erase(seed);

      while(!local_q.empty())
      {
        std::size_t current = local_q.front();
        local_q.pop();

        for(std::size_t neighbor : tet_adjacency[current])
        {
          auto it = remaining.find(neighbor);
          if(it != remaining.end())
          {
            local_q.push(neighbor);
            remaining.erase(it);
          }
        }
      }

      if(!remaining.empty())
      {
        if (verbose) {
          std::cerr << "more than one CC around vertex" << std::endl;
        }
        return false;
      }
    }

    // Check that around each edge, tetrahedra form a single connected component
    for(const auto& e_and_tets : edge_to_tets)
    {
      const std::vector<std::size_t>& inc_tets = e_and_tets.second;
      if(inc_tets.size() == 0)
        continue;

      std::set<std::size_t> remaining(inc_tets.begin(), inc_tets.end());
      std::queue<std::size_t> local_q;

      std::size_t seed = *remaining.begin();
      local_q.push(seed);
      remaining.erase(seed);

      while(!local_q.empty())
      {
        std::size_t current = local_q.front();
        local_q.pop();

        for(std::size_t neighbor : tet_adjacency[current])
        {
          auto it = remaining.find(neighbor);
          if(it != remaining.end())
          {
            local_q.push(neighbor);
            remaining.erase(it);
          }
        }
      }

      if(!remaining.empty())
      {
        if (verbose) {
          std::cerr << "more than one CC around edge" << std::endl;
        }
        return false;
      }
    }

    if (verbose) {
      std::cout << "Tetrahedron soup is a valid triangulation" << std::endl;
    }

    return true;
  }

  /** \ingroup PkgSMDS3Functions
  * builds a 3D triangulation from a soup of tetrahedra.
  *
  * @tparam TetrahedronRange a model of `Range` whose value type is
  * a `Tetrahedron_3`. The point type of the tetrahedra must be convertible
  * to `Triangulation::Point`
  * @tparam Triangulation a 3D triangulation class that has
  * a vertex base model of `SimplicialMeshVertexBase_3`
  * and a cell base model of `SimplicialMeshCellBase_3`
  *
  * @param tets the set of finite tetrahedra of a valid \cgal triangulation.
  * Each element in the range is the geometric description of the
  * corresponding cell in the triangulation.
  *
  * @returns the 3D triangulation built from \p tets
  *
  * @pre the tetrahedron soup must form a partition of the convex hull of `tets`
  *
  * @sa @ref SMDS_3/tetrahedron_soup_to_c3t3_example.cpp
  *
  */
  template<typename Triangulation, typename TetrahedronRange>
  Triangulation tetrahedron_soup_to_triangulation_3(const TetrahedronRange& tets)
  {
    using Tr    = Triangulation;
    using Point = typename Tr::Point;

    Triangulation tr;
    std::vector<Point> points;
    std::vector<std::array<int, 4> > finite_cells;
    boost::unordered_map<std::array<int, 3>, typename Tr::Cell::Surface_patch_index> border_facets;
    boost::unordered_map<Point, int> p2i;

    CGAL_assertion_code(
      typename Triangulation::Geom_traits::Orientation_3 orientation =
      tr.geom_traits().orientation_3_object();
    );

    for (typename TetrahedronRange::value_type tet : tets)
    {
      CGAL_assertion(tet.orientation() == CGAL::POSITIVE);
      std::array<int, 4> cell;

      for (int i = 0; i < 4; ++i)
      {
        const Point& pi = tet[i];
        if (p2i.find(pi) == p2i.end())
        {
          points.push_back(pi);
          int index = static_cast<int>(points.size() - 1);
          p2i.insert(std::make_pair(pi, index));
          cell[i] = index;
        }
        else
          cell[i] = p2i.at(pi);
      }

      CGAL_assertion(orientation(points[cell[0]],
        points[cell[1]], points[cell[2]], points[cell[3]]) == CGAL::POSITIVE);

      finite_cells.push_back(cell);
    }

    typename Tr::Cell::Subdomain_index default_si(1);
    CGAL::SMDS_3::build_triangulation_one_subdomain(tr, points, finite_cells, default_si, border_facets,
      /*verbose = */false, /*replace_domain_0 = */false, /*allow_non_manifold =*/false);

    CGAL_assertion(CGAL::SMDS_3::internal::is_convex(tr));

    return tr;
  }

  /** \ingroup PkgSMDS3Functions
  * builds a 3D triangulation from a soup of tetrahedra.
  *
  * @tparam PointRange a model of the concept `RandomAccessContainer`
  * whose value type is the point type.
  * The point type must be convertible to `Triangulation::Point`.
  * @tparam TetrahedronRange a model of the concept `RandomAccessContainer` whose
  * value type is a model of the concept `RandomAccessContainer` whose value type is `std::size_t`
  * @tparam Triangulation a 3D triangulation class that has
  * a vertex base model of `MeshVertexBase_3`
  * and a cell base model of `MeshCellBase_3`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param points points of the soup of tetrahedra
  * @param tets the set of finite tetrahedra of a valid \cgal triangulation.
  * Each element in the range describes a tetrahedron using the indices of the points
  * in `points`.
  * It must
  * describe a non self-intersecting set of tetrahedra, that cover the convex hull of the
  * corresponding point set. The tetrahedra must form a valid triangulation with each
  * pair of neighboring cells sharing exactly one triangle. Combinatorial validity and
  * validity of the geometric embedding are required.
  *
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{surface_facets}
  *     \cgalParamDescription{each element in the range describes a surface facet using the indices of points
  *       in `points` (indices 0 to 2), and the associated `Surface_patch_index` (index 3)}
  *     \cgalParamType{a class model of `AssociativeContainer`
  *                    whose key type is model of `RandomAccessContainer` containing `int`
  *                    and mapped type is `Tr::Cell::Surface_patch_index`}
  *     \cgalParamDefault{an empty `boost::unordered_map<std::array<int, 3>, typename Tr::Cell::Surface_patch_index>`}
  *     \cgalParamExtra{to avoid copies of large data sets, this parameter can be passed using `std::cref`}
  *   \cgalParamNEnd
  *   \cgalParamNBegin{subdomain_indices}
  *     \cgalParamDescription{each element in the range gives the
  *       `Triangulation::Cell::Subdomain_index` corresponding to the tetrahedron (cell)
  *        of same index in `tets`}
  *     \cgalParamType{a class model of `RandomAccessContainer` whose value type
  *                    is `Triangulation::Cell::Subdomain_index`}
  *     \cgalParamDefault{each finite cell of the output triangulation is
  *                    set to have `1` as `Subdomain_index`}
  *     \cgalParamExtra{to avoid copies of large data sets, this parameter can be passed using `std::cref`}
  *   \cgalParamNEnd
  *\cond SKIP_IN_MANUAL
  *   \cgalParamNBegin{allow_non_manifold}
  *     \cgalParamDescription{allows the construction of a triangulation with non-manifold edges
  *       and non manifold vertices. The triangulation is invalid if this situation is met,
  *       so it should be used only in advanced cases, and the triangulation will be hardly usable.}
  *     \cgalParamType{bool}
  *     \cgalParamDefault{false}
  *   \cgalParamNEnd
  *\endcond
  * \cgalNamedParamsEnd
  *
  * @returns the 3D triangulation built from parameters
  *
  * @pre `points` contains each point only once
  * @pre the tetrahedron soup must form a partition of the convex hull of `tets`
  * @post `is_valid()` returns `true` for the returned triangulation
  *
  * @sa \link Polygon_mesh_processing::polygon_soup_to_polygon_mesh() `CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh()` \endlink
  * @sa @ref SMDS_3/tetrahedron_soup_to_c3t3_example.cpp
  *
  */
  template<typename Triangulation,
           typename PointRange,
           typename TetrahedronRange,
           typename NamedParameters = parameters::Default_named_parameters>
  Triangulation tetrahedron_soup_to_triangulation_3(const PointRange& points,
                                                    const TetrahedronRange& tets,
                                                    const NamedParameters& np = parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::get_parameter_reference;

    using Default_facet_map
      = boost::unordered_map<std::array<int, 3>,
                             typename Triangulation::Cell::Surface_patch_index>;
    using Facet_map_ref_type
      = typename internal_np::Lookup_named_param_def<
                             internal_np::surface_facets_t,
                             NamedParameters,
                             Default_facet_map>::reference;
    using Subdomain_index = typename Triangulation::Cell::Subdomain_index;
    using Default_subdomains = std::vector<Subdomain_index>;
    using Subdomains_ref_type
      = typename internal_np::Lookup_named_param_def<
                              internal_np::subdomain_indices_t,
                              NamedParameters,
                              Default_subdomains>::reference;

    Triangulation tr;
    Default_facet_map empty_map;
    Default_subdomains subdomain_indices(tets.size(), Subdomain_index(1));

    Facet_map_ref_type facets = choose_parameter(
          get_parameter_reference(np, internal_np::surface_facets),
          empty_map);
    Subdomains_ref_type subdomains = choose_parameter(
          get_parameter_reference(np, internal_np::subdomain_indices),
          subdomain_indices);
    const bool non_manifold = choose_parameter(
          get_parameter(np, internal_np::allow_non_manifold),
          false);

    CGAL::SMDS_3::build_triangulation_with_subdomains_range(tr, points, tets, subdomains, facets,
      /*verbose = */false, /*replace_domain_0 = */false, non_manifold);

    CGAL_assertion(CGAL::SMDS_3::internal::is_convex(tr));

    return tr;
  }

} //namespace CGAL


#endif // CGAL_SMDS_3_TETRAHEDRON_SOUP_TO_C3T3_H
