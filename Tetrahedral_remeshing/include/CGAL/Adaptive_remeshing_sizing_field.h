// Copyright (c) 2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description : Defines a sizing field adapted to a triangulation
//******************************************************************************

#ifndef CGAL_TETRAHEDRAL_REMESHING_ADAPTIVE_SIZING_FIELD_H
#define CGAL_TETRAHEDRAL_REMESHING_ADAPTIVE_SIZING_FIELD_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Tetrahedral_remeshing_sizing_field.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/property_map.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/property_maps.h>

#include <vector>
#include <array>

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
#include <CGAL/property_map.h>
#include <CGAL/IO/write_ply_points.h>
#include <boost/property_map/property_map.hpp>
#endif

namespace CGAL
{

/**
 * @ingroup PkgTetrahedralRemeshingSizing
 *
 * @class Adaptive_remeshing_sizing_field
 * @tparam Tr underlying triangulation type
 *
 * A sizing field for tetrahedral remeshing,
 * that keeps the same mesh density throughout the remeshing process.
 *
 * \cgalModels{RemeshingSizingField_3}
 */
template <typename Tr>
class Adaptive_remeshing_sizing_field
  : public Tetrahedral_remeshing_sizing_field<typename Tr::Geom_traits>
{
  // Types
public:
  typedef typename Tr::Geom_traits              GT;
  typedef typename GT::FT                       FT;
  typedef typename Tr::Geom_traits::Point_3     Point_3; //Bare_point
  typedef typename Tr::Vertex::Index            Index;

private:
  typedef typename Tr::Point                    Tr_point;
  typedef typename Tr::Vertex_handle            Vertex_handle;
  typedef typename Tr::Cell_handle              Cell_handle;

  struct Point_with_info
  {
    Point_3 p;
    FT size;
    int dimension;
  };

private:
  struct Point_property_map
  {
    using Self = Point_property_map;
    using value_type = Point_3;
    using reference = value_type; //TODO : why can't that be value_type& ?
    using key_type = Point_with_info;
    using category = boost::readable_property_map_tag;

    const value_type operator[](const key_type& pwi) const { return pwi.p; }
    friend const value_type get(const Self&, const key_type& pwi) { return pwi.p; }
  };

private:
  using Kd_traits = CGAL::Search_traits_adapter<Point_with_info,
                                                Point_property_map,
                                                CGAL::Search_traits_3<GT> >;
  using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Kd_traits>;
  using Kd_tree = typename Neighbor_search::Tree;
  using Distance = typename Neighbor_search::Distance;
  using Splitter = typename Neighbor_search::Splitter;

#ifndef DOXYGEN_RUNNING
public:
  template<typename ECMap, typename FCMap, typename CellSelector>
  Adaptive_remeshing_sizing_field(const Tr& tr,
                                  const ECMap& ecmap,
                                  const FCMap& fcmap,
                                  const CellSelector& cell_selector)
    : m_kd_tree_3(points_with_info(tr, 3, ecmap, fcmap, cell_selector),
                  Splitter(),
                  Kd_traits(Point_property_map()))
    , m_kd_tree_2(points_with_info(tr, 2, ecmap, fcmap, cell_selector),
                Splitter(),
                Kd_traits(Point_property_map()))
  {
    m_kd_tree_3.build();
    m_kd_tree_2.build();

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    dump_kd_trees();
#endif
  }
#endif

private:

  template<typename ECMap, typename FCMap, typename CellSelector>
  std::vector<Point_with_info> points_with_info(const Tr& tr,
                                                const int dim,
                                                const ECMap& ecmap,
                                                const FCMap& fcmap,
                                                const CellSelector& cell_selector) const
  {
    namespace Tet_remeshing = CGAL::Tetrahedral_remeshing;
    using Facet = typename Tr::Facet;
    using Edge = typename Tr::Edge;

    auto cp = tr.geom_traits().construct_point_3_object();
    auto centroid = tr.geom_traits().construct_centroid_3_object();
    auto midpt = tr.geom_traits().construct_midpoint_3_object();

    std::vector<Point_with_info> points;
    for (const Vertex_handle v : tr.finite_vertex_handles())
    {
      if(  (dim == 3 && dim == v->in_dimension())//inside volume
        || (dim < 3  && v->in_dimension() < 3))  //on surface
      {
        const FT size = average_edge_length_around(v, tr, ecmap, fcmap, cell_selector);
        if (CGAL::is_zero(size))
          continue;
        points.push_back(Point_with_info{ cp(tr.point(v)), size, v->in_dimension() });
      }
    }

    // add internal points
    for (const Cell_handle c : tr.finite_cell_handles())
    {
      if (!get(cell_selector, c))
        continue;

      // inside cells
      if(dim == 3)
      {
        const FT size = average_edge_length_3(c, tr);
        if (CGAL::is_zero(size))
          continue;
        points.push_back(Point_with_info{ centroid(tr.tetrahedron(c)), size, 3 });
      }
      else
      {
        // on surface facets
        for (int i = 0; i < 4; ++i)
        {
          const Cell_handle cn = c->neighbor(i);
          if(  get(cell_selector, c) != get(cell_selector, cn)
            || get(fcmap, Facet(c, i))
            || c->is_facet_on_surface(i) )
          {
            const FT size = average_edge_length_2(c, i, tr);
            if (CGAL::is_zero(size))
              continue;
            points.push_back(Point_with_info{ centroid(tr.triangle(c, i)), size, 2 });
          }
        }
        // on complex edges
        for (const Edge& e : Tet_remeshing::cell_edges(c, tr))
        {
          if(get(ecmap, Tet_remeshing::make_vertex_pair(e)))
          {
            const FT size = Tet_remeshing::approximate_edge_length(e, tr);
            if (CGAL::is_zero(size))
              continue;
            points.push_back(Point_with_info{midpt(tr.segment(e)), size, 1 });
          }
        }
      }
    }
    return points;
  }

  const Kd_tree& kd_tree(const int dim) const
  {
    if (dim == 3)
      return m_kd_tree_3;
    else
      return m_kd_tree_2;
  }

  int nb_neighbors(const int dim) const
  {
    if (dim == 3)      return 30;
    else               return 6;
  }

public:
  /**
  * Returns size at point `p`, assumed to be included in the input
  * subcomplex with dimension `dim` and index `index`.
  */
#ifdef DOXYGEN_RUNNING
  FT operator()(const Point_3& p, const int& dim, const Index& index) const
#else
  FT operator()(const Point_3& p, const int& dim, const Index& ) const
#endif
  {
    Point_property_map pp_map;
    Distance dist(pp_map);

    // Find nearest vertex and local size before remeshing
    Neighbor_search search(kd_tree(dim),
                           p, //query point
                           nb_neighbors(dim), //nb nearest neighbors
                           0, //epsilon
                           true, //search nearest
                           dist);

    FT max_size = 0.;
    std::vector<Point_with_info> neighbors;

    for (const auto& neighbor : search)
    {
      [[maybe_unused]] const auto& [pi, size, dimension] = neighbor.first;

      max_size = (std::max)(max_size, size);

      if(  (dim == 3 && dimension == dim) //volume point
        || (dim < 3  && dimension < 3) )  //surface point
      {
        neighbors.push_back(neighbor.first);
      }
    }

    if (neighbors.empty())
      return max_size;

    return interpolate_on_n_vertices(p, neighbors);
  }

private:
  /**
   * Returns size at point `p`, by interpolation among neighbors
   */
  FT interpolate_on_n_vertices(
    const Point_3& p,
    const std::vector<Point_with_info>& vertices) const;

  template<typename ECMap, typename FCMap, typename CellSelector>
  FT average_edge_length_around(const Vertex_handle v, const Tr& tr,
    const ECMap& ecmap, const FCMap& fcmap, const CellSelector& cell_selector) const;

  FT average_edge_length_3(const Cell_handle c, const Tr& tr) const;
  FT average_edge_length_2(const Cell_handle c, const int i, const Tr& tr) const;

private:
  Kd_tree m_kd_tree_3; //volumes
  Kd_tree m_kd_tree_2; //surfaces

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
private:
  void dump_kd_trees() const
  {
    std::ofstream ofs("kd_trees.ply");

    using Point_Size = std::pair<Point_3, float>;
    std::vector<Point_Size> points;

    for (auto it = m_kd_tree_3.begin(); it != m_kd_tree_3.end(); ++it)
      points.push_back(std::make_pair(it->p, it->size));

    for (auto it = m_kd_tree_2.begin(); it != m_kd_tree_2.end(); ++it)
      points.push_back(std::make_pair(it->p, it->size));

    using Point_map = CGAL::First_of_pair_property_map<Point_Size>;
    using Size_map = CGAL::Second_of_pair_property_map<Point_Size>;

    CGAL::IO::write_PLY_with_properties(ofs,
      points,
      CGAL::make_ply_point_writer(Point_map()),
      std::make_pair(Size_map(), CGAL::IO::PLY_property<double>("intensity")));
  }
#endif

};//end of class Adaptive_remeshing_sizing_field

/*!
 * @ingroup PkgTetrahedralRemeshingSizing
 *
 * @brief Create an adaptive sizing field for tetrahedral remeshing
 *
 * \relates Adaptive_remeshing_sizing_field
 *
 * This method is a <em>named constructor</em>. It constructs a
 * sizing field of type `Adaptive_remeshing_sizing_field<Tr>`,
 * designed to keep the density unchanged throughout the remeshing
 * process.
 *
 * @returns an `Adaptive_remeshing_sizing_field<Tr>`
 *
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 * \param tr the input triangulation
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters"
 *           among the ones listed below. All of them must be the same as the ones
 *           given to `CGAL::tetrahedral_isotropic_remeshing()`.
 *
 *\cgalNamedParamsBegin
 *  \cgalParamNBegin{ edge_is_constrained_map }
 *    \cgalParamDescription{a property map containing the constrained-or-not
 *        status of each edge of `tr`.}
 *    \cgalParamDefault{a default property map where no edge is constrained}
 *  \cgalParamNEnd
 *  \cgalParamNBegin{ facet_is_constrained_map }
 *    \cgalParamDescription{a property map containing the constrained-or-not
 *        status of each facet of `tr`.}
 *   \cgalParamDefault{a default property map where no facet is constrained}
 *  \cgalParamNEnd
 *  \cgalParamNBegin{ cell_is_selected_map }
 *    \cgalParamDescription{a property map containing the selected
 *         - or - not status for each cell of `tr` for remeshing.}
 *    \cgalParamDefault{a default property map where all cells of the domain are selected}
 *  \cgalParamNEnd
 *\cgalNamedParamsEnd
 */
template<typename Tr, typename CGAL_NP_TEMPLATE_PARAMETERS>
Adaptive_remeshing_sizing_field<Tr>
create_adaptive_remeshing_sizing_field(const Tr& tr,
  const CGAL_NP_CLASS& np = parameters::default_values())
{
  CGAL_CHECK_AUTHORIZED_NAMED_PARAMETERS(np, edge_is_constrained_t, facet_is_constrained_t, cell_selector_t);

  using parameters::choose_parameter;
  using parameters::get_parameter;

  using V = typename Tr::Vertex_handle;
  using Edge_vv = std::pair<V, V>;
  using Facet = typename Tr::Facet;

  auto ecmap = choose_parameter(get_parameter(np, internal_np::edge_is_constrained),
    CGAL::Constant_property_map<Edge_vv, bool>(false));
  auto fcmap = choose_parameter(get_parameter(np, internal_np::facet_is_constrained),
    CGAL::Constant_property_map<Facet, bool>(false));
  auto cell_selector = choose_parameter(get_parameter(np, internal_np::cell_selector),
    CGAL::Tetrahedral_remeshing::internal::All_cells_selected<Tr>());

  return Adaptive_remeshing_sizing_field<Tr>(tr, ecmap, fcmap, cell_selector);
}


template <typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
interpolate_on_n_vertices(
  const Point_3& p,
  const std::vector<Point_with_info>& points_with_info) const
{
  // Interpolate value using values at vertices
  const auto sqd = GT().compute_squared_distance_3_object();

  FT sum_weights = 0.;
  FT sum_sizes = 0.;

  for (const auto pwi : points_with_info)
  {
    const FT size = pwi.size;
    const FT sqdist = sqd(pwi.p, p);

    if (is_zero(sqdist))
      return size;

    const FT weight = 1. / CGAL::approximate_sqrt(sqdist);
    sum_weights += weight;
    sum_sizes += weight * size;
  }

  CGAL_assertion(sum_weights > 0);
  CGAL_assertion(sum_sizes > 0);
  return sum_sizes / sum_weights;
}


template<typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
average_edge_length_3(const typename Tr::Cell_handle c, const Tr& tr) const
{
  using namespace CGAL::Tetrahedral_remeshing;

  FT sum = 0.;
  for (const typename Tr::Edge& e : cell_edges(c, tr))
  {
    sum += approximate_edge_length(e, tr);
  }
  return sum / FT(6);
}

template<typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
average_edge_length_2(const typename Tr::Cell_handle c,
                      const int i,
                      const Tr& tr) const
{
  // facet(c,i)
  namespace Tet_Remeshing = CGAL::Tetrahedral_remeshing;

  FT sum = 0.;
  short nb_surface_edges = 0;
  for (const typename Tr::Edge& e : Tet_Remeshing::facet_edges(c, i, tr))
  {
    CGAL_assertion_code(const auto vs = tr.vertices(e));
    CGAL_assertion(vs[0]->in_dimension() < 3);
    CGAL_assertion(vs[1]->in_dimension() < 3);

    sum += Tet_Remeshing::approximate_edge_length(e, tr);
    ++nb_surface_edges;
  }

  return sum / FT(nb_surface_edges);
}

template <typename Tr>
template <typename ECMap, typename FCMap, typename CellSelector>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
average_edge_length_around(const Vertex_handle v, const Tr& tr,
                           const ECMap& ecmap, const FCMap& fcmap,
                           const CellSelector& cell_selector) const
{
  using Edge = typename Tr::Edge;
  using Facet = typename Tr::Facet;
  using namespace CGAL::Tetrahedral_remeshing;

  std::vector<Edge> tmp_edges;
  tr.incident_edges(v, std::back_inserter(tmp_edges));

  std::vector<Edge> edges;
  switch (v->in_dimension())
  {
  case 3:
    std::copy_if(tmp_edges.begin(), tmp_edges.end(),
      std::back_inserter(edges),
      [&tr, &cell_selector](const Edge& e)
      {
        return is_selected(e, tr, cell_selector);
      });
    break;

  case 2:
    std::copy_if(tmp_edges.begin(), tmp_edges.end(),
      std::back_inserter(edges),
      [&fcmap, &cell_selector, &tr](const Edge& e)
      {
        auto fcirc = tr.incident_facets(e);
        const auto fend = fcirc;
        do
        {
          const Facet& f = *fcirc;
          if ( get(fcmap, f)
            || get(cell_selector, f.first) != get(cell_selector, f.first->neighbor(f.second)))
            return true;
        } while (++fcirc != fend);

        return false;
      });
    break;

  case 1:
  case 0:
    std::copy_if(tmp_edges.begin(), tmp_edges.end(),
      std::back_inserter(edges),
      [&ecmap](const Edge& e)
      {
        const auto evv = make_vertex_pair(e);
        return get(ecmap, evv);
      });
    break;

  default: // could be -1 because of far points
    break;
  }

  if (edges.empty())
    return 0;

  FT sum = 0.;
  for (const Edge& e : edges)
  {
    sum += approximate_edge_length(e, tr);
  }

  return sum / static_cast<FT>(edges.size());
}

} //namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ADAPTIVE_SIZING_FIELD_H
