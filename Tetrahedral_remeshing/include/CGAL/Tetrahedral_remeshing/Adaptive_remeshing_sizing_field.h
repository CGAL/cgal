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

#include <CGAL/Tetrahedral_remeshing/Sizing_field.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/property_map.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/property_maps.h>

#include <vector>
#include <array>


namespace CGAL
{
namespace Tetrahedral_remeshing
{

/**
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
  : public Sizing_field<typename Tr::Geom_traits>
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
    : m_kd_tree(points_with_info(tr, ecmap, fcmap, cell_selector),
                Splitter(),
                Kd_traits(Point_property_map()))
  {
    m_kd_tree.build();
  }
#endif

private:

  template<typename ECMap, typename FCMap, typename CellSelector>
  std::vector<Point_with_info> points_with_info(const Tr& tr,
                                                const ECMap& ecmap,
                                                const FCMap& fcmap,
                                                const CellSelector& cell_selector) const
  {
    auto cp = tr.geom_traits().construct_point_3_object();
    std::vector<Point_with_info> points;
    for (const Vertex_handle v : tr.finite_vertex_handles())
    {
      points.push_back(
        Point_with_info{ cp(tr.point(v)),
                         average_edge_length_around(v, tr, ecmap, fcmap, cell_selector),
                         v->in_dimension() });
    }
    return points;
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
    const int nb_neighbors = (dim == 3) ? 20 : 6;

    // Find nearest vertex and local size before remeshing
    Point_property_map pp_map;
    Distance dist(pp_map);
    Neighbor_search search(m_kd_tree,
                           p, //query point
                           nb_neighbors, //nb nearest neighbors
                           0, //epsilon
                           true, //search nearest
                           dist);

    FT sum = 0;
    for (const auto& neighbor : search)
    {
      [[maybe_unused]] const auto& [pi, size, dimension] = neighbor.first;
      // todo : should we discard points when dim != dimension?

      sum += size;
    }

    CGAL_assertion(sum > 0);
    return sum / static_cast<FT>(nb_neighbors);
  }

private:
  /**
   * Returns size at point `p`, by interpolation into tetrahedron.
   * TODO : interpolate
   */
  FT interpolate_on_four_vertices(
    const Point_3& p,
    const std::array<Point_with_info, 4>& vertices) const;

  template<typename ECMap, typename FCMap, typename CellSelector>
  FT average_edge_length_around(const Vertex_handle v, const Tr& tr,
    const ECMap& ecmap, const FCMap& fcmap, const CellSelector& cell_selector) const;

private:
  Kd_tree m_kd_tree;

};//end of class Adaptive_remeshing_sizing_field

/*!
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
interpolate_on_four_vertices(
  const Point_3& p,
  const std::array<Point_with_info, 4>& vertices) const
{
  // Interpolate value using values at vertices
  const FT& va = vertices[0].size;
  const FT& vb = vertices[1].size;
  const FT& vc = vertices[2].size;
  const FT& vd = vertices[3].size;

  const Point_3& a = vertices[0].p;
  const Point_3& b = vertices[1].p;
  const Point_3& c = vertices[2].p;
  const Point_3& d = vertices[3].p;

  const auto sqd = FT().compute_squared_distance_3_object();

  const FT wa = 1. / sqd(a, p);
  const FT wb = 1. / sqd(b, p);
  const FT wc = 1. / sqd(c, p);
  const FT wd = 1. / sqd(d, p);

  // If den is 0, then compute the average value
  if (is_zero(wa + wb + wc + wd))
    return (va + vb + vc + vd) / 4.;
  else
    return (wa * va + wb * vb + wc * vc + wd * vd) / (wa + wb + wc + wd);
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
        const auto evv = CGAL::Tetrahedral_remeshing::make_vertex_pair(e);
        return get(ecmap, evv);
      });
    break;

  default: // could be -1 because of far points
    break;
  }

  CGAL_assertion(!edges.empty());
  if (edges.empty())
    return 0;

  auto cp = tr.geom_traits().construct_point_3_object();
  FT sum = 0.;
  for (const Edge& e : edges)
  {
    sum += CGAL::approximate_sqrt(
      CGAL::squared_distance(cp(e.first->vertex(e.second)->point()),
                             cp(e.first->vertex(e.third)->point())));;
  }

  return sum / static_cast<FT>(edges.size());
}

} // end namespace Tetrahedral_remeshing

} //namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ADAPTIVE_SIZING_FIELD_H