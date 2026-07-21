// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_INTERNAL_SPLIT_LONG_EDGES_H
#define CGAL_INTERNAL_SPLIT_LONG_EDGES_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <unordered_map>
#include <functional>
#include <utility>
#include <optional>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{

template<typename C3t3>
bool positive_orientation_after_edge_split(const typename C3t3::Edge& e,
                                           const typename C3t3::Cell_handle circ,
                                           const typename C3t3::Triangulation::Geom_traits::Point_3& steiner,
                                           const C3t3&)
{
  using Point = typename C3t3::Triangulation::Geom_traits::Point_3;

  const auto v1 = e.first->vertex(e.second);
  const auto v2 = e.first->vertex(e.third);

  std::array<Point, 4> pts = {point(circ->vertex(0)->point()),
                              point(circ->vertex(1)->point()),
                              point(circ->vertex(2)->point()),
                              point(circ->vertex(3)->point())};
  // 1st half-cell
  const int i1 = circ->index(v1);
  const Point p1 = pts[i1];
  pts[i1] = steiner;
  if(CGAL::orientation(pts[0], pts[1], pts[2], pts[3]) != CGAL::POSITIVE)
    return false;

  // 2nd half-cell
  pts[i1] = p1;
  pts[circ->index(v2)] = steiner;
  if(CGAL::orientation(pts[0], pts[1], pts[2], pts[3]) != CGAL::POSITIVE)
    return false;

  return true;
}

template <typename C3t3>
std::optional<typename C3t3::Triangulation::Geom_traits::Point_3>
construct_steiner_point(const typename C3t3::Edge& e,
                        const C3t3& c3t3)
{
  using Cell_circulator = typename C3t3::Triangulation::Cell_circulator;
  using Cell_handle = typename C3t3::Triangulation::Cell_handle;
  using Point = typename C3t3::Triangulation::Geom_traits::Point_3;
  using FT = typename C3t3::Triangulation::Geom_traits::FT;

  const auto& gt = c3t3.triangulation().geom_traits();
  const auto& tr = c3t3.triangulation();
  const auto& p1 = point(e.first->vertex(e.second)->point());
  const auto& p2 = point(e.first->vertex(e.third)->point());
  const auto vec = gt.construct_vector_3_object()(p1, p2);

  const std::array<FT, 6> coeff = {0.33, 0.66,    //1/3 and 2/3
                                   0.3, 0.7,      // 0.5 +/- 0.2
                                   0.25, 0.75};   // 0.5 +/- 0.25

  std::size_t attempt_id = 0;
  while(attempt_id < coeff.size())
  {
    Point steiner = gt.construct_translated_point_3_object()(
        p1, gt.construct_scaled_vector_3_object()(vec, coeff[attempt_id]));
    ++attempt_id;

    bool steiner_successful = true;
    Cell_circulator circ = tr.incident_cells(e);
    Cell_circulator end = circ;
    do
    {
      Cell_handle c = circ;
      if(!positive_orientation_after_edge_split(e, c, steiner, c3t3))
      {
        steiner_successful = false;
        break;
      }
    } while(++circ != end);

    if(steiner_successful)
      return steiner;
  }

  return std::nullopt;
}

template<typename C3t3, typename CellSelector>
typename C3t3::Vertex_handle split_edge(const typename C3t3::Edge& e,
                                        CellSelector cell_selector,
                                        C3t3& c3t3)
{
  typedef typename C3t3::Triangulation       Tr;
  typedef typename C3t3::Subdomain_index     Subdomain_index;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename C3t3::Curve_index         Curve_index;
  typedef typename Tr::Geom_traits::Point_3 Point;
  typedef typename Tr::Facet                Facet;
  typedef typename Tr::Vertex_handle        Vertex_handle;
  typedef typename Tr::Cell_handle          Cell_handle;
  typedef typename Tr::Cell_circulator      Cell_circulator;

  Tr& tr = c3t3.triangulation();
  const Vertex_handle v1 = e.first->vertex(e.second);
  const Vertex_handle v2 = e.first->vertex(e.third);

  Point m = tr.geom_traits().construct_midpoint_3_object()
    (point(v1->point()), point(v2->point()));

  //backup subdomain info of incident cells before making changes
  short dimension = 0;
  if(c3t3.is_in_complex(e))
    dimension = 1;
  else
  {
    const std::size_t nb_patches = nb_incident_surface_patches(e, c3t3);
    if(nb_patches == 1)
      dimension = 2;
    else if(nb_patches == 0)
      dimension = 3;
    else
      CGAL_assertion(false);//e should be in complex
  }
  CGAL_assertion(dimension > 0);

  // remove complex edge before splitting
  const Curve_index curve_index = (dimension == 1) ? c3t3.curve_index(e) : Curve_index();

  struct Cell_info {
    Subdomain_index subdomain_index_;
    bool selected_;
  };
  struct Facet_info {
    Vertex_handle opp_vertex_;
    Surface_patch_index patch_index_;
  };
  boost::unordered_map<Facet, Cell_info, boost::hash<Facet>> cells_info;
  boost::unordered_map<Facet, Facet_info, boost::hash<Facet>> facets_info;

  // check orientation and collect incident cells to avoid circulating twice
  bool steiner_point_found = false;
  boost::container::small_vector<Cell_handle, 30> inc_cells;
  Cell_circulator circ = tr.incident_cells(e);
  Cell_circulator end = circ;
  do
  {
    inc_cells.push_back(circ);
    if (tr.is_infinite(circ) || steiner_point_found)
    {
      ++circ;
      continue;
    }

    const Cell_handle c = circ;
    if(!positive_orientation_after_edge_split(e, c, m, c3t3))
    {
      const std::optional<Point> steiner = construct_steiner_point(e, c3t3);
      if (steiner != std::nullopt)
      {
        m = *steiner;
        steiner_point_found = true;
      }
      else
        return Vertex_handle();
    }
    ++circ;
  }
  while (circ != end);

  if (dimension == 1)
    c3t3.remove_from_complex(e);

  for(Cell_handle c : inc_cells)
  {
    const int index_v1 = c->index(v1);
    const int index_v2 = c->index(v2);

    //keys are the opposite facets to the ones not containing e,
    //because they will not be modified
    const Subdomain_index subdomain = c3t3.subdomain_index(c);
    const bool selected = get(cell_selector, c);
    const Facet opp_facet1 = tr.mirror_facet(Facet(c, index_v1));
    const Facet opp_facet2 = tr.mirror_facet(Facet(c, index_v2));

    // volume data
    cells_info.insert(std::make_pair(opp_facet1, Cell_info{subdomain, selected}));
    cells_info.insert(std::make_pair(opp_facet2, Cell_info{subdomain, selected}));
    treat_before_delete(c, cell_selector, c3t3);

    // surface data for facets of the cells to be split
    const int findex = CGAL::Triangulation_utils_3::next_around_edge(index_v1, index_v2);
    Surface_patch_index patch = c3t3.surface_patch_index(c, findex);
    Vertex_handle opp_vertex = c->vertex(findex);
    facets_info.insert(std::make_pair(opp_facet1, Facet_info{opp_vertex, patch}));
    facets_info.insert(std::make_pair(opp_facet2, Facet_info{opp_vertex, patch}));

    if(c3t3.is_in_complex(c, findex))
      c3t3.remove_from_complex(c, findex);
  }

  // insert midpoint
  Vertex_handle new_v = tr.tds().insert_in_edge(e);
  new_v->set_point(typename Tr::Point(m));
  new_v->set_dimension(dimension);

  // update c3t3 with subdomain and surface patch indices
  std::vector<Cell_handle> new_cells;
  tr.incident_cells(new_v, std::back_inserter(new_cells));
  for (Cell_handle new_cell : new_cells)
  {
    const Facet fi(new_cell, new_cell->index(new_v));
    const Facet mfi = tr.mirror_facet(fi);

    //get subdomain info back
    CGAL_assertion(cells_info.find(mfi) != cells_info.end());
    Cell_info c_info = cells_info.at(mfi);
    treat_new_cell(new_cell, c_info.subdomain_index_,
                   cell_selector, c_info.selected_, c3t3);

    // get surface info back
    CGAL_assertion(facets_info.find(mfi) != facets_info.end());
    const Facet_info v_and_opp_patch = facets_info.at(mfi);

    // facet opposite to new_v (status wrt c3t3 is unchanged)
    new_cell->set_surface_patch_index(new_cell->index(new_v),
                                      mfi.first->surface_patch_index(mfi.second));

    // new half-facet (added or not to c3t3 depending on the stored surface patch index)
    if (Surface_patch_index() == v_and_opp_patch.patch_index_)
      new_cell->set_surface_patch_index(new_cell->index(v_and_opp_patch.opp_vertex_),
                                        Surface_patch_index());
    else
      c3t3.add_to_complex(new_cell,
                          new_cell->index(v_and_opp_patch.opp_vertex_),
                          v_and_opp_patch.patch_index_);

    // newly created internal facet
    for (int i = 0; i < 4; ++i)
    {
      const Vertex_handle vi = new_cell->vertex(i);
      if (vi == v1 || vi == v2)
      {
        new_cell->set_surface_patch_index(i, Surface_patch_index());
        break;
      }
    }

    //the 4th facet (new_v, v_and_opp_patch.first, v1 or v2)
    // will have its patch tagged from the other side, if needed
  }

  // re-insert complex sub-edges
  if (dimension == 1)
  {
    c3t3.add_to_complex(new_v, v1, curve_index);
    c3t3.add_to_complex(new_v, v2, curve_index);
  }

  set_index(new_v, c3t3);

  return new_v;
}

/**
* returns [can_be_split, is_on_boundary]
*/
template<typename C3T3, typename CellSelector>
auto can_be_split(const typename C3T3::Edge& e,
                  const C3T3& c3t3,
                  const bool protect_boundaries,
                  const CellSelector& cell_selector)
{
  struct Splittable
  {
    bool can_be_split;
    bool on_boundary;
  };

  if (is_outside(e, c3t3, cell_selector))
    return Splittable{false, false};

  const bool boundary = c3t3.is_in_complex(e)
                     || is_boundary(c3t3, e, cell_selector);

  if (protect_boundaries)
  {
    if (boundary)
      return Splittable{false, boundary};

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    if (!is_internal(e, c3t3, cell_selector))
    {
      std::cerr << "e is not inside!?" << std::endl;
      typename C3T3::Vertex_handle v1 = e.first->vertex(e.second);
      typename C3T3::Vertex_handle v2 = e.first->vertex(e.third);
      std::cerr << v1->point() << " " << v2->point() << std::endl;
    }
#endif

    CGAL_assertion(is_internal(e, c3t3, cell_selector));
    return Splittable{true, boundary};
  }
  else
  {
    return Splittable{is_selected(e, c3t3.triangulation(), cell_selector), boundary};
  }
}



template<typename C3t3,
         typename SizingFunction,
         typename CellSelector,
         typename Visitor>
class EdgeSplitOperation
    : public ElementaryOperation<C3t3,
                                 std::pair<typename C3t3::Triangulation::Geom_traits::FT,
                                           std::pair<typename C3t3::Triangulation::Vertex_handle,
                                                     typename C3t3::Triangulation::Vertex_handle>>,
                                 std::vector<std::pair<typename C3t3::Triangulation::Geom_traits::FT,
                                                       std::pair<typename C3t3::Triangulation::Vertex_handle,
                                                                 typename C3t3::Triangulation::Vertex_handle>>>>
{
public:
  using Tr = typename C3t3::Triangulation;
  using Vertex_handle = typename Tr::Vertex_handle;
  using Cell_handle = typename Tr::Cell_handle;
  using Edge = typename Tr::Edge;
  using Edge_vv = std::pair<Vertex_handle, Vertex_handle>;
  using FT = typename Tr::Geom_traits::FT;

  using Long_edge = std::pair<FT, Edge_vv>;
  using Long_edges_with_lengths = std::vector<Long_edge>;
  using BaseOperation = ElementaryOperation<C3t3, Long_edge, Long_edges_with_lengths>;
  using ElementType = typename BaseOperation::ElementType;
  static_assert(std::is_same_v<ElementType, Long_edge>, "ElementType must be Long_edge");
  using ElementSource = typename BaseOperation::ElementSource;

private:
  const SizingFunction& m_sizing;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;
  Visitor& m_visitor;

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
  mutable std::ofstream m_can_be_split_ofs;
  mutable std::ofstream m_split_failed_ofs;
  mutable std::ofstream m_midpoints_ofs;
#endif

public:
  EdgeSplitOperation(const SizingFunction& sizing,
                     const CellSelector& cell_selector,
                     const bool protect_boundaries,
                     Visitor& visitor)
      : m_sizing(sizing)
      , m_cell_selector(cell_selector)
      , m_protect_boundaries(protect_boundaries)
      , m_visitor(visitor) {}

  ElementSource get_element_source(const C3t3& c3t3) const override
  {
    Long_edges_with_lengths long_edges;
    const Tr& tr = c3t3.triangulation();

    for (Edge e : tr.finite_edges())
    {
      auto [splittable, boundary] = can_be_split(e, c3t3, m_protect_boundaries, m_cell_selector);
      if (!splittable)
        continue;

      const std::optional<FT> sqlen = is_too_long(e, boundary, m_sizing, c3t3, m_cell_selector);
      if (sqlen != std::nullopt)
        long_edges.push_back(std::make_pair(*sqlen, make_vertex_pair(e)));
    }

    // longest first; stable to match the original bimap's ordering
    std::stable_sort(long_edges.begin(), long_edges.end(),
                     [](const std::pair<FT, Edge_vv>& a, const std::pair<FT, Edge_vv>& b) {
                       return a.first > b.first;
                     });

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    {
      std::ofstream ofs("long_edges.polylines.txt");
      for (const auto& le : long_edges)
        ofs << "2 " << point(le.second.first->point())
            << " " << point(le.second.second->point()) << std::endl;
    }
    m_can_be_split_ofs.open("can_be_split_edges.polylines.txt");
    m_split_failed_ofs.open("split_failed.polylines.txt");
    m_midpoints_ofs.open("midpoints.off");
    m_midpoints_ofs << "OFF" << std::endl;
    m_midpoints_ofs << long_edges.size() << " 0 0" << std::endl;
#endif
    return long_edges;
  }

  bool execute_operation(const ElementType& element, C3t3& c3t3) override
  {
    Tr& tr = c3t3.triangulation();
    const Edge_vv& e = element.second;

    Cell_handle cell;
    int i1, i2;
    if (!tr.tds().is_edge(e.first, e.second, cell, i1, i2))
      return false;

    Edge edge(cell, i1, i2);

    // check that splittability has not changed
    if (!can_be_split(edge, c3t3, m_protect_boundaries, m_cell_selector).can_be_split)
      return false;
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    m_can_be_split_ofs << "2 " << edge.first->vertex(edge.second)->point()
                        << " " << edge.first->vertex(edge.third)->point() << std::endl;
#endif

    m_visitor.before_split(tr, edge);
    Vertex_handle vh = split_edge(edge, m_cell_selector, c3t3);

    if (vh == Vertex_handle())
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      m_split_failed_ofs << "2 " << edge.first->vertex(edge.second)->point() << " "
                         << edge.first->vertex(edge.third)->point() << std::endl;
#endif
      return false;
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    m_midpoints_ofs << vh->point() << std::endl;
#endif
    m_visitor.after_split(tr, vh);
    return true;
  }

  std::string operation_name() const override { return "Split long edges"; }
};

} // internal
} // Tetrahedral_remeshing
} // CGAL

#endif // CGAL_INTERNAL_SPLIT_LONG_EDGES_H
