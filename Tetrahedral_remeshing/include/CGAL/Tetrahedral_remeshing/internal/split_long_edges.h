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

#include <CGAL/Tetrahedral_remeshing/internal/Elementary_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <unordered_map>
#include <functional>
#include <utility>
#include <optional>
#include <array>

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

// Where a split may put the vertex it creates, as fractions of the edge.
//
// `construct_steiner_point()` walks this list when the midpoint would leave a
// cell inverted. The parallel lock zone enumerates the same values, because
// the spatial lock is keyed on POSITION and has to hold the grid cell the new
// vertex lands in -- see `Edge_split_operation::lock_zone()`. The two must
// stay in step, so there is one list.
template<typename FT>
const std::array<FT, 6>& steiner_coefficients()
{
  static const std::array<FT, 6> coeff = {0.33, 0.66,    //1/3 and 2/3
                                          0.3, 0.7,      // 0.5 +/- 0.2
                                          0.25, 0.75};   // 0.5 +/- 0.25
  return coeff;
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

  const auto& coeff = steiner_coefficients<FT>();

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
class Edge_split_operation
    : public Elementary_operation<C3t3,
                                 std::pair<typename C3t3::Triangulation::Vertex_handle,
                                           typename C3t3::Triangulation::Vertex_handle>,
                                 std::vector<std::pair<typename C3t3::Triangulation::Vertex_handle,
                                                       typename C3t3::Triangulation::Vertex_handle>>>
{
public:
  using Tr = typename C3t3::Triangulation;
  using Vertex_handle = typename Tr::Vertex_handle;
  using Cell_handle = typename Tr::Cell_handle;
  using Edge = typename Tr::Edge;
  using Edge_vv = std::pair<Vertex_handle, Vertex_handle>;
  using FT = typename Tr::Geom_traits::FT;

  // Candidates are stored as vertex pairs, captured at collection time. A raw
  // Edge (Cell_handle, i, j) would go stale: each split destroys and recycles
  // cells, so by the time the executor reaches a later candidate its Cell_handle
  // may point at a different cell. Vertices are never removed by a split, so the
  // vertex pair stays valid and is re-resolved to the current edge via is_edge().
  using Long_edges = std::vector<Edge_vv>;
  using Base_operation = Elementary_operation<C3t3, Edge_vv, Long_edges>;
  using Element_type = typename Base_operation::Element_type;
  static_assert(std::is_same_v<Element_type, Edge_vv>, "Element_type must be Edge_vv");
  using ElementSource = typename Base_operation::Element_range;

private:
  // Every position `split_edge()` may give the new vertex: the midpoint first,
  // then the Steiner fallbacks, which all lie on the edge between 0.25 and
  // 0.75 of it. Taking a point twice is cheap -- the grid cell is already
  // this thread's -- and on a grid coarser than the edge they are one cell.
  static bool lock_split_destinations(const Element_type& element, const Tr& tr)
  {
    using FT = typename Tr::Geom_traits::FT;
    const auto& gt = tr.geom_traits();
    const auto& p1 = point(element.first->point());
    const auto& p2 = point(element.second->point());

    if (!tr.try_lock_point(gt.construct_midpoint_3_object()(p1, p2)))
      return false;

    const auto vec = gt.construct_vector_3_object()(p1, p2);
    for (const FT c : steiner_coefficients<FT>())
    {
      if (!tr.try_lock_point(gt.construct_translated_point_3_object()(
              p1, gt.construct_scaled_vector_3_object()(vec, c))))
        return false;
    }
    return true;
  }

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
  Edge_split_operation(const SizingFunction& sizing,
                     const CellSelector& cell_selector,
                     const bool protect_boundaries,
                     Visitor& visitor)
      : m_sizing(sizing)
      , m_cell_selector(cell_selector)
      , m_protect_boundaries(protect_boundaries)
      , m_visitor(visitor) {}

  ElementSource get_elements(const C3t3& c3t3) const override
  {
    struct Long_edge_with_length
    {
      Edge edge;
      FT sqlength;
    };
    std::vector<Long_edge_with_length> long_edges_with_lengths;
    const Tr& tr = c3t3.triangulation();

    for (Edge e : tr.finite_edges())
    {
      auto [splittable, boundary] = can_be_split(e, c3t3, m_protect_boundaries, m_cell_selector);
      if (!splittable)
        continue;

      const std::optional<FT> sqlen = is_too_long(e, boundary, m_sizing, c3t3, m_cell_selector);
      if (sqlen != std::nullopt)
        long_edges_with_lengths.push_back(Long_edge_with_length{e, sqlen.value()});
    }

    // longest first; stable to match the original bimap's ordering
    std::stable_sort(long_edges_with_lengths.begin(), long_edges_with_lengths.end(),
                     [](const Long_edge_with_length& a, const Long_edge_with_length& b) {
                       return a.sqlength > b.sqlength;
                     });

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    {
      std::ofstream ofs("long_edges.polylines.txt");
      for (const auto& le : long_edges_with_lengths)
        ofs << "2 " << point(le.edge.first->point())
            << " " << point(le.edge.second->point()) << std::endl;
    }
    m_can_be_split_ofs.open("can_be_split_edges.polylines.txt");
    m_split_failed_ofs.open("split_failed.polylines.txt");
    m_midpoints_ofs.open("midpoints.off");
    m_midpoints_ofs << "OFF" << std::endl;
    m_midpoints_ofs << long_edges_with_lengths.size() << " 0 0" << std::endl;
#endif

    Long_edges long_edges;
    long_edges.reserve(long_edges_with_lengths.size());
    for(const auto& ef : long_edges_with_lengths)
      long_edges.push_back(make_vertex_pair(ef.edge));
    return long_edges;
  }

  bool execute_operation(const Element_type& element, C3t3& c3t3) override
  {
    Tr& tr = c3t3.triangulation();
    const Edge_vv& e = element;

    Cell_handle cell;
    int i1, i2;
    if constexpr (is_parallel)
    {
      // lock_zone() ran first and located the edge in order to lock its ring.
      // No match means it found no cell carrying the edge, i.e. this pair is
      // no longer an edge and there is nothing to split. Reusing what it found
      // also keeps the star walk off the parallel path: tds().is_edge() MARKS
      // every cell it visits, and the star reaches past this zone.
      const auto& located = last_located_edge<Vertex_handle, Cell_handle>();
      if (!located.matches(e.first, e.second))
        return false;
      cell = located.c;
      i1 = located.i0;
      i2 = located.i1;
    }
    else if (!tr.tds().is_edge(e.first, e.second, cell, i1, i2))
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

  // `lock_zone()` is called by the parallel executor alone; the sequential one
  // goes straight to `execute_operation()`, which must then locate the edge
  // itself.
  static constexpr bool is_parallel
    = std::is_convertible_v<typename Tr::Concurrency_tag, CGAL::Parallel_tag>;

  /**
  * Locks the grid cells of every position `execute_operation()` may write, and
  * nothing else.
  *
  * The zone taken by default is the union of the two endpoint stars. A split
  * destroys and recreates the RING -- the cells incident to the edge -- and
  * re-stitches the MIRROR cells across the ring's outer facets through
  * `set_neighbor()`; the cells created have only ring vertices, one ring apex
  * and the new vertex as corners. Nothing outside ring u mirror is written,
  * and, once the edge has been located here, nothing outside it is read
  * either: `can_be_split()` and `split_edge()` both work off a
  * `Cell_circulator` around the edge and the mirror facets it names. That is
  * a strict subset of the two stars, and it is reached in O(ring) rather than
  * O(star).
  *
  * That footprint is measured, not reasoned from the source. Two instruments,
  * each blind to what the other sees, agree on it: a snapshot/diff of every
  * byte of every object in a radius-3 ball, over 741 splits on two meshes,
  * found every net write at graph distance 1 and not one on a star cell that
  * is not facet-adjacent to the ring; and an instruction-level trace of 10
  * splits, with every cell and vertex of the triangulation labelled so that an
  * access outside the ball would be reported rather than missed, found no load
  * and no store anywhere outside ring u mirror.
  *
  * A cell is protected by holding its four vertices, so ring u mirror is
  * locked cell by cell.
  *
  * The two endpoints are locked first for two reasons: `v->cell()` is only
  * dependable once `v` is held, and holding `element.first` is what makes the
  * search safe to read -- every cell it visits contains that vertex, and a
  * thread writing such a cell would have to hold all four of its vertices.
  *
  * The vertex a split CREATES needs the same protection as one a collapse
  * MOVES, and for the same reason: the lock grid is keyed on position, so
  * until the destination's grid cell is held another thread can be working
  * there. Split never did this, and the both-stars zone only got away with it
  * because the new position usually falls in a grid cell some star vertex
  * already covered: asserting, for every object a split changed, that the
  * thread held what the protocol demands, the shipped zone left 180 newly
  * created cells uncovered across 300 splits -- 22 of those 300 operations --
  * and holding the destination takes that to 0. The destination is
  * the midpoint, or, when the midpoint would invert a cell, one of the
  * fractions in `steiner_coefficients()`; which one is decided inside the
  * zone, so all of them are taken.
  *
  * Returning false leaves partial locks behind; the executor releases them all
  * before retrying.
  */
  bool lock_zone(const Element_type& element, const C3t3& c3t3) const
  {
    using Cell_circulator = typename Tr::Cell_circulator;
    const Tr& tr = c3t3.triangulation();

    last_located_edge<Vertex_handle, Cell_handle>().clear();

    if (!tr.try_lock_vertex(element.first) || !tr.try_lock_vertex(element.second))
      return false;

    if (!lock_split_destinations(element, tr))
      return false;

    Cell_handle edge_cell;
    const Vertex_handle other = element.second;
    if (!tr.find_first_incident_cell_threadsafe(element.first,
          [other](const Cell_handle c) { return c->has_vertex(other); },
          edge_cell))
      return true; // no longer an edge; execute_operation() will decline it

    const int i0 = edge_cell->index(element.first);
    const int i1 = edge_cell->index(element.second);
    const Edge edge(edge_cell, i0, i1);

    Cell_circulator circ = tr.incident_cells(edge);
    const Cell_circulator done = circ;
    do
    {
      const Cell_handle c = circ;
      // the ring cell, and the two cells across its outer facets
      if (!tr.try_lock_cell(c)
       || !tr.try_lock_cell(c->neighbor(c->index(element.first)))
       || !tr.try_lock_cell(c->neighbor(c->index(element.second))))
        return false;
    }
    while (++circ != done);

    last_located_edge<Vertex_handle, Cell_handle>()
      .set(element.first, element.second, edge_cell, i0, i1);
    return true;
  }

  // longest edge first is the point of the ordering built in get_elements()
  static constexpr bool requires_ordered_processing = true;

  std::string operation_name() const override { return "Split long edges"; }
};

} // internal
} // Tetrahedral_remeshing
} // CGAL

#endif // CGAL_INTERNAL_SPLIT_LONG_EDGES_H
