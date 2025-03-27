// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_3_TRIANGULATION_HELPERS_H
#define CGAL_MESH_3_TRIANGULATION_HELPERS_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/enum.h>
#include <CGAL/functional.h>
#include <CGAL/type_traits.h>
#include <CGAL/Time_stamper.h>

#include <boost/unordered_set.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>
#include <type_traits>


namespace CGAL {

namespace Mesh_3 {

enum class Helpers_tr_type {
  DELAUNAY = 0,
  REGULAR = 1,

  PERIODIC = 0,
  NON_PERIODIC = 1,

  ALL = 2
};

template <std::size_t N>
static constexpr auto indices = std::make_index_sequence<N>{};

template<typename Tr,
         Helpers_tr_type,
         Helpers_tr_type>
struct Triangulation_helpers_utils;

template <typename Tr, typename Derived, Helpers_tr_type>
struct Triangulation_helpers_utils_decorator;

template <typename Tr, typename Derived>
struct Triangulation_helpers_utils_decorator<Tr, Derived, Helpers_tr_type::NON_PERIODIC>
{
  using EK = Exact_predicates_exact_constructions_kernel;

  template <typename Obj>
  static decltype(auto) to_exact(const Tr&, Obj&& obj) {
    using Geom_traits = typename Tr::Geom_traits;
    using To_exact = Cartesian_converter<Geom_traits, EK>;
    return To_exact()(std::forward<Obj>(obj));
  }

  template <typename Obj>
  static decltype(auto) back_from_exact(const Tr&, Obj&& obj) {
    using Geom_traits = typename Tr::Geom_traits;
    using Back_from_exact = Cartesian_converter<EK, Geom_traits>;
    return Back_from_exact()(std::forward<Obj>(obj));
  }

  static auto exact_point(const Tr& tr, typename Tr::Vertex_handle v) {
    return to_exact(tr, tr.point(v));
  }

  template <typename Exact_points>
  static auto exact_circumcenter(const Tr& tr, const Exact_points& points) {
    auto circumcenter_fct = Derived::construct_exact_circumcenter_object(tr);
    return std::apply(circumcenter_fct, points);
  }

  template <typename Exact_points>
  static auto circumcenter(const Tr& tr, const Exact_points& points) {
    auto circumcenter_fct = Derived::construct_circumcenter_object(tr);
    return std::apply(circumcenter_fct, points);
  }

  static auto dual_exact(const Tr& tr, typename Tr::Cell_handle c)
  {
    auto circumcenter_fct = Derived::construct_exact_circumcenter_object(tr);
    return back_from_exact(tr, std::apply(circumcenter_fct, exact_points(tr, c)));
  }

  static auto dual_segment_exact(const Tr& tr, typename Tr::Facet facet)
  {
    auto [c, i] = facet;
    auto n = c->neighbor(i);

    return std::make_pair(dual_exact(tr, c), dual_exact(tr, n));
  }

  static auto dual_segment(const Tr& tr, typename Tr::Facet facet)
  {
    auto [c, i] = facet;
    auto n = c->neighbor(i);

    return std::make_pair(Derived::dual(tr, c), Derived::dual(tr, n));
  }

  // for a facet on the convex hull, return its version viewed from the
  // inside of the triangulation
  static auto dual_ray_aux(const Tr& tr, typename Tr::Facet facet) {
    auto [c, i] = facet;
    auto n = c->neighbor(i);
    auto in = n->index(c);

    if(tr.is_infinite(c)) {
      CGAL_precondition(!tr.is_infinite(n));
      return std::make_pair(n, in);
    } else {
      CGAL_precondition(tr.is_infinite(n));
      return facet;
    }
  }

  static auto dual_ray_exact(const Tr& tr, typename Tr::Facet f) {
    auto cstr_exact_plane = EK().construct_plane_3_object();
    auto cstr_exact_perpendicular_line = EK().construct_perpendicular_line_3_object();
    auto cstr_exact_ray = EK().construct_ray_3_object();

    const auto facet_inside = dual_ray_aux(tr, f);
    const auto facet_exact_points = exact_points(tr, facet_inside);

    auto [cell_inside, index] = facet_inside;
    auto fourth_exact_point = exact_point(tr, cell_inside->vertex(index));
    using Point = CGAL::cpp20::remove_cvref_t<decltype(facet_exact_points[0])>;
    using CRef = std::add_lvalue_reference_t<std::add_const_t<Point>>;
    std::tuple<CRef, CRef, CRef, CRef> cell_exact_points_const_ref{facet_exact_points[0],
                                                                   facet_exact_points[1],
                                                                   facet_exact_points[2],
                                                                   fourth_exact_point};
    const auto facet_exact_center = exact_circumcenter(tr, facet_exact_points);
    const auto exact_plane =
        std::apply(cstr_exact_plane, Derived::exact_bare_points(tr, facet_exact_points));
    const auto exact_line = cstr_exact_perpendicular_line(exact_plane, facet_exact_center);

    const auto exact_ray =
        cstr_exact_ray(exact_circumcenter(tr, cell_exact_points_const_ref), exact_line);
    return back_from_exact(tr, exact_ray);
  }

  static auto dual_ray(const Tr& tr, typename Tr::Facet facet) {
    auto cstr_ray = tr.geom_traits().construct_ray_3_object();
    auto cstr_perpendicular_line = tr.geom_traits().construct_perpendicular_line_3_object();
    auto cstr_plane = tr.geom_traits().construct_plane_3_object();

    const auto facet_inside = dual_ray_aux(tr, facet);
    const auto facet_points = points(tr, facet_inside);
    const auto facet_center = circumcenter(tr, facet_points);
    const auto plane = std::apply(cstr_plane, Derived::bare_points(tr, facet_points));
    const auto line = cstr_perpendicular_line(plane, facet_center);
    const auto ray = cstr_ray(Derived::dual(tr, facet_inside.first), line);
    return ray;
  }

  template <std::size_t... index, typename Vertices>
  static auto exact_points_of_vertices(const Tr& tr, const Vertices& vertices,
                                       std::index_sequence<index...>)
  {
    return std::array{exact_point(tr, vertices[index])...};
  }

  template <typename ...Args>
  static auto exact_points(const Tr& tr, Args&&... args) {
    auto vertices = tr.vertices(std::forward<Args>(args)...);

    constexpr auto N = vertices.size();
    return exact_points_of_vertices(tr, vertices, indices<N>);
  }

  template <std::size_t... index, typename Vertices>
  static auto points_of_vertices(const Tr& tr, const Vertices& vertices,
                                std::index_sequence<index...>)
  {
    return std::array{tr.point(vertices[index])...};
  }

  template <typename ...Args>
  static auto points(const Tr& tr, Args&&... args) {
    auto vertices = tr.vertices(std::forward<Args>(args)...);

    constexpr auto N = vertices.size();
    return points_of_vertices(tr, vertices, indices<N>);
  }

  template <typename Vector>
  static void set_point(const Tr&, typename Tr::Vertex_handle v, const Vector& /*move*/,
                        const typename Tr::Point& new_position)
  {
    v->set_point(new_position);
  }

  template <typename Facet, typename Vertex_handle>
  static auto get_incident_triangle(const Tr& tr, const Facet& f, const Vertex_handle)
  {
    return tr.triangle(f);
  }

  template <typename Bare_point>
  static const Bare_point& get_closest_point(const Tr&, const Bare_point& /*p*/, const Bare_point& q)
  {
    return q;
  }

  template <typename Bare_point>
  static auto min_squared_distance(const Tr& tr, const Bare_point& p, const Bare_point& q)
  {
    return tr.geom_traits().compute_squared_distance_3_object()(p, q);
  }
};

template <typename Tr>
struct Triangulation_helpers_utils<Tr, Helpers_tr_type::NON_PERIODIC, Helpers_tr_type::REGULAR>
    : public Triangulation_helpers_utils_decorator<
          Tr,
          Triangulation_helpers_utils<Tr, Helpers_tr_type::NON_PERIODIC, Helpers_tr_type::REGULAR>,
          Helpers_tr_type::NON_PERIODIC>
{
  using Self = Triangulation_helpers_utils<Tr, Helpers_tr_type::NON_PERIODIC, Helpers_tr_type::REGULAR>;

  // returns a callable that constructs a point from a bare point
  static auto construct_triangulation_point_object(const Tr& tr)
  {
    return tr.geom_traits().construct_weighted_point_3_object();
  }

  // returns a callable that constructs a bare point from a weighted point
  static auto construct_bare_point_object(const Tr& tr)
  {
    return tr.geom_traits().construct_point_3_object();
  }

  static auto compute_weight_object(const Tr& tr)
  {
    return tr.geom_traits().compute_weight_3_object();
  }

  static auto construct_exact_bare_point_object(const Tr&)
  {
    typename Self::EK ek;
    return ek.construct_point_3_object();
  }

  static auto construct_exact_circumcenter_object(const Tr&)
  {
    typename Self::EK ek;
    return ek.construct_weighted_circumcenter_3_object();
  }

  static auto construct_circumcenter_object(const Tr& tr)
  {
    return tr.geom_traits().construct_weighted_circumcenter_3_object();
  }

  template <std::size_t... index, typename Points>
  static auto exact_bare_points_aux(const Tr& tr, const Points& points, std::index_sequence<index...>) {
    return std::array{construct_exact_bare_point_object(tr)(points[index])...};
  };

  template <typename Exact_points>
  static auto exact_bare_points(const Tr& tr, const Exact_points& points) {
    constexpr auto N = points.size();
    return exact_bare_points_aux(tr, points, indices<N>);
  }

  template <std::size_t... index, typename Points>
  static auto bare_points_aux(const Tr& tr, const Points& points, std::index_sequence<index...>) {
    return std::array{construct_bare_point_object(tr)(points[index])...};
  };

  template <typename Points>
  static auto bare_points(const Tr& tr, const Points& points) {
    constexpr auto N = points.size();
    return bare_points_aux(tr, points, indices<N>);
  }

  static auto dual(const Tr& tr, typename Tr::Cell_handle c)
  {
    return c->weighted_circumcenter(tr.geom_traits());
  }

  template <typename ...Args>
  static bool greater_or_equal_power_distance(const Tr& tr, Args&&... args)
  {
    auto compare_power_distance =
        tr.geom_traits().compare_power_distance_3_object();
    return compare_power_distance(std::forward<Args>(args)...) != CGAL::SMALLER;
  }

  template <typename ...Args>
  static auto side_of_power_sphere(const Tr& tr, Args&&... args)
  {
    return tr.side_of_power_sphere(std::forward<Args>(args)...);
  }
};

template <typename Tr>
struct Triangulation_helpers_utils<Tr, Helpers_tr_type::NON_PERIODIC, Helpers_tr_type::DELAUNAY>
    : public Triangulation_helpers_utils_decorator<
          Tr,
          Triangulation_helpers_utils<Tr, Helpers_tr_type::NON_PERIODIC, Helpers_tr_type::DELAUNAY>,
          Helpers_tr_type::NON_PERIODIC>
{
  using Self = Triangulation_helpers_utils<Tr, Helpers_tr_type::NON_PERIODIC, Helpers_tr_type::DELAUNAY>;
  static auto construct_triangulation_point_object(const Tr&)
  {
    return CGAL::Identity<>();
  }

  static auto construct_bare_point_object(const Tr&)
  {
    return CGAL::Identity<>();
  }

  static auto compute_weight_object(const Tr&)
  {
    using FT = typename Tr::Geom_traits::FT;
    return [](const auto&) { return FT(0); };
  }

  static auto construct_exact_bare_point_object(const Tr&)
  {
    return CGAL::Identity<>();
  }

  static auto construct_exact_circumcenter_object(const Tr&)
  {
    typename Self::EK ek;
    return ek.construct_circumcenter_3_object();
  }

  static auto construct_circumcenter_object(const Tr& tr)
  {
    return tr.geom_traits().construct_circumcenter_3_object();
  }

  template <typename Points>
  static decltype(auto) exact_bare_points(const Tr&, Points&& points)
  {
    return std::forward<Points>(points);
  }

  template <typename Points>
  static decltype(auto) bare_points(const Tr&, Points&& points)
  {
    return std::forward<Points>(points);
  }

  static auto dual(const Tr& tr, typename Tr::Cell_handle c)
  {
    return c->circumcenter(tr.geom_traits());
  }

  template <typename ...Args>
  static bool greater_or_equal_power_distance(const Tr& tr, Args&&... args)
  {
    auto compare_distance =
        tr.geom_traits().compare_distance_3_object();
    return compare_distance(std::forward<Args>(args)...) != CGAL::SMALLER;
  }

  template <typename ...Args>
  static auto side_of_power_sphere(const Tr& tr, Args&&... args)
  {
    return tr.side_of_sphere(std::forward<Args>(args)...);
  }
};

template <typename Tr>
class Triangulation_helpers
    : public Triangulation_helpers_utils<
          Tr,
          is_periodic_triangulation_v<Tr> ? Helpers_tr_type::PERIODIC : Helpers_tr_type::NON_PERIODIC,
          is_regular_triangulation_v<Tr> ? Helpers_tr_type::REGULAR : Helpers_tr_type::DELAUNAY>
{
  typedef typename Tr::Geom_traits              GT;

  typedef typename GT::FT                       FT;
  typedef typename GT::Vector_3                 Vector_3;

  // If `Tr` is not a triangulation that has defined Bare_point,
  // use Point_3 as defined in the traits class.
  typedef Bare_point_type_t<Tr>                 Bare_point;

  // 'Point' is either a bare point or a weighted point, depending on the triangulation.
  // Since 'Triangulation_helpers' can be templated by an unweighted triangulation,
  // this is one of the rare cases where we use 'Point'.
  typedef typename Tr::Point                    Point;

  typedef typename Tr::Vertex                   Vertex;
  typedef typename Tr::Vertex_handle            Vertex_handle;
  typedef typename Tr::Facet                    Facet;
  typedef typename Tr::Cell                     Cell;
  typedef typename Tr::Cell_handle              Cell_handle;
  typedef typename Tr::Cell_iterator            Cell_iterator;

  typedef std::vector<Cell_handle>              Cell_vector;

  /**
   * A functor to reset visited flag of each facet of a cell
   */
  struct Reset_facet_visited
  {
    void operator()(const Cell_handle& c) const
    {
      for ( int i=0; i<4 ; ++i )
        c->reset_visited(i);
    }
  };

  /**
   * A functor to get the point of a vertex vh, but replacing
   * it by m_p when vh == m_vh
   */
  struct Point_getter
  {
    /// When the requested will be about vh, the returned point will be p
    /// instead of vh->point()
    Point_getter(const Vertex_handle &vh, const Point&p)
      : m_vh(vh), m_p(p)
    {}

    const Point& operator()(const Vertex_handle &vh) const
    {
      return (vh == m_vh ? m_p : vh->point());
    }

  private:
    const Vertex_handle m_vh;
    const Point &m_p;
  };

public:
  /**
   * Returns `true` if moving `v` to `p` makes no topological
   * change in `tr`.
   */
  bool no_topological_change(Tr& tr,
                             const Vertex_handle v,
                             const Vector_3& move,
                             const Point& p,
                             Cell_vector& cells_tos) const;
  bool no_topological_change__without_set_point(
                             const Tr& tr,
                             const Vertex_handle v,
                             const Point& p,
                             Cell_vector& cells_tos) const;

  bool no_topological_change(Tr& tr,
                             const Vertex_handle v,
                             const Vector_3& move,
                             const Point& p) const;
  bool no_topological_change__without_set_point(
                             const Tr& tr,
                             const Vertex_handle v,
                             const Point& p) const;

  bool inside_protecting_balls(const Tr& tr,
                               const Vertex_handle v,
                               const Bare_point& p) const;

  /**
   * Returns the squared distance from `vh` to its closest vertex.
   *
   * \pre `vh` is not the infinite vertex
   */
  template<typename Tag> // Two versions to distinguish using 'Has_visited_for_vertex_extractor'
  FT get_sq_distance_to_closest_vertex(const Tr& tr,
                                       const Vertex_handle& vh,
                                       const Cell_vector& incident_cells,
                                       typename std::enable_if_t<Tag::value>* = nullptr) const;

  // @todo are the two versions really worth it, I can't tell the difference from a time POV...
  template<typename Tag>
  FT get_sq_distance_to_closest_vertex(const Tr& tr,
                                       const Vertex_handle& vh,
                                       const Cell_vector& incident_cells,
                                       typename std::enable_if_t<!Tag::value>* = nullptr) const;

private:
  /**
   * Returns `true` if `v` is well_oriented on each cell of `cell_tos`.
   */
  // For sequential version
  bool well_oriented(const Tr& tr,
                     const Cell_vector& cell_tos) const;
  // For parallel version
  bool well_oriented(const Tr& tr,
                     const Cell_vector& cell_tos,
                     const Point_getter& pg) const;
};

template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change(Tr& tr,
                      const Vertex_handle v0,
                      const Vector_3& move,
                      const Point& p,
                      Cell_vector& cells_tos) const
{
  // Check here that the Triangulation_helpers class is empty.
  // That cannot be check within the class definition, where the class is still incomplete.
  static_assert(std::is_empty_v<Triangulation_helpers<Tr>>,
    "Triangulation_helpers should be an empty utility class");

  if(tr.point(v0) == p)
    return true;

  // For periodic triangulations, calling this function actually takes longer than
  // just removing and re-inserting the point directly because the periodic mesh
  // triangulation's side_of_power_sphere() function is somewhat brute-forcy:
  // it has to check for 27 copies of the query point to determine correctly the side.
  // Thus, we simply return 'false' if the triangulation is a periodic triangulation.
  //
  // Note that the function was nevertheless adapted to work with periodic triangulation
  // so this hack can be disabled if one day 'side_of_power_sphere()' is improved.
  if(is_periodic_triangulation_v<Tr>)
    return false;

  typename GT::Construct_opposite_vector_3 cov =
      tr.geom_traits().construct_opposite_vector_3_object();

  bool np = true;
  const Point fp = tr.point(v0);

  // move the point
  this->set_point(tr, v0, move, p);

  if(!well_oriented(tr, cells_tos))
  {
    // Reset (restore) v0
    this->set_point(tr, v0, cov(move), fp);
    return false;
  }

  // Reset visited tags of facets
  std::for_each(cells_tos.begin(), cells_tos.end(), Reset_facet_visited());

  for ( typename Cell_vector::iterator cit = cells_tos.begin() ;
                                       cit != cells_tos.end() ;
                                       ++cit )
  {
    Cell_handle c = *cit;
    for(int j=0; j<4; j++)
    {
      // Treat each facet only once
      if(c->is_facet_visited(j))
        continue;

      // Set facet and its mirrored version as visited
      Cell_handle cj = c->neighbor(j);
      int mj = tr.mirror_index(c, j);
      c->set_facet_visited(j);
      cj->set_facet_visited(mj);

      if(tr.is_infinite(c->vertex(j)))
      {
        if(this->side_of_power_sphere(tr, c, tr.point(cj->vertex(mj)), false)
           != CGAL::ON_UNBOUNDED_SIDE)
        {
          np = false;
          break;
        }
      }
      else
      {
        if(this->side_of_power_sphere(tr, cj, tr.point(c->vertex(j)), false)
           != CGAL::ON_UNBOUNDED_SIDE)
        {
          np = false;
          break;
        }
      }
    }
  }

  // Reset (restore) v0
  this->set_point(tr, v0, cov(move), fp);

  return np;
}

template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change__without_set_point(
  const Tr& tr,
  const Vertex_handle v0,
  const Point& p,
  Cell_vector& cells_tos) const
{
  bool np = true;

  Point_getter pg(v0, p);

  if(!well_oriented(tr, cells_tos, pg))
  {
    return false;
  }

  // Reset visited tags of facets
  std::for_each(cells_tos.begin(), cells_tos.end(), Reset_facet_visited());

  for ( typename Cell_vector::iterator cit = cells_tos.begin() ;
        cit != cells_tos.end() ;
        ++cit )
  {
    Cell_handle c = *cit;
    for(int j=0; j<4; j++)
    {
      // Treat each facet only once
      if(c->is_facet_visited(j))
        continue;

      // Set facet and it's mirror's one visited
      Cell_handle cj = c->neighbor(j);
      int mj = tr.mirror_index(c, j);
      c->set_facet_visited(j);
      cj->set_facet_visited(mj);

      Vertex_handle v1 = c->vertex(j);
      typedef typename Tr::Triangulation_data_structure TDS;
      typedef typename TDS::Cell_range Cell_range;
      typedef typename TDS::Vertex_range Vertex_range;
      if(tr.is_infinite(v1))
      {
        // Build a copy of c, and replace V0 by a temporary vertex (position "p")
        typename Cell_handle::value_type c_copy (*c);
        Cell_range::Time_stamper_impl::initialize_time_stamp(&c_copy);
        int i_v0;
        typename Vertex_handle::value_type v;
        Vertex_range::Time_stamper_impl::initialize_time_stamp(&v);
        if (c_copy.has_vertex(v0, i_v0))
        {
          v.set_point(p);
          c_copy.set_vertex(i_v0, Vertex_range::s_iterator_to(v));
        }

        Cell_handle c_copy_h = Cell_range::s_iterator_to(c_copy);
        if(this->side_of_power_sphere(tr, c_copy_h, pg(cj->vertex(mj)), false)
           != CGAL::ON_UNBOUNDED_SIDE)
        {
          np = false;
          break;
        }
      }
      else
      {
        // Build a copy of cj, and replace V0 by a temporary vertex (position "p")
        typename Cell_handle::value_type cj_copy (*cj);
        Cell_range::Time_stamper_impl::initialize_time_stamp(&cj_copy);
        int i_v0;
        typename Vertex_handle::value_type v;
        Vertex_range::Time_stamper_impl::initialize_time_stamp(&v);
        if (cj_copy.has_vertex(v0, i_v0))
        {
          v.set_point(p);
          cj_copy.set_vertex(i_v0, Vertex_range::s_iterator_to(v));
        }

        Cell_handle cj_copy_h = Cell_range::s_iterator_to(cj_copy);
        if(this->side_of_power_sphere(tr, cj_copy_h, pg(v1), false)
           != CGAL::ON_UNBOUNDED_SIDE)
        {
          np = false;
          break;
        }
      }
    }
  }

  return np;
}


template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change(Tr& tr,
                      const Vertex_handle v0,
                      const Vector_3& move,
                      const Point& p) const
{
  Cell_vector cells_tos;
  cells_tos.reserve(64);
  tr.incident_cells(v0, std::back_inserter(cells_tos));
  return no_topological_change(tr, v0, move, p, cells_tos);
}

template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change__without_set_point(
                      const Tr& tr,
                      const Vertex_handle v0,
                      const Point& p) const
{
  Cell_vector cells_tos;
  cells_tos.reserve(64);
  tr.incident_cells(v0, std::back_inserter(cells_tos));
  return no_topological_change__without_set_point(tr, v0, p, cells_tos);
}


template<typename Tr>
bool
Triangulation_helpers<Tr>::
inside_protecting_balls(const Tr& tr,
                        const Vertex_handle v,
                        const Bare_point& p) const
{
  if constexpr (is_regular_triangulation_v<Tr>) {
    if(tr.number_of_vertices() == 0)
      return false;

    typename GT::Compare_weighted_squared_radius_3 cwsr =
      tr.geom_traits().compare_weighted_squared_radius_3_object();

    Cell_handle hint = (v == Vertex_handle()) ? Cell_handle() : v->cell();
    Vertex_handle nv = tr.nearest_power_vertex(p, hint);
    const Point& nvwp = tr.point(nv);

    if(cwsr(nvwp, FT(0)) == CGAL::SMALLER)
    {
      auto cp = construct_bare_point_object(tr);
      const Point& nvwp = tr.point(nv);
      // 'true' if the distance between 'p' and 'nv' is smaller or equal than the weight of 'nv'
      return (cwsr(nvwp , - this->min_squared_distance(tr, p, cp(nvwp))) != CGAL::LARGER);
    }
  } // end if constexpr (is_regular_triangulation_v<Tr>)
  return false;
}

/// Return the squared distance from vh to its closest vertex
/// if `Has_visited_for_vertex_extractor` is `true`
template<typename Tr>
template<typename Tag>
typename Triangulation_helpers<Tr>::FT
Triangulation_helpers<Tr>::
get_sq_distance_to_closest_vertex(const Tr& tr,
                                  const Vertex_handle& vh,
                                  const Cell_vector& incident_cells,
                                  typename std::enable_if_t<Tag::value>*) const
{
  CGAL_precondition(!tr.is_infinite(vh));

  typedef std::vector<Vertex_handle>              Vertex_container;

  // There is no need to use tr.min_squared_distance() here because we are computing
  // distances between 'v' and a neighboring vertex within a common cell, which means
  // that even if we are using a periodic triangulation, the distance is correctly computed.
  typename GT::Compute_squared_distance_3 csqd = tr.geom_traits().compute_squared_distance_3_object();
  auto cp = construct_bare_point_object(tr);

  Vertex_container treated_vertices;
  FT min_sq_dist = std::numeric_limits<FT>::infinity();

  for(typename Cell_vector::const_iterator cit = incident_cells.begin();
                                           cit != incident_cells.end(); ++cit)
  {
    const Cell_handle c = (*cit);
    const int k = (*cit)->index(vh);
    const Point& wpvh = tr.point(c, k);

    // For each vertex of the cell
    for(int i=1; i<4; ++i)
    {
      const int n = (k+i)&3;
      const Vertex_handle& vn = c->vertex(n);

      if(vn == Vertex_handle() ||
         tr.is_infinite(vn) ||
         vn->visited_for_vertex_extractor)
        continue;

      vn->visited_for_vertex_extractor = true;
      treated_vertices.push_back(vn);

      const Point& wpvn = tr.point(c, n);
      const FT sq_d = csqd(cp(wpvh), cp(wpvn));

      if(sq_d < min_sq_dist)
        min_sq_dist = sq_d;
    }
  }

  for(std::size_t i=0; i < treated_vertices.size(); ++i)
    treated_vertices[i]->visited_for_vertex_extractor = false;

  return min_sq_dist;
}

/// Return the squared distance from vh to its closest vertex
/// if `Has_visited_for_vertex_extractor` is `false`
template<typename Tr>
template<typename Tag>
typename Triangulation_helpers<Tr>::FT
Triangulation_helpers<Tr>::
get_sq_distance_to_closest_vertex(const Tr& tr,
                                  const Vertex_handle& vh,
                                  const Cell_vector& incident_cells,
                                  typename std::enable_if_t<!Tag::value>*) const
{
  CGAL_precondition(!tr.is_infinite(vh));

  typedef CGAL::Hash_handles_with_or_without_timestamps      Hash_fct;
  typedef boost::unordered_set<Vertex_handle, Hash_fct>      Vertex_container;
  typedef typename Vertex_container::iterator                VC_it;

  // There is no need to use tr.min_squared_distance() here because we are computing
  // distances between 'v' and a neighboring vertex within a common cell, which means
  // that even if we are using a periodic triangulation, the distance is correctly computed.
  typename GT::Compute_squared_distance_3 csqd = tr.geom_traits().compute_squared_distance_3_object();
  auto cp = construct_bare_point_object(tr);

  Vertex_container treated_vertices;
  FT min_sq_dist = std::numeric_limits<FT>::infinity();

  for(typename Cell_vector::const_iterator cit = incident_cells.begin();
                                           cit != incident_cells.end(); ++cit)
  {
    const Cell_handle c = (*cit);
    const int k = (*cit)->index(vh);
    const Point& wpvh = tr.point(c, k);

    // For each vertex of the cell
    for(int i=1; i<4; ++i)
    {
      const int n = (k+i)&3;
      const Vertex_handle& vn = c->vertex(n);

      if(vn == Vertex_handle() ||
         tr.is_infinite(vn))
        continue;

      std::pair<VC_it, bool> is_insert_successful = treated_vertices.insert(vn);
      if(! is_insert_successful.second) // vertex has already been treated
        continue;

      const Point& wpvn = tr.point(c, n);
      const FT sq_d = csqd(cp(wpvh), cp(wpvn));

      if(sq_d < min_sq_dist)
        min_sq_dist = sq_d;
    }
  }

  return min_sq_dist;
}


/// This function well_oriented is called by no_topological_change after the
/// position of the vertex has been (tentatively) modified.
template<typename Tr>
bool
Triangulation_helpers<Tr>::
well_oriented(const Tr& tr,
              const Cell_vector& cells_tos) const
{
  typedef typename Tr::Geom_traits GT;
  typename GT::Orientation_3 orientation = tr.geom_traits().orientation_3_object();
  auto cp = this->construct_bare_point_object(tr);

  typename Cell_vector::const_iterator it = cells_tos.begin();
  for( ; it != cells_tos.end() ; ++it)
  {
    Cell_handle c = *it;
    if( tr.is_infinite(c) )
    {
      int iv = c->index(tr.infinite_vertex());
      Cell_handle cj = c->neighbor(iv);
      int mj = tr.mirror_index(c, iv);

      const Point& mjwp = tr.point(cj, mj);
      const Point& fwp1 = tr.point(c, (iv+1)&3);
      const Point& fwp2 = tr.point(c, (iv+2)&3);
      const Point& fwp3 = tr.point(c, (iv+3)&3);

      if(orientation(cp(mjwp), cp(fwp1), cp(fwp2), cp(fwp3)) != CGAL::NEGATIVE)
        return false;
    }
    else
    {
      const Point& cwp0 = tr.point(c, 0);
      const Point& cwp1 = tr.point(c, 1);
      const Point& cwp2 = tr.point(c, 2);
      const Point& cwp3 = tr.point(c, 3);

      if(orientation(cp(cwp0), cp(cwp1), cp(cwp2), cp(cwp3)) != CGAL::POSITIVE)
      return false;
    }
  }
  return true;
}

/// Another version for the parallel version
/// Here, the set_point is not done before, but we use a Point_getter instance
/// to get the point of a vertex.
template<typename Tr>
bool
Triangulation_helpers<Tr>::
well_oriented(const Tr& tr,
              const Cell_vector& cells_tos,
              const Point_getter& pg) const
{
  typedef typename Tr::Geom_traits GT;
  typename GT::Orientation_3 orientation = tr.geom_traits().orientation_3_object();
  auto cp = this->construct_bare_point_object(tr);

  typename Cell_vector::const_iterator it = cells_tos.begin();
  for( ; it != cells_tos.end() ; ++it)
  {
    Cell_handle c = *it;
    if( tr.is_infinite(c) )
    {
      int iv = c->index(tr.infinite_vertex());
      Cell_handle cj = c->neighbor(iv);
      int mj = tr.mirror_index(c, iv);
      if(orientation(cp(pg(cj->vertex(mj))),
                     cp(pg(c->vertex((iv+1)&3))),
                     cp(pg(c->vertex((iv+2)&3))),
                     cp(pg(c->vertex((iv+3)&3)))) != CGAL::NEGATIVE)
        return false;
    }
    else if(orientation(cp(pg(c->vertex(0))),
                        cp(pg(c->vertex(1))),
                        cp(pg(c->vertex(2))),
                        cp(pg(c->vertex(3)))) != CGAL::POSITIVE)
      return false;
  }
  return true;
}

} // end namespace Mesh_3

} //namespace CGAL

#endif // CGAL_MESH_3_TRIANGULATION_HELPERS_H
