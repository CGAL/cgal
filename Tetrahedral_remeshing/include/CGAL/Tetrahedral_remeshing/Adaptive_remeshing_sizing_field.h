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

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <vector>
#include <array>


namespace CGAL
{
namespace Tetrahedral_remeshing
{
/**
 * @class Adaptive_remeshing_sizing_field
 * @tparam Tr a triangulation
 *
 * @todo add template parameter to take an input aabb_tree
 */
template <typename Tr>
class Adaptive_remeshing_sizing_field
  : public Sizing_field<typename Tr::Geom_traits>
{
  // Types
  typedef typename Tr::Geom_traits              GT;
  typedef typename Tr::Geom_traits::Point_3     Bare_point;
  typedef typename Tr::Point                    Tr_point;
  typedef typename GT::FT                       FT;

  typedef typename Tr::Vertex_handle            Vertex_handle;
  typedef typename Tr::Cell_handle              Cell_handle;

  struct Point_with_info
  {
    Bare_point p;
    FT size;
    int dimension;
  };

private:
  struct Point_property_map
  {
    using Self = Point_property_map;
    using value_type = Bare_point;
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

  using Triangle_vec = std::vector<typename Tr::Triangle>;
  using Triangle_iter = typename Triangle_vec::iterator;
  using Triangle_primitive = CGAL::AABB_triangle_primitive<GT, Triangle_iter>;
  using AABB_triangle_traits = CGAL::AABB_traits<GT, Triangle_primitive>;
  using AABB_triangle_tree = CGAL::AABB_tree<AABB_triangle_traits>;

public:
  /**
  * Constructor
  */
  Adaptive_remeshing_sizing_field(const Tr& tr)
    : m_gt(tr.geom_traits())
    , m_kd_tree(points_with_info(tr), Splitter(), Kd_traits(Point_property_map()))
  {
    m_kd_tree.build();
    build_aabb_trees(tr);
  }

private:
  std::vector<Point_with_info> points_with_info(const Tr& tr) const
  {
    auto cp = tr.geom_traits().construct_point_3_object();
    std::vector<Point_with_info> points;
    for (const Vertex_handle v : tr.finite_vertex_handles())
    {
      points.push_back(Point_with_info{ cp(tr.point(v)),
                                        average_circumradius_around(v, tr),
                                        v->in_dimension() });
    }
    return points;
  }

public:
  /**
  * Returns size at point `p`
  */
  template <typename Index>
  FT operator()(const Bare_point& p, const int& dim, const Index& i) const
  {
    // Find nearest vertex and local size before remeshing
    Point_property_map pp_map;
    Distance dist(pp_map);
    Neighbor_search search(m_kd_tree,
                           p, //query point
                           1, //nb nearest neighbors
                           0, //epsilon
                           true, //search nearest
                           dist);
    const auto [pi, size, dimension] = search.begin()->first;

    if (dim < 3)
      return size;

    // measure distance to input surfaces
    Bare_point closest_point = p;
    FT shortest_distance = (std::numeric_limits<FT>::max)();
    Surface_patch_index closest_patch{};
    for(std::size_t i = 0; i < m_aabb_trees.size(); ++i)
    {
      const Bare_point closest = m_aabb_trees[i].closest_point(p);
      const FT sq_dist = m_gt.compute_squared_distance_3_object()(p, closest);
      if(sq_dist < shortest_distance)
      {
        shortest_distance = sq_dist;
        closest_point = closest;
        closest_patch = m_i2p[i];
      }
    }
    shortest_distance = CGAL::approximate_sqrt(shortest_distance);
    FT div_max_distance = 1. / 10.;
      //?? (10 as a magic number for distance_field.getMaxDist();)

    return size - size / (shortest_distance * div_max_distance);
  }

private:
  /**
  * Fills sizing field, using size associated to points in `tr_`
  */
  auto build_kd_tree(const Tr& tr);

  /**
  * Fills aabb trees, using triangles in `tr_`, for projection to input surfaces
  */
  void build_aabb_trees(const Tr& tr);

  /**
   * Returns size at point `p`, by interpolation into tetrahedron.
   */
  FT interpolate_on_four_vertices(
    const Bare_point& p,
    const std::array<Point_with_info, 4>& vertices) const;

  FT sq_circumradius_length(const Cell_handle cell, const Vertex_handle v, const Tr& tr) const;
  FT average_circumradius_around(const Vertex_handle v, const Tr& tr) const;

private:
  Kd_tree m_kd_tree;

  using Surface_patch_index = typename Tr::Cell::Surface_patch_index;
  std::map<Surface_patch_index, std::size_t> m_p2i;
  std::vector<Surface_patch_index> m_i2p;
  std::vector<Triangle_vec> m_triangles;
  std::vector<AABB_triangle_tree> m_aabb_trees;

  const GT& m_gt;
};


template <typename Tr>
void
Adaptive_remeshing_sizing_field<Tr>::
build_aabb_trees(const Tr& tr)
{
  // collect patch indices, and triangles for each patch
  for (const auto& f : tr.finite_facets())
  {
    if (!f.first->is_facet_on_surface(f.second))
      continue;

    const Surface_patch_index patch = f.first->surface_patch_index(f.second);
    if (m_p2i.find(patch) == m_p2i.end())
    {
      m_p2i.insert({patch, m_aabb_trees.size()});
      m_i2p.push_back(patch);

      m_triangles.emplace_back();
      m_triangles.back().push_back(tr.triangle(f));
    }
    else
      m_triangles[m_p2i[patch]].push_back(tr.triangle(f));
  }

  CGAL_assertion(m_triangles.size() == m_i2p.size());
  CGAL_assertion(m_triangles.size() == m_p2i.size());

  for (std::size_t i = 0; i < m_triangles.size(); ++i)
  {
    m_aabb_trees.push_back(AABB_triangle_tree(m_triangles[i].begin(), m_triangles[i].end()));
    m_aabb_trees.back().build();
    m_aabb_trees.back().accelerate_distance_queries();
  }
}

template <typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
interpolate_on_four_vertices(
  const Bare_point& p,
  const std::array<Point_with_info, 4>& vertices) const
{
  // Interpolate value using values at vertices
  const FT& va = boost::get<1>(vertices[0]);
  const FT& vb = boost::get<1>(vertices[1]);
  const FT& vc = boost::get<1>(vertices[2]);
  const FT& vd = boost::get<1>(vertices[3]);

  const Bare_point& a = boost::get<0>(vertices[0]);
  const Bare_point& b = boost::get<0>(vertices[1]);
  const Bare_point& c = boost::get<0>(vertices[2]);
  const Bare_point& d = boost::get<0>(vertices[3]);

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


//template <typename Tr>
//typename Adaptive_remeshing_sizing_field<Tr>::FT
//Adaptive_remeshing_sizing_field<Tr>::
//interpolate_on_facet_vertices(const Bare_point& p, const Cell_handle& cell) const
//{
//  typename GT::Compute_area_3 area =  tr_.geom_traits().compute_area_3_object();
//
//  typename GT::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();
//  // Find infinite vertex and put it in k0
//  int k0 = 0;
//  int k1 = 1;
//  int k2 = 2;
//  int k3 = 3;
//
//  if ( tr_.is_infinite(cell->vertex(1)) )
//    std::swap(k0,k1);
//  if ( tr_.is_infinite(cell->vertex(2)) )
//    std::swap(k0,k2);
//  if ( tr_.is_infinite(cell->vertex(3)) )
//    std::swap(k0,k3);
//
//  // Interpolate value using tet vertices values
//  const FT& va = cell->vertex(k1)->meshing_info();
//  const FT& vb = cell->vertex(k2)->meshing_info();
//  const FT& vc = cell->vertex(k3)->meshing_info();
//
//  const Tr_point& wa = tr_.point(cell, k1);
//  const Tr_point& wb = tr_.point(cell, k2);
//  const Tr_point& wc = tr_.point(cell, k3);
//  const Bare_point& a = cp(wa);
//  const Bare_point& b = cp(wb);
//  const Bare_point& c = cp(wc);
//
//  const FT abp = area(a, b, p);
//  const FT acp = area(a, c, p);
//  const FT bcp = area(b, c, p);
//
//  CGAL_assertion(abp >= 0);
//  CGAL_assertion(acp >= 0);
//  CGAL_assertion(bcp >= 0);
//
//  // If area is 0, then compute the average value
//  if ( is_zero(abp+acp+bcp) )
//    return (va+vb+vc)/3.;
//
//  return ( (abp*vc + acp*vb + bcp*va ) / (abp+acp+bcp) );
//}

template <typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
sq_circumradius_length(const Cell_handle cell,
                       const Vertex_handle v,
                       const Tr& tr) const
{
  auto cp  = tr.geom_traits().construct_point_3_object();
  auto sq_distance = tr.geom_traits().compute_squared_distance_3_object();
  auto cc = tr.geom_traits().construct_circumcenter_3_object();

  const auto t = tr.tetrahedron(cell);
  const Bare_point circumcenter = cc(t[0], t[1], t[2], t[3]);
  const Tr_point& position = tr.point(cell, cell->index(v));

  return sq_distance(cp(position), circumcenter);
}

template <typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
average_circumradius_around(const Vertex_handle v, const Tr& tr) const
{
  std::vector<Cell_handle> incident_cells;
  incident_cells.reserve(64);
  tr.incident_cells(v, std::back_inserter(incident_cells));

  using SI = typename Tr::Triangulation_data_structure::Cell::Subdomain_index;

  FT sum_len(0);
  unsigned int nb = 0;

  for (Cell_handle c : incident_cells)
  {
    if (c->subdomain_index() != SI())
    {
      sum_len += CGAL::approximate_sqrt(sq_circumradius_length(c, v, tr));
      ++nb;
    }
  }

  // nb == 0 could happen if there is an isolated point.
  if (0 != nb)
  {
    return sum_len / nb;
  }
  else
  {
    // Use outside cells to compute size of point
    for (Cell_handle c : incident_cells)
    {
      if (!tr.is_infinite(c))
      {
        sum_len += CGAL::approximate_sqrt(sq_circumradius_length(c, v, tr));
        ++nb;
      }
    }

    CGAL_assertion(nb != 0);
    CGAL_assertion(sum_len != 0);
    return sum_len / nb;
  }
}

} // end namespace Tetrahedral_remeshing

} //namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ADAPTIVE_SIZING_FIELD_H
