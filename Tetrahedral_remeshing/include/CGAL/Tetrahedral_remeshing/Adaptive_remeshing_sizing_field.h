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

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <vector>
#include <array>


namespace CGAL
{
namespace Tetrahedral_remeshing
{
/**
 * @class Adaptive_remeshing_sizing_field
 * @tparam Tr a triangulation
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

private:
  using Point_and_size = boost::tuple<Bare_point, FT>;
  using Traits = CGAL::Search_traits_adapter<Point_and_size,
    CGAL::Nth_of_tuple_property_map<0, Point_and_size>,
    CGAL::Search_traits_3<GT> >;
  using K_neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits>;
  using Distance = typename K_neighbor_search::Distance;
  using Tree = typename K_neighbor_search::Tree;

public:
  /**
  * Constructor
  */
  Adaptive_remeshing_sizing_field(const Tr& tr)
    : m_gt(tr.geom_traits())
  {
    build_kd_tree(tr);
  }

  /**
  * Returns size at point `p`
  */
  template <typename Index>
  FT operator()(const Bare_point& p, const int& dim, const Index& i) const
  {
    // Find nearest vertex
    Distance tr_dist;
    K_neighbor_search search(m_kd_tree, p, 4/*nb nearest neighbors*/);
    std::array<Point_and_size, 4> vertices;
    int vi = 0;
    for (typename K_neighbor_search::iterator it = search.begin();
         it != search.end();
         ++it, ++vi)
    {
      const auto& pt_size = it->first;
      vertices[vi] = pt_size;
    }
    return interpolate_on_four_vertices(p, vertices);
  }

private:
  /**
  * Fills sizing field, using size associated to points in `tr_`
  */
  void build_kd_tree(const Tr& tr);

  /**
   * Returns size at point `p`, by interpolation into tetrahedron.
   */
  FT interpolate_on_four_vertices(
    const Bare_point& p,
    const std::array<Point_and_size, 4>& vertices) const;

  ///**
  // * Returns size at point `p`, by interpolation into facet (`cell` is assumed
  // * to be an infinite cell).
  // */
  //FT interpolate_on_facet_vertices(const Bare_point& p,
  //                                 const Cell_handle& cell) const;

  FT sq_circumradius_length(const Cell_handle cell, const Vertex_handle v, const Tr& tr) const;
  FT average_circumradius_length(const Vertex_handle v, const Tr& tr) const;

private:
  Tree m_kd_tree;
  const GT& m_gt;
};


template <typename Tr>
void
Adaptive_remeshing_sizing_field<Tr>::
build_kd_tree(const Tr& tr)
{
  auto cp = m_gt.construct_point_3_object();

  std::vector<Bare_point> points;
  std::vector<FT>         sizes;

  // Fill Kd tree with local size
  for (const Vertex_handle v : tr.finite_vertex_handles())
  {
    const Tr_point& position = tr.point(v);
    points.push_back(cp(position));
    sizes.push_back(average_circumradius_length(v, tr));
  }

  m_kd_tree.insert(boost::make_zip_iterator(boost::make_tuple(points.begin(), sizes.begin())),
                   boost::make_zip_iterator(boost::make_tuple(points.end(), sizes.end())));
}

template <typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
interpolate_on_four_vertices(
  const Bare_point& p,
  const std::array<Point_and_size, 4>& vertices) const
{
  auto volume = m_gt.compute_volume_3_object();

  // Interpolate value using tet vertices values
  const FT& va = boost::get<1>(vertices[0]);
  const FT& vb = boost::get<1>(vertices[1]);
  const FT& vc = boost::get<1>(vertices[2]);
  const FT& vd = boost::get<1>(vertices[3]);

  const Bare_point& a = boost::get<0>(vertices[0]);
  const Bare_point& b = boost::get<0>(vertices[1]);
  const Bare_point& c = boost::get<0>(vertices[2]);
  const Bare_point& d = boost::get<0>(vertices[3]);

  const FT abcp = CGAL::abs(volume(a,b,c,p));
  const FT abdp = CGAL::abs(volume(a,d,b,p));
  const FT acdp = CGAL::abs(volume(a,c,d,p));
  const FT bcdp = CGAL::abs(volume(b,d,c,p));

  // If volume is 0, then compute the average value
  if ( is_zero(abcp+abdp+acdp+bcdp) )
    return (va+vb+vc+vd) / 4.;

  return ( (abcp*vd + abdp*vc + acdp*vb + bcdp*va) / (abcp+abdp+acdp+bcdp) );
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
average_circumradius_length(const Vertex_handle v, const Tr& tr) const
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
