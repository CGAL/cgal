// Copyright (c) 2024 GeometryFactory (France).
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
// File Description : Defines a Kd-tree for the medial axis of a 3D triangulation
//******************************************************************************

#ifndef CGAL_TETRAHEDRAL_REMESHING_MEDIAL_AXIS_KD_TREE_H
#define CGAL_TETRAHEDRAL_REMESHING_MEDIAL_AXIS_KD_TREE_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <CGAL/Delaunay_triangulation_3.h>

#include <vector>
#include <array>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
/**
 * @class Medial_axis_kd_tree
 * @tparam Tr a triangulation
 *
 * @todo add template parameter to take an input aabb_tree
 */
template <typename Tr>
class Medial_axis_kd_tree
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

public:
    /**
    * Constructor
    */
  Medial_axis_kd_tree(const Tr & tr)
    : m_kd_tree(poles(tr), Splitter(), Kd_traits(Point_property_map()))
  {
    m_kd_tree.build();
  }

private:
  std::vector<Point_with_info> poles(const Tr& tr) const
  {
    CGAL::Delaunay_triangulation_3<GT> dt;

    auto tr_cp = tr.geom_traits().construct_point_3_object();
    auto dt_cp = dt.geom_traits().construct_point_3_object();

    for(auto v : tr.finite_vertex_handles())
    {
      if(v->in_dimension() < 3)
        dt.insert(dt_cp(v->point()));
    }

    for(auto f : tr.finite_facets())
    {
      if(f.first->is_facet_on_surface(f.second))
      {
        auto t = tr.triangle(f);
        dt.insert(dt_cp(CGAL::centroid(tr_cp(t[0]), tr_cp(t[1]), tr_cp(t[2]))));
        dt.insert(dt_cp(CGAL::midpoint(tr_cp(t[0]), tr_cp(t[1]))));
        dt.insert(dt_cp(CGAL::midpoint(tr_cp(t[0]), tr_cp(t[2]))));
        dt.insert(dt_cp(CGAL::midpoint(tr_cp(t[1]), tr_cp(t[2]))));
      }
    }

    std::ofstream os("dt.xyz");
    for(auto v : dt.finite_vertex_handles())
      os << dt_cp(v->point()) << std::endl;
    os.close();

    std::vector<Point_with_info> points;
    for (auto c : dt.finite_cell_handles())
    {
      auto dual_pt_in_dt = dt_cp(dt.dual(c));
      Bare_point p = tr_cp(dual_pt_in_dt);
      points.push_back(Point_with_info{p});
    }

    std::ofstream ofs("medial_axis.xyz");
    for (auto p : points)
      ofs << p.p << std::endl;
    ofs.close();

    return points;
  }

public:
  FT distance_to_medial_axis(const Bare_point& p) const
  {
    const int nb_nearest_neighbors = 5;
    Point_property_map pp_map;
    Distance dist(pp_map);
    Neighbor_search search(m_kd_tree,
                           p, //query point
                           nb_nearest_neighbors, //nb nearest neighbors
                           0, //epsilon
                           true, //search nearest
                           dist);

    typename GT::Vector_3 sump = CGAL::NULL_VECTOR;
    for (auto it = search.begin(); it != search.end(); ++it)
      sump += typename GT::Vector_3(CGAL::ORIGIN, it->first.p);
    sump = sump / nb_nearest_neighbors;

    Bare_point psump = CGAL::ORIGIN + sump;
    return CGAL::approximate_sqrt(CGAL::squared_distance(p, psump));

//    const auto pi = search.begin()->first;
//    return CGAL::approximate_sqrt(CGAL::squared_distance(p, pi.p));
  }

private:
  Kd_tree m_kd_tree;
};
//end class Medial_axis_kd_tree

} // end namespace internal
} // end namespace Tetrahedral_remeshing
} //namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_MEDIAL_AXIS_KD_TREE_H
