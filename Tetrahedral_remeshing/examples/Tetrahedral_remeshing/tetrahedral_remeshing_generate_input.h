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
// Author(s)     : Jane Tournois, Noura Faraj

#include <CGAL/Random.h>

#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <utility>
#include <cassert>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
  template<typename Tr>
  void insert_random_points_in_cube(const std::size_t& nbv, Tr& tr)
  {
    CGAL::Random rng;
    typedef typename Tr::Point Point;
    std::vector<Point> pts;
    while (pts.size() < nbv)
    {
      const double x = rng.get_double(-1., 1.);
      const double y = rng.get_double(-1., 1.);
      const double z = rng.get_double(-1., 1.);

      pts.push_back(Point(x, y, z));
    }
    tr.insert(pts.begin(), pts.end());
  }

  template<typename Plane, typename Tr>
  void insert_points_on_plane(const Plane& plane, const std::size_t& nbv, Tr& tr)
  {
    CGAL::Random rng;
    typedef typename Tr::Point Point;
    std::vector<Point> pts;
    while (pts.size() < nbv)
    {
      const double x = rng.get_double(-1., 1.);
      const double y = rng.get_double(-1., 1.);
      const double z = rng.get_double(-1., 1.);

      pts.push_back(plane.projection(Point(x, y, z)));
    }
    tr.insert(pts.begin(), pts.end());
  }

  template<typename Tr>
  void add_edge(typename Tr::Vertex_handle v1,
    typename Tr::Vertex_handle v2,
    const Tr& tr,
    std::unordered_set<std::pair<typename Tr::Vertex_handle,
                                 typename Tr::Vertex_handle>,
                       boost::hash<std::pair<typename Tr::Vertex_handle,
                                 typename Tr::Vertex_handle>>>& constraints)
  {
    typename Tr::Cell_handle c;
    int i, j;
    if (tr.is_edge(v1, v2, c, i, j))
      constraints.insert(std::make_pair(v1, v2));
  }

  template<typename Tr>
  void make_constraints_from_cube_edges(
    Tr& tr,
    std::unordered_set<std::pair<typename Tr::Vertex_handle,
                                 typename Tr::Vertex_handle>,
                       boost::hash<std::pair<typename Tr::Vertex_handle,
                                             typename Tr::Vertex_handle>>
  >& constraints)
  {
    typedef typename Tr::Point Point;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Cell_handle Cell_handle;

    const Point p0(-2., -2., -2.);
    const Point p1(-2., -2., -2.);
    const Point p2(2., -2., -2.);
    const Point p3(2., -2., 2.);
    const Point p4(-2., 2., -2.);
    const Point p5(-2., 2., 2.);
    const Point p6(2., 2., -2.);
    const Point p7(2., 2., 2.);

    typename Tr::Locate_type lt;
    int li, lj;
    Cell_handle c = tr.locate(p0, lt, li, lj);
    Vertex_handle v0 = c->vertex(li);
    c = tr.locate(p1, lt, li, lj);
    Vertex_handle v1 = c->vertex(li);

    c = tr.locate(p2, lt, li, lj);
    Vertex_handle v2 = c->vertex(li);
    c = tr.locate(p3, lt, li, lj);
    Vertex_handle v3 = c->vertex(li);

    c = tr.locate(p4, lt, li, lj);
    Vertex_handle v4 = c->vertex(li);
    c = tr.locate(p5, lt, li, lj);
    Vertex_handle v5 = c->vertex(li);

    c = tr.locate(p6, lt, li, lj);
    Vertex_handle v6 = c->vertex(li);
    c = tr.locate(p7, lt, li, lj);
    Vertex_handle v7 = c->vertex(li);

    // constrain cube edges
    add_edge(v0, v1, tr, constraints);
    add_edge(v1, v2, tr, constraints);
    add_edge(v2, v3, tr, constraints);
    add_edge(v3, v0, tr, constraints);

    add_edge(v4, v5, tr, constraints);
    add_edge(v5, v6, tr, constraints);
    add_edge(v6, v7, tr, constraints);
    add_edge(v7, v4, tr, constraints);

    add_edge(v0, v4, tr, constraints);
    add_edge(v1, v5, tr, constraints);
    add_edge(v2, v6, tr, constraints);
    add_edge(v3, v7, tr, constraints);
  }

  template<typename Tr>
  void generate_input_cube(const std::size_t& n,
    Tr& tr,
    std::unordered_set<std::pair<typename Tr::Vertex_handle,
                                 typename Tr::Vertex_handle>,
                       boost::hash<std::pair<typename Tr::Vertex_handle,
                                             typename Tr::Vertex_handle>>   >& constraints)
  {
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Point Point;
    CGAL::Random rng;

    // points in a sphere
    std::vector<Point> pts;
    while (pts.size() < n)
      pts.push_back(Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));
    tr.insert(pts.begin(), pts.end());

    // vertices of a larger cube
    Vertex_handle v0 = tr.insert(Point(-2., -2., -2.));
    Vertex_handle v1 = tr.insert(Point(-2., -2., 2.));

    Vertex_handle v2 = tr.insert(Point(2., -2., -2.));
    Vertex_handle v3 = tr.insert(Point(2., -2., 2.));

    Vertex_handle v4 = tr.insert(Point(-2., 2., -2.));
    Vertex_handle v5 = tr.insert(Point(-2., 2., 2.));

    Vertex_handle v6 = tr.insert(Point(2., 2., -2.));
    Vertex_handle v7 = tr.insert(Point(2., 2., 2.));

    // constrain cube edges
    add_edge(v0, v1, tr, constraints);
    add_edge(v1, v2, tr, constraints);
    add_edge(v2, v3, tr, constraints);
    add_edge(v3, v0, tr, constraints);

    add_edge(v4, v5, tr, constraints);
    add_edge(v5, v6, tr, constraints);
    add_edge(v6, v7, tr, constraints);
    add_edge(v7, v4, tr, constraints);

    add_edge(v0, v4, tr, constraints);
    add_edge(v1, v5, tr, constraints);
    add_edge(v2, v6, tr, constraints);
    add_edge(v3, v7, tr, constraints);

    assert(tr.is_valid(true));
  }
}
}
