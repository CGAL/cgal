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
  void generate_input_cube(const std::size_t& n, Tr& tr)
  {
    typedef typename Tr::Point Point;
    CGAL::Random rng;

    // points in a sphere
    std::vector<Point> pts;
    while (pts.size() < n)
      pts.push_back(Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));
    tr.insert(pts.begin(), pts.end());

    for (auto v : tr.finite_vertex_handles())
      v->set_dimension(3);

    // vertices of a larger cube
    const std::array<Point, 8> pc = { Point(-2., -2., -2.),
                                      Point(-2., -2., 2.),
                                      Point(2., -2., -2.),
                                      Point(2., -2., 2.),
                                      Point(-2., 2., -2.),
                                      Point(-2., 2., 2.),
                                      Point(2., 2., -2.),
                                      Point(2., 2., 2.) };

    for (const Point& p : pc)
    {
      tr.insert(p)->set_dimension(0);
    }

    assert(tr.is_valid(true));

    // set subdomain indices
    for (auto c : tr.finite_cell_handles())
      c->set_subdomain_index(1);

    // set surface patches
    for (typename Tr::Facet f : tr.finite_facets())
    {
      typename Tr::Facet mf = tr.mirror_facet(f);
      if(tr.is_infinite(f.first) || tr.is_infinite(mf.first))
      {
        f.first->set_surface_patch_index(f.second, 2);
        mf.first->set_surface_patch_index(mf.second, 2);
      }
    }
  }
}
}
