// Copyright (c) 2021 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_CANVAS_CANVAS_SEEDS_H
#define CGAL_CANVAS_CANVAS_SEEDS_H

#include <CGAL/Canvas/Metric.h>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <istream>
#include <string>
#include <vector>

namespace CGAL {
namespace Canvas {

template<typename Canvas>
class Canvas_seeds
{
public:
  using Geom_traits = typename Canvas::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;

  using Metric = Metric_base<Geom_traits>;

protected:
  // the seeds and their corresponding metrics
  std::vector<Point_3> m_seeds;
  std::vector<Metric> m_seeds_metrics;

  const Canvas& m_canvas;

public:
  Canvas_seeds(const Canvas& canvas)
    :
      m_seeds(),
      m_seeds_metrics(),
      m_canvas(canvas)
  { }

public:
  const Point_3& operator[](const int i) const { return m_seeds[i]; }
  Point_3& operator[] (const int i) { return m_seeds[i]; }
  std::size_t size() const { return m_seeds.size(); }

  const std::vector<Metric>& metrics() const { return m_seeds_metrics; }

  std::size_t insert_new_seed(const FT x, const FT y, const FT z)
  {
#if (VERBOSITY > 5)
    std::cout << "added new seed: " << x << " " << y << " " << z;
    std::cout << " (" << m_seeds.size() << ")" << std::endl;
#endif

    m_seeds.emplace_back(x, y, z);
    m_seeds_metrics.push_back(m_canvas.metric_field()->compute_metric(m_seeds.back()));

    return m_seeds.size();
  }

  std::size_t read_seeds(const char* seeds_str)
  {
    std::ifstream in(seeds_str);
    if(!in)
    {
      std::cerr << "Warning: could not open seed file" << std::endl;
      return 0.;
    }

    std::string word;
    std::size_t useless, nv, dim;
    FT r_x, r_y, r_z;

    in >> word >> useless; // MeshVersionFormatted i
    in >> word >> dim; // Dimension d
    in >> word >> nv;
    std::cout << "seeds nv: " << nv << std::endl;
    CGAL_assertion(dim == 3);

    std::size_t old_size = m_seeds.size();
    m_seeds.reserve(old_size + nv);
    m_seeds_metrics.reserve(old_size + nv);

    for(std::size_t i=0; i<nv; ++i)
    {
      in >> r_x >> r_y >> r_z >> useless;
      insert_new_seed(r_x, r_y, r_z);
    }

#if (VERBOSITY > 0)
    std::cout << "m_seeds: " << m_seeds.size() << std::endl;
#endif

    return (m_seeds.size() - old_size);
  }
};

} // namespace Canvas
} // namespace CGAL

#endif // CGAL_CANVAS_CANVAS_SEEDS_H
