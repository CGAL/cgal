// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_VISIBILITY_H
#define CGAL_KSR_3_VISIBILITY_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/Random.h>
#include <CGAL/assertions.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/point_generators_3.h>

// Internal includes.
#include <CGAL/KSR/enum.h>
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR/property_map.h>
#include <CGAL/KSR_3/Data_structure.h>

namespace CGAL {
namespace KSR_3 {

  template<
  typename PointMap_3,
  typename GeomTraits>
  class Visibility {

  public:
    using Point_map_3 = PointMap_3;
    using Kernel      = GeomTraits;

    using FT = typename Kernel::FT;
    using Point_2 = typename Kernel::Point_2;
    using Point_3 = typename Kernel::Point_3;
    using Indices = std::vector<std::size_t>;

    using Data_structure = KSR_3::Data_structure<Kernel>;
    using Volume_cell    = typename Data_structure::Volume_cell;

    using IK         = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IPoint_3   = typename IK::Point_3;
    using Delaunay_3 = CGAL::Delaunay_triangulation_3<IK>;
    using Generator  = CGAL::Random_points_in_tetrahedron_3<IPoint_3>;
    using Converter  = CGAL::Cartesian_converter<Kernel, IK>;

    using Visibility_label = KSR::Visibility_label;

    Visibility(
      const Data_structure& data,
      const std::vector<std::size_t>& roof_points,
      const Point_map_3& point_map_3) :
    m_data(data),
    m_roof_points(roof_points),
    m_point_map_3(point_map_3),
    m_num_samples(100),
    m_random(0)
    { }

    void compute(std::vector<Volume_cell>& volumes) {

      CGAL_assertion(volumes.size() > 0);
      if (volumes.size() == 0) return;
      for (auto& volume : volumes) {
        estimate_volume_label(volume);
      }
      CGAL_assertion_msg(false, "TODO: FINISH VISIBILITY!");
    }

  private:
    const Data_structure& m_data;
    const std::vector<std::size_t>& m_roof_points;
    const Point_map_3& m_point_map_3;

    const Converter m_converter;
    const std::size_t m_num_samples;
    Random m_random;

    void estimate_volume_label(Volume_cell& volume) {

      const auto stats = estimate_in_out_values(volume);
      CGAL_assertion(stats.first  >= FT(0) && stats.first  <= FT(1));
      CGAL_assertion(stats.second >= FT(0) && stats.second <= FT(1));
      CGAL_assertion(
        CGAL::abs(stats.first + stats.second - FT(1)) < KSR::tolerance<FT>());

      if (stats.first >= FT(1) / FT(2)) {
        volume.visibility = Visibility_label::INSIDE;
      } else {
        volume.visibility = Visibility_label::OUTSIDE;
      }
      volume.inside  = stats.first;
      volume.outside = stats.second;

      // std::cout << "visibility in/out: " <<
      //   volume.inside << "/" << volume.outside << std::endl;
    }

    const std::pair<FT, FT> estimate_in_out_values(
      const Volume_cell& volume) {

      std::size_t in = 0, out = 0;
      std::vector<Point_3> samples;
      create_samples(volume, samples);
      compute_stats(samples, in, out);
      if (in == 0 && out == 0) {
        in = 1; out = 1;
      }

      const FT tmp_in  = static_cast<FT>(in);
      const FT tmp_out = static_cast<FT>(out);
      const FT sum = tmp_in + tmp_out;
      CGAL_assertion(sum > FT(0));

      const FT final_in  = tmp_in  / sum;
      const FT final_out = tmp_out / sum;

      return std::make_pair(final_in, final_out);
    }

    void create_samples(
      const Volume_cell& polyhedron,
      std::vector<Point_3>& samples) {

      const auto& pvertices = polyhedron.pvertices;
      Delaunay_3 delaunay_3;
      for (const auto& pvertex : pvertices) {
        CGAL_assertion(m_data.has_ivertex(pvertex));
        const auto ivertex = m_data.ivertex(pvertex);
        delaunay_3.insert(m_converter(m_data.point_3(ivertex)));
      }

      std::vector<IPoint_3> points;
      for (auto cit = delaunay_3.finite_cells_begin();
      cit != delaunay_3.finite_cells_end(); ++cit) {
        const auto& tet = delaunay_3.tetrahedron(cit);

        const FT volume = CGAL::abs(tet.volume());
        if (volume < KSR::tolerance<FT>()) {
          continue;
        }

        Generator generator(tet, m_random);
        std::copy_n(generator, m_num_samples, std::back_inserter(points));
      }

      samples.clear();
      samples.reserve(points.size());
      for (const auto& point : points) {
        samples.push_back(Point_3(
          static_cast<FT>(point.x()),
          static_cast<FT>(point.y()),
          static_cast<FT>(point.z())));
      }
      CGAL_assertion(samples.size() == points.size());
      // std::cout << "num samples: " << samples.size() << std::endl;
    }

    void compute_stats(
      const std::vector<Point_3>& samples,
      std::size_t& in, std::size_t& out) {

      for (const auto& sample : samples) {
        handle_sample_point(sample, in, out);
      }
    }

    void handle_sample_point(
      const Point_3& query,
      std::size_t& in, std::size_t& out) {

      in  += 1;
      out += 1;
    }
  };

} // KSR_3
} // CGAL

#endif // CGAL_KSR_3_VISIBILITY_H
