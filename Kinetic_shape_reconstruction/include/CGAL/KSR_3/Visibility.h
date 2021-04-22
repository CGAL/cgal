// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
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
  typename  GeomTraits,
  typename  PointMap_3,
  typename VectorMap_3>
  class Visibility {

  public:
    using Kernel       = GeomTraits;
    using Point_map_3  = PointMap_3;
    using Vector_map_3 = VectorMap_3;

    using FT       = typename Kernel::FT;
    using Point_3  = typename Kernel::Point_3;
    using Vector_3 = typename Kernel::Vector_3;
    using Indices  = std::vector<std::size_t>;

    using Data_structure = KSR_3::Data_structure<Kernel>;
    using PFace          = typename Data_structure::PFace;
    using Volume_cell    = typename Data_structure::Volume_cell;

    using IK         = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IPoint_3   = typename IK::Point_3;
    using Delaunay_3 = CGAL::Delaunay_triangulation_3<IK>;
    using Generator  = CGAL::Random_points_in_tetrahedron_3<IPoint_3>;
    using Converter  = CGAL::Cartesian_converter<Kernel, IK>;

    using Visibility_label = KSR::Visibility_label;

    Visibility(
      const Data_structure& data,
      const std::map<PFace, Indices>& pface_points,
      const Point_map_3& point_map_3,
      const Vector_map_3& normal_map_3) :
    m_data(data),
    m_pface_points(pface_points),
    m_point_map_3(point_map_3),
    m_normal_map_3(normal_map_3),
    m_num_samples(100),
    m_random(0) {
      CGAL_assertion(m_pface_points.size() > 0);
    }

    void compute(std::vector<Volume_cell>& volumes) const {

      CGAL_assertion(volumes.size() > 0);
      if (volumes.size() == 0) return;
      for (auto& volume : volumes) {
        estimate_volume_label(volume);
      }
    }

  private:
    const Data_structure& m_data;
    const std::map<PFace, Indices>& m_pface_points;
    const Point_map_3& m_point_map_3;
    const Vector_map_3& m_normal_map_3;
    const std::size_t m_num_samples;
    Random m_random;

    void estimate_volume_label(Volume_cell& volume) const {

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
      const Volume_cell& volume) const {

      std::size_t in = 0, out = 0;
      std::vector<Point_3> samples;
      create_samples(volume, samples);
      compute_stats( volume, samples, in, out);
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
      const Volume_cell& volume,
      std::vector<Point_3>& samples) const {

      samples.push_back(volume.centroid);
      if (true) return;

      // If we need more samples, we use Delaunay.
      const Converter converter;
      const auto& pvertices = volume.pvertices;
      Delaunay_3 delaunay_3;
      for (const auto& pvertex : pvertices) {
        CGAL_assertion(m_data.has_ivertex(pvertex));
        const auto ivertex = m_data.ivertex(pvertex);
        delaunay_3.insert(converter(m_data.point_3(ivertex)));
      }

      std::vector<IPoint_3> points;
      for (auto cit = delaunay_3.finite_cells_begin();
      cit != delaunay_3.finite_cells_end(); ++cit) {
        const auto& tet = delaunay_3.tetrahedron(cit);

        const FT volume_size = CGAL::abs(tet.volume());
        if (volume_size < KSR::tolerance<FT>()) {
          continue;
        }

        // Generator generator(tet, m_random);
        // std::copy_n(generator, m_num_samples, std::back_inserter(points));
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
      const Volume_cell& volume,
      const std::vector<Point_3>& samples,
      std::size_t& in, std::size_t& out) const {

      CGAL_assertion(samples.size() >= 1);
      for (const auto& sample : samples) {
        const bool success = handle_sample_point(volume, sample, in, out);
        if (!success) return;
      }
    }

    const bool handle_sample_point(
      const Volume_cell& volume,
      const Point_3& query,
      std::size_t& in, std::size_t& out) const {

      bool found = false;
      const auto& pfaces = volume.pfaces;
      for (const auto& pface : pfaces) {
        CGAL_assertion(m_pface_points.find(pface) != m_pface_points.end());
        const auto& indices = m_pface_points.at(pface);
        if (indices.size() == 0) continue;
        found = true;

        for (const std::size_t index : indices) {
          const auto& point  = get(m_point_map_3 , index);
          const auto& normal = get(m_normal_map_3, index);

          const Vector_3 vec(point, query);
          const FT dot_product = vec * normal;
          if (dot_product < FT(0)) {
            in  += 1;
          } else {
            out += 1;
          }
        }
      }

      if (!found) {
        out += 1; return false;
      }
      return true;
    }
  };

} // KSR_3
} // CGAL

#endif // CGAL_KSR_3_VISIBILITY_H
