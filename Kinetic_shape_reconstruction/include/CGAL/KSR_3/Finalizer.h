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
// Author(s)     : Simon Giraudot, Dmitry Anisimov

#ifndef CGAL_KSR_3_FINALIZER_H
#define CGAL_KSR_3_FINALIZER_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR/parameters.h>
#include <CGAL/KSR/conversions.h>

#include <CGAL/KSR_3/Data_structure.h>

namespace CGAL {
namespace KSR_3 {

template<typename GeomTraits>
class Finalizer {

public:
  using Kernel = GeomTraits;

private:
  using FT          = typename Kernel::FT;
  using Point_2     = typename Kernel::Point_2;
  using Point_3     = typename Kernel::Point_3;
  using Vector_2    = typename Kernel::Vector_2;
  using Vector_3    = typename Kernel::Vector_3;
  using Segment_3   = typename Kernel::Segment_3;
  using Line_3      = typename Kernel::Line_3;
  using Plane_3     = typename Kernel::Plane_3;
  using Direction_2 = typename Kernel::Direction_2;

  using Data_structure = KSR_3::Data_structure<Kernel>;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  using PVertex = typename Data_structure::PVertex;
  using PEdge   = typename Data_structure::PEdge;
  using PFace   = typename Data_structure::PFace;

  using Support_plane      = typename Data_structure::Support_plane;
  using Intersection_graph = typename Data_structure::Intersection_graph;
  using Volume_cell        = typename Data_structure::Volume_cell;

  using Mesh           = typename Data_structure::Mesh;
  using Vertex_index   = typename Data_structure::Vertex_index;
  using Face_index     = typename Data_structure::Face_index;
  using Edge_index     = typename Data_structure::Edge_index;
  using Halfedge_index = typename Data_structure::Halfedge_index;

  struct Vertex_info {
    bool tagged;
    PVertex pvertex;
    IVertex ivertex;
    Vertex_info() :
    tagged(false),
    pvertex(Data_structure::null_pvertex()),
    ivertex(Data_structure::null_ivertex())
    { }
  };

  struct Face_info {
    std::size_t index;
    std::size_t input;
    Face_info() :
    index(KSR::uninitialized()),
    input(KSR::uninitialized())
    { }
  };

  using Parameters     = KSR::Parameters_3<FT>;
  using Kinetic_traits = KSR::Kinetic_traits_3<Kernel>;

public:
  Finalizer(Data_structure& data, const Parameters& parameters) :
  m_data(data), m_parameters(parameters), m_kinetic_traits(parameters.use_hybrid_mode)
  { }

  void create_polyhedra() {

    std::cout.precision(20);
    CGAL_assertion(m_data.check_bbox());
    CGAL_assertion(m_data.check_interior());
    CGAL_assertion(m_data.check_vertices());
    CGAL_assertion(m_data.check_edges());
    create_volumes();
    CGAL_assertion(m_data.check_faces());
  }

  void clear() {
    // to be done
  }

private:
  Data_structure& m_data;
  const Parameters& m_parameters;
  Kinetic_traits m_kinetic_traits;

  /*******************************
  **     EXTRACTING VOLUMES     **
  ********************************/

  void create_volumes() {

    // Initialize an empty volume map.
    auto& volumes = m_data.volumes();
    auto& map_volumes = m_data.pface_neighbors();
    auto& volume_level_map = m_data.volume_level_map();

    volumes.clear();
    std::map<int, Point_3> centroids;
    map_volumes.clear();
    for (std::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      const auto pfaces = m_data.pfaces(i);
      for (const auto pface : pfaces)
        map_volumes[pface] = std::make_pair(-1, -1);
    }

    // First, traverse only boundary volumes.
    bool is_found_new_volume = false;
    std::size_t volume_size = 0;
    int num_volumes  = 0;
    int volume_index = 0;
    int volume_level = 0;
    for (std::size_t i = 0; i < 6; ++i) {
      const auto pfaces = m_data.pfaces(i);
      for (const auto pface : pfaces) {
        CGAL_assertion(pface.first < 6);
        std::tie(is_found_new_volume, volume_size) = traverse_boundary_volume(
          pface, volume_index, num_volumes, map_volumes, centroids);
        if (is_found_new_volume) {
          CGAL_assertion(m_data.check_volume(volume_index, volume_size, map_volumes));
          ++volume_index;
        }
      }
    }
    if (m_parameters.verbose) {
      std::cout << "* found boundary volumes: "<< volume_index << std::endl;
    }
    num_volumes = volume_index;
    CGAL_assertion(num_volumes > 0);
    volume_level_map[volume_level] = static_cast<std::size_t>(num_volumes);
    ++volume_level;

    // Then traverse all other volumes if any.
    std::vector<PFace> other_pfaces;
    for (std::size_t i = 6; i < m_data.number_of_support_planes(); ++i) {
      const auto pfaces = m_data.pfaces(i);
      for (const auto pface : pfaces) {
        CGAL_assertion(pface.first >= 6);
        other_pfaces.push_back(pface);
      }
    }

    std::sort(other_pfaces.begin(), other_pfaces.end(),
      [&](const PFace& pface1, const PFace& pface2) -> bool {
        const auto pedges1 = m_data.pedges_of_pface(pface1);
        const auto pedges2 = m_data.pedges_of_pface(pface2);
        return pedges1.size() > pedges2.size();
      }
    );

    bool quit = true;
    do {
      quit = true;
      const int before = volume_index;
      for (const auto& other_pface : other_pfaces) {
        std::tie(is_found_new_volume, volume_size) = traverse_interior_volume(
          other_pface, volume_index, num_volumes, map_volumes, centroids);
        if (is_found_new_volume) {
          quit = false;
          CGAL_assertion(m_data.check_volume(volume_index, volume_size, map_volumes));
          ++volume_index;
        }
      }
      const int after = volume_index;
      if (m_parameters.verbose) {
        std::cout << "* found interior volumes: "<< after - before << std::endl;
      }
      num_volumes = volume_index;
      CGAL_assertion(after >= before);
      if (after > before) {
        volume_level_map[volume_level] = static_cast<std::size_t>(after - before);
        ++volume_level;
      }

    } while (!quit);

    // Now, set final volumes and their neighbors.
    for (const auto& item : map_volumes) {
      const auto& pface = item.first;
      const auto& pair  = item.second;

      if (pair.first == -1) {
        dump_pface(m_data, pface, "face-debug");
        std::cout << "DEBUG face: " << m_data.str(pface) << " "   << std::endl;
        std::cout << "DEBUG  map: " << pair.first << " : " << pair.second << std::endl;
      }

      CGAL_assertion(pair.first != -1);
      if (volumes.size() <= static_cast<std::size_t>(pair.first))
        volumes.resize(pair.first + 1);
      volumes[pair.first].add_pface(pface, pair.second);

      if (pface.first < 6 && pair.second == -1) continue;
      CGAL_assertion(pair.second != -1);
      if (volumes.size() <= static_cast<std::size_t>(pair.second))
        volumes.resize(pair.second + 1);
      volumes[pair.second].add_pface(pface, pair.first);
    }
    for (auto& volume : volumes)
      create_cell_pvertices(volume);

    if (m_parameters.verbose) {
      std::cout << "* created volumes: " << volumes.size() << std::endl;
      if (m_parameters.export_all) dump_volumes(m_data, "final");
      for (std::size_t i = 0; i < volumes.size(); ++i) {
        const auto& volume = volumes[i];
        CGAL_assertion(volume.pfaces.size() > 3);
        if (m_parameters.debug) {
          std::cout <<
          " VOLUME "     << std::to_string(i) << ": "
          " pvertices: " << volume.pvertices.size() <<
          " pfaces: "    << volume.pfaces.size()    << std::endl;
        }
      }
    }

    CGAL_assertion(volumes.size() == centroids.size());
    for (std::size_t i = 0; i < volumes.size(); ++i) {
      auto& volume = volumes[i];
      volume.set_index(i);
      volume.set_centroid(centroids.at(i));
    }
  }

  const std::pair<bool, std::size_t> traverse_boundary_volume(
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::map<int, Point_3>& centroids) const {

    CGAL_assertion(num_volumes  == 0);
    CGAL_assertion(volume_index >= 0);
    if (pface.first >= 6) return std::make_pair(false, 0);
    CGAL_assertion(pface.first < 6);
    const auto& pair = map_volumes.at(pface);
    CGAL_assertion(pair.second == -1);
    if (pair.first != -1) return std::make_pair(false, 0);
    CGAL_assertion(pair.first == -1);

    std::deque<PFace> queue;
    queue.push_front(pface);

    Point_3 volume_centroid;
    std::size_t volume_size = 0;

    while (!queue.empty()) {
      // print_queue(volume_index, queue);
      const auto query = queue.front();
      queue.pop_front();
      propagate_pface(
        false, query, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    }

    if (m_parameters.debug) {
      std::cout << "- FOUND VOLUME " << volume_index << ", (SIZE/BARYCENTER): "
      << volume_size << " / " << volume_centroid << std::endl;
    }

    centroids[volume_index] = volume_centroid;
    return std::make_pair(true, volume_size);
  }

  const std::pair<bool, std::size_t> traverse_interior_volume(
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::map<int, Point_3>& centroids) const {

    CGAL_assertion(volume_index > 0);
    CGAL_assertion(volume_index >= num_volumes);

    if (pface.first < 6) return std::make_pair(false, 0);
    CGAL_assertion(pface.first >= 6);
    const auto& pair = map_volumes.at(pface);
    if (pair.second != -1) {
      CGAL_assertion(pair.first != -1);
      return std::make_pair(false, 0);
    }
    CGAL_assertion(pair.second == -1);
    if (pair.first == -1) {
      CGAL_assertion(pair.second == -1);
      return std::make_pair(false, 0);
    }
    CGAL_assertion(pair.first != -1);
    if (pair.first >= num_volumes) return std::make_pair(false, 0);
    CGAL_assertion(pair.first < num_volumes);

    std::deque<PFace> queue;
    queue.push_front(pface);

    Point_3 volume_centroid;
    std::size_t volume_size = 0;

    while (!queue.empty()) {
      // print_queue(volume_index, queue);
      const auto query = queue.front();
      queue.pop_front();
      propagate_pface(
        false, query, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    }
    if (m_parameters.debug) {
      std::cout << "- FOUND VOLUME " << volume_index << ", (SIZE/BARYCENTER): "
      << volume_size << " / " << volume_centroid << std::endl;
    }
    centroids[volume_index] = volume_centroid;
    return std::make_pair(true, volume_size);
  }

  void print_queue(
    const int volume_index,
    const std::deque<PFace>& queue) const {

    // if (volume_index != -1) return;
    std::cout << "QUEUE: " << std::endl;
    for (const auto& pface : queue) {
      std::cout << volume_index << " "
      << pface.first << " " << pface.second << std::endl;
    }
  }

  void propagate_pface(
    const bool verbose,
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<int, Point_3>& centroids,
    std::size_t& volume_size,
    Point_3& volume_centroid,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::deque<PFace>& queue) const {

    const bool is_boundary = is_boundary_pface(
      pface, volume_index, num_volumes, map_volumes);
    if (is_boundary) {
      propagate_boundary_pface(
        verbose, pface, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    } else {
      propagate_interior_pface(
        verbose, pface, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    }
  }

  bool is_boundary_pface(
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<PFace, std::pair<int, int> >& map_volumes) const {

    CGAL_assertion(volume_index >= 0);
    if (pface.first < 6) return true;
    CGAL_assertion(pface.first >= 6);
    if (num_volumes == 0) return false;
    CGAL_assertion(num_volumes  > 0);
    CGAL_assertion(volume_index > 0);
    CGAL_assertion(volume_index >= num_volumes);

    const auto& pair = map_volumes.at(pface);
    if (pair.first == -1) {
      CGAL_assertion(pair.second == -1);
      return false;
    }
    CGAL_assertion(pair.first != -1);
    if (pair.first < num_volumes) return true;
    CGAL_assertion(pair.first >= num_volumes);
    return false;
  }

  void propagate_boundary_pface(
    const bool verbose,
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<int, Point_3>& centroids,
    std::size_t& volume_size,
    Point_3& volume_centroid,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::deque<PFace>& queue) const {

    auto& pair = map_volumes.at(pface);
    if (pair.first >= num_volumes) return;
    CGAL_assertion(pair.first < num_volumes);
    if (pair.second != -1) {
      CGAL_assertion(pair.first != -1);
      return;
    }
    CGAL_assertion(pair.second == -1);

    if (pair.first == -1) {
      pair.first  = volume_index;
    } else {
      pair.second = volume_index;
    }

    Point_3 centroid = m_data.centroid_of_pface(pface);
    if (num_volumes > 0) {
      // std::cout << "SHIFTING CENTROID" << std::endl;

      CGAL_assertion(pair.first < num_volumes);
      CGAL_assertion(centroids.find(pair.first) != centroids.end());
      const auto& other_centroid = centroids.at(pair.first);
      const auto plane = m_data.plane_of_pface(pface);
      auto vec1 = plane.orthogonal_vector();
      vec1 = KSR::normalize(vec1);
      auto vec2 = Vector_3(centroid, other_centroid);
      vec2 = KSR::normalize(vec2);

      // TODO: CAN WE AVOID THIS VALUE?
      const FT tol = KSR::tolerance<FT>();
      const FT dot_product = vec1 * vec2;

      if (dot_product < FT(0)) {
        centroid += tol * vec1;
      } else {
        centroid -= tol * vec1;
      }
      volume_centroid = CGAL::barycenter(
        volume_centroid, static_cast<FT>(volume_size), centroid, FT(1));

    } else {
      volume_centroid = CGAL::barycenter(
        volume_centroid, static_cast<FT>(volume_size), centroid, FT(1));
    }

    // std::cout << "volume centroid: " << volume_centroid << std::endl;
    ++volume_size;

    if (verbose) {
      // std::cout << "BND PFACE MAP: (" <<
      // pair.first << ", " << pair.second << ")" << std::endl;
      std::cout << "DUMPING BND PFACE: " <<
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second) << std::endl;

      dump_pface(m_data, pface, "bnd-pface-" +
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second));
    }

    std::vector<PFace> nfaces, bnd_nfaces, int_nfaces, all_nfaces;
    const auto pedges = m_data.pedges_of_pface(pface);
    for (const auto pedge : pedges) {
      CGAL_assertion(m_data.has_iedge(pedge));
      m_data.incident_faces(m_data.iedge(pedge), nfaces);

      split_pfaces(
        pface, volume_index, num_volumes, map_volumes, nfaces,
        bnd_nfaces, int_nfaces, all_nfaces);

      if (num_volumes == 0) {
        CGAL_assertion(bnd_nfaces.size() == 1);
        CGAL_assertion(int_nfaces.size() == 0 || int_nfaces.size() == 1);
      }

      if (int_nfaces.size() == 1) {
        queue.push_back(int_nfaces[0]);
        continue;
      }

      if (int_nfaces.size() == 0 && bnd_nfaces.size() == 1) {
        queue.push_front(bnd_nfaces[0]);
        continue;
      }

      if (all_nfaces.size() == 0) {
        dump_info(m_data, pface, pedge, nfaces);
        std::cout << "DEBUG: num nfaces: " << nfaces.size() << std::endl;
      }
      CGAL_assertion(all_nfaces.size() > 0);

      const auto found_nface = find_using_2d_directions(
      volume_index, volume_centroid, pface, pedge, all_nfaces);
      if (found_nface == m_data.null_pface()) continue;

      if (is_boundary_pface(
        found_nface, volume_index, num_volumes, map_volumes)) {
        queue.push_front(found_nface);
      } else {
        queue.push_back(found_nface);
      }
    }
  }

  void propagate_interior_pface(
    const bool verbose,
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<int, Point_3>& /* centroids */,
    std::size_t& volume_size,
    Point_3& volume_centroid,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::deque<PFace>& queue) const {

    CGAL_assertion(num_volumes >= 0);
    auto& pair = map_volumes.at(pface);
    if (pair.first != -1 && pair.second != -1) return;
    CGAL_assertion(pair.second == -1);
    if (pair.first == volume_index) return;
    CGAL_assertion(pair.first != volume_index);
    if (pair.first != -1) {
      pair.second = volume_index;
    } else {
      pair.first  = volume_index;
    }

    if (verbose) {
      std::cout << "pface: " << m_data.str(pface) << std::endl;
      std::cout << "pair: " <<
      std::to_string(pair.first) << "/" << std::to_string(pair.second) << std::endl;
    }

    const Point_3 centroid = m_data.centroid_of_pface(pface);
    volume_centroid = CGAL::barycenter(
      volume_centroid, static_cast<FT>(volume_size), centroid, FT(1));
    // std::cout << "volume centroid: " << volume_centroid << std::endl;
    ++volume_size;

    if (verbose) {
      // std::cout << "INT PFACE MAP: (" <<
      // pair.first << ", " << pair.second << ")" << std::endl;
      std::cout << "DUMPING INT PFACE: " <<
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second) << std::endl;
      dump_pface(m_data, pface, "int-pface-" +
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second));
    }

    std::vector<PFace> nfaces, bnd_nfaces, int_nfaces, all_nfaces;
    const auto pedges = m_data.pedges_of_pface(pface);
    for (const auto pedge : pedges) {
      CGAL_assertion(m_data.has_iedge(pedge));
      m_data.incident_faces(m_data.iedge(pedge), nfaces);
      split_pfaces(
        pface, volume_index, num_volumes, map_volumes, nfaces,
        bnd_nfaces, int_nfaces, all_nfaces);

      if (all_nfaces.size() == 0) {
        dump_info(m_data, pface, pedge, nfaces);
        std::cout << "DEBUG: num nfaces: " << nfaces.size() << std::endl;
      }
      CGAL_assertion(all_nfaces.size() > 0);

      const auto found_nface = find_using_2d_directions(
      volume_index, volume_centroid, pface, pedge, all_nfaces);
      if (found_nface == m_data.null_pface()) continue;

      if (is_boundary_pface(
        found_nface, volume_index, num_volumes, map_volumes)) {
        queue.push_front(found_nface);
      } else {
        queue.push_back(found_nface);
      }
    }
  }

  void split_pfaces(
    const PFace& current,
    const int volume_index,
    const int num_volumes,
    const std::map<PFace, std::pair<int, int> >& map_volumes,
    const std::vector<PFace>& pfaces,
    std::vector<PFace>& bnd_pfaces,
    std::vector<PFace>& int_pfaces,
    std::vector<PFace>& all_pfaces) const {

    bnd_pfaces.clear();
    int_pfaces.clear();
    all_pfaces.clear();
    for (const auto& pface : pfaces) {
      if (pface == current) continue;
      CGAL_assertion(pface != current);
      all_pfaces.push_back(pface);

      const auto& pair = map_volumes.at(pface);
      if (num_volumes > 0 && pair.first != -1) {
        if (pair.first < num_volumes && pair.second != -1) {
          if (pair.second < num_volumes) {
            continue;
          }
          CGAL_assertion(pair.second >= num_volumes);
        }
      }
      if (is_boundary_pface(
        pface, volume_index, num_volumes, map_volumes))  {
        bnd_pfaces.push_back(pface);
      } else {
        int_pfaces.push_back(pface);
      }
    }
  }

  const PFace find_using_2d_directions(
    const int /* volume_index */,
    const Point_3& volume_centroid,
    const PFace& pface,
    const PEdge& pedge,
    const std::vector<PFace>& nfaces) const {

    CGAL_assertion(nfaces.size() > 0);
    if (nfaces.size() == 1) return nfaces[0];
    const bool debug = false;
      // ( volume_index == 31 &&
      //   pface.first == 8 &&
      //   static_cast<std::size_t>(pface.second) == 7);

    if (debug) {
      dump_info(m_data, pface, pedge, nfaces);
    }
    CGAL_assertion(nfaces.size() > 1);

    Point_3 center = m_data.centroid_of_pface(pface);
    const Segment_3 segment = m_data.segment_3(pedge);
    const Line_3 line(segment.source(), segment.target());
    Point_3 midp = CGAL::midpoint(segment.source(), segment.target());
    // std::cout << "midp: " << midp << std::endl;
    Vector_3 norm(segment.source(), segment.target());
    norm = KSR::normalize(norm);
    const Plane_3 plane(midp, norm);

    std::vector<Point_3> points;
    points.reserve(nfaces.size() + 2);

    points.push_back(midp);
    points.push_back(center);
    for (const auto& nface : nfaces) {
      center = m_data.centroid_of_pface(nface);
      points.push_back(center);
    }
    CGAL_assertion(points.size() >= 3);

    for (auto& point : points) {
      point = plane.projection(point);
    }

    if (debug) {
      dump_frame(points, "volumes/directions-init");
    }

    const FT cx = volume_centroid.x();
    const FT cy = volume_centroid.y();
    const FT cz = volume_centroid.z();
    FT d = (norm.x() * cx + norm.y() * cy + norm.z() * cz + plane.d());

    // std::cout << "1 d: " << d << std::endl;
    // std::cout << "1 norm: " << norm << std::endl;
    const Plane_3 tr_plane(midp + norm * d, norm);
    Point_3 inter;
    const bool is_intersection_found =
      m_kinetic_traits.intersection(line, tr_plane, inter);
    if (!is_intersection_found) {
      std::cout << "d = " << d << std::endl;
    }
    CGAL_assertion(is_intersection_found);
    // std::cout << "inter: " << inter << std::endl;

    d = KSR::distance(midp, inter);
    norm = Vector_3(midp, inter);
    // std::cout << "2 d: " << d << std::endl;
    // std::cout << "2 norm: " << norm << std::endl;

    if (d != FT(0)) {
      CGAL_assertion(norm != Vector_3(FT(0), FT(0), FT(0)));
      norm = KSR::normalize(norm);
      for (auto& point : points) {
        point += norm * d;
      }
    }

    if (debug) {
      auto extended = points;
      extended.push_back(volume_centroid);
      dump_frame(extended, "volumes/directions");
    }

    std::vector< std::pair<Direction_2, PFace> > dir_edges;
    dir_edges.reserve(nfaces.size() + 1);

    const Point_2 proj_0 = plane.to_2d(points[0]);
    for (std::size_t i = 1; i < points.size(); ++i) {
      const Point_2 proj_i = plane.to_2d(points[i]);
      const Vector_2 vec(proj_0, proj_i);
      if (i == 1) {
        dir_edges.push_back(std::make_pair(Direction_2(vec), pface));
      } else {
        dir_edges.push_back(std::make_pair(Direction_2(vec), nfaces[i - 2]));
      }
    }
    CGAL_assertion(dir_edges.size() == nfaces.size() + 1);

    const Point_2 proj_vc = plane.to_2d(volume_centroid);
    const Vector_2 vec(proj_0, proj_vc);
    const Direction_2 ref_dir(vec);

    std::sort(dir_edges.begin(), dir_edges.end(), [&](
      const std::pair<Direction_2, PFace>& p,
      const std::pair<Direction_2, PFace>& q) -> bool {
        return p.first < q.first;
      }
    );

    const std::size_t n = dir_edges.size();
    for (std::size_t i = 0; i < n; ++i) {
      if (dir_edges[i].second == pface) {

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        const auto& dir_prev = dir_edges[im].first;
        const auto& dir_curr = dir_edges[i].first;
        const auto& dir_next = dir_edges[ip].first;

        if (debug) {
          dump_pface(m_data, dir_edges[im].second, "prev");
          dump_pface(m_data, dir_edges[ip].second, "next");
        }

        if (ref_dir.counterclockwise_in_between(dir_prev, dir_curr)) {
          if (debug) {
            std::cout << "found prev" << std::endl;
            exit(EXIT_SUCCESS);
          }
          return dir_edges[im].second;
        } else if (ref_dir.counterclockwise_in_between(dir_curr, dir_next)) {
          if (debug) {
            std::cout << "found next" << std::endl;
            exit(EXIT_SUCCESS);
          }
          return dir_edges[ip].second;
        } else {
          // return m_data.null_pface();
          std::cout << "ERROR: WRONG ORIENTATIONS!" << std::endl;
          dump_info(m_data, pface, pedge, nfaces);
          dump_frame(points, "volumes/directions-init");
          auto extended = points;
          extended.push_back(volume_centroid);
          dump_frame(extended, "volumes/directions");
          CGAL_assertion_msg(false, "ERROR: WRONG ORIENTATIONS!");
        }
      }
    }

    CGAL_assertion_msg(false, "ERROR: NEXT PFACE IS NOT FOUND!");
    return m_data.null_pface();
  }

  void create_cell_pvertices(Volume_cell& cell) {
    cell.pvertices.clear();
    for (const auto& pface : cell.pfaces) {
      for (const auto pvertex : m_data.pvertices_of_pface(pface)) {
        cell.pvertices.insert(pvertex);
      }
    }
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_FINALIZER_H
