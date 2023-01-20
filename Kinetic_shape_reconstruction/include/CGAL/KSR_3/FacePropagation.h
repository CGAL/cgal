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
// Author(s)     : Sven Oesau, Simon Giraudot, Dmitry Anisimov

#ifndef CGAL_KSR_3_FACEPROPAGATION_H
#define CGAL_KSR_3_FACEPROPAGATION_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR/parameters.h>
#include <CGAL/KSR/debug.h>

#include <CGAL/KSR_3/Data_structure.h>

namespace CGAL {
namespace KSR_3 {

#ifdef DOXYGEN_RUNNING
#else

template<typename GeomTraits, typename Intersection_Kernel>
class FacePropagation {

public:
  using Kernel = GeomTraits;

private:
  using FT          = typename Kernel::FT;
  using Point_2     = typename Kernel::Point_2;
  using Vector_2    = typename Kernel::Vector_2;
  using Segment_2   = typename Kernel::Segment_2;
  using Direction_2 = typename Kernel::Direction_2;
  using Line_2      = typename Kernel::Line_2;

  using Data_structure = KSR_3::Data_structure<Kernel, Intersection_Kernel>;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;
  using IFace   = typename Data_structure::IFace;

  using PVertex = typename Data_structure::PVertex;
  using PEdge   = typename Data_structure::PEdge;
  using PFace   = typename Data_structure::PFace;

  using Bbox_2     = CGAL::Bbox_2;
  using Face_index = typename Data_structure::Face_index;

  using Parameters     = KSR::Parameters_3<FT>;
  using Kinetic_traits = KSR::Kinetic_traits_3<Kernel>;

  using FaceEvent      = typename Data_structure::Support_plane::FaceEvent;

  struct FaceEventOrder {
    bool operator()(const FaceEvent &a, const FaceEvent &b) {
      return a.time > b.time;
    }
  };

public:
  FacePropagation(Data_structure& data, const Parameters& parameters) :
  m_data(data), m_parameters(parameters), m_kinetic_traits(),
  m_min_time(-FT(1)), m_max_time(-FT(1))
  { }

  const std::pair<std::size_t, std::size_t> propagate() {
    std::size_t num_queue_calls = 0;
    std::size_t num_events = 0;

    initialize_queue();

    while (!m_face_queue.empty()) {
      num_events = run(num_events);

      ++num_queue_calls;
    }

    return std::make_pair(num_queue_calls, num_events);
  }

  void clear() {
    m_queue.clear();
    m_min_time = -FT(1);
    m_max_time = -FT(1);
  }

private:
  Data_structure& m_data;
  const Parameters& m_parameters;
  Kinetic_traits m_kinetic_traits;

  FT m_min_time;
  FT m_max_time;

  std::priority_queue<typename Data_structure::Support_plane::FaceEvent, std::vector<FaceEvent>, FaceEventOrder> m_face_queue;

  /*******************************
  **       IDENTIFY EVENTS      **
  ********************************/

  void initialize_queue() {
    //m_face_queue.clear();

    if (m_parameters.debug) {
      std::cout << "initializing queue" << std::endl;
    }

    m_data.fill_event_queue(m_face_queue);
  }

  /*******************************
  **          RUNNING           **
  ********************************/

  std::size_t run(
    const std::size_t initial_iteration) {

    if (m_parameters.debug) {
      std::cout << "* unstacking queue, current size: " << m_face_queue.size() << std::endl;
    }

    std::size_t iteration = initial_iteration;
    while (!m_face_queue.empty()) {
      // m_queue.print();

      const FaceEvent event = m_face_queue.top();
      m_face_queue.pop();
      const FT current_time = event.time;

      ++iteration;

      apply(event);
    }
    return iteration;
  }

  /*******************************
  **        HANDLE EVENTS       **
  ********************************/

  void apply(const FaceEvent& event) {
    //std::cout << "support plane: " << event.support_plane << " edge: " << event.crossed_edge << " t: " << event.time << std::endl;
    if (m_data.igraph().face(event.face).part_of_partition) {
      //std::cout << " face already crossed, skipping event" << std::endl;
      return;
    }

    std::size_t line = m_data.line_idx(event.crossed_edge);
    if (!m_data.support_plane(event.support_plane).has_crossed_line(line)) {
      // Check intersection against kinetic intervals from other support planes
      int crossing = 0;
      auto kis = m_data.igraph().kinetic_intervals(event.crossed_edge);
      for (auto ki = kis.first; ki != kis.second; ki++) {
        if (ki->first == event.support_plane)
          continue;

        for (std::size_t i = 0; i < ki->second.size(); i++) {
          // Exactly on one
          if (ki->second[i].first == event.intersection_bary) {
            if (ki->second[i].second < event.time)
              crossing++;

            break;
          }

          // Within an interval
          if (ki->second[i].first > event.intersection_bary && ki->second[i - 1].first < event.intersection_bary) {
            FT interval_pos = (event.intersection_bary - ki->second[i - 1].first) / (ki->second[i].first - ki->second[i - 1].first);
            FT interval_time = interval_pos * (ki->second[i].second - ki->second[i - 1].second) + ki->second[i - 1].second;

            if (event.time > interval_time)
              crossing++;

            break;
          }
        }
      }

      // Check if the k value is sufficient for crossing the edge.
      int& k = m_data.support_plane(event.support_plane).k();
      if (k <= crossing)
        return;

      // The edge can be crossed.
      // Adjust k value
      k -= crossing;

      m_data.support_plane(event.support_plane).set_crossed_line(line);
    }

    // Associate IFace to mesh.
    PFace f = m_data.add_iface_to_mesh(event.support_plane, event.face);

    // Calculate events for new border edges.
    // Iterate inside of this face, check if each opposite edge is border and on bbox and then calculate intersection times.
    std::vector<IEdge> border;
    m_data.support_plane(event.support_plane).get_border(m_data.igraph(), f.second, border);

    for (IEdge edge : border) {
      FaceEvent fe;
      FT t = m_data.calculate_edge_intersection_time(event.support_plane, edge, fe);
      if (t > 0)
        m_face_queue.push(fe);
    }
  }
};

#endif //DOXYGEN_RUNNING

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_FACEPROPAGATION_H
