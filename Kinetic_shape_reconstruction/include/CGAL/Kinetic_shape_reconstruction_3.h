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

#ifndef CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
#define CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR_3/Event.h>
#include <CGAL/KSR_3/Event_queue.h>
#include <CGAL/KSR_3/Data_structure.h>
#include <CGAL/KSR_3/Initializer.h>

namespace CGAL {

template<typename GeomTraits>
class Kinetic_shape_reconstruction_3 {

public:
  using Kernel = GeomTraits;

private:
  using FT        = typename Kernel::FT;
  using Point_2   = typename Kernel::Point_2;
  using Vector_2  = typename Kernel::Vector_2;
  using Segment_2 = typename Kernel::Segment_2;

  using Data_structure = KSR_3::Data_structure<Kernel>;

  using PVertex = typename Data_structure::PVertex;
  using PFace   = typename Data_structure::PFace;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  using Event       = KSR_3::Event<Data_structure>;
  using Event_queue = KSR_3::Event_queue<Data_structure>;

  using EK = CGAL::Exact_predicates_exact_constructions_kernel;
  using Initializer = KSR_3::Initializer<EK>;

private:
  Data_structure m_data;
  Event_queue m_queue;
  FT m_min_time;
  FT m_max_time;
  const bool m_verbose;
  Initializer m_initializer;

public:
  Kinetic_shape_reconstruction_3(const bool verbose = true) :
  m_min_time(-FT(1)),
  m_max_time(-FT(1)),
  m_verbose(verbose),
  m_initializer(m_verbose)
  { }

  template<
  typename InputRange,
  typename PolygonMap>
  const bool partition(
    const InputRange& input_range,
    const PolygonMap polygon_map,
    const unsigned int k = 1,
    const double enlarge_bbox_ratio = 1.1,
    const bool reorient = false) {

    if (m_verbose) std::cout.precision(20);
    if (input_range.size() == 0) {
      CGAL_warning_msg(input_range.size() != 0,
      "WARNING: YOUR INPUT IS EMPTY. RETURN WITH NO CHANGE!");
      return false;
    }

    if (k == 0) {
      CGAL_warning_msg(k != 0,
      "WARNING: YOU SET K TO 0. THE VALID VALUES ARE {1,2,...}. RETURN WITH NO CHANGE!");
      return false;
    }

    if (enlarge_bbox_ratio < 1.0) {
      CGAL_warning_msg(enlarge_bbox_ratio >= 1.0,
      "WARNING: YOU SET ENLARGE_BBOX_RATIO < 1.0. THE VALID RANGE IS [1.0, +INF). RETURN WITH NO CHANGE!");
      return false;
    }

    const FT time_step = static_cast<FT>(m_initializer.initialize(
      input_range, polygon_map, k, enlarge_bbox_ratio, reorient));
    m_initializer.convert(m_data);

    // if (m_verbose) {
    //   std::cout << std::endl << "POLYGON SPLITTER SUCCESS!" << std::endl << std::endl;
    //   exit(EXIT_SUCCESS);
    // }

    if (m_verbose) {
      std::cout << std::endl << "--- RUNNING THE QUEUE:" << std::endl;
      std::cout << "propagation started ..." << std::endl;
    }
    std::size_t num_iterations = 0;
    m_min_time = FT(0);
    m_max_time = time_step;
    CGAL_assertion(m_min_time >= FT(0) && m_max_time >= m_min_time);
    while (initialize_queue()) {

      run(k);
      m_min_time = m_max_time;
      m_max_time += time_step;
      m_data.check_integrity();
      ++num_iterations;

      // if (m_verbose) {
      //   std::cout << ".";
      //   if (num_iterations == 50) {
      //     std::cout << std::endl;
      //   }
      // }

      // if (num_iterations > 100) {
      //   CGAL_assertion_msg(false, "WHY SO MANY ITERATIONS?");
      // }
    }
    if (m_verbose) {
      std::cout << "... propagation finished" << std::endl;
    }

    if (m_verbose) {
      std::cout << std::endl << "--- FINALIZING KSR:" << std::endl;
      std::cout << "* checking final mesh integrity ...";
    }
    m_data.check_integrity();
    if (m_verbose) {
      dump(m_data, "iter_1000-final-result");
      std::cout << " done" << std::endl;
    }

    if (m_verbose) {
      std::cout << "* getting volumes:" << std::endl;
    }
    m_data.create_polyhedrons();
    return true;
  }

  template<typename OutputIterator>
  OutputIterator output_partition_edges_to_segment_soup(
    OutputIterator edges) const {

    CGAL_assertion_msg(false, "TODO: IMPLEMENT OUTPUT PARTITION EDGES!");
    return edges;
  }

  template<typename VertexOutputIterator, typename FaceOutputIterator>
  void output_partition_faces_to_polygon_soup(
    VertexOutputIterator vertices, FaceOutputIterator faces,
    const bool with_bbox = false) const {

    CGAL_assertion_msg(false, "TODO: IMPLEMENT OUTPUT PARTITION FACES!");
  }

  template<typename PolyhedronOutputIterator>
  PolyhedronOutputIterator output_partition_polyhedrons(
    PolyhedronOutputIterator polyhedrons) const {

    CGAL_assertion_msg(false, "TODO: IMPLEMENT OUTPUT PARTITION POLYHEDRONS!");
    return polyhedrons;
  }

  template<typename InputRange, typename PointMap, typename VectorMap>
  void reconstruct(
    const InputRange& input_range,
    const PointMap point_map,
    const VectorMap normal_map) {

    CGAL_assertion_msg(false, "TODO: ADD RECONSTRUCTION!");
  }

private:
  const bool initialize_queue() {
    std::cout << "Initializing queue for events in [" <<
    m_min_time << ";" << m_max_time << "]" << std::endl;

    m_data.update_positions(m_max_time);
    bool still_running = false;

    KSR::vector<IEdge> iedges;
    KSR::vector<Segment_2> segments_2;
    KSR::vector<CGAL::Bbox_2> segment_bboxes;
    for (KSR::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      init_search_structures(i, iedges, segments_2, segment_bboxes);
      for (const PVertex pvertex : m_data.pvertices(i)) {
        if (compute_events_of_vertex(pvertex, iedges, segments_2, segment_bboxes)) {
          still_running = true;
        }
      }
    }
    m_data.update_positions(m_min_time);
    return still_running;
  }

  void init_search_structures(
    const KSR::size_t i,
    KSR::vector<IEdge>& iedges,
    KSR::vector<Segment_2>& segments_2,
    KSR::vector<CGAL::Bbox_2>& segment_bboxes) {

    iedges.clear();
    segments_2.clear();
    segment_bboxes.clear();

    // To get random access, copy in vector (suboptimal to do this
    // all the time, maybe this should be done once and for all and
    // replace the set).
    iedges.reserve(m_data.iedges(i).size());
    std::copy(m_data.iedges(i).begin(), m_data.iedges(i).end(), std::back_inserter(iedges));

    // Precompute segments and bboxes.
    segments_2.reserve(iedges.size());
    segment_bboxes.reserve(iedges.size());
    for (const IEdge& iedge : iedges) {
      segments_2.push_back(m_data.segment_2(i, iedge));
      segment_bboxes.push_back(segments_2.back().bbox());
    }
  }

  bool compute_events_of_vertex (const PVertex& pvertex,
                                 const KSR::vector<IEdge>& iedges,
                                 const KSR::vector<Segment_2>& segments_2,
                                 const KSR::vector<CGAL::Bbox_2>& segment_bboxes)
  {
    std::cout.precision(20);
    if (m_data.is_frozen(pvertex))
      return false;

    Segment_2 sv (m_data.point_2 (pvertex, m_min_time),
                  m_data.point_2 (pvertex, m_max_time));
    CGAL::Bbox_2 sv_bbox = sv.bbox();

    if (m_data.has_iedge(pvertex)) // constrained vertex
    {
      // const auto cutime = m_data.current_time();
      // m_data.update_positions(m_min_time);
      // std::cout << "Computing events for the constrained pvertex: " << m_data.str(pvertex) << ": " << m_data.point_3(pvertex) << std::endl;
      // m_data.update_positions(cutime);

      // Test left and right vertices on mesh face.
      PVertex prev;
      PVertex next;
      std::tie (prev, next) = m_data.prev_and_next (pvertex);

      for (const PVertex& pother : { prev, next })
      {
        if (pother == Data_structure::null_pvertex()
            || !m_data.is_active(pother)
            || m_data.has_iedge (pother))
          continue;

        Segment_2 so (m_data.point_2 (pother, m_min_time),
                      m_data.point_2 (pother, m_max_time));
        CGAL::Bbox_2 so_bbox = so.bbox();

        if (!do_overlap (sv_bbox, so_bbox))
          continue;

        Point_2 point;
        if (!KSR::intersection(sv, so, point))
          continue;

        FT dist = CGAL::approximate_sqrt (CGAL::squared_distance (sv.source(), point));
        FT time = dist / m_data.speed(pvertex);

        m_queue.push (Event (true, pvertex, pother, m_min_time + time));

        // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
        // std::cout << "pother: " << m_data.point_3(pother) << std::endl;
      }

      // Test end-vertices of intersection edge.
      IEdge iedge = m_data.iedge(pvertex);
      for (const IVertex& ivertex : { m_data.source(iedge), m_data.target(iedge) })
      {
        if (!m_data.is_active(ivertex))
          continue;
        Point_2 pi = m_data.to_2d (pvertex.first, ivertex);
        // std::cout << m_data.str(pvertex) << std::endl;
        // std::cout << m_data.point_3(ivertex) << std::endl;
        // std::cout << "2 " << m_data.to_3d(pvertex.first, sv.source()) <<
        // " " << m_data.to_3d(pvertex.first, sv.target()) << std::endl;
        // std::cout << "2 " << m_data.to_3d(pvertex.first, sv.source()) <<
        // " " << m_data.to_3d(pvertex.first, pi) << std::endl;
        // std::cout << sv.to_vector() << std::endl;
        // std::cout << Vector_2 (sv.source(), pi) << std::endl;
        // std::cout << sv.to_vector() * Vector_2 (sv.source(), pi) << std::endl;
        if (sv.to_vector() * Vector_2 (sv.source(), pi) < 0)
          continue;

        FT dist = CGAL::approximate_sqrt(CGAL::squared_distance (sv.source(), pi));
        FT time = dist / m_data.speed(pvertex);

        if (time < m_max_time - m_min_time) {
          m_queue.push (Event (true, pvertex, ivertex, m_min_time + time));

          // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
          // std::cout << "ivertex: " << m_data.point_3(ivertex) << std::endl;
        }
      }
    }
    else // unconstrained vertex
    {
      PVertex prev = m_data.prev(pvertex);
      PVertex next = m_data.next(pvertex);

      // Test all intersection edges.
      for (std::size_t j = 0; j < iedges.size(); ++ j)
      {
        const IEdge& iedge = iedges[j];

        if (m_data.iedge(prev) == iedge ||
            m_data.iedge(next) == iedge)
          continue;
        if (!m_data.is_active(iedge))
          continue;

        if (!CGAL::do_overlap (sv_bbox, segment_bboxes[j]))
          continue;

        Point_2 point;
        if (!KSR::intersection (sv, segments_2[j], point))
          continue;

        FT dist = CGAL::approximate_sqrt (CGAL::squared_distance (m_data.point_2 (pvertex, m_min_time), point));
        FT time = dist / m_data.speed(pvertex);

        m_queue.push (Event (false, pvertex, iedge, m_min_time + time));

        // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
        // std::cout << "iedge: " << m_data.segment_3(iedge) << std::endl;
      }
    }
    return true;
  }

  const bool are_parallel(
    const Segment_2& seg1, const Segment_2& seg2) {

    const FT tol = FT(1) / FT(100000);
    FT m1 = FT(100000), m2 = FT(100000);

    const FT d1 = (seg1.target().x() - seg1.source().x());
    const FT d2 = (seg2.target().x() - seg2.source().x());

    if (CGAL::abs(d1) > tol)
      m1 = (seg1.target().y() - seg1.source().y()) / d1;
    if (CGAL::abs(d2) > tol)
      m2 = (seg2.target().y() - seg2.source().y()) / d2;

    // return CGAL::parallel(seg1, seg2); // exact version

    if (CGAL::abs(m1 - m2) < tol)
      return true;
    return false;
  }

  void run(const unsigned int k)
  {
    std::cout << "Unstacking queue size: " << m_queue.size() << std::endl;

    KSR::size_t iterations = 0;

    static int iter = 0;

    while (!m_queue.empty())
    {
      Event ev = m_queue.pop();

      FT current_time = ev.time();

      if (iter < 10)
      {
        dump (m_data, "iter_0" + std::to_string(iter));
        dump_event (m_data, ev, "iter_0" + std::to_string(iter));
      }
      else
      {
        dump (m_data, "iter_" + std::to_string(iter));
        dump_event (m_data, ev, "iter_" + std::to_string(iter));
      }

      m_data.update_positions (current_time);

      std::cout << "* APPLYING " << iter << ": " << ev << std::endl << std::endl;

      ++ iter;

      // if (iter == 380) {
      //   exit(EXIT_FAILURE);
      // }

      apply(k, ev);
      m_data.check_integrity();

      // m_data.update_positions (0.5 * (current_time + m_queue.next().time()));
      // dump (m_data, "after_" + std::to_string(iter - 1));
      // m_data.update_positions (current_time);
      ++ iterations;
    }
  }

  void apply (
    const unsigned int k,
    const Event& ev)
  {
    PVertex pvertex = ev.pvertex();

    if (ev.is_pvertex_to_pvertex())
    {
      PVertex pother = ev.pother();

      remove_events (pvertex);
      remove_events (pother);

      CGAL_assertion (m_data.has_iedge(pvertex));

      if (m_data.has_iedge(pother)) // two constrained vertices meet
      {
        CGAL_assertion_msg(false, "TODO: ADD CASE TWO CONSTRAINED PVERTICES MEET!");
      }
      else // one constrained vertex meets a free vertex
      {
        if (m_data.transfer_vertex(pvertex, pother)) {

          if (m_data.has_iedge(pvertex))
            remove_events(m_data.iedge(pvertex), pvertex.first); // should we remove it here?
          if (m_data.has_iedge(pother))
            remove_events(m_data.iedge(pother), pother.first); // should we remove it here?
          compute_events_of_vertices(ev.time(), std::array<PVertex,2>{pvertex, pother});

          PVertex prev, next;
          std::tie(prev, next) = m_data.border_prev_and_next(pvertex);

          PVertex pthird = prev;
          if (pthird == pother)
            pthird = next;
          else
            CGAL_assertion(next == pother);

          // remove_events(pthird);
          if (m_data.has_iedge(pthird))
            remove_events(m_data.iedge(pthird), pthird.first); // should we remove it here?
          compute_events_of_vertices(ev.time(), std::array<PVertex,1>{pthird});

        } else {

          if (m_data.has_iedge(pvertex))
            remove_events(m_data.iedge(pvertex), pvertex.first); // should we remove it here?
          compute_events_of_vertices(ev.time(), std::array<PVertex,1>{pvertex});
        }
      }
    }
    else if (ev.is_pvertex_to_iedge())
    {
      PVertex prev = m_data.prev(pvertex);
      PVertex next = m_data.next(pvertex);
      IEdge iedge = ev.iedge();
      PFace pface = m_data.pface_of_pvertex (pvertex);

      Segment_2 seg_edge = m_data.segment_2 (pvertex.first, iedge);

      bool done = false;
      for (const PVertex& pother : { prev, next })
      {
        Segment_2 seg (m_data.point_2(pother, ev.time()),
                       m_data.point_2(pvertex, ev.time()));
        CGAL_assertion (seg.squared_length() != 0);

        bool both_are_free = true;
        if (
          m_data.iedge(pvertex) != m_data.null_iedge() ||
          m_data.iedge(pother)  != m_data.null_iedge()) {
          both_are_free = false;
        }

        if (both_are_free && are_parallel(seg, seg_edge)) {

          remove_events(pvertex);
          remove_events(pother);

          bool collision, bbox_reached;
          std::tie(collision, bbox_reached) = m_data.collision_occured(pvertex, iedge);
          // std::tie(collision, bbox_reached) = m_data.is_occupied(pvertex, iedge);
          std::cout << "collision/bbox: " << collision << "/" << bbox_reached << std::endl;

          bool collision_other, bbox_reached_other;
          std::tie(collision_other, bbox_reached_other) = m_data.collision_occured(pother, iedge);
          // std::tie(collision_other, bbox_reached_other) = m_data.is_occupied(pother, iedge);
          std::cout << "other/bbox: " << collision_other << "/" << bbox_reached_other << std::endl;

          std::cout << "k intersections: " << m_data.k(pface) << std::endl;
          bool stop = false;
          if (bbox_reached) {

            CGAL_assertion(bbox_reached_other); // can we have a case with only one box side reached?
            std::cout << "pv po k bbox" << std::endl;
            stop = true;

          } else if (bbox_reached_other) {

            CGAL_assertion(bbox_reached); // can we have a case with only one box side reached?
            std::cout << "pv po k bbox" << std::endl;
            stop = true;

          } else if ((collision || collision_other) && m_data.k(pface) == 1) {

            std::cout << "pv po k stop" << std::endl;
            stop = true;

          } else if ((collision || collision_other) && m_data.k(pface) > 1) {

            std::cout << "pv po k continue" << std::endl;
            m_data.k(pface)--;

          } else {

            std::cout << "pv po continue" << std::endl;
            CGAL_assertion(m_data.iedge(pvertex) == m_data.iedge(pother));
            if (m_data.is_occupied(pvertex, iedge).first) {
              CGAL_assertion_msg(false, "TODO: TWO PVERTICES SNEAK ON THE OTHER SIDE EVEN WHEN WE HAVE A POLYGON!");
            }
          }
          CGAL_assertion(m_data.k(pface) >= 1);

          if (stop) // polygon stops
          {
            m_data.crop_polygon(pvertex, pother, iedge);
            remove_events(iedge, pvertex.first);
            compute_events_of_vertices(ev.time(), std::array<PVertex,2>{pvertex, pother});
          }
          else // polygon continues beyond the edge
          {
            PVertex pv0, pv1;
            std::tie(pv0, pv1) = m_data.propagate_polygon(m_data.k(pface), pvertex, pother, iedge);
            remove_events(iedge, pvertex.first);
            compute_events_of_vertices(ev.time(), std::array<PVertex, 4>{pvertex, pother, pv0, pv1});
          }

          done = true;
          break;
        }
      }

      if (!done) {

        remove_events(pvertex);

        bool collision, bbox_reached;
        std::tie(collision, bbox_reached) = m_data.collision_occured(pvertex, iedge);
        // std::tie(collision, bbox_reached) = m_data.is_occupied(pvertex, iedge);
        std::cout << "collision/bbox: " << collision << "/" << bbox_reached << std::endl;

        std::cout << "k intersections: " << m_data.k(pface) << std::endl;
        bool stop = false;
        if (bbox_reached) {

          std::cout << "pv k bbox" << std::endl;
          stop = true;

        } else if (collision && m_data.k(pface) == 1) {

          std::cout << "pv k stop" << std::endl;
          stop = true;

        } else if (collision && m_data.k(pface) > 1) {

          std::cout << "pv k continue" << std::endl;
          m_data.k(pface)--;

        } else {
          std::cout << "pv continue" << std::endl;
        }
        CGAL_assertion(m_data.k(pface) >= 1);

        if (stop) // polygon stops
        {
          const PVertex pvnew = m_data.crop_polygon(pvertex, iedge);
          remove_events(iedge, pvertex.first);
          compute_events_of_vertices(ev.time(), std::array<PVertex,2>{pvertex, pvnew});
        }
        else // polygon continues beyond the edge
        {
          const std::array<PVertex, 3> pvnew = m_data.propagate_polygon(m_data.k(pface), pvertex, iedge);
          remove_events(iedge, pvertex.first);
          compute_events_of_vertices(ev.time(), pvnew);
        }
      }
    }
    else if (ev.is_pvertex_to_ivertex())
    {
      // First, let's gather all vertices that will get merged.
      std::vector<PVertex> pvertices
        = m_data.pvertices_around_ivertex (ev.pvertex(), ev.ivertex());

      for (auto& pv: pvertices)
        std::cerr << m_data.point_3(pv) << std::endl;
      std::cerr << std::endl;

      std::cerr << "Found " << pvertices.size() << " pvertices ready to be merged" << std::endl;

      // Remove associated events.
      // for (const PVertex& pvertex : pvertices)
      //   remove_events (pvertex);
      for (std::size_t i = 1; i < pvertices.size() - 1; ++i)
        remove_events (pvertices[i]);

      // Merge them and get the newly created vertices.
      // std::cout << "came from: " << m_data.segment_3(m_data.iedge(ev.pvertex())) << std::endl;
      std::vector<IEdge> crossed;
      std::vector<PVertex> new_pvertices
        = m_data.merge_pvertices_on_ivertex(m_min_time, m_max_time, ev.pvertex(), pvertices, ev.ivertex(), crossed);

      // Remove all events of the crossed iedges.
      for (const auto& iedge : crossed)
        remove_events(iedge, pvertex.first);

      // And compute new events.
      CGAL_assertion(new_pvertices.size() > 0);
      compute_events_of_vertices (ev.time(), new_pvertices);
    }
    else
    {
      CGAL_assertion_msg (false, "ERROR: INVALID EVENT!");
    }
  }

  void remove_events (
    const IEdge& iedge,
    const KSR::size_t support_plane_idx) {

    m_queue.erase_vertex_events(iedge, support_plane_idx);
    // std::cout << "erasing events for iedge " << m_data.str(iedge) << std::endl;
    // std::cout << m_data.segment_3(iedge) << std::endl;
  }

  void remove_events (const PVertex& pvertex) {
    m_queue.erase_vertex_events (pvertex);
    // std::cout << "erasing events for pvertex " << m_data.str(pvertex) << " : " << m_data.point_3(pvertex) << std::endl;
  }

  template <typename PVertexRange>
  void compute_events_of_vertices (
    const FT last_event_time, const PVertexRange& pvertices) {

    m_min_time = m_data.current_time();
    m_data.update_positions(m_max_time);

    KSR::vector<IEdge> iedges;
    KSR::vector<Segment_2> segments_2;
    KSR::vector<CGAL::Bbox_2> segment_bboxes;
    init_search_structures(pvertices.front().first, iedges, segments_2, segment_bboxes);

    for (const PVertex& pvertex : pvertices)
      m_data.deactivate(pvertex);

    for (const PVertex& pvertex : pvertices) {
      m_data.set_last_event_time(pvertex, last_event_time);
      compute_events_of_vertex(pvertex, iedges, segments_2, segment_bboxes);
    }

    for (const PVertex& pvertex : pvertices)
      m_data.activate(pvertex);

    m_data.update_positions(m_min_time);
  }

};

} // namespace CGAL

#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
