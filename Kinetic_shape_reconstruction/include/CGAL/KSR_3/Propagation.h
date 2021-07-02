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

#ifndef CGAL_KSR_3_PROPAGATION_H
#define CGAL_KSR_3_PROPAGATION_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR_3/Event.h>
#include <CGAL/KSR_3/Event_queue.h>
#include <CGAL/KSR_3/Data_structure.h>

namespace CGAL {
namespace KSR_3 {

template<typename GeomTraits>
class Propagation {

public:
  using Kernel = GeomTraits;

private:
  using FT          = typename Kernel::FT;
  using Point_2     = typename Kernel::Point_2;
  using Vector_2    = typename Kernel::Vector_2;
  using Segment_2   = typename Kernel::Segment_2;
  using Direction_2 = typename Kernel::Direction_2;
  using Line_2      = typename Kernel::Line_2;

  using Data_structure = KSR_3::Data_structure<Kernel>;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  using PVertex = typename Data_structure::PVertex;
  using PEdge   = typename Data_structure::PEdge;
  using PFace   = typename Data_structure::PFace;

  using Event       = KSR_3::Event<Data_structure>;
  using Event_queue = KSR_3::Event_queue<Data_structure>;

  using Bbox_2     = CGAL::Bbox_2;
  using Face_index = typename Data_structure::Face_index;

public:
  Propagation(
    const bool verbose, const bool dprint, const bool debug, Data_structure& data) :
  m_verbose(verbose), m_export(dprint), m_debug(debug), m_data(data),
  m_queue(m_debug), m_min_time(-FT(1)), m_max_time(-FT(1))
  { }

  const std::pair<std::size_t, std::size_t> propagate(const FT time_step) {

    std::size_t num_queue_calls = 0;
    m_min_time = FT(0);
    m_max_time = time_step;
    CGAL_assertion(m_min_time >= FT(0) && m_max_time >= m_min_time);
    std::size_t num_events = 0;
    while (initialize_queue()) {

      num_events = run(num_events);
      m_min_time = m_max_time;
      m_max_time += time_step;
      CGAL_assertion(m_data.check_integrity());
      ++num_queue_calls;

      if (m_verbose && !m_debug) {
        if ((num_queue_calls % 50) == 0) {
          std::cout << ".................................................." << std::endl;
        }
      }

      if (num_queue_calls > 1000000) {
        CGAL_assertion_msg(false, "DEBUG ERROR: WHY SO MANY ITERATIONS?");
        break;
      }
    }
    return std::make_pair(num_queue_calls, num_events);
  }

  void clear() {
    m_queue.clear();
    m_min_time = -FT(1);
    m_max_time = -FT(1);
  }

private:
  const bool m_verbose;
  const bool m_export;
  const bool m_debug;
  Data_structure& m_data;
  Event_queue m_queue;
  FT m_min_time;
  FT m_max_time;

  /*******************************
  **       IDENTIFY EVENTS      **
  ********************************/

  bool initialize_queue() {

    if (m_debug) {
      std::cout << "* initializing queue for events in [" <<
      m_min_time << ";" << m_max_time << "]" << std::endl;
    }

    m_data.update_positions(m_max_time);
    bool still_running = false;
    for (std::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      const auto& iedges   = m_data.iedges(i);
      const auto& segments = m_data.isegments(i);
      const auto& bboxes   = m_data.ibboxes(i);
      for (const auto pvertex : m_data.pvertices(i)) {
        if (compute_events_of_pvertex(pvertex, iedges, segments, bboxes)) {
          still_running = true;
        }
      }
    }
    m_data.update_positions(m_min_time);
    return still_running;
  }

  bool compute_events_of_pvertex(
    const PVertex& pvertex,
    const std::vector<IEdge>& iedges,
    const std::vector<Segment_2>& segments,
    const std::vector<Bbox_2>& bboxes) {

    CGAL_assertion(iedges.size() > 0);
    CGAL_assertion(iedges.size() == segments.size());
    CGAL_assertion(iedges.size() == bboxes.size());

    std::cout.precision(20);
    if (m_data.is_frozen(pvertex)) {
      return false;
    }

    CGAL_assertion(
      CGAL::abs(m_max_time - m_min_time) >= KSR::tolerance<FT>());
    const auto pv_min = m_data.point_2(pvertex, m_min_time);
    const auto pv_max = m_data.point_2(pvertex, m_max_time);
    const Segment_2 pv_segment(pv_min, pv_max);
    const auto pv_bbox = pv_segment.bbox();

    if (m_data.has_iedge(pvertex)) {
      compute_events_of_constrained_pvertex(
        pvertex, pv_segment, pv_bbox);
    } else {
      compute_events_of_unconstrained_pvertex(
        pvertex, pv_segment, pv_bbox, iedges, segments, bboxes);
    }
    m_queue.finalize_pushing();
    return true;
  }

  void compute_events_of_constrained_pvertex(
    const PVertex& pvertex, const Segment_2& pv_segment, const Bbox_2& pv_bbox) {

    // const bool is_event_found =
    // try_pvertices_to_ivertex_event(pvertex, pv_segment, pv_bbox);
    // if (!is_event_found) return;

    try_pvertex_to_pvertex_constrained_event(pvertex, pv_segment, pv_bbox);
    try_pvertex_to_ivertex_constrained_event(pvertex, pv_segment);
  }

  bool try_pvertices_to_ivertex_event(
    const PVertex& pvertex, const Segment_2& pv_segment, const Bbox_2& pv_bbox) {
    bool is_event_found = false;

    PVertex prev, next;
    std::tie(prev, next) = m_data.prev_and_next(pvertex);
    for (const auto& pother : { prev, next }) {
      if (pother == m_data.null_pvertex()
          || !m_data.is_active(pother)
          ||  m_data.has_iedge(pother)) {
        continue;
      }

      const Segment_2 po_segment(
        m_data.point_2(pother, m_min_time),
        m_data.point_2(pother, m_max_time));
      const auto po_bbox = po_segment.bbox();

      if (!do_overlap(pv_bbox, po_bbox)) {
        continue;
      }

      Point_2 inter;
      if (!KSR::intersection(pv_segment, po_segment, inter)) {
        continue;
      }

      CGAL_assertion(m_data.has_iedge(pvertex));
      const auto iedge = m_data.iedge(pvertex);

      const auto isource = m_data.source(iedge);
      const auto itarget = m_data.target(iedge);

      const auto source = m_data.point_2(pvertex.first, isource);
      const auto target = m_data.point_2(pvertex.first, itarget);

      const FT ptol = KSR::point_tolerance<FT>();
      const FT dist1 = KSR::distance(inter, source);
      const FT dist2 = KSR::distance(inter, target);

      // std::cout << "ptol: " << ptol << std::endl;
      // std::cout << "dist 1: " << dist1 << std::endl;
      // std::cout << "dist 2: " << dist2 << std::endl;

      Point_2 ipoint;
      IVertex ivertex = m_data.null_ivertex();
      if (dist1 < ptol) {
        CGAL_assertion(dist2 >= ptol);
        ipoint = source; ivertex = isource;
      } else if (dist2 < ptol) {
        CGAL_assertion(dist1 >= ptol);
        ipoint = target; ivertex = itarget;
      }

      if (ivertex != m_data.null_ivertex()) {
        CGAL_assertion(ipoint != Point_2());

        const auto& pinit = pv_segment.source();
        const FT distance = KSR::distance(pinit, ipoint);
        const FT time = distance / m_data.speed(pvertex);

        // Should I break here?
        is_event_found = true;
        CGAL_assertion(time < m_max_time - m_min_time);
        m_queue.push(Event(true, pvertex, pother, ivertex, m_min_time + time));
      }
    }

    CGAL_assertion_msg(false, "TODO: TRY PVERTICES TO IVERTEX EVENT!");
    return is_event_found;
  }

  void try_pvertex_to_pvertex_constrained_event(
    const PVertex& pvertex, const Segment_2& pv_segment, const Bbox_2& pv_bbox) {

    // std::cout << "min time: " << m_min_time << std::endl;
    // std::cout << "max time: " << m_max_time << std::endl;
    // std::cout << "cur time: " << m_data.current_time() << std::endl;

    // std::cout << "pvertex: " << m_data.str(pvertex) << std::endl;
    // std::cout << "direction: " << m_data.direction(pvertex) << std::endl;
    // std::cout << "p: " << m_data.point_3(pvertex, m_min_time) << std::endl;
    // std::cout << "q: " << m_data.point_3(pvertex, m_max_time) << std::endl;

    // std::cout << "pv segment: " <<
    //   m_data.to_3d(pvertex.first, source_p) << " " <<
    //   m_data.to_3d(pvertex.first, target_p) << std::endl;

    PVertex prev, next;
    std::tie(prev, next) = m_data.prev_and_next(pvertex);
    for (const auto& pother : { prev, next }) {
      if (pother == m_data.null_pvertex()
          || !m_data.is_active(pother)
          ||  m_data.has_iedge(pother)) {
        continue;
      }

      CGAL_assertion_msg(KSR::distance(
        pv_segment.source(), pv_segment.target()) >= KSR::point_tolerance<FT>(),
      "TODO: ZERO LENGTH PV_SEGMENT FOUND!");

      const Segment_2 po_segment(
        m_data.point_2(pother, m_min_time),
        m_data.point_2(pother, m_max_time));
      CGAL_assertion_msg(KSR::distance(
        po_segment.source(), po_segment.target()) >= KSR::point_tolerance<FT>(),
      "TODO: ZERO LENGTH PO_SEGMENT FOUND!");

      const auto po_bbox = po_segment.bbox();
      if (!do_overlap(pv_bbox, po_bbox)) {
        continue;
      }

      Point_2 inter;
      if (!KSR::intersection(pv_segment, po_segment, inter)) {
        continue;
      }

      const auto& pinit = pv_segment.source();
      const FT distance = KSR::distance(pinit, inter);
      const FT time = distance / m_data.speed(pvertex);

      // Constrained pvertex to another pvertex event.
      CGAL_assertion(time < m_max_time - m_min_time);
      m_queue.push(Event(true, pvertex, pother, m_min_time + time));

      // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
      // std::cout << "pother: "  << m_data.point_3(pother)  << std::endl;
    }
  }

  void try_pvertex_to_ivertex_constrained_event(
    const PVertex& pvertex, const Segment_2& pv_segment) {

    CGAL_assertion(m_data.has_iedge(pvertex));
    const auto iedge = m_data.iedge(pvertex);
    const auto& source_p = pv_segment.source();
    const auto& target_p = pv_segment.target();
    const FT ptol = KSR::point_tolerance<FT>();

    // std::cout << "min time: " << m_min_time << std::endl;
    // std::cout << "max time: " << m_max_time << std::endl;
    // std::cout << "cur time: " << m_data.current_time() << std::endl;

    // std::cout << "pvertex: " << m_data.str(pvertex) << std::endl;
    // std::cout << "direction: " << m_data.direction(pvertex) << std::endl;
    // std::cout << "p: " << m_data.point_3(pvertex, m_min_time) << std::endl;
    // std::cout << "q: " << m_data.point_3(pvertex, m_max_time) << std::endl;

    // std::cout << "pv segment: " <<
    //   m_data.to_3d(pvertex.first, source_p) << " " <<
    //   m_data.to_3d(pvertex.first, target_p) << std::endl;

    CGAL_assertion_msg(KSR::distance(
      m_data.point_3(m_data.source(iedge)),
      m_data.point_3(m_data.target(iedge))) >= ptol,
    "TODO: ZERO-LENGTH IEDGE FOUND!");

    for (const auto& ivertex : { m_data.source(iedge), m_data.target(iedge) }) {
      if (!m_data.is_active(ivertex)) {
        continue;
      }

      const Point_2 ipoint = m_data.to_2d(pvertex.first, ivertex);
      // std::cout << "po segment: " <<
      //   m_data.to_3d(pvertex.first, source_p) << " " <<
      //   m_data.to_3d(pvertex.first, ipoint) << std::endl;

      const FT distance = KSR::distance(source_p, ipoint);
      if (distance < ptol) {

        const auto overtex = m_data.opposite(iedge, ivertex);
        const Point_2 opoint = m_data.to_2d(pvertex.first, overtex);
        CGAL_assertion_msg(KSR::distance(source_p, opoint) >= ptol,
          "TODO: ZERO-LENGTH VECTOR FOUND!");

        // Here, in the dot product, we can have maximum 1 zero-length vector.
        const Vector_2 vec1(source_p, target_p);
        const Vector_2 vec2(source_p, opoint);
        const FT dot_product = vec1 * vec2;
        if (dot_product >= FT(0)) continue;

      } else {

        // Here, in the dot product, we can have maximum 1 zero-length vector.
        CGAL_assertion(distance >= ptol);
        const Vector_2 vec1(source_p, target_p);
        const Vector_2 vec2(source_p, ipoint);
        const FT dot_product = vec1 * vec2;
        if (dot_product < FT(0)) continue; // opposite directions
      }

      // Constrained pvertex to ivertex event.
      // std::cout << "before" << std::endl;
      const FT time = distance / m_data.speed(pvertex);
      if (time < m_max_time - m_min_time) {

        // std::cout << "after" << std::endl;
        CGAL_assertion(time < m_max_time - m_min_time);
        m_queue.push(Event(true, pvertex, ivertex, m_min_time + time));

        // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
        // std::cout << "ivertex: " << m_data.point_3(ivertex) << std::endl;
      }
    }
  }

  void compute_events_of_unconstrained_pvertex(
    const PVertex& pvertex,
    const Segment_2& pv_segment,
    const Bbox_2& pv_bbox,
    const std::vector<IEdge>& iedges,
    const std::vector<Segment_2>& segments,
    const std::vector<Bbox_2>& bboxes) {

    try_pvertex_to_iedge_unconstrained_event(
      pvertex, pv_segment, pv_bbox, iedges, segments, bboxes);
  }

  void try_pvertex_to_iedge_unconstrained_event(
    const PVertex& pvertex,
    const Segment_2& pv_segment,
    const Bbox_2& pv_bbox,
    const std::vector<IEdge>& iedges,
    const std::vector<Segment_2>& segments,
    const std::vector<Bbox_2>& bboxes) {

    // std::cout << "min time: " << m_min_time << std::endl;
    // std::cout << "max time: " << m_max_time << std::endl;
    // std::cout << "cur time: " << m_data.current_time() << std::endl;

    // std::cout << "pvertex: " << m_data.str(pvertex) << std::endl;
    // std::cout << "direction: " << m_data.direction(pvertex) << std::endl;
    // std::cout << "p: " << m_data.point_3(pvertex, m_min_time) << std::endl;
    // std::cout << "q: " << m_data.point_3(pvertex, m_max_time) << std::endl;

    // std::cout << "pv segment: " <<
    //   m_data.to_3d(pvertex.first, source_p) << " " <<
    //   m_data.to_3d(pvertex.first, target_p) << std::endl;

    const auto prev = m_data.prev(pvertex);
    const auto next = m_data.next(pvertex);
    for (std::size_t i = 0; i < iedges.size(); ++i) {
      const auto& iedge = iedges[i];

      if (m_data.iedge(prev) == iedge ||
          m_data.iedge(next) == iedge) {
        continue;
      }

      if (!m_data.is_active(iedge)) {
        continue;
      }

      CGAL_assertion_msg(KSR::distance(
        pv_segment.source(), pv_segment.target()) >= KSR::point_tolerance<FT>(),
      "TODO: ZERO LENGTH PV_SEGMENT FOUND!");

      CGAL_assertion_msg(KSR::distance(
        segments[i].source(), segments[i].target()) >= KSR::point_tolerance<FT>(),
      "TODO: ZERO LENGTH PI_SEGMENT FOUND!");

      if (!CGAL::do_overlap(pv_bbox, bboxes[i])) {
        continue;
      }

      Point_2 inter;
      if (!KSR::intersection(pv_segment, segments[i], inter)) {
        continue;
      }

      // Try to add unconstrained pvertex to ivertex event.
      const auto& pinit = pv_segment.source();
      // const bool is_event_found = try_pvertex_to_ivertex_unconstrained_event(
      //   pvertex, iedge, inter, pinit);

      // Otherwise we add unconstrained pvertex to iedge event.
      // if (!is_event_found) {

      const FT distance = KSR::distance(pinit, inter);
      const FT time = distance / m_data.speed(pvertex);
      CGAL_assertion(time < m_max_time - m_min_time);
      m_queue.push(Event(false, pvertex, iedge, m_min_time + time));

      // }

      // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
      // std::cout << "iedge: "   << m_data.segment_3(iedge) << std::endl;
    }
  }

  bool try_pvertex_to_ivertex_unconstrained_event(
    const PVertex& pvertex, const IEdge& iedge,
    const Point_2& inter, const Point_2& pinit) {

    bool is_event_found = false;
    const auto isource = m_data.source(iedge);
    const auto itarget = m_data.target(iedge);

    const auto source = m_data.point_2(pvertex.first, isource);
    const auto target = m_data.point_2(pvertex.first, itarget);

    const FT tol = KSR::tolerance<FT>();
    const FT dist1 = KSR::distance(inter, source);
    const FT dist2 = KSR::distance(inter, target);

    // std::cout << "tol: " << tol << std::endl;
    // std::cout << "dist 1: " << dist1 << std::endl;
    // std::cout << "dist 2: " << dist2 << std::endl;

    Point_2 ipoint;
    IVertex ivertex = m_data.null_ivertex();
    if (dist1 < tol) {
      CGAL_assertion(dist2 >= tol);
      ipoint = source; ivertex = isource;
    } else if (dist2 < tol) {
      CGAL_assertion(dist1 >= tol);
      ipoint = target; ivertex = itarget;
    }

    if (ivertex != m_data.null_ivertex()) {
      CGAL_assertion(ipoint != Point_2());
      const FT distance = KSR::distance(pinit, ipoint);
      const FT time = distance / m_data.speed(pvertex);
      CGAL_assertion(time < m_max_time - m_min_time);
      m_queue.push(Event(false, pvertex, ivertex, m_min_time + time));
      is_event_found = true;
    }

    CGAL_assertion_msg(false, "TODO: ADD PVERTEX TO IVERTEX UNCONSTRAINED EVENT!");
    return is_event_found;
  }

  /*******************************
  **          RUNNING           **
  ********************************/

  std::size_t run(
    const std::size_t initial_iteration) {

    if (m_debug) {
      std::cout << "* unstacking queue, current size: " << m_queue.size() << std::endl;
    }

    std::size_t iteration = initial_iteration;
    while (!m_queue.empty()) {
      // m_queue.print();

      const Event event = m_queue.pop();
      const FT current_time = event.time();

      // const std::size_t sp_debug_idx = 17;
      if (m_export /* && event.pvertex().first == sp_debug_idx */) {
        if (iteration < 10) {
          dump(m_data, "iter-0" + std::to_string(iteration));
          // dump_2d_surface_mesh(m_data, sp_debug_idx, "iter-" + std::to_string(iteration) +
          //   "-surface-mesh-" + std::to_string(sp_debug_idx));
          dump_event(m_data, event, "iter-0" + std::to_string(iteration));
        } else {
          dump(m_data, "iter-" + std::to_string(iteration));
          // dump_2d_surface_mesh(m_data, sp_debug_idx, "iter-" + std::to_string(iteration) +
          //   "-surface-mesh-" + std::to_string(sp_debug_idx));
          dump_event(m_data, event, "iter-" + std::to_string(iteration));
        }
      }

      m_data.update_positions(current_time);
      if (m_debug) {
        std::cout << std::endl << "* APPLYING " << iteration << ": " << event << std::endl;
      }
      ++iteration;

      // if (iteration == 54) {
      //   exit(EXIT_FAILURE);
      // }

      apply(event);
      CGAL_assertion(m_data.check_integrity());
    }
    return iteration;
  }

  void apply(const Event& event) {

    const auto pvertex = event.pvertex();
    if (event.is_pvertices_to_ivertex()) {

      const auto pother  = event.pother();
      const auto ivertex = event.ivertex();
      apply_event_pvertices_meet_ivertex(pvertex, pother, ivertex, event);

    } else if (event.is_pvertex_to_pvertex()) {
      const auto pother = event.pother();

      remove_events(pvertex);
      remove_events(pother);

      if (m_data.has_iedge(pvertex)) {
        CGAL_assertion(m_data.has_iedge(pvertex));
        if (m_data.has_iedge(pother)) {
          apply_event_two_constrained_pvertices_meet(pvertex, pother, event);
        } else {
          apply_event_constrained_pvertex_meets_free_pvertex(pvertex, pother, event);
        }
      } else {
        CGAL_assertion(!m_data.has_iedge(pvertex));
        if (!m_data.has_iedge(pother)) {
          apply_event_two_unconstrained_pvertices_meet(pvertex, pother, event);
        } else {
          CGAL_assertion_msg(false, "ERROR: THIS EVENT SHOULD NOT EVER HAPPEN!");
          apply_event_constrained_pvertex_meets_free_pvertex(pother, pvertex, event);
        }
      }
    } else if (event.is_pvertex_to_iedge()) {

      const auto iedge = event.iedge();
      if (m_data.has_iedge(pvertex)) {
        apply_event_constrained_pvertex_meets_iedge(pvertex, iedge, event);
      } else {
        const bool is_event_happend = apply_event_unconstrained_pedge_meets_iedge(
          pvertex, iedge, event);
        if (!is_event_happend) {
          apply_event_unconstrained_pvertex_meets_iedge(pvertex, iedge, event);
        }
      }
    } else if (event.is_pvertex_to_ivertex()) {

      const auto ivertex = event.ivertex();
      if (m_data.has_iedge(pvertex)) {
        apply_event_constrained_pvertex_meets_ivertex(pvertex, ivertex, event);
      } else {
        apply_event_unconstrained_pvertex_meets_ivertex(pvertex, ivertex, event);
      }
    } else {
      CGAL_assertion_msg(false, "ERROR: INVALID EVENT FOUND!");
    }
  }

  /*******************************
  **        HANDLE EVENTS       **
  ********************************/

  // INVALID EVENTS!
  void apply_event_two_unconstrained_pvertices_meet(
    const PVertex& /* pvertex */,
    const PVertex& /* pother  */,
    const Event&   /* event   */) {

    // if (m_debug) {
    //   std::cout << "WARNING: SKIPPING TWO UNCONSTRAINED PVERTICES MEET EVENT!" << std::endl;
    //   std::cout << "WARNING: THIS EVENT IS DUST!" << std::endl;
    //   m_queue.print();
    //   CGAL_assertion_msg(
    //     false, "TODO: CHECK, SEE CONSTR. PV. MEETS IEDGE FOR MORE DETAILS!");
    // }
    // return; // skip

    CGAL_assertion_msg(false,
    "ERROR: TWO UNCONSTRAINED PVERTICES MEET! DO WE HAVE A CONCAVE POLYGON?");
  }

  void apply_event_two_constrained_pvertices_meet(
    const PVertex& /* pvertex */,
    const PVertex& /* pother  */,
    const Event&   /* event   */) {

    // if (m_debug) {
    //   std::cout << "WARNING: SKIPPING TWO CONSTRAINED PVERTICES MEET EVENT!" << std::endl;
    //   std::cout << "WARNING: THIS EVENT IS DUST!" << std::endl;
    //   m_queue.print();
    //   CGAL_assertion_msg(
    //     false, "TODO: CHECK, SEE CONSTR. PV. MEETS IEDGE FOR MORE DETAILS!");
    // }
    // return; // skip

    CGAL_assertion_msg(false,
    "ERROR: TWO CONSTRAINED PVERTICES MEET! CAN IT HAPPEN?");
  }

  void apply_event_constrained_pvertex_meets_iedge(
    const PVertex& /* pvertex */,
    const IEdge&   /* iedge   */,
    const Event&   /* event   */) {

    // if (m_debug) {
    //   std::cout << "WARNING: SKIPPING CONSTRAINED PVERTEX MEETS IEDGE EVENT!" << std::endl;
    //   std::cout << "WARNING: THIS EVENT IS DUST!" << std::endl;
    //   // m_queue.print();
    // }

    // In this case what happens is:
    // We push multiple events between pvertex and iedges.
    // Several of these events with the closest time are handled,
    // however several can stay because they happen along iedges, which
    // are not direct neighbors of the handled events and so they stay.
    // Here, however, these events are useless and can be safely ignored.
    // This is a solution that is off for now. I found a better one.
    // return;

    CGAL_assertion_msg(false,
    "ERROR: CONSTRAINED PVERTEX MEETS IEDGE! WHAT IS WRONG?");
  }

  void apply_event_pvertices_meet_ivertex(
    const PVertex& pvertex, const PVertex& pother,
    const IVertex& /* ivertex */, const Event& /* event */) {

    CGAL_assertion( m_data.has_iedge(pvertex));
    CGAL_assertion(!m_data.has_iedge(pother));

    // if (m_debug) {
    //   std::cout << "WARNING: SKIPPING PVERTICES MEET IVERTEX EVENT!" << std::endl;
    //   std::cout << "WARNING: THIS EVENT IS DUST!" << std::endl;
    //   m_queue.print();
    //   CGAL_assertion_msg(
    //     false, "TODO: CHECK, SEE CONSTR. PV. MEETS IEDGE FOR MORE DETAILS!");
    // }
    // return; // skip

    CGAL_assertion_msg(false,
    "ERROR: PVERTICES MEET IVERTEX! IT SHOULD NOT EVER HAPPEN!");
  }

  void apply_event_unconstrained_pvertex_meets_ivertex(
    const PVertex& pvertex, const IVertex& /* ivertex */, const Event& /* event */) {

    CGAL_assertion(!m_data.has_iedge(pvertex));
    CGAL_assertion( m_data.has_one_pface(pvertex));

    // if (m_debug) {
    //   std::cout << "WARNING: SKIPPING UNCONSTRAINED PVERTEX MEETS IVERTEX EVENT!" << std::endl;
    //   std::cout << "WARNING: THIS EVENT IS DUST!" << std::endl;
    //   m_queue.print();
    //   CGAL_assertion_msg(
    //     false, "TODO: CHECK, SEE CONSTR. PV. MEETS IEDGE FOR MORE DETAILS!");
    // }
    // return; // skip

    CGAL_assertion_msg(false,
    "ERROR: UNCONSTRAINED PVERTEX MEETS IVERTEX! IT SHOULD NOT EVER HAPPEN!");

    // apply_event_pvertex_meets_ivertex(pvertex, ivertex, event);
  }

  // VALID EVENTS!
  void apply_event_pvertex_meets_ivertex(
    const PVertex& pvertex, const IVertex& ivertex, const Event& event) {

    // First, let's gather all pvertices that will get merged.
    const std::vector<PVertex> crossed_pvertices =
      m_data.pvertices_around_ivertex(pvertex, ivertex);

    // Remove associated events.
    CGAL_assertion(crossed_pvertices.size() >= 3);
    for (std::size_t i = 1; i < crossed_pvertices.size() - 1; ++i) {
      remove_events(crossed_pvertices[i]);
    }

    // Merge them and get the newly created pvertices.
    CGAL_assertion(!m_data.has_ivertex(pvertex));
    std::vector< std::pair<IEdge, bool> > crossed_iedges;
    const std::vector<PVertex> pvertices =
      merge_pvertices_on_ivertex(
        m_min_time, m_max_time, ivertex, event.pvertex(),
        crossed_pvertices, crossed_iedges);

    // Remove all events of the crossed iedges.
    CGAL_assertion(crossed_iedges.size() >= 1);
    for (const auto& crossed_iedge : crossed_iedges) {
      // TODO: SHOULD I LEAVE THIS CHECK? WILL IT MAKE THE CODE FASTER?
      // if (crossed_iedges[ip].second) {
        // bla bla
      // }
      const auto& iedge = crossed_iedge.first;
      remove_events(iedge, pvertex.first);
    }

    // In general, pvertices in this container are newly created that is
    // they are either cropped or propagated. However, in parallel cases,
    // they may be the ones, which are prev or next or, in other words, either
    // first or last in the crossed_pvertices above. The first and last there
    // are skipped and their events are not removed, so we remove them here,
    // to be certain.
    for (const auto& pvertex : pvertices) {
      if (pvertex == m_data.null_pvertex()) continue;
      remove_events(pvertex);
    }

    // And compute new events.
    CGAL_assertion(pvertices.size() > 0);
    compute_events_of_pvertices(event.time(), pvertices);
    // CGAL_assertion_msg(false, "TODO: PVERTEX MEETS IVERTEX!");
  }

  void apply_event_unconstrained_pvertex_meets_iedge(
    const PVertex& pvertex, const IEdge& iedge, const Event& event) {

    CGAL_assertion(!m_data.has_iedge(pvertex));
    CGAL_assertion( m_data.has_one_pface(pvertex));

    remove_events(pvertex);
    const bool stop = check_stop_condition(pvertex, iedge);

    if (stop) { // polygon stops
      const PVertex pother =
        crop_pvertex_along_iedge(pvertex, iedge);
      const std::array<PVertex, 2> pvertices = {pvertex, pother};
      remove_events(iedge, pvertex.first);
      compute_events_of_pvertices(event.time(), pvertices);
    } else { // polygon continues beyond the iedge
      const std::array<PVertex, 3> pvertices =
        propagate_pvertex_beyond_iedge(pvertex, iedge);
      remove_events(iedge, pvertex.first);
      compute_events_of_pvertices(event.time(), pvertices);
    }
    CGAL_assertion(m_data.has_iedge(pvertex));
    // CGAL_assertion_msg(false, "TODO: UNCONSTRAINED PVERTEX MEETS IEDGE!");
  }

  bool apply_event_unconstrained_pedge_meets_iedge(
    const PVertex& pvertex, const IEdge& iedge, const Event& event) {

    bool is_event_happend = false;
    const auto prev = m_data.prev(pvertex);
    const auto next = m_data.next(pvertex);
    const auto isegment = m_data.segment_2(pvertex.first, iedge);

    for (const auto& pother : { prev, next }) {
      const Segment_2 segment(
        m_data.point_2(pother , event.time()),
        m_data.point_2(pvertex, event.time()));
      CGAL_assertion(segment.squared_length() != FT(0));

      bool both_are_free = true;
      if (m_data.has_iedge(pvertex) || m_data.has_iedge(pother)) {
        both_are_free = false;
      }

      if (both_are_free && KSR::are_parallel(segment, isegment)) {
        CGAL_assertion(!m_data.has_iedge(pother));
        CGAL_assertion(!m_data.has_iedge(pvertex));

        CGAL_assertion(m_data.has_one_pface(pother));
        CGAL_assertion(m_data.has_one_pface(pvertex));

        remove_events(pother);
        remove_events(pvertex);

        const bool stop = check_stop_condition(pvertex, pother, iedge);

        if (stop) { // polygon stops
          crop_pedge_along_iedge(pvertex, pother, iedge);
          const auto pvertices = std::array<PVertex, 2>{pvertex, pother};
          remove_events(iedge, pvertex.first);
          compute_events_of_pvertices(event.time(), pvertices);
        } else { // polygon continues beyond the edge
          PVertex pv0, pv1;
          std::tie(pv0, pv1) = propagate_pedge_beyond_iedge(pvertex, pother, iedge);
          const auto pvertices = std::array<PVertex, 4>{pvertex, pother, pv0, pv1};
          remove_events(iedge, pvertex.first);
          compute_events_of_pvertices(event.time(), pvertices);
        }

        CGAL_assertion(m_data.has_iedge(pother));
        CGAL_assertion(m_data.has_iedge(pvertex));
        CGAL_assertion(m_data.iedge(pvertex) == m_data.iedge(pother));
        is_event_happend = true;
        // CGAL_assertion_msg(false, "TODO: UNCONSTRAINED PEDGE MEETS IEDGE!");
        break;
      }
    }
    return is_event_happend;
  }

  void apply_event_constrained_pvertex_meets_ivertex(
    const PVertex& pvertex, const IVertex& ivertex, const Event& event) {

    CGAL_assertion(m_data.has_iedge(pvertex));
    apply_event_pvertex_meets_ivertex(pvertex, ivertex, event);
    // CGAL_assertion_msg(false, "TODO: CONSTRAINED PVERTEX MEETS IVERTEX!");
  }

  void apply_event_constrained_pvertex_meets_free_pvertex(
    const PVertex& pvertex, const PVertex& pother, const Event& event) {

    CGAL_assertion( m_data.has_iedge(pvertex));
    CGAL_assertion(!m_data.has_iedge(pother));

    if (transfer_pvertex_via_iedge(pvertex, pother)) {

      // Check the first two pvertices.
      if (m_data.has_iedge(pvertex)) {
        remove_events(m_data.iedge(pvertex), pvertex.first);
      }
      if (m_data.has_iedge(pother)) {
        remove_events(m_data.iedge(pother) , pother.first);
      }
      const auto pvertices1 = std::array<PVertex, 2>{pvertex, pother};
      compute_events_of_pvertices(event.time(), pvertices1);

      // Check the last pvertex.
      PVertex prev, next;
      std::tie(prev, next) = m_data.border_prev_and_next(pvertex);
      PVertex pthird = prev;
      if (pthird == pother) {
        pthird = next;
      } else { CGAL_assertion(pother == next); }

      if (m_data.has_iedge(pthird)) {
        remove_events(m_data.iedge(pthird), pthird.first);
      }
      const auto pvertices2 = std::array<PVertex, 1>{pthird};
      compute_events_of_pvertices(event.time(), pvertices2);

    } else {

      if (m_data.has_iedge(pvertex)) {
        remove_events(m_data.iedge(pvertex), pvertex.first);
      }
      const auto pvertices = std::array<PVertex, 1>{pvertex};
      compute_events_of_pvertices(event.time(), pvertices);
    }
    // CGAL_assertion_msg(false, "TODO: CONSTRAINED PVERTEX MEETS FREE PVERTEX!");
  }

  // STOP CONDITIONS!
  bool check_stop_condition(
    const PVertex& pvertex, const IEdge& iedge) {
    return check_pvertex_meets_iedge_global_k(pvertex, iedge);
  }

  // GLOBAL STOP CONDITIONS!
  bool check_pvertex_meets_iedge_global_k(
    const PVertex& pvertex, const IEdge& iedge) {

    if (m_debug) {
      std::cout << "- k intersections before: " << m_data.k(pvertex.first) << std::endl;
    }

    bool is_occupied_iedge, is_bbox_reached;
    std::tie(is_occupied_iedge, is_bbox_reached) = m_data.collision_occured(pvertex, iedge);
    const bool is_limit_line = m_data.update_limit_lines_and_k(pvertex, iedge, is_occupied_iedge);

    if (m_debug) {
      std::cout << "- bbox: " << is_bbox_reached  << "; " <<
      " limit: "    << is_limit_line << "; " <<
      " occupied: " << is_occupied_iedge << std::endl;
    }

    bool stop = false;
    if (is_bbox_reached) {
      if (m_debug) std::cout << "- bbox, stop" << std::endl;
      stop = true;
    } else if (is_limit_line) {
      if (m_debug) std::cout << "- limit, stop" << std::endl;
      stop = true;
    } else {
      if (m_debug) std::cout << "- free, any k, continue" << std::endl;
      stop = false;
    }
    CGAL_assertion(m_data.k(pvertex.first) >= 1);
    if (m_debug) {
      std::cout << "- k intersections after: " << m_data.k(pvertex.first) << std::endl;
    }
    // CGAL_assertion_msg(false, "TODO: CHECK PVERTEX MEETS IVERTEX GLOBAL!");
    return stop;
  }

  bool check_stop_condition(
    const PVertex& pvertex, const PVertex& pother, const IEdge& iedge) {
    return check_pedge_meets_iedge_global_k(pvertex, pother, iedge);
  }

  bool check_pedge_meets_iedge_global_k(
    const PVertex& pvertex, const PVertex& pother, const IEdge& iedge) {

    if (m_debug) {
      std::cout << "- k intersections before: " << m_data.k(pvertex.first) << std::endl;
    }

    bool is_occupied_iedge_1, is_bbox_reached_1;
    std::tie(is_occupied_iedge_1, is_bbox_reached_1) = m_data.collision_occured(pvertex, iedge);
    bool is_occupied_iedge_2, is_bbox_reached_2;
    std::tie(is_occupied_iedge_2, is_bbox_reached_2) = m_data.collision_occured(pother, iedge);

    const bool is_limit_line_1 = m_data.update_limit_lines_and_k(pvertex, iedge, is_occupied_iedge_1);
    const bool is_limit_line_2 = m_data.update_limit_lines_and_k(pother , iedge, is_occupied_iedge_2);

    if (m_debug) {
      std::cout << "- bbox1: " << is_bbox_reached_1  << "; " <<
      " limit1: "    << is_limit_line_1 << "; " <<
      " occupied1: " << is_occupied_iedge_1 << std::endl;
      std::cout << "- bbox2: " << is_bbox_reached_2  << "; " <<
      " limit2: "    << is_limit_line_2 << "; " <<
      " occupied2: " << is_occupied_iedge_2 << std::endl;
    }
    CGAL_assertion(is_limit_line_1 == is_limit_line_2);
    CGAL_assertion(is_bbox_reached_1 == is_bbox_reached_2);

    bool stop = false;
    if (is_bbox_reached_1 || is_bbox_reached_2) {
      if (m_debug) std::cout << "- bbox, stop" << std::endl;
      stop = true;
    } else if (is_limit_line_1 || is_limit_line_2) {
      if (m_debug) std::cout << "- limit, stop" << std::endl;
      stop = true;
    } else {
      if (m_debug) std::cout << "- free, any k, continue" << std::endl;
      CGAL_assertion(!m_data.is_sneaking_pedge(pvertex, pother, iedge));
      stop = false;
    }

    CGAL_assertion(m_data.k(pvertex.first) >= 1);
    if (m_debug) {
      std::cout << "- k intersections after: " << m_data.k(pvertex.first) << std::endl;
    }
    // CGAL_assertion_msg(false, "TODO: CHECK PEDGE MEETS IEDGE GLOBAL!");
    return stop;
  }

  // RECOMPUTE EVENTS!
  template<typename PVertexRange>
  void compute_events_of_pvertices(
    const FT last_event_time, const PVertexRange& pvertices) {

    m_min_time = m_data.current_time();
    m_data.update_positions(m_max_time);

    const auto& pfront = pvertices.front();
    CGAL_assertion(pfront != m_data.null_pvertex());
    const auto& iedges   = m_data.iedges(pfront.first);
    const auto& segments = m_data.isegments(pfront.first);
    const auto& bboxes   = m_data.ibboxes(pfront.first);

    for (const auto& pvertex : pvertices) {
      if (pvertex == m_data.null_pvertex()) continue;
      m_data.deactivate(pvertex);
    }
    for (const auto& pvertex : pvertices) {
      if (pvertex == m_data.null_pvertex()) continue;
      m_data.set_last_event_time(pvertex, last_event_time);
      compute_events_of_pvertex(pvertex, iedges, segments, bboxes);
    }
    for (const auto& pvertex : pvertices) {
      if (pvertex == m_data.null_pvertex()) continue;
      m_data.activate(pvertex);
    }
    m_data.update_positions(m_min_time);
  }

  // REMOVE EVENTS!
  // Remove events associated with the given iedge.
  void remove_events(const IEdge& iedge, const std::size_t support_plane_idx) {
    CGAL_assertion(iedge != m_data.null_iedge());
    m_queue.erase_vertex_events(iedge, support_plane_idx);
    // std::cout << "erasing events for iedge: " << m_data.str(iedge) << std::endl;
    // std::cout << "iedge: " << m_data.segment_3(iedge) << std::endl;
  }

  // Remove events associated with the given pvertex.
  void remove_events(const PVertex& pvertex) {
    CGAL_assertion(pvertex != m_data.null_pvertex());
    m_queue.erase_vertex_events(pvertex);
    // std::cout << "erasing events for pvertex: " << m_data.str(pvertex) << std::endl;
    // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
  }

  /*******************************
  **    OPERATIONS ON POLYGONS  **
  ********************************/

  PVertex crop_pvertex_along_iedge(
    const PVertex& pvertex, const IEdge& iedge) {

    if (m_debug) {
      std::cout.precision(20);
      std::cout << "** cropping " << m_data.str(pvertex) << " along " << m_data.str(iedge) << std::endl;
      std::cout << "- pvertex: " << m_data.point_3(pvertex) << std::endl;
      std::cout << "- iedge: "   << m_data.segment_3(iedge) << std::endl;
    }

    CGAL_assertion_msg(KSR::distance(
      m_data.point_2(pvertex.first, m_data.source(iedge)),
      m_data.point_2(pvertex.first, m_data.target(iedge))) >= KSR::point_tolerance<FT>(),
    "TODO: PVERTEX -> IEDGE, HANDLE ZERO-LENGTH IEDGE!");

    const PVertex prev(pvertex.first, m_data.support_plane(pvertex).prev(pvertex.second));
    const PVertex next(pvertex.first, m_data.support_plane(pvertex).next(pvertex.second));

    Point_2 future_point_a, future_point_b;
    Vector_2 future_direction_a, future_direction_b;
    bool is_parallel_a = false, is_parallel_b = false;
    std::tie(is_parallel_a, is_parallel_b) =
    m_data.compute_future_points_and_directions(
      pvertex, IVertex(), iedge,
      future_point_a, future_point_b,
      future_direction_a, future_direction_b);
    CGAL_assertion(future_direction_a != Vector_2());
    CGAL_assertion(future_direction_b != Vector_2());
    if (is_parallel_a || is_parallel_b) {
      if (m_debug) std::cout << "- pvertex to iedge, parallel case" << std::endl;
      // CGAL_assertion_msg(!is_parallel_a && !is_parallel_b,
      // "TODO: PVERTEX -> IEDGE, HANDLE CASE WITH PARALLEL LINES!");
    }

    const PEdge pedge(pvertex.first, m_data.support_plane(pvertex).split_vertex(pvertex.second));
    CGAL_assertion(m_data.source(pedge) == pvertex || m_data.target(pedge) == pvertex);
    const PVertex pother = m_data.opposite(pedge, pvertex);
    if (m_debug) {
      std::cout << "- new pedge: " << m_data.str(pedge)  << " between "
        << m_data.str(pvertex) << " and " << m_data.str(pother) << std::endl;
    }

    m_data.connect(pedge, iedge);
    m_data.connect(pvertex, iedge);
    m_data.connect(pother, iedge);

    m_data.support_plane(pvertex).set_point(pvertex.second, future_point_a);
    m_data.support_plane(pother).set_point(pother.second, future_point_b);
    m_data.direction(pvertex) = future_direction_a;
    m_data.direction(pother) = future_direction_b;

    if (m_debug) std::cout << "- new pvertices: " <<
    m_data.str(pother) << ": " << m_data.point_3(pother) << std::endl;

    // CGAL_assertion_msg(false, "TODO: CROP PVERTEX ALONG IEDGE!");
    return pother;
  }

  std::array<PVertex, 3> propagate_pvertex_beyond_iedge(
    const PVertex& pvertex, const IEdge& iedge) {

    if (m_debug) {
      std::cout.precision(20);
      std::cout << "** propagating " << m_data.str(pvertex) << " beyond " << m_data.str(iedge) << std::endl;
      std::cout << "- pvertex: " << m_data.point_3(pvertex) << std::endl;
      std::cout << "- iedge: "   << m_data.segment_3(iedge) << std::endl;
    }

    const Point_2 original_point = m_data.point_2(pvertex, FT(0));
    const Vector_2 original_direction = m_data.direction(pvertex);
    const PVertex pother = crop_pvertex_along_iedge(pvertex, iedge);

    const PVertex propagated = m_data.add_pvertex(pvertex.first, original_point);
    m_data.direction(propagated) = original_direction;

    if (m_debug) {
      std::cout << "- propagated: " << m_data.str(propagated) << ": " << m_data.point_3(propagated) << std::endl;
    }

    std::array<PVertex, 3> pvertices;
    pvertices[0] = pvertex;
    pvertices[1] = pother;
    pvertices[2] = propagated;

    const PFace new_pface = m_data.add_pface(pvertices);
    CGAL_assertion(new_pface != m_data.null_pface());
    CGAL_assertion(new_pface.second != Face_index());
    if (m_debug) {
      std::cout << "- new pface " << m_data.str(new_pface) << ": " <<
      m_data.centroid_of_pface(new_pface) << std::endl;
    }

    // CGAL_assertion_msg(false, "TODO: PROPAGATE PVERTEX BEYOND IEDGE!");
    return pvertices;
  }

  void crop_pedge_along_iedge(
    const PVertex& pvertex, const PVertex& pother, const IEdge& iedge) {

    if (m_debug) {
      std::cout.precision(20);
      std::cout << "** cropping pedge [" << m_data.str(pvertex) << "-" << m_data.str(pother)
      << "] along " << m_data.str(iedge) << std::endl;
      std::cout << "- pvertex: " << m_data.point_3(pvertex) << std::endl;
      std::cout << "- pother: "  << m_data.point_3(pother)  << std::endl;
      std::cout << "- iedge: "   << m_data.segment_3(iedge) << std::endl;
    }

    CGAL_assertion(pvertex.first == pother.first);
    CGAL_assertion_msg(KSR::distance(
      m_data.point_2(pvertex.first, m_data.source(iedge)),
      m_data.point_2(pvertex.first, m_data.target(iedge))) >= KSR::point_tolerance<FT>(),
    "TODO: PEDGE -> IEDGE, HANDLE ZERO-LENGTH IEDGE!");
    Point_2 future_point; Vector_2 future_direction;

    { // cropping pvertex ...
      const PVertex prev(pvertex.first, m_data.support_plane(pvertex).prev(pvertex.second));
      const PVertex next(pvertex.first, m_data.support_plane(pvertex).next(pvertex.second));

      if (m_debug) {
        std::cout << "- prev pv: " << m_data.point_3(prev) << std::endl;
        std::cout << "- next pv: " << m_data.point_3(next) << std::endl;
      }

      PVertex pthird = m_data.null_pvertex();
      if (pother == prev) {
        pthird = next;
      } else {
        CGAL_assertion(pother == next);
        pthird = prev;
      }
      CGAL_assertion(pthird != m_data.null_pvertex());

      if (m_debug) {
        std::cout << "- pthird pv: " << m_data.point_3(pthird) << std::endl;
      }

      const bool is_parallel = m_data.compute_future_point_and_direction(
        0, IVertex(), pvertex, pthird, iedge, future_point, future_direction);
      CGAL_assertion(future_direction != Vector_2());
      if (is_parallel) {
        if (m_debug) std::cout << "- pedge to iedge 1, parallel case" << std::endl;
        // CGAL_assertion_msg(!is_parallel,
        // "TODO: PEDGE -> IEDGE 1, HANDLE CASE WITH PARALLEL LINES!");
      }

      m_data.direction(pvertex) = future_direction;
      m_data.support_plane(pvertex).set_point(pvertex.second, future_point);
      m_data.connect(pvertex, iedge);
    }

    { // cropping pother ...
      const PVertex prev(pother.first, m_data.support_plane(pother).prev(pother.second));
      const PVertex next(pother.first, m_data.support_plane(pother).next(pother.second));

      if (m_debug) {
        std::cout << "- prev po: " << m_data.point_3(prev) << std::endl;
        std::cout << "- next po: " << m_data.point_3(next) << std::endl;
      }

      PVertex pthird = m_data.null_pvertex();
      if (pvertex == prev) {
        pthird = next;
      } else {
        CGAL_assertion(pvertex == next);
        pthird = prev;
      }
      CGAL_assertion(pthird != m_data.null_pvertex());

      if (m_debug) {
        std::cout << "- pthird po: " << m_data.point_3(pthird) << std::endl;
      }

      const bool is_parallel = m_data.compute_future_point_and_direction(
        0, IVertex(), pother, pthird, iedge, future_point, future_direction);
      CGAL_assertion(future_direction != Vector_2());
      if (is_parallel) {
        if (m_debug) std::cout << "- pedge to iedge 2, parallel case" << std::endl;
        // CGAL_assertion_msg(!is_parallel,
        // "TODO: PEDGE -> IEDGE 2, HANDLE CASE WITH PARALLEL LINES!");
      }

      m_data.direction(pother) = future_direction;
      m_data.support_plane(pother).set_point(pother.second, future_point);
      m_data.connect(pother, iedge);
    }

    const PEdge pedge(pvertex.first, m_data.support_plane(pvertex).edge(pvertex.second, pother.second));
    m_data.connect(pedge, iedge);

    // CGAL_assertion_msg(false, "TODO: CROP PEDGE ALONG IEDGE!");
  }

  std::pair<PVertex, PVertex> propagate_pedge_beyond_iedge(
    const PVertex& pvertex, const PVertex& pother, const IEdge& iedge) {

    if (m_debug) {
      std::cout.precision(20);
      std::cout << "** propagating pedge [" << m_data.str(pvertex) << "-" << m_data.str(pother)
      << "] beyond " << m_data.str(iedge) << std::endl;
      std::cout << "- pvertex: " << m_data.point_3(pvertex) << std::endl;
      std::cout << "- pother: "  << m_data.point_3(pother)  << std::endl;
      std::cout << "- iedge: "   << m_data.segment_3(iedge) << std::endl;
    }

    const Point_2 original_point_1 = m_data.point_2(pvertex, FT(0));
    const Point_2 original_point_2 = m_data.point_2(pother, FT(0));

    const Vector_2 original_direction_1 = m_data.direction(pvertex);
    const Vector_2 original_direction_2 = m_data.direction(pother);

    crop_pedge_along_iedge(pvertex, pother, iedge);

    const PVertex propagated_1 = m_data.add_pvertex(pvertex.first, original_point_1);
    m_data.direction(propagated_1) = original_direction_1;

    const PVertex propagated_2 = m_data.add_pvertex(pother.first, original_point_2);
    m_data.direction(propagated_2) = original_direction_2;

    if (m_debug) {
      std::cout << "- propagated 1: " << m_data.str(propagated_1) << ": " <<
      m_data.point_3(propagated_1) << std::endl;
      std::cout << "- propagated 2: " << m_data.str(propagated_2) << ": " <<
      m_data.point_3(propagated_2) << std::endl;
    }

    std::array<PVertex, 4> pvertices;
    pvertices[0] = pvertex;
    pvertices[1] = pother;
    pvertices[2] = propagated_2;
    pvertices[3] = propagated_1;

    const PFace new_pface = m_data.add_pface(pvertices);
    CGAL_assertion(new_pface != m_data.null_pface());
    CGAL_assertion(new_pface.second != Face_index());
    if (m_debug) {
      std::cout << "- new pface " << m_data.str(new_pface) << ": " << m_data.centroid_of_pface(new_pface) << std::endl;
    }

    // CGAL_assertion_msg(false, "TODO: PROPAGATE PEDGE BEYOND IEDGE!");
    return std::make_pair(propagated_2, propagated_1);
  }

  bool transfer_pvertex_via_iedge(
    const PVertex& pvertex, const PVertex& pother) {

    if (m_debug) {
      std::cout.precision(20);
      CGAL_assertion(m_data.has_iedge(pvertex));
      std::cout << "** transfering " << m_data.str(pother) << " through " << m_data.str(pvertex) << " via "
        << m_data.str(m_data.iedge(pvertex)) << std::endl;
      std::cout << "- pvertex: " << m_data.point_3(pvertex) << std::endl;
      std::cout << "- pother: "  << m_data.point_3(pother)  << std::endl;
    }
    CGAL_assertion(pvertex.first == pother.first);

    // Is pvertex adjacent to one or two pfaces?
    PFace source_pface, target_pface;
    std::tie(source_pface, target_pface) = m_data.pfaces_of_pvertex(pvertex);
    const auto common_pface = m_data.pface_of_pvertex(pother);
    if (common_pface == target_pface) {
      if (m_debug) std::cout << "- swap pfaces" << std::endl;
      std::swap(source_pface, target_pface);
    }
    CGAL_assertion(common_pface == source_pface);

    if (m_debug) {
      std::cout << "- initial pfaces: " << std::endl;
      if (source_pface != m_data.null_pface()) {
        std::cout << "source " << m_data.str(source_pface) << ": " <<
        m_data.centroid_of_pface(source_pface) << std::endl;
      }
      if (target_pface != m_data.null_pface()) {
        std::cout << "target " << m_data.str(target_pface) << ": " <<
        m_data.centroid_of_pface(target_pface) << std::endl;
      }
    }

    // Get pthird.
    PVertex pthird = m_data.next(pother);
    if (pthird == pvertex) pthird = m_data.prev(pother);
    if (m_debug) std::cout << "- pthird: " << m_data.point_3(pthird) << std::endl;

    // Get future point and direction.
    CGAL_assertion(m_data.has_iedge(pvertex));
    const auto iedge = m_data.iedge(pvertex);
    const auto source_p = m_data.point_2(pvertex.first, m_data.source(iedge));
    const auto target_p = m_data.point_2(pvertex.first, m_data.target(iedge));
    CGAL_assertion_msg(KSR::distance(source_p, target_p) >= KSR::point_tolerance<FT>(),
    "TODO: TRANSFER PVERTEX, HANDLE ZERO-LENGTH IEDGE!");
    const Line_2 iedge_line(source_p, target_p);

    Point_2 future_point;
    Vector_2 future_direction;
    const bool is_parallel =
    m_data.compute_future_point_and_direction(
      0, IVertex(), pother, pthird, iedge, future_point, future_direction);
    CGAL_assertion(future_direction != Vector_2());
    if (is_parallel) {
      if (m_debug) std::cout << "- transfer pvertex, parallel case" << std::endl;
      // CGAL_assertion_msg(!is_parallel,
      // "TODO: TRANSFER PVERTEX, HANDLE CASE WITH PARALLEL LINES!");
    }

    if (target_pface == m_data.null_pface()) { // in case we have 1 pface

      m_data.support_plane(pvertex).set_point(pvertex.second, future_point);
      m_data.direction(pvertex) = future_direction;
      const auto he = m_data.mesh(pvertex).halfedge(pother.second, pvertex.second);
      CGAL::Euler::join_vertex(he, m_data.mesh(pvertex));

      // CGAL_assertion_msg(false,
      // "TODO: TRANSFER PVERTEX 1, ADD NEW FUTURE POINTS AND DIRECTIONS!");

    } else { // in case we have both pfaces

      m_data.disconnect_iedge(pvertex);
      PEdge pedge = m_data.null_pedge();
      for (const auto edge : m_data.pedges_around_pvertex(pvertex)) {
        if (m_data.iedge(edge) == iedge) {
          pedge = edge; break;
        }
      }
      CGAL_assertion(pedge != m_data.null_pedge());

      auto he = m_data.mesh(pedge).halfedge(pedge.second);
      if (m_data.mesh(pedge).face(he) != common_pface.second) {
        he = m_data.mesh(pedge).opposite(he);
      }
      CGAL_assertion(m_data.mesh(pedge).face(he) == common_pface.second);

      if (m_data.mesh(pedge).target(he) == pvertex.second) {
        // if (m_debug) std::cout << "- shifting target" << std::endl;
        CGAL::Euler::shift_target(he, m_data.mesh(pedge));
      } else {
        CGAL_assertion(m_data.mesh(pedge).source(he) == pvertex.second);
        // if (m_debug) std::cout << "- shifting source" << std::endl;
        CGAL::Euler::shift_source(he, m_data.mesh(pedge));
      }

      const auto pother_p = m_data.point_2(pother);
      const Point_2 pinit = iedge_line.projection(pother_p);
      m_data.direction(pvertex) = m_data.direction(pother);
      const auto fp = pinit - m_data.direction(pother) * m_data.current_time();
      m_data.support_plane(pvertex).set_point(pvertex.second, fp);

      m_data.support_plane(pother).set_point(pother.second, future_point);
      m_data.direction(pother) = future_direction;
      m_data.connect(pother, iedge);

      // CGAL_assertion_msg(false,
      // "TODO: TRANSFER PVERTEX 2, ADD NEW FUTURE POINTS AND DIRECTIONS!");
    }

    if (m_debug) {
      std::cout << "- new pfaces: " << std::endl;
      if (source_pface != m_data.null_pface()) {
        std::cout << "source " << m_data.str(source_pface) << ": " <<
        m_data.centroid_of_pface(source_pface) << std::endl;
      }
      if (target_pface != m_data.null_pface()) {
        std::cout << "target " << m_data.str(target_pface) << ": " <<
        m_data.centroid_of_pface(target_pface) << std::endl;
      }
    }

    // CGAL_assertion_msg(false, "TODO: TRANSFER PVERTEX VIA IEDGE!");
    return (target_pface != m_data.null_pface());
  }

  std::vector<PVertex> merge_pvertices_on_ivertex(
    const FT min_time, const FT max_time, const IVertex& ivertex,
    const PVertex& event_pvertex, const std::vector<PVertex>& pvertices,
    std::vector< std::pair<IEdge, bool> >& crossed_iedges) {

    if (m_debug) {
      std::cout.precision(20);
      std::cout << "** merging " << m_data.str(pvertices[1]) << " on " << m_data.str(ivertex) << std::endl;
      std::cout << "- pvertex: " << m_data.point_3(pvertices[1]) << std::endl;
      std::cout << "- ivertex: " << m_data.point_3(ivertex) << std::endl;
    }

    CGAL_assertion(pvertices.size() >= 3);
    const std::size_t sp_idx = pvertices.front().first;
    const PVertex prev = pvertices.front();
    const PVertex next = pvertices.back();
    const PVertex pvertex = pvertices[1];

    if (m_debug) {
      const auto iedge = m_data.iedge(pvertex);
      if (iedge != m_data.null_iedge()) {
        std::cout << "- start from: " << m_data.str(iedge) << " " <<
        m_data.segment_3(iedge) << std::endl;
      } else {
        std::cout << "- start from: unconstrained setting" << std::endl;
      }
    }

    // Copy front/back to remember position/direction.
    PVertex front, back;
    if (pvertices.size() < 3) {
      CGAL_assertion_msg(false, "ERROR: INVALID CONNECTIVITY CASE!");
    } else if (pvertices.size() == 3 || pvertices.size() == 4) {
      std::tie(front, back) = m_data.front_and_back_34(pvertex);
    } else if (pvertices.size() >= 5) {
      std::tie(front, back) = m_data.front_and_back_5(
        pvertices[1], pvertices[pvertices.size() - 2]);
    } else {
      CGAL_assertion_msg(false, "ERROR: INVALID CONNECTIVITY CASE!");
    }

    if (m_debug) {
      std::cout << "- found neighbors: " << std::endl <<
      "prev = " << m_data.point_3(prev)  << std::endl <<
      "fron = " << m_data.point_3(front) << std::endl <<
      "back = " << m_data.point_3(back)  << std::endl <<
      "next = " << m_data.point_3(next)  << std::endl;
    }

    // Should we use here event_pvertex or pvertex?
    // If we use pvertex, we miss important iedges!
    std::vector<IEdge> fiedges, biedges;
    m_data.get_iedges_front_back(event_pvertex, pvertices, fiedges, biedges);
    std::pair<Point_2, Vector_2> query_pvertex = std::make_pair(
      m_data.point_2(event_pvertex, FT(0)), m_data.direction(event_pvertex));

    // Freeze pvertices.
    const Point_2 ipoint = m_data.point_2(sp_idx, ivertex);
    for (std::size_t i = 1; i < pvertices.size() - 1; ++i) {
      const PVertex& curr = pvertices[i];
      m_data.support_plane(curr).direction(curr.second) = CGAL::NULL_VECTOR;
      m_data.support_plane(curr).set_point(curr.second, ipoint);
    }
    m_data.connect(pvertex, ivertex);
    if (m_debug) {
      std::cout << "- frozen pvertex: " << m_data.str(pvertex) << " : " << m_data.point_3(pvertex) << std::endl;
    }

    // Join pvertices.
    for (std::size_t i = 2; i < pvertices.size() - 1; ++i) {
      const auto he = m_data.mesh(sp_idx).halfedge(pvertices[i].second, pvertex.second);
      m_data.disconnect_ivertex(pvertices[i]);
      CGAL::Euler::join_vertex(he, m_data.mesh(sp_idx));
    }

    // Get all connected iedges.
    std::vector< std::pair<IEdge, Direction_2> > iedges;
    m_data.get_and_sort_all_connected_iedges(sp_idx, ivertex, iedges);
    CGAL_assertion(iedges.size() > 0);

    // Get sub-event type.
    bool back_constrained = false;
    if (
      (m_data.iedge(next) != m_data.null_iedge() &&
      (m_data.source(m_data.iedge(next)) == ivertex || m_data.target(m_data.iedge(next)) == ivertex)) ||
      (m_data.ivertex(next) != m_data.null_ivertex() && m_data.is_iedge(m_data.ivertex(next), ivertex))) {
      back_constrained = true;
    }

    bool front_constrained = false;
    if (
      (m_data.iedge(prev) != m_data.null_iedge() &&
      (m_data.source(m_data.iedge(prev)) == ivertex || m_data.target(m_data.iedge(prev)) == ivertex)) ||
      (m_data.ivertex(prev) != m_data.null_ivertex() && m_data.is_iedge(m_data.ivertex(prev), ivertex))) {
      front_constrained = true;
    }

    if (back_constrained && !front_constrained) {
      if (m_debug) std::cout << "- reverse iedges" << std::endl;
      std::reverse(iedges.begin(), iedges.end());
    }

    if (m_debug) {
      std::cout << "- initial iedges: " << iedges.size() << std::endl;
      for (const auto& iedge : iedges) {
        std::cout << m_data.str(iedge.first) << ": " <<
        m_data.segment_3(iedge.first) << std::endl;
      }
    }

    // Handle sub-events.
    crossed_iedges.clear();
    std::vector<PVertex> new_pvertices;

    if (back_constrained && front_constrained) {
      apply_closing_case(pvertex);
    } else if (back_constrained) {
      apply_back_border_case(
        min_time, max_time, query_pvertex,
        pvertex, ivertex, back, prev, fiedges,
        iedges, crossed_iedges, new_pvertices);
    } else if (front_constrained) {
      apply_front_border_case(
        min_time, max_time,
        pvertex, ivertex, front, next, biedges,
        iedges, crossed_iedges, new_pvertices);
    } else {
      apply_open_case(
        min_time, max_time,
        pvertex, ivertex, front, back, prev, next,
        fiedges, biedges, iedges, crossed_iedges, new_pvertices);
    }

    m_data.support_plane(sp_idx).remove_vertex(front.second);
    m_data.support_plane(sp_idx).remove_vertex(back.second);

    // Push also the remaining pvertex so that its events are recomputed.
    new_pvertices.push_back(pvertex);
    if (m_data.iedge(pvertex) != m_data.null_iedge()) {
      crossed_iedges.push_back(std::make_pair(m_data.iedge(pvertex), true));
    }

    if (m_debug) {
      std::size_t num_new_pvertices = 0;
      for (const auto& new_pvertex : new_pvertices) {
        if (new_pvertex != m_data.null_pvertex()) ++num_new_pvertices;
      }
      std::cout << "- number of new pvertices: "  << num_new_pvertices     << std::endl;
      std::cout << "- number of crossed iedges: " << crossed_iedges.size() << std::endl;
    }

    // CGAL_assertion_msg(false, "TODO: MERGE PVERTICES ON IVERTEX!");
    return new_pvertices;
  }

  void apply_closing_case(const PVertex& pvertex) const {

    if (m_debug) {
      std::cout.precision(20);
      std::cout << "*** CLOSING CASE" << std::endl;
    }
    CGAL_assertion(m_data.has_complete_graph(pvertex));

    // CGAL_assertion_msg(false, "TODO: CLOSING CASE!");
  }

  void apply_back_border_case(
    const FT min_time, const FT max_time,
    const std::pair<Point_2, Vector_2>& event_pvertex,
    const PVertex& pvertex, const IVertex& ivertex,
    const PVertex& back, const PVertex& prev,
    const std::vector<IEdge>& fiedges,
    const std::vector< std::pair<IEdge, Direction_2> >& iedges,
    std::vector< std::pair<IEdge, bool> >& crossed_iedges,
    std::vector<PVertex>& new_pvertices) {

    if (m_debug) {
      std::cout.precision(20);
      std::cout << "*** BACK BORDER CASE" << std::endl;
    }

    // We use this modification in order to avoid collinear directions.
    const FT tol = KSR::tolerance<FT>();
    CGAL_assertion(m_data.has_iedge(pvertex));
    const std::size_t other_side_limit = m_data.line_idx(pvertex);
    const FT prev_time = m_data.last_event_time(prev);
    const FT curr_time = m_data.current_time();
    CGAL_assertion(prev_time >= FT(0));
    CGAL_assertion(curr_time >= FT(0));

    // std::cout << "minn time: " <<  min_time << std::endl;
    // std::cout << "curr time: " << curr_time << std::endl;
    // std::cout << "maxx time: " <<  max_time << std::endl;
    // std::cout << "lrev time: " << prev_time << std::endl;

    const FT prev_diff = CGAL::abs(curr_time - prev_time);
    // CGAL_assertion(prev_diff >= tol);
    // if (prev_diff < tol) {
    //   std::cout << "TODO: BACK, EVENTS ARE HAPPENNING AT THE SAME TIME!" << std::endl;
    //   exit(EXIT_FAILURE);
    // }

    FT ntime = max_time;
    if (prev_diff < tol) {
      ntime = m_queue.get_next_time(min_time, max_time, curr_time);
      // std::cout << "next time: " << ntime << std::endl;
    }

    Point_2 shifted_prev;
    const auto pp_curr = m_data.point_2(prev, curr_time);
    if (prev_diff < tol) {
      if (m_debug) std::cout << "- back, same time events, prev" << std::endl;
      CGAL_assertion(CGAL::abs(ntime - curr_time) >= tol);
      const auto pp_futr = m_data.point_2(prev, ntime);
      const auto dirp = Vector_2(pp_curr, pp_futr);

      // Should we reverse fiedges to satisfy the order?
      CGAL_assertion_msg(fiedges.size() <= 2,
      "TODO: BACK, CAN WE HAVE MORE THAN 2 FIEDGES?");

      bool found_iedge = false;
      for (const auto& pair : iedges) {
        const auto& iedge = pair.first;
        CGAL_assertion(iedge != m_data.null_iedge());
        // std::cout << "iedge: " << m_data.str(iedge) << ", " << m_data.segment_3(iedge) << std::endl;
        // std::cout << "fiedge: " << (fiedges.size() > 0) << std::endl;
        // std::cout << "fiedge: " << m_data.segment_3(fiedges.back()) << std::endl;
        if (fiedges.size() > 0 && iedge == fiedges.back()) {
          if (m_debug) std::cout << "- found same time iedge, prev" << std::endl;
          found_iedge = true; break;
        }
      }

      if (found_iedge) {
        shifted_prev = pp_curr + dirp / FT(2);
        if (m_debug) std::cout << "- excluding iedge, prev" << std::endl;
        // CGAL_assertion_msg(false, "TODO: CHECK BACK PREV CASE 1!");
      } else {
        shifted_prev = pp_curr - dirp / FT(2);
        if (m_debug) std::cout << "- including iedge, prev" << std::endl;
        // CGAL_assertion_msg(false, "TODO: CHECK BACK PREV CASE 2!");
      }
    } else {
      const auto pp_last = m_data.point_2(prev, prev_time);
      const auto dirp = Vector_2(pp_last, pp_curr);
      shifted_prev = pp_curr - dirp / FT(10);
      if (m_debug) std::cout << "- including iedge, prev" << std::endl;
    }

    if (m_debug) {
      std::cout << "- shifting prev: " << m_data.to_3d(pvertex.first, shifted_prev) << std::endl;
    }

    const auto ipoint = m_data.point_2(pvertex.first, ivertex);
    const Direction_2 ref_direction_prev(shifted_prev - ipoint);

    // Find the first iedge.
    std::size_t first_idx = std::size_t(-1);
    const std::size_t n = iedges.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;

      const auto& i_dir  = iedges[i].second;
      const auto& ip_dir = iedges[ip].second;
      CGAL_assertion(iedges[i].first != iedges[ip].first);
      if (ref_direction_prev.counterclockwise_in_between(ip_dir, i_dir)) {
        first_idx = ip; break;
      }
    }
    CGAL_assertion(first_idx != std::size_t(-1));
    // std::cout << "- curr: " << m_data.segment_3(iedges[first_idx].first) << std::endl;

    // Find all crossed iedges.
    crossed_iedges.clear();
    CGAL_assertion(crossed_iedges.size() == 0);
    std::size_t iedge_idx = first_idx;
    std::size_t iteration = 0;
    while (true) {
      const auto& iedge = iedges[iedge_idx].first;
      // std::cout << "- next: " << m_data.segment_3(iedge) << std::endl;

      const bool is_bbox_reached  = ( m_data.collision_occured(pvertex, iedge)   ).second;
      const bool is_limit_reached = ( m_data.line_idx(iedge) == other_side_limit );
      if (m_debug) {
        std::cout << "- bbox: " << is_bbox_reached << "; limit: " << is_limit_reached << std::endl;
      }

      crossed_iedges.push_back(std::make_pair(iedge, false));
      if (is_bbox_reached || is_limit_reached) {
        break;
      }

      iedge_idx = (iedge_idx + 1) % n;
      if (iteration >= iedges.size()) {
        CGAL_assertion_msg(false, "ERROR: BACK, WHY SO MANY ITERATIONS?");
      } ++iteration;
    }

    CGAL_assertion(crossed_iedges.size() > 0);
    if (m_debug) {
      std::cout << "- crossed " << crossed_iedges.size() << " iedges: " << std::endl;
      for (const auto& crossed_iedge : crossed_iedges) {
        std::cout << m_data.str(crossed_iedge.first) << ": " <<
        m_data.segment_3(crossed_iedge.first) << std::endl;
      }
    }

    // Compute future points and directions.
    Point_2 future_point; Vector_2 future_direction;
    IEdge prev_iedge = m_data.null_iedge();
    const auto iedge_0 = crossed_iedges[0].first;
    CGAL_assertion_msg(KSR::distance(
      m_data.point_2(pvertex.first, m_data.source(iedge_0)),
      m_data.point_2(pvertex.first, m_data.target(iedge_0))) >= KSR::point_tolerance<FT>(),
    "TODO: BACK, HANDLE ZERO-LENGTH IEDGE!");

    { // future point and direction
      bool is_parallel = false;
      if (KSR::distance(m_data.point_2(back), m_data.point_2(prev)) < KSR::point_tolerance<FT>()) {
        // is_parallel = m_data.compute_future_point_and_direction(
        //   0, back, prev, iedge_0, future_point, future_direction); // does not work!
        if (m_debug) std::cout << "- back = prev, equal points case" << std::endl;
        is_parallel = m_data.compute_future_point_and_direction(
          0, ivertex, event_pvertex, prev, iedge_0, future_point, future_direction);
        // CGAL_assertion_msg(false, "TODO: BACK, FIX CASE WITH EQUAL BACK AND PREV!");
      } else {
        if (m_debug) std::cout << "- back, prev, not equal points case" << std::endl;
        is_parallel = m_data.compute_future_point_and_direction(
          0, ivertex, back, prev, iedge_0, future_point, future_direction);
      }
      if (is_parallel) {
        if (m_data.is_intersecting_iedge(min_time, max_time, prev, iedge_0)) {
          prev_iedge = iedge_0;
        }
      }
    }

    // Crop the pvertex.
    new_pvertices.clear();
    new_pvertices.resize(crossed_iedges.size(), m_data.null_pvertex());

    { // crop
      PVertex cropped = m_data.null_pvertex();
      if (prev_iedge == iedge_0) {
        if (m_debug) std::cout << "- back, prev, parallel case" << std::endl;

        // In case, we are parallel, we update the future point and direction.
        cropped = prev;
        const auto pprev = ( m_data.border_prev_and_next(prev) ).first;
        m_data.compute_future_point_and_direction(
          0, ivertex, prev, pprev, prev_iedge, future_point, future_direction);

      } else {
        if (m_debug) std::cout << "- back, prev, standard case" << std::endl;
        cropped = PVertex(pvertex.first, m_data.support_plane(pvertex).split_edge(pvertex.second, prev.second));
      }

      CGAL_assertion(cropped != m_data.null_pvertex());
      CGAL_assertion(cropped.first == pvertex.first);
      CGAL_assertion(cropped != pvertex);

      const PEdge pedge(pvertex.first, m_data.support_plane(pvertex).edge(pvertex.second, cropped.second));
      new_pvertices[0] = cropped;

      m_data.connect(pedge, iedge_0);
      m_data.connect(cropped, iedge_0);

      CGAL_assertion(future_direction != Vector_2());
      m_data.support_plane(cropped).set_point(cropped.second, future_point);
      m_data.direction(cropped) = future_direction;
      if (m_debug) std::cout << "- cropped: " <<
        m_data.str(cropped) << ", " << m_data.point_3(cropped) << std::endl;
      CGAL_assertion(m_data.is_correctly_oriented(
        cropped.first, future_direction, ivertex, iedge_0));
    }

    // Create new pfaces if any.
    m_data.add_pfaces(
      pvertex, ivertex, back, prev, false, true, true,
      crossed_iedges, new_pvertices);

    // CGAL_assertion_msg(false, "TODO: BACK BORDER CASE!");
  }

  void apply_front_border_case(
    const FT min_time, const FT max_time,
    const PVertex& pvertex, const IVertex& ivertex,
    const PVertex& front, const PVertex& next,
    const std::vector<IEdge>& biedges,
    const std::vector< std::pair<IEdge, Direction_2> >& iedges,
    std::vector< std::pair<IEdge, bool> >& crossed_iedges,
    std::vector<PVertex>& new_pvertices) {

    if (m_debug) {
      std::cout.precision(20);
      std::cout << "*** FRONT BORDER CASE" << std::endl;
    }

    // We use this modification in order to avoid collinear directions.
    const FT tol = KSR::tolerance<FT>();
    CGAL_assertion(m_data.has_iedge(pvertex));
    const std::size_t other_side_limit = m_data.line_idx(pvertex);
    const FT next_time = m_data.last_event_time(next);
    const FT curr_time = m_data.current_time();
    CGAL_assertion(next_time >= FT(0));
    CGAL_assertion(curr_time >= FT(0));

    // std::cout << "minn time: " <<  min_time << std::endl;
    // std::cout << "curr time: " << curr_time << std::endl;
    // std::cout << "maxx time: " <<  max_time << std::endl;
    // std::cout << "lext time: " << next_time << std::endl;

    const FT next_diff = CGAL::abs(curr_time - next_time);
    // CGAL_assertion(next_diff >= tol);
    // if (next_diff < tol) {
    //   std::cout << "TODO: FRONT, EVENTS ARE HAPPENNING AT THE SAME TIME!" << std::endl;
    //   exit(EXIT_FAILURE);
    // }

    FT ntime = max_time;
    if (next_diff < tol) {
      ntime = m_queue.get_next_time(min_time, max_time, curr_time);
      // std::cout << "next time: " << ntime << std::endl;
    }

    Point_2 shifted_next;
    const auto pn_curr = m_data.point_2(next, curr_time);
    if (next_diff < tol) {
      if (m_debug) std::cout << "- front, same time events, next" << std::endl;
      CGAL_assertion(CGAL::abs(ntime - curr_time) >= tol);
      const auto pn_futr = m_data.point_2(next, ntime);
      const auto dirn = Vector_2(pn_curr, pn_futr);

      CGAL_assertion_msg(biedges.size() <= 2,
      "TODO: FRONT, CAN WE HAVE MORE THAN 2 BIEDGES?");

      bool found_iedge = false;
      for (const auto& pair : iedges) {
        const auto& iedge = pair.first;
        CGAL_assertion(iedge != m_data.null_iedge());
        // std::cout << "iedge: " << m_data.str(iedge) << ", " << m_data.segment_3(iedge) << std::endl;
        // std::cout << "biedge: " << (biedges.size() > 0) << std::endl;
        // std::cout << "biedge: " << m_data.segment_3(biedges.front()) << std::endl;
        if (biedges.size() > 0 && iedge == biedges.front()) {
          if (m_debug) std::cout << "- found same time iedge, next" << std::endl;
          found_iedge = true; break;
        }
      }

      if (found_iedge) {
        shifted_next = pn_curr + dirn / FT(2);
        if (m_debug) std::cout << "- excluding iedge, next" << std::endl;
        // CGAL_assertion_msg(false, "TODO: CHECK FRONT NEXT CASE 1!");
      } else {
        shifted_next = pn_curr - dirn / FT(2);
        if (m_debug) std::cout << "- including iedge, next" << std::endl;
        // CGAL_assertion_msg(false, "TODO: CHECK FRONT NEXT CASE 2!");
      }
    } else {
      const auto pn_last = m_data.point_2(next, next_time);
      const auto dirn = Vector_2(pn_last, pn_curr);
      shifted_next = pn_curr - dirn / FT(10);
      if (m_debug) std::cout << "- including iedge, next" << std::endl;
    }

    if (m_debug) {
      std::cout << "- shifting next: " << m_data.to_3d(pvertex.first, shifted_next) << std::endl;
    }

    const auto ipoint = m_data.point_2(pvertex.first, ivertex);
    const Direction_2 ref_direction_next(shifted_next - ipoint);

    // Find the first iedge.
    std::size_t first_idx = std::size_t(-1);
    const std::size_t n = iedges.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;

      const auto& i_dir  = iedges[i].second;
      const auto& ip_dir = iedges[ip].second;
      CGAL_assertion(iedges[i].first != iedges[ip].first);
      if (ref_direction_next.counterclockwise_in_between(i_dir, ip_dir)) {
        first_idx = ip; break;
      }
    }
    CGAL_assertion(first_idx != std::size_t(-1));
    // std::cout << "- curr: " << m_data.segment_3(iedges[first_idx].first) << std::endl;

    // Find all crossed iedges.
    crossed_iedges.clear();
    CGAL_assertion(crossed_iedges.size() == 0);
    std::size_t iedge_idx = first_idx;
    std::size_t iteration = 0;
    while (true) {
      const auto& iedge = iedges[iedge_idx].first;
      // std::cout << "- next: " << m_data.segment_3(iedge) << std::endl;

      const bool is_bbox_reached  = ( m_data.collision_occured(pvertex, iedge)   ).second;
      const bool is_limit_reached = ( m_data.line_idx(iedge) == other_side_limit );
      if (m_debug) {
        std::cout << "- bbox: " << is_bbox_reached << "; limit: " << is_limit_reached << std::endl;
      }

      crossed_iedges.push_back(std::make_pair(iedge, false));
      if (is_bbox_reached || is_limit_reached) {
        break;
      }

      iedge_idx = (iedge_idx + 1) % n;
      if (iteration >= iedges.size()) {
        CGAL_assertion_msg(false, "ERROR: FRONT, WHY SO MANY ITERATIONS?");
      } ++iteration;
    }

    CGAL_assertion(crossed_iedges.size() > 0);
    if (m_debug) {
      std::cout << "- crossed " << crossed_iedges.size() << " iedges: " << std::endl;
      for (const auto& crossed_iedge : crossed_iedges) {
        std::cout << m_data.str(crossed_iedge.first) << ": " <<
        m_data.segment_3(crossed_iedge.first) << std::endl;
      }
    }

    // Compute future points and directions.
    Point_2 future_point; Vector_2 future_direction;
    IEdge next_iedge = m_data.null_iedge();
    const auto iedge_0 = crossed_iedges[0].first;
    CGAL_assertion_msg(KSR::distance(
      m_data.point_2(pvertex.first, m_data.source(iedge_0)),
      m_data.point_2(pvertex.first, m_data.target(iedge_0))) >= KSR::point_tolerance<FT>(),
    "TODO: FRONT, HANDLE ZERO-LENGTH IEDGE!");

    { // future point and direction
    bool is_parallel = false;
      if (KSR::distance(m_data.point_2(front), m_data.point_2(next)) < KSR::point_tolerance<FT>()) {
        if (m_debug) std::cout << "- front = next, equal points case" << std::endl;
        CGAL_assertion_msg(false,
        "TODO: FRONT, FIX CASE WITH EQUAL FRONT AND NEXT! SEE BACK CASE FOR REFERENCE!");
      } else {
        if (m_debug) std::cout << "- front, next, not equal points case" << std::endl;
        is_parallel = m_data.compute_future_point_and_direction(
          0, ivertex, front, next, iedge_0, future_point, future_direction);
      }
      if (is_parallel) {
        if (m_data.is_intersecting_iedge(min_time, max_time, next, iedge_0)) {
          next_iedge = iedge_0;
        }
      }
    }

    // Crop the pvertex.
    new_pvertices.clear();
    new_pvertices.resize(crossed_iedges.size(), m_data.null_pvertex());

    { // crop
      PVertex cropped = m_data.null_pvertex();
      if (next_iedge == iedge_0) {
        if (m_debug) std::cout << "- front, next, parallel case" << std::endl;

        // In case, we are parallel, we update the future point and direction.
        cropped = next;
        const auto nnext = ( m_data.border_prev_and_next(next) ).second;
        m_data.compute_future_point_and_direction(
          0, ivertex, next, nnext, next_iedge, future_point, future_direction);

      } else {
        if (m_debug) std::cout << "- front, next, standard case" << std::endl;
        cropped = PVertex(pvertex.first, m_data.support_plane(pvertex).split_edge(pvertex.second, next.second));
      }

      CGAL_assertion(cropped != m_data.null_pvertex());
      CGAL_assertion(cropped.first == pvertex.first);
      CGAL_assertion(cropped != pvertex);

      const PEdge pedge(pvertex.first, m_data.support_plane(pvertex).edge(pvertex.second, cropped.second));
      new_pvertices[0] = cropped;

      m_data.connect(pedge, iedge_0);
      m_data.connect(cropped, iedge_0);

      CGAL_assertion(future_direction != Vector_2());
      m_data.support_plane(cropped).set_point(cropped.second, future_point);
      m_data.direction(cropped) = future_direction;
      if (m_debug) std::cout << "- cropped: " <<
        m_data.str(cropped) << ", " << m_data.point_3(cropped) << std::endl;
      CGAL_assertion(m_data.is_correctly_oriented(
        cropped.first, future_direction, ivertex, iedge_0));
    }

    // Create new pfaces if any.
    m_data.add_pfaces(
      pvertex, ivertex, front, next, false, false, true,
      crossed_iedges, new_pvertices);

    // CGAL_assertion_msg(false, "TODO: FRONT BORDER CASE!");
  }

  void apply_open_case(
    const FT min_time, const FT max_time,
    const PVertex& pvertex, const IVertex& ivertex,
    const PVertex& /* front */, const PVertex& /* back */,
    const PVertex& prev , const PVertex& next,
    const std::vector<IEdge>& fiedges,
    const std::vector<IEdge>& biedges,
    const std::vector< std::pair<IEdge, Direction_2> >& iedges,
    std::vector< std::pair<IEdge, bool> >& crossed_iedges,
    std::vector<PVertex>& new_pvertices) {

    if (m_debug) {
      std::cout.precision(20);
      std::cout << "*** OPEN CASE" << std::endl;
    }

    // We use this modification in order to avoid collinear directions.
    const FT prev_time = m_data.last_event_time(prev);
    const FT curr_time = m_data.current_time();
    const FT next_time = m_data.last_event_time(next);

    // std::cout << "minn time: " <<  min_time << std::endl;
    // std::cout << "curr time: " << curr_time << std::endl;
    // std::cout << "maxx time: " <<  max_time << std::endl;

    // std::cout << "lrev time: " << prev_time << std::endl;
    // std::cout << "lext time: " << next_time << std::endl;

    const FT tol = KSR::tolerance<FT>();
    CGAL_assertion(prev_time >= FT(0));
    CGAL_assertion(curr_time >= FT(0));
    CGAL_assertion(next_time >= FT(0));
    const FT prev_diff = CGAL::abs(curr_time - prev_time);
    const FT next_diff = CGAL::abs(curr_time - next_time);

    // CGAL_assertion(prev_diff >= tol);
    // CGAL_assertion(next_diff >= tol);
    // if (prev_diff < tol || next_diff < tol) {
    //   std::cout << "TODO: OPEN, EVENTS ARE HAPPENNING AT THE SAME TIME!" << std::endl;
    //   exit(EXIT_FAILURE);
    // }

    FT ntime = max_time;
    if (prev_diff < tol || next_diff < tol) {
      ntime = m_queue.get_next_time(min_time, max_time, curr_time);
      // std::cout << "next time: " << ntime << std::endl;
    }

    Point_2 shifted_prev;
    const auto pp_curr = m_data.point_2(prev, curr_time);
    if (prev_diff < tol) {
      if (m_debug) std::cout << "- open, same time events, prev" << std::endl;
      CGAL_assertion(CGAL::abs(ntime - curr_time) >= tol);
      const auto pp_futr = m_data.point_2(prev, ntime);
      const auto dirp = Vector_2(pp_curr, pp_futr);

      CGAL_assertion_msg(fiedges.size() <= 2,
      "TODO: OPEN PREV, CAN WE HAVE MORE THAN 2 FIEDGES?");

      bool found_iedge = false;
      for (const auto& pair : iedges) {
        const auto& iedge = pair.first;
        CGAL_assertion(iedge != m_data.null_iedge());
        // std::cout << "iedge: " << m_data.str(iedge) << ", " << m_data.segment_3(iedge) << std::endl;
        // std::cout << "fiedge: " << (fiedges.size() > 0) << std::endl;
        // std::cout << "fiedge: " << m_data.segment_3(fiedges.back()) << std::endl;
        if (fiedges.size() > 0 && iedge == fiedges.back()) {
          if (m_debug) std::cout << "- found same time iedge, prev" << std::endl;
          found_iedge = true; break;
        }
      }

      if (found_iedge) {
        shifted_prev = pp_curr + dirp / FT(2);
        if (m_debug) std::cout << "- excluding iedge, prev" << std::endl;
        // CGAL_assertion_msg(false, "TODO: CHECK OPEN PREV CASE 1!");
      } else {
        shifted_prev = pp_curr - dirp / FT(2);
        if (m_debug) std::cout << "- including iedge, prev" << std::endl;
        CGAL_assertion_msg(false, "TODO: CHECK OPEN PREV CASE 2!");
      }
    } else {
      const auto pp_last = m_data.point_2(prev, prev_time);
      const auto dirp = Vector_2(pp_last, pp_curr);
      shifted_prev = pp_curr - dirp / FT(10);
      if (m_debug) std::cout << "- including iedge, prev" << std::endl;
    }

    Point_2 shifted_next;
    const auto pn_curr = m_data.point_2(next, curr_time);
    if (next_diff < tol) {
      if (m_debug) std::cout << "- open, same time events, next" << std::endl;
      CGAL_assertion(CGAL::abs(ntime - curr_time) >= tol);
      const auto pn_futr = m_data.point_2(next, ntime);
      const auto dirn = Vector_2(pn_curr, pn_futr);

      CGAL_assertion_msg(biedges.size() <= 2,
      "TODO: OPEN NEXT, CAN WE HAVE MORE THAN 2 BIEDGES?");

      bool found_iedge = false;
      for (const auto& pair : iedges) {
        const auto& iedge = pair.first;
        CGAL_assertion(iedge != m_data.null_iedge());
        // std::cout << "iedge: " << m_data.str(iedge) << ", " << m_data.segment_3(iedge) << std::endl;
        // std::cout << "biedge: " << (biedges.size() > 0) << std::endl;
        // std::cout << "biedge: " << m_data.segment_3(biedges.front()) << std::endl;
        if (biedges.size() > 0 && iedge == biedges.front()) {
          if (m_debug) std::cout << "- found same time iedge, next" << std::endl;
          found_iedge = true; break;
        }
      }

      if (found_iedge) {
        shifted_next = pn_curr + dirn / FT(2);
        if (m_debug) std::cout << "- excluding iedge, next" << std::endl;
        // CGAL_assertion_msg(false, "TODO: CHECK OPEN NEXT CASE 1!");
      } else {
        shifted_next = pn_curr - dirn / FT(2);
        if (m_debug) std::cout << "- including iedge, next" << std::endl;
        CGAL_assertion_msg(false, "TODO: CHECK OPEN NEXT CASE 2!");
      }
    } else {
      const auto pn_last = m_data.point_2(next, next_time);
      const auto dirn = Vector_2(pn_last, pn_curr);
      shifted_next = pn_curr - dirn / FT(10);
      if (m_debug) std::cout << "- including iedge, next" << std::endl;
    }

    if (m_debug) {
      std::cout << "- shifting prev: " << m_data.to_3d(pvertex.first, shifted_prev) << std::endl;
      std::cout << "- shifting next: " << m_data.to_3d(pvertex.first, shifted_next) << std::endl;
    }

    const auto ipoint = m_data.point_2(pvertex.first, ivertex);
    const Direction_2 ref_direction_prev(shifted_prev - ipoint);
    const Direction_2 ref_direction_next(shifted_next - ipoint);

    // Find the first iedge.
    std::size_t first_idx = std::size_t(-1);
    const std::size_t n = iedges.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;

      const auto& i_dir  = iedges[i].second;
      const auto& ip_dir = iedges[ip].second;
      CGAL_assertion(iedges[i].first != iedges[ip].first);
      if (ref_direction_next.counterclockwise_in_between(i_dir, ip_dir)) {
        first_idx = ip; break;
      }
    }
    CGAL_assertion(first_idx != std::size_t(-1));
    // std::cout << "- curr: " << m_data.segment_3(iedges[first_idx].first) << std::endl;

    // Find all crossed iedges.
    crossed_iedges.clear();
    CGAL_assertion(crossed_iedges.size() == 0);
    std::size_t iedge_idx = first_idx;
    std::size_t iteration = 0;
    while (true) {
      const auto& iedge = iedges[iedge_idx].first;
      // std::cout << "- next: " << m_data.segment_3(iedge) << std::endl;

      if (iteration == iedges.size()) {
        CGAL_assertion_msg(iedges.size() == 2,
        "ERROR: OPEN, CAN WE HAVE THIS CASE IN THE CONSTRAINED SETTING?");
        break;
      }

      const auto& ref_direction = iedges[iedge_idx].second;
      if (!ref_direction.counterclockwise_in_between(
        ref_direction_next, ref_direction_prev)) {
        break;
      }

      crossed_iedges.push_back(std::make_pair(iedge, false));
      iedge_idx = (iedge_idx + 1) % n;
      if (iteration >= iedges.size()) {
        CGAL_assertion_msg(false, "ERROR: OPEN, WHY SO MANY ITERATIONS?");
      } ++iteration;
    }

    CGAL_assertion(crossed_iedges.size() > 0);
    if (m_debug) {
      std::cout << "- crossed " << crossed_iedges.size() << " iedges: " << std::endl;
      for (const auto& crossed_iedge : crossed_iedges) {
        std::cout << m_data.str(crossed_iedge.first) << ": " <<
        m_data.segment_3(crossed_iedge.first) << std::endl;
      }
    }

    if (crossed_iedges.size() == 1) {

      CGAL_assertion_msg(KSR::distance(
        m_data.point_2(pvertex.first, m_data.source(crossed_iedges[0].first)),
        m_data.point_2(pvertex.first, m_data.target(crossed_iedges[0].first))) >= KSR::point_tolerance<FT>(),
      "TODO: OPEN, 1 EDGE CASE, HANDLE ZERO-LENGTH IEDGE!");

      new_pvertices.clear();
      new_pvertices.resize(crossed_iedges.size(), m_data.null_pvertex());

      if (m_debug) std::cout << "- open, 1 edge case" << std::endl;

      Point_2 future_point;
      Vector_2 future_direction;
      bool is_parallel = false;
      if (KSR::distance(m_data.point_2(prev), m_data.point_2(next)) < KSR::point_tolerance<FT>()) {
        if (m_debug) std::cout << "- prev = next, equal points case" << std::endl;
        CGAL_assertion_msg(false,
        "TODO: OPEN, 1 EDGE CASE, FIX CASE WITH EQUAL PREV AND NEXT! SEE BACK CASE FOR REFERENCE!");
      } else {
        if (m_debug) std::cout << "- prev, next, not equal points case" << std::endl;
        is_parallel = m_data.compute_future_point_and_direction(
          pvertex, ivertex, prev, next,
          crossed_iedges[0].first, future_point, future_direction);
      }

      if (is_parallel) {
        if (m_debug) std::cout << "- parallel case" << std::endl;
        CGAL_assertion_msg(!is_parallel, "TODO: OPEN, 1 EDGE CASE, ADD PARALLEL CASE!");
      } else {
        if (m_debug) std::cout << "- standard case" << std::endl;
      }

      const auto cropped1 = PVertex(pvertex.first, m_data.support_plane(pvertex).split_edge(pvertex.second, next.second));
      const auto cropped2 = PVertex(pvertex.first, m_data.support_plane(pvertex).split_edge(pvertex.second, prev.second));
      m_data.add_pface(std::array<PVertex, 3>{pvertex, cropped1, cropped2});
      const auto he = m_data.mesh(pvertex).halfedge(cropped1.second, cropped2.second);
      const auto ei = m_data.mesh(pvertex).edge(he);
      CGAL::Euler::collapse_edge(ei, m_data.mesh(pvertex));
      const auto cropped = cropped2;

      CGAL_assertion(cropped != m_data.null_pvertex());
      CGAL_assertion(cropped.first == pvertex.first);
      CGAL_assertion(cropped != pvertex);

      const PEdge pedge(pvertex.first, m_data.support_plane(pvertex).edge(pvertex.second, cropped.second));
      new_pvertices[0] = cropped;

      m_data.connect(pedge, crossed_iedges[0].first);
      m_data.connect(cropped, crossed_iedges[0].first);

      CGAL_assertion(future_direction != Vector_2());
      m_data.support_plane(cropped).set_point(cropped.second, future_point);
      m_data.direction(cropped) = future_direction;
      if (m_debug) std::cout << "- cropped: " <<
        m_data.str(cropped) << ", " << m_data.point_3(cropped) << std::endl;
      CGAL_assertion(m_data.is_correctly_oriented(
        cropped.first, future_direction, ivertex, crossed_iedges[0].first));

      // CGAL_assertion_msg(false, "TODO: OPEN, HANDLE 1 EDGE CASE!");
      return;
    }

    // Compute future points and directions.
    CGAL_assertion(crossed_iedges.size() >= 2);
    std::vector<Point_2> future_points(2);
    std::vector<Vector_2> future_directions(2);
    IEdge prev_iedge = m_data.null_iedge();
    IEdge next_iedge = m_data.null_iedge();

    { // first future point and direction
      CGAL_assertion_msg(KSR::distance(
        m_data.point_2(pvertex.first, m_data.source(crossed_iedges.front().first)),
        m_data.point_2(pvertex.first, m_data.target(crossed_iedges.front().first))) >= KSR::point_tolerance<FT>(),
      "TODO: OPEN, FRONT, HANDLE ZERO-LENGTH IEDGE!");

      if (m_debug) std::cout << "- getting future point and direction, front" << std::endl;
      bool is_parallel = false;
      if (KSR::distance(m_data.point_2(prev), m_data.point_2(next)) < KSR::point_tolerance<FT>()) {
        if (m_debug) std::cout << "- prev = next, equal points case" << std::endl;
        CGAL_assertion_msg(false,
        "TODO: OPEN, FRONT, FIX CASE WITH EQUAL PREV AND NEXT! SEE BACK CASE FOR REFERENCE!");
      } else {
        if (m_debug) std::cout << "- prev, next, not equal points case" << std::endl;
        is_parallel = m_data.compute_future_point_and_direction(
          pvertex, ivertex, prev, next,
          crossed_iedges.front().first, future_points.front(), future_directions.front());
      }
      if (is_parallel) {
        if (m_data.is_intersecting_iedge(min_time, max_time, prev, crossed_iedges.front().first)) {
          prev_iedge = crossed_iedges.front().first;
        }
        if (m_data.is_intersecting_iedge(min_time, max_time, next, crossed_iedges.front().first)) {
          next_iedge = crossed_iedges.front().first;
        }
      }
    }

    // second future point and direction
    {
      CGAL_assertion_msg(KSR::distance(
        m_data.point_2(pvertex.first, m_data.source(crossed_iedges.back().first)),
        m_data.point_2(pvertex.first, m_data.target(crossed_iedges.back().first))) >= KSR::point_tolerance<FT>(),
      "TODO: OPEN, BACK, HANDLE ZERO-LENGTH IEDGE!");

      if (m_debug) std::cout << "- getting future point and direction, back" << std::endl;
      bool is_parallel = false;
      if (KSR::distance(m_data.point_2(prev), m_data.point_2(next)) < KSR::point_tolerance<FT>()) {
        if (m_debug) std::cout << "- prev = next, equal points case" << std::endl;
        CGAL_assertion_msg(false,
        "TODO: OPEN, BACK, FIX CASE WITH EQUAL PREV AND NEXT! SEE BACK CASE FOR REFERENCE!");
      } else {
        if (m_debug) std::cout << "- prev, next, not equal points case" << std::endl;
        is_parallel = m_data.compute_future_point_and_direction(
          pvertex, ivertex, prev, next,
          crossed_iedges.back().first, future_points.back(), future_directions.back());
      }
      if (is_parallel) {
        if (m_data.is_intersecting_iedge(min_time, max_time, prev, crossed_iedges.back().first)) {
          prev_iedge = crossed_iedges.back().first;
        }
        if (m_data.is_intersecting_iedge(min_time, max_time, next, crossed_iedges.back().first)) {
          next_iedge = crossed_iedges.back().first;
        }
      }
    }

    // Crop the pvertex.
    new_pvertices.clear();
    new_pvertices.resize(crossed_iedges.size(), m_data.null_pvertex());

    { // first crop
      PVertex cropped = m_data.null_pvertex();
      if (next_iedge == crossed_iedges.front().first) {
        if (m_debug) std::cout << "- open, next, parallel case" << std::endl;

        // In case, we are parallel, we update the future point and direction.
        cropped = next;
        const auto nnext = ( m_data.border_prev_and_next(next) ).second;
        m_data.compute_future_point_and_direction(
          0, ivertex, next, nnext, next_iedge, future_points.front(), future_directions.front());

      } else {
        if (m_debug) std::cout << "- open, next, standard case" << std::endl;
        cropped = PVertex(pvertex.first, m_data.support_plane(pvertex).split_edge(pvertex.second, next.second));
      }

      CGAL_assertion(cropped != m_data.null_pvertex());
      CGAL_assertion(cropped.first == pvertex.first);
      CGAL_assertion(cropped != pvertex);

      const PEdge pedge(pvertex.first, m_data.support_plane(pvertex).edge(pvertex.second, cropped.second));
      new_pvertices.front() = cropped;

      m_data.connect(pedge, crossed_iedges.front().first);
      m_data.connect(cropped, crossed_iedges.front().first);

      CGAL_assertion(future_directions.front() != Vector_2());
      m_data.support_plane(cropped).set_point(cropped.second, future_points.front());
      m_data.direction(cropped) = future_directions.front();
      if (m_debug) std::cout << "- cropped 1: " <<
        m_data.str(cropped) << ", " << m_data.point_3(cropped) << std::endl;
      CGAL_assertion(m_data.is_correctly_oriented(
        cropped.first, future_directions.front(), ivertex, crossed_iedges.front().first));
    }

    { // second crop
      PVertex cropped = m_data.null_pvertex();
      if (prev_iedge == crossed_iedges.back().first) {
        if (m_debug) std::cout << "- open, prev, parallel case" << std::endl;

        // In case, we are parallel, we update the future point and direction.
        cropped = prev;
        const auto pprev = ( m_data.border_prev_and_next(prev) ).first;
        m_data.compute_future_point_and_direction(
          0, ivertex, prev, pprev, prev_iedge, future_points.back(), future_directions.back());

      } else {
        if (m_debug) std::cout << "- open, prev, standard case" << std::endl;
        cropped = PVertex(pvertex.first, m_data.support_plane(pvertex).split_edge(pvertex.second, prev.second));
      }

      CGAL_assertion(cropped != m_data.null_pvertex());
      CGAL_assertion(cropped.first == pvertex.first);
      CGAL_assertion(cropped != pvertex);

      const PEdge pedge(pvertex.first, m_data.support_plane(pvertex).edge(pvertex.second, cropped.second));
      new_pvertices.back() = cropped;

      m_data.connect(pedge, crossed_iedges.back().first);
      m_data.connect(cropped, crossed_iedges.back().first);

      CGAL_assertion(future_directions.back() != Vector_2());
      m_data.support_plane(cropped).set_point(cropped.second, future_points.back());
      m_data.direction(cropped) = future_directions.back();
      if (m_debug) std::cout << "- cropped 2: " <<
        m_data.str(cropped) << ", " << m_data.point_3(cropped) << std::endl;
      CGAL_assertion(m_data.is_correctly_oriented(
        cropped.first, future_directions.back(), ivertex, crossed_iedges.back().first));
    }

    // Create new pfaces if any.
    m_data.add_pfaces(
      pvertex, ivertex, prev, next, true, false, true,
      crossed_iedges, new_pvertices);

    // CGAL_assertion_msg(false, "TODO: OPEN CASE!");
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_PROPAGATION_H
