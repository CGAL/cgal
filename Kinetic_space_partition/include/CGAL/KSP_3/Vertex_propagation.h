// Copyright (c) 2026 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Florent Lafarge, Dmitry Anisimov, Simon Giraudot

#ifndef CGAL_KSP_3_VERTEXPROPAGATION_H
#define CGAL_KSP_3_VERTEXPROPAGATION_H

#include <CGAL/license/Kinetic_space_partition.h>

// Internal includes.
#include <CGAL/KSP/utils.h>
#include <CGAL/KSP/debug.h>
#include <CGAL/KSP/parameters.h>

#include <CGAL/KSP_3/Data_structure.h>

#include <CGAL/number_utils.h>

namespace CGAL {
namespace KSP_3 {
namespace internal {

#ifdef DOXYGEN_RUNNING
#else

template<typename GeomTraits, typename IntersectionKernel>
class Vertex_propagation {

public:
  using Kernel = GeomTraits;
  using Intersection_kernel = IntersectionKernel;

  std::size_t handled_events = 0;
  std::size_t a_events = 0;
  std::size_t b_events = 0;
  std::size_t c1_events = 0;
  std::size_t c2_events = 0;
  std::size_t d1_events = 0;
  std::size_t d2_events = 0;
  std::size_t t_events = 0;

  inline static int sp_filter = -1;

private:
  using FT = typename Kernel::FT;
  using IkFT = typename Intersection_kernel::FT;
  using IkPoint_2 = typename Intersection_kernel::Point_2;
  using IkPoint_3 = typename Intersection_kernel::Point_3;
  using IkVector_2 = typename Intersection_kernel::Vector_2;
  using IkVector_3 = typename Intersection_kernel::Vector_3;
  using IkRay_2 = typename Intersection_kernel::Ray_2;
  using IkLine_2 = typename Intersection_kernel::Line_2;
  using IkLine_3 = typename Intersection_kernel::Line_3;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using Vector_2 = typename Kernel::Vector_2;
  using Segment_2 = typename Kernel::Segment_2;
  using Direction_2 = typename Kernel::Direction_2;
  using Line_2 = typename Kernel::Line_2;

  using Data_structure = CGAL::KSP_3::internal::Data_structure<Kernel, Intersection_kernel>;
  using Support_plane = typename Data_structure::Support_plane;
  using Vertex = typename Support_plane::Vertex;
  using Cached_event = typename Support_plane::Cached_event;

  using IVertex = typename Data_structure::IVertex;
  using IEdge = typename Data_structure::IEdge;
  using IFace = typename Data_structure::IFace;

  using PVertex = typename Data_structure::PVertex;
  using PEdge = typename Data_structure::PEdge;
  using PFace = typename Data_structure::PFace;

  using Bbox_2 = CGAL::Bbox_2;
  using Face_index = typename Data_structure::Face_index;

  using Parameters = KSP::internal::Parameters_3<FT>;

  using Event = typename Data_structure::Support_plane::Vertex_event;

  using From_exact = CGAL::Cartesian_converter<Intersection_kernel, Kernel>;

  template<typename Event>
  struct Event_order {
    bool operator()(const Event& a, const Event& b) {
      return a.time > b.time;
    }
  };

public:
  Vertex_propagation(Data_structure& data, const Parameters& parameters) :
    m_data(data), m_parameters(parameters), vts(data.kinetic_vertices()),
    m_min_time(-FT(1)), m_max_time(-FT(1))
  { }

  std::size_t propagate(std::size_t k) {
    std::size_t num_events = 0;

    for (Vertex &v : vts) {
      if (v.moving)
        v.k = k;
    }

    initialize_queue();

    while (!m_event_queue.empty())
      run(num_events);

    return num_events;
  }

  void clear() {
    m_event_queue.clear();
    m_min_time = -FT(1);
    m_max_time = -FT(1);
  }

private:
  Data_structure& m_data;
  const Parameters& m_parameters;
  std::vector<Vertex> &vts;

  FT m_min_time;
  FT m_max_time;

  std::size_t skipped = 0;
  std::size_t case_a = 0;
  std::size_t case_b = 0;
  std::size_t case_c1 = 0;
  std::size_t case_c2 = 0;
  std::size_t case_d = 0;

  std::size_t counter = 0;

  std::priority_queue<Event, std::vector<Event>, Event_order<Event>> m_event_queue;
  std::map<std::pair<IkFT, int>, std::list<Event> > queue;

  /*******************************
  **       IDENTIFY EVENTS      **
  ********************************/

  void add_next_events(std::size_t v_idx, IkFT time) {
    Vertex &v = vts[v_idx];
    CGAL_assertion(!v.cached_events.empty());
    // Just adding the first event may not be sufficient as several events may happen at the same time.
    std::size_t first = 0;
    std::size_t idx = 0;

    for (auto it = v.cached_events.begin(); it != v.cached_events.end();)
      if (it->t <= time)
        it = v.cached_events.erase(it);
      else ++it;

    CGAL_assertion(!v.cached_events.empty());
    IkFT fmin = v.cached_events.begin()->t;

    for (auto it = ++v.cached_events.begin();it !=v.cached_events.end(); it++, idx++) {
      if (it->t < fmin) {
        first = idx;
        fmin = it->t;
      }
    }

    idx = 0;
    for (auto it = v.cached_events.begin(); it != v.cached_events.end(); it++, idx++) {
      if (idx == first || it->t == fmin) {
        Event ve(*it, v_idx);
        if (v.constraints.empty())
          ve.crossed_edge = it->e;
        else
          ve.destination = it->d;

        m_event_queue.push(std::move(ve));

        v.queued_events++;
      }
    }
  }

  void check_vertices_on_edge(IEdge edge) {
    for (auto& [sp_other, ov] : m_data.igraph().edge(edge).vertices) {
      Vertex& v1 = vts[ov.first];
      Vertex& v2 = vts[ov.second];
      if (!v1.moving && !v2.moving) {
        CGAL_assertion(std::size_t(-1) != v2.ivertex);
        CGAL_assertion(v1.ivertex != std::size_t(-1));
        CGAL_assertion(v1.ivertex != v2.ivertex);
        IEdge e = v1.constraint_edge;
        std::size_t line_idx = m_data.igraph().edge(e).line;
        if (std::find(v2.constraints.begin(), v2.constraints.end(), line_idx) == v2.constraints.end())
          e = v1.other_constraint_edge;
        line_idx = m_data.igraph().edge(e).line;
        CGAL_assertion(std::find(v2.constraints.begin(), v2.constraints.end(), line_idx) != v2.constraints.end());
        CGAL_assertion(line_idx != std::size_t(-1));
        IkFT u1 = m_data.get_u_on_line(v1.ivertex, line_idx);
        IkFT u2 = m_data.get_u_on_line(v2.ivertex, line_idx);
        CGAL_assertion(u1 != u2);
      }
      else if (v1.moving && v2.moving) {
        CGAL_assertion(v1.v != IkVector_2(0, 0));
        CGAL_assertion(v2.v != IkVector_2(0, 0));
        CGAL_assertion(CGAL::is_negative(v1.v * v2.v)); // pointing in different directions
        CGAL_assertion(CGAL::is_zero(v1.v.x() * v2.v.y() - v1.v.y() * v2.v.x())); // collinear
      }
      else {
        if (v1.moving) {
          CGAL_assertion(v1.v != IkVector_2(0, 0));
          if (v1.p0 == v2.p0) {
            CGAL_assertion(v1.t_init == v2.t_init);
            CGAL_assertion(v1.t_init != 0);
          }
          else
            CGAL_assertion(CGAL::is_negative((v2.p0 - v1.p0) * v1.v));
        }
        else {
          CGAL_assertion(v2.moving);
          CGAL_assertion(v2.v != IkVector_2(0, 0));
          if (v1.p0 == v2.p0) {
            CGAL_assertion(v1.t_init == v2.t_init);
            CGAL_assertion(v1.t_init != 0);
          }
          else
            CGAL_assertion(CGAL::is_negative((v1.p0 - v2.p0) * v2.v));
        }
      }
    }
  }

  void check_edge_is_turning(std::size_t v1_idx, std::size_t v2_idx) {
    IkFT t_min = (vts[v1_idx].t_init < vts[v2_idx].t_init) ? vts[v2_idx].t_init : vts[v1_idx].t_init;
    IkPoint_2 v1 = vts[v1_idx].p0 + (t_min - vts[v1_idx].t_init) * vts[v1_idx].v;
    IkPoint_2 v2 = vts[v2_idx].p0 + (t_min - vts[v2_idx].t_init) * vts[v2_idx].v;
    IkPoint_2 v1_t = v1 + vts[v1_idx].v;
    IkPoint_2 v2_t = v2 + vts[v2_idx].v;
    IkVector_2 e1 = v2 - v1;
    IkVector_2 e2 = v2_t - v1_t;
    CGAL_assertion((e1.x * e2.y() - e1.y() * e2.x()) == 0); // Collinear
    CGAL_assertion((e1 * e2) > 0); // not inverted
  }

  template<typename IndexRange>
  void export_face(std::string filename, const IndexRange &poly, IFace face, Support_plane &sp) {
    return;
    const std::vector<IkPoint_2> &pts = m_data.igraph().face(face).pts;
    if (!pts.empty()) {
      std::ofstream out(filename + "_face_" + std::to_string(face) + ".polylines.txt");
      out << (pts.size() + 1);
      From_exact from_exact;
      for (const IkPoint_2& p : pts) {
        Point_3 p3d = from_exact(sp.to_3d(p));
        out << " " << CGAL::to_double(p3d.x()) << " " << CGAL::to_double(p3d.y()) << " " << CGAL::to_double(p3d.z());
      }

      Point_3 p3d = from_exact(sp.to_3d(pts[0]));
      out << " " << CGAL::to_double(p3d.x()) << " " << CGAL::to_double(p3d.y()) << " " << CGAL::to_double(p3d.z());
      out << std::endl;
      out.close();
    }

    if (!poly.empty()) {
      std::ofstream out(filename + "_poly.polylines.txt");
      out << (poly.size() + 1);
      for (const auto& i : poly) {
        IkPoint_3 p3d = sp.to_3d(vts[i].p0);
        out << " " << CGAL::to_double(p3d.x()) << " " << CGAL::to_double(p3d.y()) << " " << CGAL::to_double(p3d.z());
      }
      IkPoint_3 p3d = sp.to_3d(vts[*poly.begin()].p0);
      out << " " << CGAL::to_double(p3d.x()) << " " << CGAL::to_double(p3d.y()) << " " << CGAL::to_double(p3d.z()) << std::endl;
      out.close();

      std::ofstream out2(filename + "_poly_front.polylines.txt");
      out2 << (poly.size() + 1);

      for (const auto& i : poly) {
        IkPoint_3 p3d = sp.to_3d(vts[i].p0 + (0.1 * vts[i].v));
        out2 << " " << CGAL::to_double(p3d.x()) << " " << CGAL::to_double(p3d.y()) << " " << CGAL::to_double(p3d.z());
      }
      p3d = sp.to_3d(vts[*poly.begin()].p0 + (0.1 * vts[*poly.begin()].v));
      out2 << " " << CGAL::to_double(p3d.x()) << " " << CGAL::to_double(p3d.y()) << " " << CGAL::to_double(p3d.z()) << std::endl;

      out2.close();
    }
  }

  IkFT intersection_time(const Vertex& v, const IkLine_2& l, bool& collinear) const {
    IkFT a = l.a(), b = l.b(), c = l.c();

    IkFT t_den = a * v.v.x() + b * v.v.y();
    if (t_den == 0) {
      // The direction of the vector is collinear to the line
      collinear = true;
      return -1;
    }

    collinear = false;

    IkFT t_num = -(a * v.p0.x() + b * v.p0.y() + c);

    return t_num / t_den + v.t_init;
  }

  IkFT intersection_time(const Vertex& v, const IkPoint_2& p) const {
    if (CGAL::abs(v.v.x()) > CGAL::abs(v.v.y()))
      return (p.x() - v.p0.x()) / v.v.x() + v.t_init;
    else
      return (p.y() - v.p0.y()) / v.v.y() + v.t_init;
  }

  void check_vertex_collisions(std::size_t v_idx) {
    std::vector<Vertex> &vts = m_data.kinetic_vertices();
    Vertex &v = vts[v_idx];
    CGAL_assertion(!v.constraints.empty());

    int dim = (CGAL::abs(v.v[0]) < CGAL::abs(v.v[1])) ? 1 : 0;

    for (auto others : m_data.igraph().edge(v.constraint_edge).vertices) {
      if (others.first == v.sp_idx)
        continue;

      CGAL_assertion(others.second.first != std::size_t(-1));
      CGAL_assertion(others.second.second != std::size_t(-1));

      Vertex& v1 = vts[others.second.first];
      Vertex& v2 = vts[others.second.second];

      // optimization possible: provide pair of vertices, instead of checking each alone
      // check below if vertex is in front or behind benefits.

      // Which other vertex has opposing direction to this vertex?
      Support_plane& sp = m_data.support_plane(v.sp_idx);
      Support_plane& sp_other = m_data.support_plane(others.first);

      if (v1.moving && v_idx < others.second.first) {
        IkVector_3 v1v3d = sp_other.exact_plane().base1() * v1.v[0] + sp_other.exact_plane().base2() * v1.v[1];
        IkVector_2 v1v = IkVector_2(v1v3d * sp.exact_plane().base1(), v1v3d * sp.exact_plane().base2());
        IkFT f1 = v.v * v1v;
        double d1 = CGAL::to_double(f1);
        CGAL_assertion_code(
        IkVector_3 v2v3d = sp_other.exact_plane().base1() * v2.v[0] + sp_other.exact_plane().base2() * v2.v[1];
        IkVector_2 v2v = IkVector_2(v2v3d * sp.exact_plane().base1(), v2v3d * sp.exact_plane().base2());
        IkFT f2 = v.v * v2v;
        double d2 = CGAL::to_double(f2););
        if (v.v * v1v < 0) {
          CGAL_assertion(!v2.moving || v.v * v2v > 0); // Only one can have opposing direction
          IkPoint_2 p1 = sp.to_2d(sp_other.to_3d(v1.p0));

          IkFT t_diff = v.t_init - v1.t_init;
          IkFT t;
          if (t_diff <= 0)
            t = (p1[dim] - v.p0[dim] + v.v[dim] * t_diff) / (v.v[dim] - v1v[dim]);
          else
            t = (p1[dim] + v1v[dim] * t_diff - v.p0[dim]) / (v.v[dim] - v1v[dim]);

          if (t < 0)
            continue;

          Event ve;
          ve.vertex = v_idx;
          ve.other = others.second.first;
          ve.time = t;

          m_event_queue.push(std::move(ve));
          v.queued_events++;

          Event ve_other;
          ve_other.vertex = others.second.first;
          ve_other.time = t;
          ve_other.other = v_idx;

          m_event_queue.push(std::move(ve_other));
          v1.queued_events++;
        }
      }

      if (v2.moving && v_idx < others.second.second) {
        IkVector_3 v2v3d = sp_other.exact_plane().base1() * v2.v[0] + sp_other.exact_plane().base2() * v2.v[1];
        IkVector_2 v2v = IkVector_2(v2v3d * sp.exact_plane().base1(), v2v3d * sp.exact_plane().base2());
        IkFT f2 = v.v * v2v;
        double d2 = CGAL::to_double(f2);
        CGAL_assertion_code(
        IkVector_3 v1v3d = sp_other.exact_plane().base1() * v1.v[0] + sp_other.exact_plane().base2() * v1.v[1];
        IkVector_2 v1v = IkVector_2(v1v3d * sp.exact_plane().base1(), v1v3d * sp.exact_plane().base2());
        IkFT f1 = v.v * v1v;
        double d1 = CGAL::to_double(f1););
        if (v.v * v2v < 0) {
          CGAL_assertion(!v1.moving || v.v * v1v > 0); // Only one can have opposing direction
          IkPoint_2 p2 = sp.to_2d(sp_other.to_3d(v2.p0));
          IkFT t;

          IkFT t_diff = v.t_init - v2.t_init;
          if (t_diff <= 0)
            t = (p2[dim] - v.p0[dim] + v.v[dim] * t_diff) / (v.v[dim] - v2v[dim]);
          else
            t = (p2[dim] + v2v[dim] * t_diff - v.p0[dim]) / (v.v[dim] - v2v[dim]);
          if (t < 0)
            continue;

          Event ve;
          ve.vertex = v_idx;
          ve.other = others.second.second;
          ve.time = t;

          m_event_queue.push(std::move(ve));
          v.queued_events++;

          Event ve_other;
          ve_other.vertex = others.second.second;
          ve_other.time = t;
          ve_other.other = v_idx;

          m_event_queue.push(std::move(ve_other));
          v2.queued_events++;
        }
      }
    }
  }

  std::size_t has_adjacent_with_same_target(std::size_t v_idx) {
    const Vertex &v = vts[v_idx];
    const std::list<std::size_t>& poly = m_data.support_plane(v.sp_idx).data().polygons[v.face];

    const std::list<std::size_t>::const_iterator it = std::find(poly.begin(), poly.end(), v_idx);
    std::list<std::size_t>::const_iterator pit = (it == poly.begin()) ? std::prev(poly.end()) : std::prev(it);
    if (vts[*pit].itarget == v.itarget)
      return *pit;

    // check if next has same target
    pit = (std::next(it) == poly.end()) ? poly.begin() : std::next(it);
    if (vts[*pit].itarget == v.itarget)
      return *pit;

    return -1;
  }

  IEdge get_other_edge(IVertex former, IVertex current, std::size_t line_idx) {
    auto& v_on_constraint_line = m_data.igraph().vertices_on_line(line_idx);

    auto it = v_on_constraint_line.begin();
    for (;it != v_on_constraint_line.end();it++)
      if (it->second == current)
        break;

    CGAL_assertion(it != v_on_constraint_line.end());

    if (it != v_on_constraint_line.begin() && std::prev(it)->second != former)
      return m_data.igraph().edge(current, std::prev(it)->second);
    else
      return m_data.igraph().edge(current, std::next(it)->second);
  }

  IEdge get_next_edge(IEdge edge, IVertex current, std::size_t line_idx) {
    IEdge prolongation_edge;
    IVertex former = m_data.igraph().source(edge);
    if (former == current)
      former = m_data.igraph().target(edge);

    return get_other_edge(former, current, line_idx);
  }

  void calculate_events(std::size_t v_idx, IkFT time) {
    Vertex &v = vts[v_idx];
    CGAL_assertion(v.moving);
    Support_plane& sp = m_data.support_plane(v.sp_idx);
    if (v.constraints.empty()) {
      std::list<IEdge>& border = sp.data().borders[v.face];
      for (const IEdge& e : border) {
        std::size_t line_idx = m_data.igraph().edge(e).line;
        // Already calculated?
        if (v.known_intersections.find(line_idx) != v.known_intersections.end())
          continue;

        Cached_event ce;
        bool collinear;
        auto it = sp.data().lines.find(line_idx);
        if (it == sp.data().lines.end()) {
          const auto& vertices = m_data.igraph().vertices_on_line(line_idx);
          it = sp.data().lines.insert(std::make_pair(line_idx, IkLine_2(sp.to_2d(m_data.igraph().point_3(vertices[0].second)), sp.to_2d(m_data.igraph().point_3(vertices.back().second))))).first;
        }

        IkLine_2& l = it->second;
        ce.t = intersection_time(v, l, collinear);
        if (collinear || ce.t < time) {
          v.known_intersections[line_idx].t = -1;
          continue;
        }

        ce.p = v.p0 + (ce.t - v.t_init) * v.v;

        auto l_3d = m_data.igraph().line(line_idx);

        ce.u = (sp.to_3d(ce.p) - l_3d.point()) * l_3d.to_vector();
        ce.e = m_data.igraph().locate_edge_on_line(line_idx, ce.u);
        //ve.face = m_data.igraph().get_other_face(sp_idx, e, face);

        if (ce.e != m_data.igraph().null_iedge()) {
          v.cached_events.push_back(ce);
          v.known_intersections[line_idx].t = ce.t;
        }
        else
          v.known_intersections[line_idx].t = -1;
      }

      add_next_events(v_idx, time);
    }
    else {
      CGAL_assertion(v.constraints.size() == 1);
      if (v.other == std::size_t(-1) || v_idx < v.other) {
        std::size_t adj = has_adjacent_with_same_target(v_idx);

        if ((adj == -1 || vts[adj].cached_events.empty()) && v.other != -1)
          adj = has_adjacent_with_same_target(v.other);

        /// Does it have a constrained neighbor with the same target (= C1 event)
        if (adj != -1 && !vts[adj].cached_events.empty()) {
          auto it = vts[adj].cached_events.begin();
          while (it != vts[adj].cached_events.end()) {
            if (it->d == v.itarget)
              break;
            ++it;
          }
          if (it != vts[adj].cached_events.end()) {
            Event ve;
            ve.vertex = v_idx;
            ve.time = it->t;
            ve.destination = v.itarget;
            ve.p = it->p;

            if (v.other != std::size_t(-1) && vts[v.other].moving) {
              Event ve_other;
              ve_other.vertex = v.other;
              ve_other.destination = ve.destination;
              ve_other.time = ve.time;
              ve_other.p = ve.p;

              m_event_queue.push(std::move(ve_other));
              vts[v.other].queued_events++;
            }

            m_event_queue.push(std::move(ve));
            v.queued_events++;
            return;
          }
        }

        Cached_event ce;
        Event ve;
        ve.vertex = v_idx;
        CGAL_assertion(v.itarget != -1);
        IkPoint_2 target = sp.data().exact_plane.to_2d(m_data.igraph().point_3(v.itarget));
        IkFT t1 = intersection_time(v, target);
        CGAL_assertion(t1 > time);

        ve.time = t1;
        ce.p = target;
        ce.t = t1;
        ce.d = v.itarget;
        v.cached_events.push_back(ce);
        ve.destination = v.itarget;
        ve.p = target;

        if (v.other != std::size_t(-1) && vts[v.other].moving) {
          Event ve_other;
          ve_other.vertex = v.other;
          ve_other.destination = ve.destination;
          ve_other.time = ve.time;
          ve_other.p = ve.p;

          Cached_event ce;
          ce.d = ve.destination;
          ce.t = ve.time;
          ce.p = ve.p;
          vts[v.other].cached_events.push_back(ce);

          m_event_queue.push(std::move(ve_other));
          vts[v.other].queued_events++;
        }

        m_event_queue.push(std::move(ve));
        v.queued_events++;
      }
    }
  }

  void initialize_queue() {
    for (std::size_t sp_idx = 6; sp_idx < m_data.number_of_support_planes(); sp_idx++) {
      Support_plane& sp = m_data.support_plane(sp_idx);
      std::vector<Vertex> &v = sp.data().vertices;
      // for each polygon
      for (auto [face, poly] : sp.data().polygons) {
      //  extract border edges

        for (std::size_t i : poly) {
          if (!v[i].moving)
            continue;

          calculate_events(i, 0);
        }
      }
    }
  }

  /*******************************
  **          RUNNING           **
  ********************************/

  std::size_t run(
    const std::size_t initial_iteration) {

    std::size_t iteration = initial_iteration;
    while (!m_event_queue.empty()) {
      std::list<Event> events;

      do {
        const Event e = m_event_queue.top();
        m_event_queue.pop();
        if (!vts[e.vertex].moving) {
          skipped++;
          continue;
        }
        events.push_back(e);
      } while (!m_event_queue.empty() && (events.empty() || m_event_queue.top().time == events.front().time));

      if (events.empty())
        return iteration;

      ++iteration;

      apply(events);
    }

    if (t_events != 0)
      std::cout << "\n" << t_events << std::endl;

    return iteration;
  }

  void export_timestamp(const std::string &filename, IkFT time) {
    From_exact from_exact;
    for (std::size_t i = 6;i<m_data.number_of_support_planes();i++) {
      std::ofstream out(filename + "_sp_" + std::to_string(i) + ".ply");

      std::vector<Point_3> points;
      std::vector<std::vector<std::size_t>> active, still;
      std::vector<std::size_t> active_face;
      std::vector<std::size_t> still_face;
      Support_plane& sp = m_data.support_plane(i);
      for (auto [face, poly] : sp.data().polygons) {
        bool is_active = false;
        for (std::size_t idx : poly)
          if (vts[idx].moving) {
            is_active = true;
            break;
          }

        if (is_active) {
          active_face.push_back(face);
          active.resize(active.size() + 1);
        }
        else {
          still_face.push_back(face);
          still.resize(still.size() + 1);
        }

        for (std::size_t idx : poly) {
          Point_3 p;
          if (vts[idx].moving)
            p = from_exact(sp.to_3d(vts[idx].p0 + (time - vts[idx].t_init) * vts[idx].v));
          else {
            p = from_exact(sp.to_3d(vts[idx].p0));
            CGAL_assertion(vts[idx].v == IkVector_2(0, 0));
          }
          points.push_back(p);
          if (is_active)
            active.back().push_back(points.size() - 1);
          else
            still.back().push_back(points.size() - 1);
        }
      }

      out << "ply\n"
        << "format ascii 1.0\n"
        << "element vertex " << points.size() << "\n"
        << "property double x\n"
        << "property double y\n"
        << "property double z\n"
        << "element face " << (active.size() + still.size()) << "\n"
        << "property list uchar int vertex_indices\n"
        << "property uchar red\n"
        << "property uchar green\n"
        << "property uchar blue\n"
        << "property ushort id\n"
        << "end_header\n";

      for (const Point_3 &p : points)
        out << p << "\n";

      int idx = 0;
      for (const std::vector<std::size_t>& poly : active) {
        CGAL::Random rand(static_cast<unsigned int>(active_face[idx]));
        out << poly.size();
        for (std::size_t idx : poly)
          out << " " << idx;

        out << " " << rand.get_int(128, 192) << " " << rand.get_int(128, 192) << " " << rand.get_int(128, 192) << " " << active_face[idx] << std::endl;
        idx++;
      }

      idx = 0;
      for (const std::vector<std::size_t>& poly : still) {
        CGAL::Random rand(static_cast<unsigned int>(still_face[idx]));
        out << poly.size();
        for (std::size_t idx : poly)
          out << " " << idx;

        out << " " << rand.get_int(32, 92) << " " << rand.get_int(32, 92) << " " << rand.get_int(32, 92) << " " << still_face[idx] << std::endl;
        idx++;
      }

      out.close();
    }
  }

  void export_event(const std::string &filename, const Event &e, bool in_3d = true) {
    IEdge edge = e.crossed_edge;
    Vertex &v = vts[e.vertex];
    Support_plane &sp = m_data.support_plane(v.sp_idx);

    if (e.crossed_edge == m_data.igraph().null_iedge()) {
      for (IEdge border : m_data.igraph().face(v.face).edges) {
        if (border != v.constraint_edge && (m_data.igraph().target(border) == e.destination || m_data.igraph().source(border) == e.destination)) {
          edge = border;
          break;
        }
      }
      CGAL_assertion(edge != m_data.igraph().null_iedge());
    }

    std::list<std::size_t> &poly = sp.data().polygons[v.face];

    std::size_t line_idx = m_data.igraph().edge(edge).line;
    auto l_it = sp.data().lines.find(line_idx);
    if (l_it == sp.data().lines.end()) {
      const std::vector<std::pair<IkFT, IVertex>> &vol = m_data.igraph().vertices_on_line(line_idx);
      l_it = sp.data().lines.insert(std::make_pair(line_idx, IkLine_2(sp.to_2d(m_data.igraph().point_3(vol[0].second)), sp.to_2d(m_data.igraph().point_3(vol.back().second))))).first;
    }
    IkLine_2 l = l_it->second;
    std::ofstream out(filename + "line_" + std::to_string(line_idx) + ".polylines.txt");
    FT fmin = FLT_MAX, fmax = -FLT_MAX;
    From_exact from_exact;
    Vector_2 dir = from_exact(l.to_vector() * (FT(1.0) / CGAL::approximate_sqrt(l.to_vector().squared_length())));
    Vector_3 dir3d = sp.exact_plane().base1() * dir.x() + sp.exact_plane().base2() * dir.y();

    for (const auto& i : poly) {
      const auto& v = vts[i];
      FT t = dir * from_exact(((v.p0 + (e.time - v.t_init) * v.v) - l.point()));
      fmin = (std::min)(fmin, t);
      fmax = (std::max)(fmax, t);
    }
    FT fspan = fmax - fmin;
    fmin -= fspan * 0.2;
    fmax += fspan * 0.2;

    if (in_3d) {
      Point_3 p1 = sp.to_3d(from_exact(l.point()) + fmin * dir);
      out << 2 << " " << p1.x() << " " << p1.y() << " " << p1.z();
      Point_3 p2 = sp.to_3d(from_exact(l.point()) + fmax * dir);
      out << " " << p2.x() << " " << p2.y() << " " << p2.z();
    }
    else {
      out << 2 << " " << CGAL::to_double(l.point().x() + fmin * dir.x()) << " " << CGAL::to_double(l.point().y() + fmin * dir.y()) << " 0";
      out << " " << CGAL::to_double(l.point().x() + fmax * dir.x()) << " " << CGAL::to_double(l.point().y() + fmax * dir.y()) << " 0";
    }
    out << std::endl;
    out.close();

    {
      std::ofstream out(filename + "_edge_s" + std::to_string(m_data.igraph().source(edge)) + "_t" + std::to_string(m_data.igraph().target(edge))  + ".polylines.txt");
      out << 2;

      if (in_3d) {
        Point_3 s = from_exact(m_data.igraph().point_3(m_data.igraph().source(edge)));
        out << " " << s.x() << " " << s.y() << " " << s.z();
        Point_3 t = from_exact(m_data.igraph().point_3(m_data.igraph().target(edge)));
        out << " " << t.x() << " " << t.y() << " " << t.z() << std::endl;
      }
      else {
        auto s = sp.to_2d(m_data.igraph().point_3(m_data.igraph().source(edge)));
        out << " " << CGAL::to_double(s.x()) << " " << CGAL::to_double(s.y()) << " 0";
        auto t = sp.to_2d(m_data.igraph().point_3(m_data.igraph().target(edge)));
        out << " " << CGAL::to_double(t.x()) << " " << CGAL::to_double(t.y()) << " 0" << std::endl;
      }
      out.close();
    }

    if (e.destination != IVertex(-1))
    {
      std::ofstream out(filename + "_destination.xyz");
      if (in_3d) {
        Point_3 s = from_exact(m_data.igraph().point_3(e.destination));
        out << s.x() << " " << s.y() << " " << s.z() << std::endl;
      }
      else {
        auto s = sp.to_2d(m_data.igraph().point_3(e.destination));
        out << " " << CGAL::to_double(s.x()) << " " << CGAL::to_double(s.y()) << " 0" << std::endl;
      }
      out.close();
    }

    {
      std::ofstream out(filename + "_vertex.polylines.txt");
      Vertex& v = vts[e.vertex];
      out << "2";
      if (in_3d) {
        Point_3 p = from_exact(sp.to_3d(v.p0));
        out << " " << p.x() << " " << p.y() << " " << p.z();
        Point_3 p2 = from_exact(sp.to_3d(v.p0 + ((e.time - v.t_init) * v.v)));
        out << " " << p2.x() << " " << p2.y() << " " << p2.z() << std::endl;
      }
      else {
        out << " " << CGAL::to_double(v.p0.x()) << " " << CGAL::to_double(v.p0.y()) << " 0";
        out << " " << CGAL::to_double(v.p0.x() + ((e.time - v.t_init) * v.v.x())) << " " << CGAL::to_double(v.p0.y() + ((e.time - v.t_init) * v.v.y())) << " 0" << std::endl;
      }
      out.close();
    }

    {
      IkFT f_creation_time = 0;
      for (const auto& i : poly)
        if (vts[i].t_init > f_creation_time)
          f_creation_time = vts[i].t_init;

      std::ofstream out(filename + "_poly.polylines.txt");
      out << (poly.size() + 1);
      for (const auto& i : poly) {
        const auto& v = vts[i];
        if (in_3d) {
          Point_3 p = from_exact(sp.to_3d(v.p0 + (f_creation_time - v.t_init) * v.v));
          out << " " << p.x() << " " << p.y() << " " << p.z();
        }
        else
          out << " " << CGAL::to_double((v.p0 + (f_creation_time - v.t_init) * v.v).x()) << " " << CGAL::to_double((v.p0 + (f_creation_time - v.t_init) * v.v).y()) << " 0";
      }
      if (in_3d) {
        Point_3 p = from_exact(sp.to_3d(vts[poly.front()].p0 + (f_creation_time - vts[poly.front()].t_init) * vts[poly.front()].v));
        out << " " << p.x() << " " << p.y() << " " << p.z() << std::endl;
      }
      else
        out << " " << CGAL::to_double((vts[poly.front()].p0 + (f_creation_time - vts[poly.front()].t_init) * vts[poly.front()].v).x()) << " " << CGAL::to_double((vts[poly.front()].p0 + (f_creation_time - vts[poly.front()].t_init) * vts[poly.front()].v).y()) << " 0" << std::endl;
      out.close();
    }

    {
      std::ofstream out(filename + "_face_" + std::to_string(v.face) + "_vertices.polylines.txt");
      out << (m_data.igraph().face(v.face).vertices.size() + 1);
      for (IVertex vert : m_data.igraph().face(v.face).vertices) {
        if (in_3d) {
          auto s = from_exact(m_data.igraph().point_3(vert));
          out << " " << s.x() << " " << s.y() << " " << s.z();
        }
        else {
          auto s = sp.to_2d(m_data.igraph().point_3(vert));
          out << " " << CGAL::to_double(s.x()) << " " << CGAL::to_double(s.y()) << " 0";
        }
      }
      IVertex vert = m_data.igraph().face(v.face).vertices.front();
      if (in_3d) {
        auto s = from_exact(m_data.igraph().point_3(vert));
        out << " " << s.x() << " " << s.y() << " " << s.z() << std::endl;
      }
      else {
        auto s = sp.to_2d(m_data.igraph().point_3(vert));
        out << " " << CGAL::to_double(s.x()) << " " << CGAL::to_double(s.y()) << " 0" << std::endl;
      }
      out.close();
    }

    {
      std::ofstream out(filename + "_face_" + std::to_string(v.face) + "_border.polylines.txt");
      for (IEdge edge : sp.data().borders[v.face]) {
        if (in_3d) {
          auto s = from_exact(m_data.igraph().point_3(m_data.igraph().source(edge)));
          out << "2 " << s.x() << " " << s.y() << " " << s.z();
          auto t = from_exact(m_data.igraph().point_3(m_data.igraph().target(edge)));
          out << " " << t.x() << " " << t.y() << " " << t.z() << std::endl;
        }
        else {
          auto s = sp.to_2d(m_data.igraph().point_3(m_data.igraph().source(edge)));
          out << "2 " << CGAL::to_double(s.x()) << " " << CGAL::to_double(s.y()) << " 0";
          auto t = sp.to_2d(m_data.igraph().point_3(m_data.igraph().target(edge)));
          out << " " << CGAL::to_double(t.x()) << " " << CGAL::to_double(t.y()) << " 0" << std::endl;
        }
      }
      out.close();
    }

    {
      std::ofstream out(filename + "_poly_front.polylines.txt");
      out << (poly.size() + 1);
      for (const auto& i : poly) {
        const auto& v = vts[i];
        if (in_3d) {
          Point_3 p = from_exact(sp.to_3d(v.p0 + ((e.time - v.t_init) * v.v)));
          out << " " << p.x() << " " << p.y() << " " << p.z();
        }
        else
          out << " " << CGAL::to_double(v.p0.x() + ((e.time - v.t_init) * v.v.x())) << " " << CGAL::to_double(v.p0.y() + ((e.time - v.t_init) * v.v.y())) << " 0";
      }
      const auto& v = vts[poly.front()];
      if (in_3d) {
        Point_3 p = from_exact(sp.to_3d(v.p0 + ((e.time - v.t_init) * v.v)));
        out << " " << p.x() << " " << p.y() << " " << p.z() << std::endl;
      }
      else
        out << " " << CGAL::to_double(v.p0.x() + ((e.time - v.t_init) * v.v.x())) << " " << CGAL::to_double(v.p0.y() + ((e.time - v.t_init) * v.v.y())) << " 0" << std::endl;
    }

    {
      const auto& v = vts[e.vertex];
      out << "2";
      if (in_3d) {
        Point_3 p = from_exact(sp.to_3d(v.p0));
        out << " " << p.x() << " " << p.y() << " " << p.z();
        Point_3 p2 = from_exact(sp.to_3d(v.p0 + ((e.time - v.t_init) * v.v)));
        out << " " << p2.x() << " " << p2.y() << " " << p2.z() << std::endl;
      }
      else {
        out << " " << CGAL::to_double(v.p0.x()) << " " << CGAL::to_double(v.p0.y()) << " 0";
        out << " " << CGAL::to_double(v.p0.x() + ((e.time - v.t_init) * v.v.x())) << " " << CGAL::to_double(v.p0.y() + ((e.time - v.t_init) * v.v.y())) << " 0" << std::endl;
      }
    }

    out.close();
  }

  /*******************************
  **        HANDLE EVENTS       **
  ********************************/
  bool can_cross(const Event& e, const IEdge& edge, bool decrement = true) {
    if (m_data.igraph().iedge_is_on_bbox(edge))
      return false;

    Vertex& v = vts[e.vertex];
    Support_plane& sp = m_data.support_plane(v.sp_idx);

    // Does a polygon already exist on the face on other side?
    IFace other_face = m_data.igraph().get_other_face(v.sp_idx, edge, v.face);
    if (sp.data().polygons.find(other_face) != sp.data().polygons.end())
      return false;

    const IkLine_3& l = m_data.igraph().line(m_data.igraph().edge(edge).line);
    std::vector<std::size_t> crossing;
    for (auto& [sp_other, ov] : m_data.igraph().edge(edge).vertices) {
      if (sp_other == v.sp_idx)
        continue;

      Support_plane &sp = m_data.support_plane(sp_other);

      Vertex& v1 = vts[ov.first];
      Vertex& v2 = vts[ov.second];

      if (e.destination != std::size_t(-1)) {
        if (vts[ov.first].ivertex == e.destination || vts[ov.second].ivertex == e.destination)
          crossing.push_back(sp_other);
      }
      else {
        IkFT u_low = l.to_vector() * (sp.to_3d(v1.p0 + (e.time - v1.t_init) * v1.v) - l.point());
        IkFT u_high = l.to_vector() * (sp.to_3d(v2.p0 + (e.time - v2.t_init) * v2.v) - l.point());
        if (u_high < u_low)
          std::swap(u_high, u_low);

        if (u_low <= e.u && e.u <= u_high)
          crossing.push_back(sp_other);
      }
    }

    if (!decrement)
      return crossing.empty();

    if (crossing.size() <= v.k) {
      v.k -= crossing.size();
      return true;
    }
    else {
      return false;
    }
  }

  bool could_cross(std::size_t sp_idx, IFace face, const IVertex &on, const IkFT &t, const IEdge& edge) {
    CGAL_assertion(on != std::size_t(-1));
    if (m_data.igraph().iedge_is_on_bbox(edge))
      return false;

    Support_plane& sp = m_data.support_plane(sp_idx);

    // Does a polygon already exist on the face on other side?
    IFace other_face = m_data.igraph().get_other_face(sp_idx, edge, face);
    if (sp.data().polygons.find(other_face) != sp.data().polygons.end())
      return false;

    IkFT u = m_data.get_u_on_line(on, m_data.igraph().edge(edge).line);

    const IkLine_3& l = m_data.igraph().line(m_data.igraph().edge(edge).line);
    std::vector<std::size_t> crossing;
    for (auto& [sp_other, ov] : m_data.igraph().edge(edge).vertices) {
      if (sp_other == sp_idx)
        continue;

      Support_plane sp = m_data.support_plane(sp_other);

      Vertex& v1 = vts[ov.first];
      Vertex& v2 = vts[ov.second];

      if (v1.ivertex == on)
        return false;

      if (v2.ivertex == on)
        return false;

      IkFT u_low = l.to_vector() * (sp.to_3d(v1.p0 + (t - v1.t_init) * v1.v) - l.point());
      IkFT u_high = l.to_vector() * (sp.to_3d(v2.p0 + (t - v2.t_init) * v2.v) - l.point());
      if (u_high < u_low)
        std::swap(u_high, u_low);

      if (u_low <= u && u <= u_high)
        return false;
    }

    return true;
  }

  bool could_cross(std::size_t sp_idx, IFace face, const IkFT &u, const IkFT& t, const IEdge &edge) {
    if (m_data.igraph().iedge_is_on_bbox(edge))
      return false;

    Support_plane& sp = m_data.support_plane(sp_idx);

    // Does a polygon already exist on the face on other side?
    IFace other_face = m_data.igraph().get_other_face(sp_idx, edge, face);
    if (sp.data().polygons.find(other_face) != sp.data().polygons.end())
      return false;

    const IkLine_3& l = m_data.igraph().line(m_data.igraph().edge(edge).line);
    std::vector<std::size_t> crossing;
    for (auto& [sp_other, ov] : m_data.igraph().edge(edge).vertices) {
      if (sp_other == sp_idx)
        continue;

      Support_plane sp = m_data.support_plane(sp_other);

      Vertex& v1 = vts[ov.first];
      Vertex& v2 = vts[ov.second];

      IkFT u_low = l.to_vector() * (sp.to_3d(v1.p0 + (t - v1.t_init) * v1.v) - l.point());
      IkFT u_high = l.to_vector() * (sp.to_3d(v2.p0 + (t - v2.t_init) * v2.v) - l.point());
      if (u_high < u_low)
        std::swap(u_high, u_low);

      if (u_low <= u && u <= u_high)
        crossing.push_back(sp_other);
    }

    return crossing.empty();
  }

  bool expands_laterally(const Event &e, std::size_t hit_line_idx, IEdge hit_edge) {
    CGAL_assertion(vts[e.vertex].moving);
    CGAL_assertion(!vts[e.vertex].constraints.empty());
    CGAL_assertion(vts[e.vertex].constraint_edge != m_data.igraph().null_iedge());
    IFace q2_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, hit_edge, vts[e.vertex].face);
    IFace q4_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, vts[e.vertex].constraint_edge, vts[e.vertex].face);
    Support_plane &sp = m_data.support_plane(vts[e.vertex].sp_idx);

    IEdge prolongation_edge = get_next_edge(vts[e.vertex].constraint_edge, e.destination, *vts[e.vertex].constraints.begin());

    // check if could cross
    if (!could_cross(vts[e.vertex].sp_idx, q2_face, e.destination, e.time, prolongation_edge))
      return false;

    IFace q3_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, prolongation_edge, q2_face);
    if (q4_face == std::size_t(-1) || q3_face == std::size_t(-1))
      return false;

    auto it = sp.data().polygons.find(q4_face);
    if (it == sp.data().polygons.end())
      return true;

    std::list<std::size_t> &poly = it->second;
    std::size_t other = vts[e.vertex].other;
    std::list<std::size_t>::iterator vit;
    if (other == std::size_t(-1)) {
      for (auto pit = poly.begin(); pit != poly.end(); ++pit)
        if (vts[*pit].constraint_edge == vts[e.vertex].constraint_edge && vts[*pit].itarget == vts[e.vertex].itarget) {
          vit = pit;
          other = *vit;
          break;
        }
    }
    else vit = std::find(poly.begin(), poly.end(), other);

    if (other == std::size_t(-1) || !vts[*vit].moving)
      return true;

    std::size_t adj = *((vit == poly.begin()) ? std::prev(poly.end()) : std::prev(vit));
    if (vts[adj].constraint_edge == vts[e.vertex].constraint_edge || vts[adj].other_constraint_edge == vts[e.vertex].constraint_edge)
      adj = *((std::next(vit) == poly.end()) ? poly.begin() : std::next(vit));
    CGAL_assertion(vts[adj].constraint_edge != vts[e.vertex].constraint_edge);

    if (!vts[adj].moving)
      return true;

    IkPoint_2 p_other = vts[*vit].p0 + (e.time - vts[*vit].t_init) * vts[*vit].v;
    IkPoint_2 p_adj = vts[adj].p0 + (e.time - vts[adj].t_init) * vts[adj].v;
    IkVector_2 adj_edge = p_other - p_adj;
    const IkLine_2 &line = sp.data().lines[hit_line_idx];

    if (!CGAL::is_zero(adj_edge.x() * line.to_vector().y() - adj_edge.y() * line.to_vector().x()))
      return true;

    IkFT t = intersection_time(vts[adj], sp.data().exact_plane.to_2d(m_data.igraph().point_3(vts[e.vertex].itarget)));

    IEdge q4_q3_edge = get_next_edge(hit_edge, e.destination, hit_line_idx);

    return !could_cross(vts[e.vertex].sp_idx, q4_face, vts[e.vertex].itarget, t, q4_q3_edge);
  }

  void case_A(const std::vector<Event>& events) {
    a_events++;
    const Event& e = events.front();
    if (events.size() > 1)
      t_events++;

    //export_event("e_A_" + std::to_string(handled_events) + "_" + std::to_string(vts[e.vertex].sp_idx) + "_", e);

    Support_plane& sp = m_data.support_plane(vts[e.vertex].sp_idx);
    sp.ea++;
    std::list<std::size_t>& poly = sp.data().polygons[vts[e.vertex].face];
    std::list<std::size_t>::iterator it = std::find(poly.begin(), poly.end(), e.vertex);
    std::list<std::size_t>::iterator it_next = ((std::next(it) == poly.end()) ? poly.begin() : std::next(it));
    std::list<std::size_t>::iterator it_prev = (it == poly.begin()) ? std::prev(poly.end()) : std::prev(it);
    CGAL_assertion(vts[e.vertex].moving);
    CGAL_assertion(vts[*it_next].moving);
    CGAL_assertion(vts[*it_prev].moving);

    // Create two constrained vertices moving into opposite directions on intersection line
    std::size_t line_idx = m_data.igraph().edge(e.crossed_edge).line;
    IkVector_2 dir_prev = sp.calculate_edge_speed(vts[e.vertex], vts[*it_prev], e.p, sp.data().lines[line_idx], e.time);
    IkVector_2 dir_next = sp.calculate_edge_speed(vts[e.vertex], vts[*it_next], e.p, sp.data().lines[line_idx], e.time);

    CGAL_assertion(dir_prev * dir_next < 0);

    vts.emplace_back(vts[e.vertex].sp_idx, e.p, dir_prev, e.time);
    vts.back().constraints.insert(line_idx);
    vts.back().face = vts[e.vertex].face;
    vts.back().itarget = m_data.get_target_ivertex(vts.back(), e.crossed_edge);
    vts.back().constraint_edge = e.crossed_edge;
    const std::size_t prev = vts.size() - 1;

    vts.emplace_back(vts[e.vertex].sp_idx, e.p, dir_next, e.time);
    vts.back().constraints.insert(line_idx);
    vts.back().face = vts[e.vertex].face;
    vts.back().itarget = (vts[prev].itarget == m_data.igraph().source(e.crossed_edge)) ? m_data.igraph().target(e.crossed_edge) : m_data.igraph().source(e.crossed_edge);
    vts.back().constraint_edge = e.crossed_edge;
    const std::size_t next = vts.size() - 1;

    poly.insert(it, prev);
    *it = next;

    std::vector<std::size_t> new_moving_vertices;
    new_moving_vertices.push_back(prev);
    new_moving_vertices.push_back(next);

    if (can_cross(e, e.crossed_edge)) {
      // Entering a new face
      IFace other_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, e.crossed_edge, vts[e.vertex].face);
      m_data.init_border(other_face);
      std::list<std::size_t>& new_poly = sp.data().polygons[other_face];
      CGAL_assertion(new_poly.empty());
      sp.active_polygons++;

      // Create paired vertices
      vts.emplace_back(vts[e.vertex].sp_idx, e.p, dir_prev, e.time);
      vts.back().constraints.insert(line_idx);
      vts.back().face = other_face;
      vts.back().itarget = vts[prev].itarget;
      vts.back().constraint_edge = e.crossed_edge;
      vts.back().other = prev;
      vts[prev].other = vts.size() - 1;
      new_moving_vertices.push_back(vts.size() - 1);
      new_poly.push_back(vts.size() - 1);

      vts.emplace_back(vts[e.vertex].sp_idx, e.p, dir_next, e.time);
      vts.back().constraints.insert(line_idx);
      vts.back().face = other_face;
      vts.back().itarget = vts[next].itarget;
      vts.back().constraint_edge = e.crossed_edge;
      vts.back().other = next;
      vts[next].other = vts.size() - 1;
      new_moving_vertices.push_back(vts.size() - 1);
      new_poly.push_back(vts.size() - 1);

      // prolongation vertex
      vts.emplace_back(vts[e.vertex].sp_idx, e.p, vts[e.vertex].v, e.time);
      vts.back().face = other_face;
      vts.back().k = vts[e.vertex].k;
      std::swap(vts.back().known_intersections, vts[e.vertex].known_intersections);
      std::swap(vts.back().cached_events, vts[e.vertex].cached_events);
      new_moving_vertices.push_back(vts.size() - 1);
      new_poly.push_back(vts.size() - 1);
    }

    // insert border
    m_data.igraph().edge(e.crossed_edge).vertices.insert(std::make_pair(vts[e.vertex].sp_idx, std::pair<std::size_t, std::size_t>(prev, next)));
    //CGAL_assertion(inserted);

    vts[e.vertex].stop(e);

    for (std::size_t v : new_moving_vertices)
      calculate_events(v, e.time);

    bool still_moving = false;
    for (std::size_t i : poly)
      if (vts[i].moving) {
        still_moving = true;
        break;
      }
    if (!still_moving)
      sp.active_polygons--;
    CGAL_assertion_code(check_vertices_on_edge(e.crossed_edge);)
  }

  void case_B(const std::vector<Event>& events, std::size_t adj) {
    b_events++;
    const Event& e = events.front();
    if (events.size() > 1)
      t_events++;
    CGAL_assertion(adj != std::size_t(-1));
    Support_plane& sp = m_data.support_plane(vts[e.vertex].sp_idx);
    sp.eb++;

    //export_event("e_B_" + std::to_string(handled_events) + "_" + std::to_string(vts[e.vertex].sp_idx) + "_", e);

    std::size_t line_idx = m_data.igraph().edge(vts[adj].constraint_edge).line;

    // Get edge on other side of v (not connected to adj) to calculate new speed
    std::list<std::size_t>& poly = sp.data().polygons[vts[e.vertex].face];
    std::list<std::size_t>::iterator it = std::find(poly.begin(), poly.end(), e.vertex);
    std::list<std::size_t>::iterator it_prev = ((std::next(it) == poly.end()) ? poly.begin() : std::next(it));
    if (*it_prev == adj)
      it_prev = (it == poly.begin()) ? std::prev(poly.end()) : std::prev(it);

    CGAL_assertion(*it != adj);

    IkVector_2 dir(0, 0);

    dir = sp.calculate_edge_speed(vts[e.vertex], vts[*it_prev], e.p, sp.data().lines[line_idx], e.time);
    IkVector_2 u_ref = vts[adj].v - vts[e.vertex].v;

    std::vector<std::size_t> new_moving_vertices;

    vts.emplace_back(vts[e.vertex].sp_idx, e.p, dir, e.time);
    vts.back().face = vts[e.vertex].face;
    vts.back().constraints = vts[adj].constraints;
    vts.back().itarget = vts[adj].itarget;
    vts.back().constraint_edge = vts[adj].constraint_edge;

    new_moving_vertices.push_back(vts.size() - 1);

    // Replace former vts[adj] in polygon and replace v in polygon

    std::list<std::size_t>::iterator it_adj = ((std::next(it) == poly.end()) ? poly.begin() : std::next(it));
    std::list<std::size_t>::iterator it_before_adj = ((std::next(it_adj) == poly.end()) ? poly.begin() : std::next(it_adj));

    if (*it_adj != adj) {
      it_adj = (it == poly.begin()) ? std::prev(poly.end()) : std::prev(it);
      it_before_adj = (it_adj == poly.begin()) ? std::prev(poly.end()) : std::prev(it_adj);
    }

    CGAL_assertion(*it_adj == adj);

    *it = vts.size() - 1;
    poly.erase(it_adj);

    // update border
    auto border_it = m_data.igraph().edge(vts[adj].constraint_edge).vertices.find(vts[e.vertex].sp_idx);
    CGAL_assertion(border_it != m_data.igraph().edge(vts[adj].constraint_edge).vertices.end());
    border_it->second.first = *it_before_adj;
    border_it->second.second = *it;

    std::size_t other = vts[adj].other;
    std::size_t other_face;
    CGAL_assertion(other == std::size_t(-1) || vts[other].other == adj);
    if (other == std::size_t(-1)) {
      // Check if there is a vertex on the other side of vts[adj];
      other_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, vts[adj].constraint_edge, vts[e.vertex].face);
      typename std::unordered_map<IFace, std::list<std::size_t> >::iterator it;
      if (other_face != std::size_t(-1) && (it = sp.data().polygons.find(other_face)) != sp.data().polygons.end()) {
        std::list<std::size_t>& other_poly = it->second;
        CGAL_assertion(!other_poly.empty());
        auto it = other_poly.begin();
        for (; it != other_poly.end(); it++)
          if (vts[*it].constraint_edge == vts[adj].constraint_edge && vts[*it].moving)
            if (vts[*it].p0 + (e.time - vts[*it].t_init) * vts[*it].v == e.p)
              break;
        if (it != other_poly.end())
          other = *it;
      }
    }
    else other_face = vts[other].face;

    if (other != std::size_t(-1)) {
      // Insert new vertex into other polygon if required
      std::list<std::size_t>& other_poly = sp.data().polygons[other_face];
      CGAL_assertion(!other_poly.empty());
      it = std::find(other_poly.begin(), other_poly.end(), other);

      CGAL_assertion(it != other_poly.end());

      it_adj = (std::next(it) == other_poly.end()) ? other_poly.begin() : std::next(it);
      it_prev = (it == other_poly.begin()) ? std::prev(other_poly.end()) : std::prev(it);
      bool before = false;
      if (!vts[*it_adj].moving || (!vts[*it_adj].constraints.empty() && *vts[*it_adj].constraints.begin() == line_idx)) {
        CGAL_assertion(vts[*it_prev].constraints.empty() || *vts[*it_prev].constraints.begin() != line_idx);
        it_adj = it_prev;
      }
      else {
        CGAL_assertion(!vts[*it_prev].moving || *vts[*it_prev].constraints.begin() == line_idx);
        before = true;
      }

      CGAL_assertion(*it_adj != *it);
      CGAL_assertion(vts[*it].moving);
      CGAL_assertion(vts[*it_adj].moving);

      IkPoint_2 a_t = vts[*it].p0 + (1 + e.time - vts[*it].t_init) * vts[*it].v;
      IkPoint_2 b_t = vts[*it_adj].p0 + (1 + e.time - vts[*it_adj].t_init) * vts[*it_adj].v;

      IkVector_2 u_adj = b_t - a_t;
      if (CGAL::is_zero(u_adj.x() * u_ref.y() - u_adj.y() * u_ref.x())) {
        vts.emplace_back(vts[e.vertex].sp_idx, e.p, dir, e.time);
        vts.back().face = vts[*it].face;
        vts.back().constraints = vts[*it].constraints;
        vts.back().itarget = vts[*it].itarget;
        vts.back().constraint_edge = vts[*it].constraint_edge;
        vts.back().other = vts.size() - 2;
        vts[vts.size() - 2].other = vts.size() - 1;
        new_moving_vertices.push_back(vts.size() - 1);

        // Disable replaced vertex
        vts[*it].stop(e);

        if (before) {
          other_poly.insert(it_adj, vts.size() - 1);
          other_poly.erase(it);
          it = std::prev(it_adj);
        }
        else {
          other_poly.insert(it, vts.size() - 1);
          it = std::prev(other_poly.erase(it));
        }

        vts.emplace_back(vts[e.vertex].sp_idx, e.p, vts[e.vertex].v, e.time);
        vts.back().k = vts[e.vertex].k;
        vts.back().face = vts[*it].face;
        std::swap(vts.back().known_intersections, vts[e.vertex].known_intersections);
        std::swap(vts.back().cached_events, vts[e.vertex].cached_events);
        new_moving_vertices.push_back(vts.size() - 1);

        if (before)
          other_poly.insert(it_adj, vts.size() - 1);
        else
          other_poly.insert(it, vts.size() - 1);
      }
    }

    vts[e.vertex].stop(e);
    vts[adj].stop(e);

    for (std::size_t v : new_moving_vertices)
      calculate_events(v, e.time);


    CGAL_assertion_code(check_vertices_on_edge(vts[adj].constraint_edge);)
  }

  void case_C1(const std::vector<Event>& events) {
    c1_events++;
    const Event &e = events.front();
    Vertex& v1 = vts[e.vertex];
    Vertex& v2 = vts[events[1].vertex];

    CGAL_assertion(v1.moving);
    CGAL_assertion(v2.moving);
    CGAL_assertion(v1.constraints.size() == 1);
    CGAL_assertion(v2.constraints.size() == 1);
    CGAL_assertion(v1.constraint_edge != m_data.igraph().null_iedge());
    CGAL_assertion(v2.constraint_edge != m_data.igraph().null_iedge());
    //export_event("e_C1_" + std::to_string(handled_events) + "_" + std::to_string(v1.sp_idx) + "_", e);

    v1.stop(e);
    v2.stop(e);

    v1.constraints.insert(*v2.constraints.begin());
    v1.other_constraint_edge = v2.constraint_edge;
    v2.constraints = v1.constraints;
    v2.other_constraint_edge = v1.constraint_edge;
    v1.ivertex = v2.ivertex = e.destination;

    Support_plane& sp = m_data.support_plane(v1.sp_idx);
    sp.ec1++;

    bool still_moving = false;
    std::list<std::size_t>& poly = sp.data().polygons[v1.face];
    for (std::size_t i : poly)
      if (vts[i].moving) {
        still_moving = true;
        break;
      }

    if (!still_moving)
      sp.active_polygons--;
    else {
      // update border
      std::list<std::size_t>::iterator it = std::find(poly.begin(), poly.end(), e.vertex);
      bool ccw = true;
      std::list<std::size_t>::iterator it_adj = ((std::next(it) == poly.end()) ? poly.begin() : std::next(it));
      if (*it_adj == events[1].vertex) {
        ccw = false;
        it_adj = (it == poly.begin()) ? std::prev(poly.end()) : std::prev(it);
      }

      std::list<IEdge>& border = sp.data().borders[v1.face];
      if (!vts[*it_adj].moving) {
        auto eit = std::find(border.begin(), border.end(), v1.constraint_edge);
        CGAL_assertion(eit != border.end());
        border.erase(eit);

        // if closed, update vertices also
        auto border_v1 = m_data.igraph().edge(v1.constraint_edge).vertices.find(v1.sp_idx);
        CGAL_assertion(border_v1 != m_data.igraph().edge(v1.constraint_edge).vertices.end());
        border_v1->second.first = *it_adj;
        border_v1->second.second = e.vertex;
      }

      // Get adjacent vertex on the other side
      if (ccw) {
        it_adj = (it == poly.begin()) ? std::prev(poly.end()) : std::prev(it);
        CGAL_assertion(*it_adj == events[1].vertex);
        it_adj = (it_adj == poly.begin()) ? std::prev(poly.end()) : std::prev(it_adj);
      }
      else {
        it_adj = ((std::next(it) == poly.end()) ? poly.begin() : std::next(it));
        CGAL_assertion(*it_adj == events[1].vertex);
        it_adj = ((std::next(it_adj) == poly.end()) ? poly.begin() : std::next(it_adj));
      }

      if (!vts[*it_adj].moving) {
        auto eit = std::find(border.begin(), border.end(), v2.constraint_edge);
        CGAL_assertion(eit != border.end());
        border.erase(eit);

        // if closed, update vertices also
        auto border_v2 = m_data.igraph().edge(v2.constraint_edge).vertices.find(v2.sp_idx);
        CGAL_assertion(border_v2 != m_data.igraph().edge(v2.constraint_edge).vertices.end());
        border_v2->second.first = *it_adj;
        border_v2->second.second = events[1].vertex;
      }
    }
  }

  void case_C2(const std::vector<Event>& events) {
    c2_events++;
    const Event& e = events.front();

    if (events.size() > 2)
      t_events++;

    Support_plane &sp = m_data.support_plane(vts[e.vertex].sp_idx);
    sp.ec2++;
    std::list<std::size_t> &poly = sp.data().polygons[vts[e.vertex].face];

    IEdge hit_edge = m_data.igraph().null_iedge();
    for (IEdge edge : sp.data().borders[vts[e.vertex].face]) {
      if (edge != vts[e.vertex].constraint_edge && (m_data.igraph().target(edge) == e.destination || m_data.igraph().source(edge) == e.destination)) {
        hit_edge = edge;
        break;
      }
    }
    CGAL_assertion(hit_edge != m_data.igraph().null_iedge());
    std::size_t hit_line_idx = m_data.igraph().edge(hit_edge).line;

    CGAL_assertion(vts[e.vertex].moving);
    CGAL_assertion(vts[e.vertex].constraints.size() == 1);

    std::list<std::size_t>::iterator it = std::find(poly.begin(), poly.end(), e.vertex);
    std::list<std::size_t>::iterator prev = (it == poly.begin()) ? std::prev(poly.end()) : std::prev(it);
    std::list<std::size_t>::iterator next = (std::next(it) == poly.end()) ? poly.begin() : std::next(it);

    // Get the adjacent edge that is not constrained by the same line as the event vertex
    auto adjacent = (vts[*prev].constraints.find(m_data.igraph().edge(vts[e.vertex].constraint_edge).line) != vts[*prev].constraints.end()) ? next : prev;
    IkVector_2 dir(0, 0);

    dir = sp.calculate_edge_speed(vts[e.vertex], vts[*adjacent], e.p, sp.data().lines[hit_line_idx], e.time);
    IVertex target = (m_data.igraph().source(hit_edge) == e.destination ? m_data.igraph().target(hit_edge) : m_data.igraph().source(hit_edge));
    IkPoint_2 rico_target = sp.to_2d(m_data.igraph().point_3(target));
    IkVector_2 dir_check = rico_target - e.p;

    // Ricochet vertex has k = 0
    vts.emplace_back(vts[e.vertex].sp_idx, e.p, dir, e.time);
    std::size_t ricochet = vts.size() - 1;
    vts[ricochet].face = vts[e.vertex].face;
    vts[ricochet].constraints.insert(hit_line_idx);
    vts[ricochet].itarget = target;
    vts[ricochet].constraint_edge = m_data.igraph().edge(e.destination, vts[ricochet].itarget);
    vts[e.vertex].other_constraint_edge = vts[ricochet].constraint_edge;

    CGAL_assertion(vts[ricochet].itarget == m_data.get_target_ivertex(vts.back(), vts[ricochet].constraint_edge));

    // insert in polygon
    if (adjacent == prev)
      poly.insert(it, ricochet);
    else
      poly.insert(next, ricochet);

    // Add vertices to newly covered edge
    auto pair = m_data.igraph().edge(hit_edge).vertices.insert(std::make_pair(vts[e.vertex].sp_idx, std::pair<std::size_t, std::size_t>(e.vertex, ricochet)));

    // vertices in already covered edge property don't need adjustment as only indices are stored

    // Expanding to other side of edge is still necessary (check if there is a polygon in that face)
    // May also be an initial vertex on line, thus checking k for crossing is required

    std::vector<std::size_t> new_moving_vertices;
    new_moving_vertices.push_back(ricochet);

    std::size_t corner = -1;
    std::size_t prolongation = -1;

    if (can_cross(e, hit_edge)) {
      // Initialize new face
      IFace other_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, hit_edge, vts[e.vertex].face);
      m_data.init_border(other_face);
      std::list<std::size_t>& new_poly = sp.data().polygons[other_face];
      CGAL_assertion(new_poly.empty());
      sp.active_polygons++;

      // stationary vertex at corner
      vts.emplace_back(vts[e.vertex].sp_idx, e.p, IkVector_2(0, 0), e.time);
      vts.back().moving = false;
      vts.back().constraints = vts[e.vertex].constraints;
      vts.back().constraints.insert(hit_line_idx);
      vts.back().constraint_edge = vts[ricochet].constraint_edge;
      vts.back().face = other_face;
      vts.back().ivertex = e.destination;

      corner = vts.size() - 1;

      // paired vertex for ricochet vertex
      vts.emplace_back(vts[e.vertex].sp_idx, e.p, vts[ricochet].v, e.time);
      vts.back().face = other_face;
      vts.back().constraints.insert(hit_line_idx);
      vts.back().itarget = vts[ricochet].itarget;
      vts.back().constraint_edge = vts[ricochet].constraint_edge;
      vts.back().other = ricochet;
      vts[ricochet].other = vts.size() - 1;
      // No need to insert as ricochet vertex has been already added

      // prolongation vertex
      vts.emplace_back(vts[e.vertex].sp_idx, e.p, vts[e.vertex].v, e.time);
      vts.back().constraints = vts[e.vertex].constraints;
      vts.back().face = other_face;
      vts.back().k = vts[e.vertex].k;
      vts.back().itarget = m_data.get_target_ivertex(vts.back(), e.destination, m_data.igraph().edge(vts[e.vertex].constraint_edge).line);
      vts.back().constraint_edge = m_data.igraph().edge(e.destination, vts.back().itarget);
      vts[corner].other_constraint_edge = vts.back().constraint_edge;
      for (Cached_event &ce : vts[e.vertex].cached_events)
        if (e.destination != ce.d)
        vts.back().cached_events.push_back(std::move(ce));

      std::swap(vts.back().known_intersections, vts[e.vertex].known_intersections);
      vts[e.vertex].cached_events.clear();

      CGAL_assertion_code(
        {
        bool found = false;
      for (IEdge e : m_data.igraph().face(other_face).edges)
        if (e == vts.back().constraint_edge) {
          found = true;
          break;
        }
      CGAL_assertion(found);
        });

      m_data.igraph().edge(vts.back().constraint_edge).vertices.insert(std::make_pair(vts[e.vertex].sp_idx, std::pair<std::size_t, std::size_t>(vts.size() - 3, vts.size() - 1)));

      prolongation = vts.size() - 1;

      new_poly.push_back(corner);
      new_poly.push_back(vts.size() - 2);
      new_poly.push_back(vts.size() - 1);

      new_moving_vertices.push_back(vts.size() - 2);
      new_moving_vertices.push_back(vts.size() - 1);

      // Check whether a polygon needs to be created on other side of edge
      // Edge for prolongation vertex is needed anyway
      std::size_t q3_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, vts.back().constraint_edge, other_face);
      if (q3_face != -1 && sp.data().polygons.find(q3_face) == sp.data().polygons.end() && could_cross(vts[e.vertex].sp_idx, other_face, e.destination, e.time, vts[prolongation].constraint_edge)) {
        // Create new polygon with prolongation vertex and ivertex
        m_data.init_border(q3_face);
        std::list<std::size_t>& poly_q3 = sp.data().polygons[q3_face];
        CGAL_assertion(poly_q3.empty());
        sp.active_polygons++;

        // Insert vertices into q3 like in the reference implementation
        std::list<std::size_t> polygon(sp.data().exact_vertices.size());
        std::vector<Vertex> vertices;
        std::iota(polygon.begin(), polygon.end(), 0);

        vertices.reserve(sp.data().exact_vertices.size());

        IkVector_2 shift = (e.p - sp.data().ikcentroid);

        // get initial polygon for support plane
        for (std::size_t v = 0; v < sp.data().exact_vertices.size(); v++)
          vertices.emplace_back(vts[e.vertex].sp_idx, sp.data().exact_vertices[v] + shift, sp.data().original_rays[v].to_vector());

        std::list<std::size_t> other;
        const IkLine_2& constraint_line = sp.data().lines[*vts[e.vertex].constraints.begin()];
        const IkLine_2& hit_line = sp.data().lines[hit_line_idx];
        m_data.split_polygon(vertices, polygon, *vts[e.vertex].constraints.begin(), other);
        // Check which one is on the good side
        Oriented_side s = constraint_line.oriented_side(e.p + 5 * vts[ricochet].v);
        CGAL_assertion(s != CGAL::ON_ORIENTED_BOUNDARY);
        if (s == CGAL::ON_POSITIVE_SIDE) {
          other.clear();
          m_data.split_polygon(vertices, polygon, hit_line_idx, other);
        }
        else {
          polygon.clear();
          m_data.split_polygon(vertices, other, hit_line_idx, polygon);
          std::swap(other, polygon);
        }

        s = hit_line.oriented_side(e.p + 5 * vts[e.vertex].v);
        CGAL_assertion(s != CGAL::ON_ORIENTED_BOUNDARY);
        if (s == CGAL::ON_POSITIVE_SIDE)
          std::swap(other, polygon);

        std::size_t q3_stationary = std::size_t(-1);
        std::size_t q3_prolongation = std::size_t(-1);
        std::size_t opposite_ricochet = std::size_t(-1);
        for (std::size_t &i : polygon) {
          if (vertices[i].constraints.empty()) {
            continue;
          }
          vts.emplace_back(vertices[i]);
          i = vts.size() - 1;
          poly_q3.push_back(i);
          vts.back().t_init = e.time;
          vts.back().p0 = e.p;
          vts.back().face = q3_face;

          if (!vts.back().moving) {
            CGAL_assertion(q3_stationary == std::size_t(-1));
            CGAL_assertion(vts.back().constraints.size() == 2);
            vts.back().constraint_edge = vts[prolongation].constraint_edge;
            vts.back().other_constraint_edge = hit_edge;
            vts.back().ivertex = e.destination;
            q3_stationary = i;
          }
          else {
            new_moving_vertices.push_back(i);
            if (!vts.back().constraints.empty()) {
              if (*vts.back().constraints.begin() == hit_line_idx) {
                CGAL_assertion(opposite_ricochet == std::size_t(-1));
                opposite_ricochet = i;
                vts.back().itarget = m_data.get_target_ivertex(vts.back(), e.destination, hit_line_idx);
                vts.back().constraint_edge = m_data.igraph().edge(e.destination, vts.back().itarget);
                vts.back().other = std::size_t(-1);
              }
              else if (*vts.back().constraints.begin() == *vts[e.vertex].constraints.begin()) {
                CGAL_assertion(q3_prolongation == std::size_t(-1));
                q3_prolongation = i;
                //vts.back().k = vts[e.vertex].k;
                vts.back().itarget = vts[prolongation].itarget;
                vts.back().constraint_edge = vts[prolongation].constraint_edge;
                vts.back().other = std::size_t(-1);
              }
              CGAL_assertion_code(else CGAL_assertion(false);)
            }
          }
        }
        vts[q3_stationary].other_constraint_edge = vts[opposite_ricochet].constraint_edge;
        CGAL_assertion(!poly_q3.empty());
        CGAL_assertion_code(sp.check_edge_is_turning(vts[vts.size() - 2], vts[vts.size() - 1]));

        CGAL_assertion_code(
          {
        bool found = false;
        for (IEdge e : m_data.igraph().face(q3_face).edges)
          if (e == vts[q3_prolongation].constraint_edge) {
            found = true;
            break;
          }
        CGAL_assertion(found);
          });

        CGAL_assertion_code(
          {
        bool found = false;
        for (IEdge e : m_data.igraph().face(q3_face).edges)
          if (e == vts[opposite_ricochet].constraint_edge) {
            found = true;
            break;
          }
        CGAL_assertion(found);
          });

        m_data.igraph().edge(vts[opposite_ricochet].constraint_edge).vertices.insert(std::make_pair(vts[e.vertex].sp_idx, std::pair<std::size_t, std::size_t>(q3_stationary, opposite_ricochet)));

        std::size_t q4_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, vts[e.vertex].constraint_edge, vts[e.vertex].face);
        if (q4_face != -1 && sp.data().polygons.find(q4_face) == sp.data().polygons.end() && could_cross(vts[e.vertex].sp_idx, q3_face, e.destination, e.time, vts[opposite_ricochet].constraint_edge)) {
          // Create new polygon with prolongation vertex and ivertex
          m_data.init_border(q4_face);
          std::list<std::size_t>& poly_q4 = sp.data().polygons[q4_face];
          CGAL_assertion(poly_q4.empty());
          sp.active_polygons++;

          std::size_t q4_stationary = std::size_t(-1);
          std::size_t q4_inverse = std::size_t(-1);
          std::size_t q4_opposite_ricochet = std::size_t(-1);
          for (std::size_t& i : other) {
            if (vertices[i].constraints.empty())
              continue;
            vts.emplace_back(vertices[i]);
            i = vts.size() - 1;
            poly_q4.push_back(i);
            vts.back().t_init = e.time;
            vts.back().p0 = e.p;
            vts.back().face = q4_face;

            if (!vts.back().moving) {
              CGAL_assertion(q4_stationary == std::size_t(-1));
              CGAL_assertion(vts.back().constraints.size() == 2);
              vts.back().constraint_edge = vts[opposite_ricochet].constraint_edge;
              vts.back().other_constraint_edge = vts[e.vertex].constraint_edge;
              vts.back().ivertex = e.destination;
              q4_stationary = i;
            }
            else {
              new_moving_vertices.push_back(i);
              if (!vts.back().constraints.empty()) {
                if (*vts.back().constraints.begin() == hit_line_idx) {
                  CGAL_assertion(q4_opposite_ricochet == std::size_t(-1));
                  q4_opposite_ricochet = i;
                  vts.back().itarget = vts[opposite_ricochet].itarget;
                  vts.back().constraint_edge = vts[opposite_ricochet].constraint_edge;
                  vts.back().other = opposite_ricochet;
                  vts[opposite_ricochet].other = q4_opposite_ricochet;
                }
                else if (*vts.back().constraints.begin() == *vts[e.vertex].constraints.begin()) {
                  CGAL_assertion(q4_inverse == std::size_t(-1));
                  q4_inverse = i;
                  vts.back().itarget = (m_data.igraph().source(vts[e.vertex].constraint_edge) == e.destination) ? m_data.igraph().target(vts[e.vertex].constraint_edge) : m_data.igraph().source(vts[e.vertex].constraint_edge);
                  vts.back().constraint_edge = vts[e.vertex].constraint_edge;
                  vts.back().other = std::size_t(-1);
                }
                CGAL_assertion_code(else CGAL_assertion(false);)
              }
            }
          }
          CGAL_assertion(!poly_q4.empty());
          CGAL_assertion_code(sp.check_edge_is_turning(vts[vts.size() - 2], vts[vts.size() - 1]));

          CGAL_assertion_code(
            {
            bool found = false;
          for (IEdge e : m_data.igraph().face(q4_face).edges)
            if (e == vts[q4_opposite_ricochet].constraint_edge) {
              found = true;
              break;
            }
          CGAL_assertion(found);
            });

          CGAL_assertion_code(
            {
            bool found = false;
          for (IEdge e : m_data.igraph().face(q4_face).edges)
            if (e == vts[q4_inverse].constraint_edge) {
              found = true;
              break;
            }
          CGAL_assertion(found);
            });
        }
      }
    }

    // - add new vertex
    // traverse if possible
    // other quadrants (may actually be more if more than 3 planes intersect in the same point)
    // stop vertex

    vts[e.vertex].stop(e);
    vts[e.vertex].constraints.insert(hit_line_idx);
    vts[e.vertex].ivertex = e.destination;

    // Update edges & borders
    // Check if v.constraint_edge is fully covered now and remove from border
    auto vts_it = m_data.igraph().edge(vts[e.vertex].constraint_edge).vertices.find(vts[e.vertex].sp_idx);
    CGAL_assertion(vts_it != m_data.igraph().edge(vts[e.vertex].constraint_edge).vertices.end());
    if (!vts[vts_it->second.first].moving && !vts[vts_it->second.second].moving) {
      auto it = std::find(sp.data().borders[vts[e.vertex].face].begin(), sp.data().borders[vts[e.vertex].face].end(), vts[e.vertex].constraint_edge);
      CGAL_assertion(it != sp.data().borders[vts[e.vertex].face].end());
      sp.data().borders[vts[e.vertex].face].erase(it);
    }

    for (std::size_t v : new_moving_vertices)
      calculate_events(v, e.time);

    CGAL_assertion_code(check_vertices_on_edge(hit_edge));
    CGAL_assertion_code(check_vertices_on_edge(vts[e.vertex].constraint_edge));
  }

  void propagate_laterally(const Event& e, const Event& e2, IEdge hit_edge, std::vector<std::size_t> &new_moving_vertices) {
    CGAL_assertion(vts[e.vertex].moving);
    CGAL_assertion(!vts[e.vertex].constraints.empty());

    IFace q2_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, hit_edge, vts[e.vertex].face);
    IFace q4_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, vts[e.vertex].constraint_edge, vts[e.vertex].face);
    IEdge prolongation_edge = get_next_edge(vts[e.vertex].constraint_edge, e.destination, *vts[e.vertex].constraints.begin());

    Support_plane &sp = m_data.support_plane(vts[e.vertex].sp_idx);

    // Insert vertices into q3 like in the reference implementation
    std::list<std::size_t> polygon(sp.data().exact_vertices.size());
    std::vector<Vertex> vertices;
    std::iota(polygon.begin(), polygon.end(), 0);

    vertices.reserve(sp.data().exact_vertices.size());

    IkVector_2 shift = (e.p - sp.data().ikcentroid);

    // get initial polygon for support plane
    for (std::size_t v = 0; v < sp.data().exact_vertices.size(); v++)
      vertices.emplace_back(vts[e.vertex].sp_idx, sp.data().exact_vertices[v] + shift, sp.data().original_rays[v].to_vector());

    std::list<std::size_t> other;
    std::size_t hit_line_idx = m_data.igraph().edge(hit_edge).line;
    const IkLine_2& constraint_line = sp.data().lines[*vts[e.vertex].constraints.begin()];
    const IkLine_2& hit_line = sp.data().lines[hit_line_idx];
    m_data.split_polygon(vertices, polygon, *vts[e.vertex].constraints.begin(), other);
    // Check which one is on the good side
    Oriented_side s = constraint_line.oriented_side(e2.p);
    CGAL_assertion(s != CGAL::ON_ORIENTED_BOUNDARY);
    if (s == CGAL::ON_POSITIVE_SIDE) {
      other.clear();
      m_data.split_polygon(vertices, polygon, hit_line_idx, other);
    }
    else {
      polygon.clear();
      m_data.split_polygon(vertices, other, hit_line_idx, polygon);
      std::swap(other, polygon);
    }

    s = hit_line.oriented_side(e.p + 5 * vts[e.vertex].v);
    CGAL_assertion(s != CGAL::ON_ORIENTED_BOUNDARY);
    if (s == CGAL::ON_POSITIVE_SIDE)
      std::swap(other, polygon);

    std::size_t q3_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, prolongation_edge, q2_face);
    CGAL_assertion(q3_face != std::size_t(-1));
    CGAL_assertion(sp.data().polygons.find(q3_face) == sp.data().polygons.end());

    // Create new polygon with prolongation vertex and ivertex
    m_data.init_border(q3_face);
    std::list<std::size_t>& poly_q3 = sp.data().polygons[q3_face];
    CGAL_assertion(poly_q3.empty());
    sp.active_polygons++;

    std::size_t q3_stationary = std::size_t(-1);
    std::size_t q3_prolongation = std::size_t(-1);
    std::size_t opposite_ricochet = std::size_t(-1);
    for (std::size_t& i : polygon) {
      if (vertices[i].constraints.empty()) {
        continue;
      }
      vts.emplace_back(vertices[i]);
      i = vts.size() - 1;
      poly_q3.push_back(i);
      vts.back().t_init = e.time;
      vts.back().p0 = e.p;
      vts.back().face = q3_face;

      if (!vts.back().moving) {
        CGAL_assertion(q3_stationary == std::size_t(-1));
        CGAL_assertion(vts.back().constraints.size() == 2);
        vts.back().constraint_edge = prolongation_edge;
        vts.back().ivertex = e.destination;
        q3_stationary = i;
      }
      else {
        new_moving_vertices.push_back(i);
        if (!vts.back().constraints.empty()) {
          if (*vts.back().constraints.begin() == hit_line_idx) {
            CGAL_assertion(opposite_ricochet == std::size_t(-1));
            opposite_ricochet = i;
            vts.back().constraint_edge = get_next_edge(hit_edge, e.destination, hit_line_idx);
            vts.back().itarget = m_data.igraph().other(vts.back().constraint_edge, e.destination);
            vts.back().other = std::size_t(-1);
          }
          else if (*vts.back().constraints.begin() == *vts[e.vertex].constraints.begin()) {
            CGAL_assertion(q3_prolongation == std::size_t(-1));
            q3_prolongation = i;
            //vts.back().k = vts[e.vertex].k;
            vts.back().itarget = m_data.igraph().other(prolongation_edge, e.destination);
            vts.back().constraint_edge = prolongation_edge;
            vts.back().other = std::size_t(-1);
          }
          CGAL_assertion_code(else CGAL_assertion(false);)
        }
      }
    }
    vts[q3_stationary].other_constraint_edge = vts[opposite_ricochet].constraint_edge;
    CGAL_assertion(!poly_q3.empty());
    CGAL_assertion_code(sp.check_edge_is_turning(vts[vts.size() - 2], vts[vts.size() - 1]));

    CGAL_assertion_code(
      {
    bool found = false;
    for (IEdge e : m_data.igraph().face(q3_face).edges)
      if (e == vts[q3_prolongation].constraint_edge) {
        found = true;
        break;
      }
    CGAL_assertion(found);
      });

    CGAL_assertion_code(
      {
    bool found = false;
    for (IEdge e : m_data.igraph().face(q3_face).edges)
      if (e == vts[opposite_ricochet].constraint_edge) {
        found = true;
        break;
      }
    CGAL_assertion(found);
      });

    m_data.igraph().edge(vts[opposite_ricochet].constraint_edge).vertices.insert(std::make_pair(vts[e.vertex].sp_idx, std::pair<std::size_t, std::size_t>(q3_stationary, opposite_ricochet)));


    auto it = sp.data().polygons.find(q4_face);
    if (it != sp.data().polygons.end())
      return;

    m_data.init_border(q4_face);
    std::list<std::size_t>& poly_q4 = sp.data().polygons[q4_face];
    CGAL_assertion(poly_q4.empty());
    sp.active_polygons++;

    std::size_t q4_stationary = std::size_t(-1);
    std::size_t q4_inverse = std::size_t(-1);
    std::size_t q4_opposite_ricochet = std::size_t(-1);
    for (std::size_t& i : other) {
      if (vertices[i].constraints.empty())
        continue;
      vts.emplace_back(vertices[i]);
      i = vts.size() - 1;
      poly_q4.push_back(i);
      vts.back().t_init = e.time;
      vts.back().p0 = e.p;
      vts.back().face = q4_face;

      if (!vts.back().moving) {
        CGAL_assertion(q4_stationary == std::size_t(-1));
        CGAL_assertion(vts.back().constraints.size() == 2);
        vts.back().constraint_edge = vts[opposite_ricochet].constraint_edge;
        vts.back().other_constraint_edge = vts[e.vertex].constraint_edge;
        vts.back().ivertex = e.destination;
        q4_stationary = i;
      }
      else {
        new_moving_vertices.push_back(i);
        if (!vts.back().constraints.empty()) {
          if (*vts.back().constraints.begin() == hit_line_idx) {
            CGAL_assertion(q4_opposite_ricochet == std::size_t(-1));
            q4_opposite_ricochet = i;
            vts.back().itarget = vts[opposite_ricochet].itarget;
            vts.back().constraint_edge = vts[opposite_ricochet].constraint_edge;
            vts.back().other = opposite_ricochet;
            vts[opposite_ricochet].other = q4_opposite_ricochet;
          }
          else if (*vts.back().constraints.begin() == *vts[e.vertex].constraints.begin()) {
            CGAL_assertion(q4_inverse == std::size_t(-1));
            q4_inverse = i;
            vts.back().itarget = (m_data.igraph().source(vts[e.vertex].constraint_edge) == e.destination) ? m_data.igraph().target(vts[e.vertex].constraint_edge) : m_data.igraph().source(vts[e.vertex].constraint_edge);
            vts.back().constraint_edge = vts[e.vertex].constraint_edge;
            vts.back().other = std::size_t(-1);
          }
          CGAL_assertion_code(else CGAL_assertion(false);)
        }
      }
    }
    CGAL_assertion(!poly_q4.empty());
    CGAL_assertion_code(sp.check_edge_is_turning(vts[vts.size() - 2], vts[vts.size() - 1]));

    CGAL_assertion_code(
      {
      bool found = false;
    for (IEdge e : m_data.igraph().face(q4_face).edges)
      if (e == vts[q4_opposite_ricochet].constraint_edge) {
        found = true;
        break;
      }
    CGAL_assertion(found);
      });

    CGAL_assertion_code(
      {
      bool found = false;
    for (IEdge e : m_data.igraph().face(q4_face).edges)
      if (e == vts[q4_inverse].constraint_edge) {
        found = true;
        break;
      }
    CGAL_assertion(found);
      });
  }

  void case_D1(const std::vector<Event>& events) {
    d1_events++;
    CGAL_assertion(events.size() == 2);

    const Event& e1 = events.front();
    const Event& e2 = *++events.begin();

    CGAL_assertion(vts[e1.vertex].moving);
    CGAL_assertion(vts[e2.vertex].moving);
    CGAL_assertion(e1.crossed_edge == m_data.igraph().null_iedge() || e2.crossed_edge == m_data.igraph().null_iedge() || e1.crossed_edge == e2.crossed_edge);
    //export_event("e_D1_" + std::to_string(handled_events) + "_" + std::to_string(vts[e1.vertex].sp_idx) + "_", e1);

    IEdge hit_edge = (e1.crossed_edge == m_data.igraph().null_iedge()) ? e2.crossed_edge : e1.crossed_edge;
    // Check if both vertices are constrained, in this case,
    if (hit_edge == m_data.igraph().null_iedge())
      hit_edge = m_data.igraph().edge(e1.destination, e2.destination);
    CGAL_assertion(hit_edge != m_data.igraph().null_iedge());
    std::size_t hit_line_idx = m_data.igraph().edge(hit_edge).line;

    bool propagates_front = can_cross(e1, hit_edge) || can_cross(e2, hit_edge);
    bool propagates_e1_side = false, propagates_e2_side = false;
    if (propagates_front && !vts[e1.vertex].constraints.empty())
      propagates_e1_side = expands_laterally(e1, hit_line_idx, hit_edge); // Check lateral propagation
    if (propagates_front && !vts[e2.vertex].constraints.empty())
      propagates_e2_side = expands_laterally(e2, hit_line_idx, hit_edge); // Check lateral propagation

    // Local propagation
    Support_plane& sp = m_data.support_plane(vts[e1.vertex].sp_idx);
    sp.ed1++;
    std::list<std::size_t>& poly = sp.data().polygons[vts[e1.vertex].face];

    bool cw = true;
    std::size_t adj_e1 = -1, adj_e2 = -1;
    std::size_t inserted_adj_e1 = -1, inserted_adj_e2 = -1;
    std::list<std::size_t>::iterator e1_v_it = std::find(poly.begin(), poly.end(), e1.vertex);
    std::list<std::size_t>::iterator e2_v_it = ((std::next(e1_v_it) == poly.end()) ? poly.begin() : std::next(e1_v_it));
    std::list<std::size_t>::iterator e1_adj_it, e2_adj_it;
    if (*e2_v_it == e2.vertex) {
      e1_adj_it = (e1_v_it == poly.begin()) ? std::prev(poly.end()) : std::prev(e1_v_it);
      e2_adj_it = (std::next(e2_v_it) == poly.end()) ? poly.begin() : std::next(e2_v_it);
    }
    else {
      cw = false;
      e2_v_it = (e1_v_it == poly.begin()) ? std::prev(poly.end()) : std::prev(e1_v_it);
      CGAL_assertion(*e2_v_it == e2.vertex);

      e1_adj_it = (std::next(e1_v_it) == poly.end()) ? poly.begin() : std::next(e1_v_it);
      e2_adj_it = (e2_v_it == poly.begin()) ? std::prev(poly.end()) : std::prev(e2_v_it);
    }
    adj_e1 = *e1_adj_it;
    adj_e2 = *e2_adj_it;
    CGAL_assertion(adj_e1 != -1);
    CGAL_assertion(adj_e2 != -1);

    std::vector<std::size_t> new_moving_vertices;

    if (vts[e1.vertex].constraints.empty()) {
      CGAL_assertion(vts[adj_e1].moving);

      IkVector_2 dir_adj_e1 = sp.calculate_edge_speed(vts[e1.vertex], vts[adj_e1], e1.p, sp.data().lines[hit_line_idx], e1.time);

      vts.emplace_back(vts[e1.vertex].sp_idx, e1.p, dir_adj_e1, e1.time);
      vts.back().constraints.insert(hit_line_idx);
      vts.back().face = vts[e1.vertex].face;
      vts.back().itarget = m_data.get_target_ivertex(vts.back(), hit_edge);
      vts.back().constraint_edge = hit_edge;

      inserted_adj_e1 = vts.size() - 1;
      new_moving_vertices.push_back(inserted_adj_e1);
      if (cw)
        poly.insert(e1_v_it, inserted_adj_e1);
      else
        poly.insert(e1_adj_it, inserted_adj_e1);
    }

    if (vts[e2.vertex].constraints.empty()) {
      CGAL_assertion(vts[adj_e2].moving);

      IkVector_2 dir_adj_e2 = sp.calculate_edge_speed(vts[e2.vertex], vts[adj_e2], e2.p, sp.data().lines[hit_line_idx], e2.time);

      vts.emplace_back(vts[e2.vertex].sp_idx, e2.p, dir_adj_e2, e2.time);
      vts.back().constraints.insert(hit_line_idx);
      vts.back().face = vts[e2.vertex].face;
      if (inserted_adj_e1 == -1)
        vts.back().itarget = m_data.get_target_ivertex(vts.back(), hit_edge);
      else
        vts.back().itarget = (vts[inserted_adj_e1].itarget == m_data.igraph().source(hit_edge)) ? m_data.igraph().target(hit_edge) : m_data.igraph().source(hit_edge);
      vts.back().constraint_edge = hit_edge;

      inserted_adj_e2 = vts.size() - 1;
      new_moving_vertices.push_back(inserted_adj_e2);
      if (cw)
        poly.insert(e2_adj_it, inserted_adj_e2);
      else
        poly.insert(e2_v_it, inserted_adj_e2);
    }

    if (propagates_front) {
      IFace other_face = m_data.igraph().get_other_face(vts[e1.vertex].sp_idx, hit_edge, vts[e1.vertex].face);
      m_data.init_border(other_face);
      std::list<std::size_t>& new_poly = sp.data().polygons[other_face];
      CGAL_assertion(new_poly.empty());
      sp.active_polygons++;

      // Create adjacent of e1.vertex prolongation
      if (!vts[e1.vertex].constraints.empty()) {
        IEdge prolongation_edge_e1 = get_next_edge(vts[e1.vertex].constraint_edge, e1.destination, *vts[e1.vertex].constraints.begin());
        vts.emplace_back(vts[e1.vertex].sp_idx, e1.p, IkVector_2(0, 0), e1.time);
        vts.back().moving = false;
        vts.back().constraints = vts[e1.vertex].constraints;
        vts.back().constraints.insert(hit_line_idx);
        vts.back().face = other_face;
        vts.back().ivertex = vts[e1.vertex].itarget;
        vts.back().constraint_edge = hit_edge;
        vts.back().other_constraint_edge = prolongation_edge_e1;
        new_poly.push_back(vts.size() - 1);
      }
      else {
        CGAL_assertion(inserted_adj_e1 != -1);
        vts.emplace_back(vts[e1.vertex].sp_idx, e1.p, vts[inserted_adj_e1].v, e1.time);
        vts.back().constraints = vts[inserted_adj_e1].constraints;
        vts.back().face = other_face;
        vts.back().itarget = vts[inserted_adj_e1].itarget;
        vts.back().constraint_edge = hit_edge;
        vts.back().other = inserted_adj_e1;
        vts[inserted_adj_e1].other = vts.size() - 1;
        new_moving_vertices.push_back(vts.size() - 1);
        new_poly.push_back(vts.size() - 1);
      }

      // prolongation e1.vertex
      vts.emplace_back(vts[e1.vertex].sp_idx, e1.p, vts[e1.vertex].v, e1.time);
      vts.back().constraints = vts[e1.vertex].constraints;
      vts.back().face = other_face;
      vts.back().k = vts[e1.vertex].k;
      if (!vts[e1.vertex].constraints.empty()) {
        vts.back().constraint_edge = vts[vts.size() - 2].other_constraint_edge;
        vts.back().itarget = m_data.igraph().other(vts.back().constraint_edge, e1.destination);
      }
      std::swap(vts.back().known_intersections, vts[e1.vertex].known_intersections);
      std::swap(vts.back().cached_events, vts[e1.vertex].cached_events);
      new_moving_vertices.push_back(vts.size() - 1);
      new_poly.push_back(vts.size() - 1);

      if (!vts[e1.vertex].constraints.empty())
        m_data.igraph().edge(vts.back().constraint_edge).vertices.insert(std::make_pair(vts[e1.vertex].sp_idx, std::pair<std::size_t, std::size_t>(vts.size() - 2, vts.size() - 1)));

      // prolongation e2.vertex
      vts.emplace_back(vts[e2.vertex].sp_idx, e2.p, vts[e2.vertex].v, e2.time);
      vts.back().constraints = vts[e2.vertex].constraints;
      vts.back().face = other_face;
      vts.back().k = vts[e2.vertex].k;
      if (!vts[e2.vertex].constraints.empty()) {
        vts.back().constraint_edge = get_next_edge(vts[e2.vertex].constraint_edge, e2.destination, *vts[e2.vertex].constraints.begin());
        vts.back().itarget = m_data.igraph().other(vts.back().constraint_edge, e2.destination);
      }
      std::swap(vts.back().known_intersections, vts[e2.vertex].known_intersections);
      std::swap(vts.back().cached_events, vts[e2.vertex].cached_events);
      new_moving_vertices.push_back(vts.size() - 1);
      new_poly.push_back(vts.size() - 1);

      // Create adjacent of e2.vertex prolongation
      if (!vts[e2.vertex].constraints.empty()) {
        vts.emplace_back(vts[e2.vertex].sp_idx, e2.p, IkVector_2(0, 0), e2.time);
        vts.back().moving = false;
        vts.back().constraints = vts[e2.vertex].constraints;
        vts.back().constraints.insert(hit_line_idx);
        vts.back().face = other_face;
        vts.back().ivertex = vts[e2.vertex].itarget;
        vts.back().constraint_edge = hit_edge;
        vts.back().other_constraint_edge = vts[vts.size() - 2].constraint_edge;
        new_poly.push_back(vts.size() - 1);
        m_data.igraph().edge(vts.back().other_constraint_edge).vertices.insert(std::make_pair(vts[e1.vertex].sp_idx, std::pair<std::size_t, std::size_t>(vts.size() - 2, vts.size() - 1)));
      }
      else {
        CGAL_assertion(inserted_adj_e2 != -1);
        vts.emplace_back(vts[e2.vertex].sp_idx, e2.p, vts[inserted_adj_e2].v, e2.time);
        vts.back().constraints = vts[inserted_adj_e2].constraints;
        vts.back().face = other_face;
        vts.back().itarget = vts[inserted_adj_e2].itarget;
        vts.back().constraint_edge = hit_edge;
        vts.back().other = inserted_adj_e2;
        vts[inserted_adj_e2].other = vts.size() - 1;
        new_moving_vertices.push_back(vts.size() - 1);
        new_poly.push_back(vts.size() - 1);
      }
    }

    // propagate laterally
    if (!vts[e1.vertex].constraints.empty() && propagates_e1_side) {
      propagate_laterally(e1, e2, hit_edge, new_moving_vertices);
    }
    // propagate laterally
    if (!vts[e2.vertex].constraints.empty() && propagates_e2_side) {
      propagate_laterally(e2, e1, hit_edge, new_moving_vertices);
    }

    // insert new boundary on edge
    if (!vts[e1.vertex].constraints.empty())
      inserted_adj_e1 = e1.vertex;
    CGAL_assertion(inserted_adj_e1 != -1);

    if (!vts[e2.vertex].constraints.empty())
      inserted_adj_e2 = e2.vertex;
    CGAL_assertion(inserted_adj_e2 != -1);

    m_data.igraph().edge(hit_edge).vertices.insert(std::make_pair(vts[e1.vertex].sp_idx, std::pair<std::size_t, std::size_t>(inserted_adj_e1, inserted_adj_e2)));

    // stop event vertices
    if (vts[e1.vertex].constraints.empty())
      poly.erase(e1_v_it);

    if (vts[e2.vertex].constraints.empty())
      poly.erase(e2_v_it);

    vts[e1.vertex].stop(e1);
    vts[e2.vertex].stop(e2);
    vts[e1.vertex].constraints.insert(hit_line_idx);
    vts[e1.vertex].other_constraint_edge = hit_edge;
    vts[e2.vertex].constraints.insert(hit_line_idx);
    vts[e2.vertex].other_constraint_edge = hit_edge;

    for (std::size_t v : new_moving_vertices)
      calculate_events(v, e1.time);
  }

  void case_D2(const std::vector<Event>& events) {
    d2_events++;
    CGAL_assertion(events.size() == 4);

    t_events++;
    std::set<std::size_t> vertices;
    std::cout << "[" << vts[events[0].vertex].sp_idx << "] " << events[0].vertex << " " << events[0].time << " D2" << std::flush;

    std::size_t ec_idx = std::size_t(-1), ec2_idx = std::size_t(-1), e1_idx = std::size_t(-1), e2_idx = std::size_t(-1);
    for (std::size_t i = 0;i<events.size();i++)
      if (!vertices.insert(events[i].vertex).second) {
        ec_idx = i;
        break;
      }

    for (std::size_t i = 0; i < events.size(); i++)
      if (events[i].vertex != events[ec_idx].vertex) {
        if (e1_idx == std::size_t(-1))
          e1_idx = i;
        else
          e2_idx = i;
      }
      else if (i != ec_idx)
        ec2_idx = i;
    CGAL_assertion(ec_idx != std::size_t(-1));
    CGAL_assertion(e1_idx != std::size_t(-1));
    CGAL_assertion(e2_idx != std::size_t(-1));
    const Event &ec = events[ec_idx], &ec2 = events[ec2_idx], &e1 = events[e1_idx], &e2 = events[e2_idx];
    CGAL_assertion(ec.crossed_edge != m_data.igraph().null_iedge());
    CGAL_assertion(events[ec2_idx].crossed_edge != m_data.igraph().null_iedge());

    std::cout << "[" << vts[events[0].vertex].sp_idx << "] " << events[0].vertex << " " << events[0].time << " D2" << std::flush;
    //export_event("e_D2_" + std::to_string(handled_events) + "_" + std::to_string(vts[events[0].vertex].sp_idx) + "_", e);

    Support_plane& sp = m_data.support_plane(vts[ec.vertex].sp_idx);
    std::list<std::size_t>& poly = sp.data().polygons[vts[e1.vertex].face];
    sp.ed2++;

    IEdge e1_hit_edge = e1.crossed_edge;
    IEdge e2_hit_edge = e2.crossed_edge;

    if (e1_hit_edge == m_data.igraph().null_iedge()) {
      CGAL_assertion(e1.destination != std::size_t(-1));
      if (m_data.igraph().source(ec.crossed_edge) == e1.destination || m_data.igraph().target(ec.crossed_edge) == e1.destination)
        e1_hit_edge = ec.crossed_edge;
      else {
        CGAL_assertion(m_data.igraph().source(ec2.crossed_edge) == e1.destination || m_data.igraph().target(ec2.crossed_edge) == e1.destination);
        e1_hit_edge = ec2.crossed_edge;
      }
    }
    CGAL_assertion(e1_hit_edge != m_data.igraph().null_iedge());
    if (e2_hit_edge == m_data.igraph().null_iedge()) {
      if (e1_hit_edge == ec.crossed_edge)
        e2_hit_edge = ec2.crossed_edge;
      else
        e2_hit_edge = ec.crossed_edge;
    }
    CGAL_assertion(e2_hit_edge != m_data.igraph().null_iedge());
    CGAL_assertion(e1_hit_edge != e2_hit_edge);

    std::size_t e1_hit_line_idx = m_data.igraph().edge(e1_hit_edge).line;
    std::size_t e2_hit_line_idx = m_data.igraph().edge(e2_hit_edge).line;

    bool e1_propagates_frontally = can_cross(ec, e1_hit_edge) || can_cross(e1, e1_hit_edge);
    bool e2_propagates_frontally = can_cross(ec, e2_hit_edge) || can_cross(e2, e2_hit_edge);
    //bool propagates_e1_side = false, propagates_e2_side = false;

    //if (e1_propagates_frontally && !vts[e1.vertex].constraints.empty())
    //  propagates_e1_side = expands_laterally(e1, e1_hit_line_idx, e1_hit_edge); // Check lateral propagation
    //if (e2_propagates_frontally && !vts[e2.vertex].constraints.empty())
    //  propagates_e2_side = expands_laterally(e2, e2_hit_line_idx, e2_hit_edge); // Check lateral propagation

    if (e1_propagates_frontally || e2_propagates_frontally) {
      std::cout << "D2 frontal propagation case not implemented" << std::endl;
    }

    // Create lateral vertices in case e1.vertex or e2.vertex are not constrained
    std::vector<std::size_t> new_moving_vertices;

    // get adjacent vertices
    std::list<std::size_t>::iterator e1_adj, e2_adj, e2_it;
    std::list<std::size_t>::iterator e1_it = std::find(poly.begin(), poly.end(), e1.vertex);
    std::list<std::size_t>::iterator it_next = ((std::next(e1_it) == poly.end()) ? poly.begin() : std::next(e1_it));
    std::list<std::size_t>::iterator it_prev = (e1_it == poly.begin()) ? std::prev(poly.end()) : std::prev(e1_it);
    bool e1_prev_adj;
    if (*it_next == ec.vertex) {
      e1_prev_adj = true;
      e1_adj = it_prev;
      e2_it = ((std::next(it_next) == poly.end()) ? poly.begin() : std::next(it_next));
      it_next = ((std::next(e2_it) == poly.end()) ? poly.begin() : std::next(e2_it));
      e2_adj = it_next;
    }
    else {
      e1_prev_adj = false;
      e1_adj = it_next;
      e2_it = (it_prev == poly.begin()) ? std::prev(poly.end()) : std::prev(it_prev);
      it_prev = (e2_it == poly.begin()) ? std::prev(poly.end()) : std::prev(e2_it);
      e2_adj = it_prev;
    }

    if (vts[e1.vertex].constraints.empty()) {
      IkVector_2 dir = sp.calculate_edge_speed(vts[e1.vertex], vts[*e1_adj], e1.p, sp.data().lines[e1_hit_line_idx], e1.time);

      vts.emplace_back(vts[e1.vertex].sp_idx, e1.p, dir, e1.time);
      vts.back().constraints.insert(e1_hit_line_idx);
      vts.back().face = vts[e1.vertex].face;
      vts.back().itarget = m_data.opposite(e1_hit_edge, ec.destination);
      vts.back().constraint_edge = e1_hit_edge;

      // insert in poly
      if (e1_prev_adj)
        poly.insert(e1_it, vts.size() - 1);
      else
        poly.insert(e1_adj, vts.size() - 1);

      // insert in edge vertices
      m_data.igraph().edge(e1_hit_edge).vertices.insert(std::make_pair(vts[e1.vertex].sp_idx, std::pair<std::size_t, std::size_t>(ec.vertex, vts.size() - 1)));

      // insert in new vertices
      new_moving_vertices.push_back(vts.size() - 1);
    }

    if (vts[e2.vertex].constraints.empty()) {
      IkVector_2 dir = sp.calculate_edge_speed(vts[e2.vertex], vts[*e2_adj], e2.p, sp.data().lines[e2_hit_line_idx], e2.time);

      vts.emplace_back(vts[e2.vertex].sp_idx, e2.p, dir, e2.time);
      vts.back().constraints.insert(e2_hit_line_idx);
      vts.back().face = vts[e2.vertex].face;
      vts.back().itarget = m_data.opposite(e2_hit_edge, ec.destination);
      vts.back().constraint_edge = e2_hit_edge;

      // insert in poly
      if (e1_prev_adj)
        poly.insert(e2_adj, vts.size() - 1);
      else
        poly.insert(e2_it, vts.size() - 1);

      // insert in edge vertices
      m_data.igraph().edge(e2_hit_edge).vertices.insert(std::make_pair(vts[e2.vertex].sp_idx, std::pair<std::size_t, std::size_t>(ec.vertex, vts.size() - 1)));

      // insert in new vertices
      new_moving_vertices.push_back(vts.size() - 1);
    }

    vts[ec.vertex].stop(ec);
    vts[e1.vertex].stop(e1);
    vts[e2.vertex].stop(e2);

    // stop 3 vertices
    for (std::size_t v : new_moving_vertices)
      calculate_events(v, e1.time);
  }

  std::vector<Event> get_connected_events(std::list<Event> &events, std::vector<std::size_t> &vertices) {
    std::vector<Event> connected;
    if (events.size() == 1) {
      connected.push_back(*events.begin());
      vertices.push_back(connected.front().vertex);
      events.clear();
      return connected;
    }

    connected.push_back(events.front());
    events.pop_front();

    std::set<std::size_t> unique_v;
    unique_v.insert(connected.front().vertex);

    auto it = events.begin();
    while (it != events.end()) {
      if (vts[it->vertex].face == vts[connected[0].vertex].face) {
        // check whether they hit the same edge (or a vertex of it) or if they all hit the same vertex
        unique_v.insert(it->vertex);
        connected.push_back(*it);
        it = events.erase(it);
        // I also need to check whether they are hitting the same edge)
      }

      if (it != events.end())
        it++;
    }

    vertices.clear();
    std::copy(unique_v.begin(), unique_v.end(), std::back_inserter(vertices));

    return connected;
  }

  std::size_t is_adjacent_vertex_constrained(std::size_t v_idx, std::size_t line_idx) {
    Support_plane &sp = m_data.support_plane(vts[v_idx].sp_idx);
    const std::list<std::size_t> &poly = sp.data().polygons[vts[v_idx].face];

    for (auto it = poly.begin();it != poly.end();++it)
      if (*it == v_idx) {
        std::list<std::size_t>::const_iterator prev = (it == poly.begin()) ? std::prev(poly.end()) : std::prev(it);

        if (vts[*prev].constraints.find(line_idx) != vts[*prev].constraints.end())
          return *prev;

        std::list<std::size_t>::const_iterator next = (std::next(it) == poly.end()) ? poly.begin() : std::next(it);
        if (vts[*next].constraints.find(line_idx) != vts[*next].constraints.end())
          return *next;
        else return std::size_t(-1);
      }
    return std::size_t(-1);
  }

  void apply(std::list<Event>& events) {

    std::list<Event> bkp = events;

    /*
    int total_active_polys = 0;
    int total_active_vertices = 0;

    std::size_t idx = 0;
    for (const Vertex &v : vts) {
      CGAL_assertion(v.face != std::size_t(-1));
      CGAL_assertion(v.sp_idx != std::size_t(-1));
      if (v.moving) {
        total_active_vertices++;
        CGAL_assertion(v.v != IkVector_2(0, 0));
        CGAL_assertion(v.queued_events != 0);
        if (v.constraints.empty()) {
          CGAL_assertion(v.other == std::size_t(-1));
          CGAL_assertion(v.itarget == IVertex(-1));
          CGAL_assertion(v.constraint_edge == IEdge(IVertex(-1), IVertex(-1), nullptr));
        }
        else {
          CGAL_assertion(v.constraints.size() == 1);
          CGAL_assertion(v.constraint_edge != IEdge(IVertex(-1), IVertex(-1), nullptr));
          CGAL_assertion(v.itarget != IVertex(-1));
        }
      }
      idx++;
    }
    for (std::size_t i = 6;i<m_data.number_of_support_planes();i++) {
      total_active_polys += m_data.support_plane(i).active_polygons;
    }*/

    //export_timestamp(std::to_string(handled_events) + "_" + std::to_string(CGAL::to_double(events.front().time)), events.front().time);

    while(!events.empty()) {
      //export_timestamp(std::to_string(handled_events) + "_" + std::to_string(CGAL::to_double(events.front().time)), events.front().time);

      std::vector<std::size_t> vertices;
      std::vector<Event> connected = get_connected_events(events, vertices);
      std::size_t adj;

      for (const Event& e : connected)
        vts[e.vertex].queued_events--;
      switch (vertices.size()) {
      case 1:
        if (!vts[vertices[0]].constraints.empty())
          case_C2(connected); // Single constrained vertex hits corner
        else if ((adj = is_adjacent_vertex_constrained(vertices[0], m_data.igraph().edge(connected[0].crossed_edge).line)) != std::size_t(-1))
          case_B(connected, adj); // Unconstrained vertex hits line with neighbor that is constrained by that line
        else
          case_A(connected); // Single unconstrained vertex hits line
        break;
      case 2:
        if (!vts[vertices[0]].constraints.empty() && !vts[vertices[1]].constraints.empty() && connected[0].destination == connected[1].destination)
          case_C1(connected); // Two constrained vertices meet in corner
        else
          case_D1(connected); // Edge with one constrained vertex hits line
        break;
      case 3:
        // check two_edges_intersect_two_lines
        case_D2(connected); // Edge with two unconstrained vertices hits line (like A, but with neighbor vertex hitting in the same time)
        break;
      default:
        std::cout << "Error: Unknown event type!" << std::endl;
      }
      handled_events++;
    }
  }
};

#endif //DOXYGEN_RUNNING

} // namespace internal
} // namespace KSP_3
} // namespace CGAL

#endif // CGAL_KSP_3_VERTEXPROPAGATION_H
