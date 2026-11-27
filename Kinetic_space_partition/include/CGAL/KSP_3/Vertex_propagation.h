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

  // Strict total order on the events. Simultaneous events are popped from the priority queue in the order given by
  // the tie-breaks below. Without them, the order of simultaneous events would depend on the heap implementation of
  // the standard library and on the insertion order, and the results would differ between platforms.
  template<typename Event>
  struct Event_order {
    static std::size_t edge_order(const IEdge& e) {
      using Edge_property = typename Data_structure::Intersection_graph::Edge_property;
      const Edge_property* p = static_cast<const Edge_property*>(e.get_property());
      return (p == nullptr) ? std::size_t(-1) : p->order;
    }

    bool operator()(const Event& a, const Event& b) const {
      if (a.time != b.time)
        return a.time > b.time;
      if (a.vertex != b.vertex)
        return a.vertex > b.vertex;
      const std::size_t ea = edge_order(a.crossed_edge), eb = edge_order(b.crossed_edge);
      if (ea != eb)
        return ea > eb;
      return a.destination > b.destination;
    }
  };

public:
  Vertex_propagation(Data_structure& data, const Parameters& parameters)
    : m_data(data), m_parameters(parameters), vts(data.kinetic_vertices()),
    m_min_time(-FT(1)), m_max_time(-FT(1))
  { }

  std::size_t propagate(std::size_t k) {
    std::size_t num_events = 0;

    for (Vertex &v : vts) {
      if (v.moving)
        v.k = static_cast<int>(k);
    }

    filter_touching_pieces();

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


  // The vertices bounding a polygon on an iedge are registered per support plane and face, as both faces adjacent to an
  // iedge in a support plane may hold a polygon.
  std::pair<std::size_t, std::size_t> edge_key(std::size_t v) const {
    return std::make_pair(vts[v].sp_idx, std::size_t(vts[v].face));
  }

  // Registers the pair (a, b) for the face of a on edge e. Returns false if the face was already registered.
  bool register_on_edge(IEdge e, std::size_t a, std::size_t b) {
    CGAL_assertion(a != std::size_t(-1) && b != std::size_t(-1));
    CGAL_assertion(vts[a].sp_idx == vts[b].sp_idx);
    CGAL_assertion(vts[a].face == vts[b].face);
    return m_data.igraph().edge(e).vertices.insert(std::make_pair(edge_key(a), std::make_pair(a, b))).second;
  }

  // Entry of the face of vertex v on edge e.
  auto edge_vertices(IEdge e, std::size_t v) {
    return m_data.igraph().edge(e).vertices.find(edge_key(v));
  }

  // Replaces vertex old_v by new_v in the entry of the face of new_v on edge e (if registered).
  void replace_on_edge(IEdge e, std::size_t old_v, std::size_t new_v) {
    auto it = edge_vertices(e, new_v);
    if (it == m_data.igraph().edge(e).vertices.end())
      return;
    if (it->second.first == old_v)
      it->second.first = new_v;
    if (it->second.second == old_v)
      it->second.second = new_v;
  }

  // Zero-area pieces are not added to the mesh by the initializer (no border yet), see filter_touching_pieces().
  bool is_degenerate_piece(Support_plane& sp, IFace f) {
    CGAL_assertion(sp.data().polygons.find(f) != sp.data().polygons.end());
    const bool degenerate = sp.data().borders.find(f) == sp.data().borders.end();
    CGAL_assertion(degenerate == m_data.is_zero_area(sp.data().polygons[f]));
    return degenerate;
  }

  // Does the polygon registered with ov on an iedge cross the iedge at t = 0, i.e., does it have pieces with an area
  // on both sides? A polygon only touching the line (zero-area piece on one side) does not block other polygons at t = 0.
  bool crosses_at_init(std::size_t sp_other, const std::pair<std::size_t, std::size_t>& ov) {
    Support_plane& sp = m_data.support_plane(sp_other);
    for (std::size_t v : { ov.first, ov.second }) {
      if (vts[v].other == std::size_t(-1))
        return false;
      if (sp.data().polygons.find(vts[v].face) == sp.data().polygons.end() || is_degenerate_piece(sp, vts[v].face))
        return false;
      if (sp.data().polygons.find(vts[vts[v].other].face) == sp.data().polygons.end() || is_degenerate_piece(sp, vts[vts[v].other].face))
        return false;
    }
    return true;
  }

  // A polygon touching a line at t = 0 is split at t = 0+ and has a zero-area piece beyond the line. Such a piece
  // only exists if the polygon is allowed to cross the line: like a crossing during the propagation, the number of
  // polygons crossing the touched part of the line at t = 0 is counted against k if the piece is reached from a piece
  // with an area. A piece reached from another zero-area piece (around the end of a line, e.g., the diagonal piece of
  // a corner touching two lines) requires the line to be free, like in case_C2 and case_D2.
  // The decision only depends on the initial state, thus not on the order of the events at t = 0.
  void filter_touching_pieces() {
    for (std::size_t sp_idx = 6; sp_idx < m_data.number_of_support_planes(); sp_idx++) {
      Support_plane& sp = m_data.support_plane(sp_idx);

      std::set<IFace> degenerate, alive;
      std::deque<IFace> queue;
      for (const auto& [face, poly] : sp.data().polygons) {
        if (is_degenerate_piece(sp, face))
          degenerate.insert(face);
        else {
          alive.insert(face);
          queue.push_back(face);
        }
      }
      if (degenerate.empty())
        continue;

      while (!queue.empty()) {
        IFace from = queue.front();
        queue.pop_front();
        const bool from_degenerate = degenerate.count(from) != 0;
        for (IEdge e : m_data.igraph().face(from).edges) {
          if (m_data.igraph().iedge_is_on_bbox(e))
            continue;
          IFace to = m_data.igraph().get_other_face(sp_idx, e, from);
          if (to == IFace(-1) || !degenerate.count(to) || alive.count(to))
            continue;

          // Part of the line touched by the piece.
          const std::size_t line_idx = m_data.igraph().edge(e).line;
          const IkLine_3& l = m_data.igraph().line(line_idx);
          bool first = true;
          IkFT lo = 0, hi = 0;
          for (std::size_t v : sp.data().polygons[to]) {
            IkPoint_3 p = sp.to_3d(vts[v].p0);
            if (!l.has_on(p))
              continue;
            IkFT u = l.to_vector() * (p - l.point());
            if (first || u < lo) lo = u;
            if (first || u > hi) hi = u;
            first = false;
          }
          if (first)
            continue;

          // Each support plane is counted once, even if its polygon is registered for both faces along the iedge.
          std::set<std::size_t> counted;
          for (const auto& [key, ov] : m_data.igraph().edge(e).vertices) {
            const std::size_t sp_other = key.first;
            if (sp_other == sp_idx || counted.count(sp_other) || !crosses_at_init(sp_other, ov))
              continue;
            Support_plane& osp = m_data.support_plane(sp_other);
            IkFT u1 = l.to_vector() * (osp.to_3d(vts[ov.first].p0) - l.point());
            IkFT u2 = l.to_vector() * (osp.to_3d(vts[ov.second].p0) - l.point());
            if (u2 < u1)
              std::swap(u1, u2);
            if (!(u2 < lo || hi < u1))
              counted.insert(sp_other);
          }
          const std::size_t count = counted.size();

          std::list<std::size_t>& piece = sp.data().polygons[to];
          int k = vts[piece.front()].k;
          for (std::size_t v : piece)
            k = (std::max)(k, vts[v].k);
          const bool allowed = from_degenerate ? (count == 0) : (signed(count) <= k);
          if (!allowed)
            continue;

          for (std::size_t v : piece)
            if (vts[v].moving)
              vts[v].k = k - static_cast<int>(count);
          alive.insert(to);
          queue.push_back(to);
        }
      }

      // Remove the pieces that cannot be reached.
      std::set<std::size_t> removed_vertices;
      for (IFace f : degenerate)
        if (!alive.count(f))
          for (std::size_t v : sp.data().polygons[f])
            removed_vertices.insert(v);

      if (removed_vertices.empty())
        continue;

      std::set<IEdge> edges;
      for (IFace f : degenerate)
        if (!alive.count(f)) {
          for (IEdge e : m_data.igraph().face(f).edges)
            edges.insert(e);
          for (std::size_t v : sp.data().polygons[f]) {
            if (vts[v].constraint_edge != m_data.igraph().null_iedge())
              edges.insert(vts[v].constraint_edge);
            if (vts[v].other_constraint_edge != m_data.igraph().null_iedge())
              edges.insert(vts[v].other_constraint_edge);
          }
        }

      // Registrations on the iedges referring to removed vertices are replaced by the vertices of a remaining piece
      // that are constrained to the iedge.
      for (IEdge e : edges) {
        auto& registered = m_data.igraph().edge(e).vertices;
        bool erased = false;
        for (auto it = registered.begin(); it != registered.end();) {
          if (it->first.first == sp_idx && (removed_vertices.count(it->second.first) || removed_vertices.count(it->second.second))) {
            it = registered.erase(it);
            erased = true;
          }
          else
            ++it;
        }
        if (!erased)
          continue;
        const std::size_t line_idx = m_data.igraph().edge(e).line;
        const auto& faces = m_data.igraph().edge(e).faces;
        auto fit = faces.find(sp_idx);
        if (fit == faces.end())
          continue;
        for (IFace f : { fit->second.first, fit->second.second }) {
          if (f == IFace(-1) || !alive.count(f) || registered.count(std::make_pair(sp_idx, std::size_t(f))))
            continue;
          std::vector<std::size_t> on_edge;
          for (std::size_t v : sp.data().polygons[f])
            if (vts[v].constraints.count(line_idx) && (vts[v].constraint_edge == e || vts[v].other_constraint_edge == e))
              on_edge.push_back(v);
          if (on_edge.empty())
            continue;
          CGAL_assertion(on_edge.size() == 2);
          registered.insert(std::make_pair(std::make_pair(sp_idx, std::size_t(f)), std::make_pair(on_edge.front(), on_edge.back())));
        }
      }

      for (std::size_t v : removed_vertices) {
        if (vts[v].other != std::size_t(-1) && !removed_vertices.count(vts[v].other))
          vts[vts[v].other].other = std::size_t(-1);
        vts[v].other = std::size_t(-1);
        vts[v].moving = false;
      }

      for (IFace f : degenerate)
        if (!alive.count(f))
          sp.data().polygons.erase(f);
    }

    // The kept pieces are added to the mesh (deferred by the initializer).
    for (std::size_t sp_idx = 6; sp_idx < m_data.number_of_support_planes(); sp_idx++) {
      Support_plane& sp = m_data.support_plane(sp_idx);
      for (const auto& [face, poly] : sp.data().polygons) {
        if (sp.data().borders.find(face) != sp.data().borders.end())
          continue;
        m_data.init_border(face, true);
      }
    }
  }


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
    for (auto& [key, ov] : m_data.igraph().edge(edge).vertices) {
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
          }
          else
            CGAL_assertion(CGAL::is_negative((v2.p0 - v1.p0) * v1.v));
        }
        else {
          CGAL_assertion(v2.moving);
          CGAL_assertion(v2.v != IkVector_2(0, 0));
          if (v1.p0 == v2.p0) {
            CGAL_assertion(v1.t_init == v2.t_init);
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
      if (others.first.first == v.sp_idx)
        continue;

      CGAL_assertion(others.second.first != std::size_t(-1));
      CGAL_assertion(others.second.second != std::size_t(-1));

      Vertex& v1 = vts[others.second.first];
      Vertex& v2 = vts[others.second.second];

      // optimization possible: provide pair of vertices, instead of checking each alone
      // check below if vertex is in front or behind benefits.

      // Which other vertex has opposing direction to this vertex?
      Support_plane& sp = m_data.support_plane(v.sp_idx);
      Support_plane& sp_other = m_data.support_plane(others.first.first);

      if (v1.moving && v_idx < others.second.first) {
        IkVector_3 v1v3d = sp_other.exact_plane().base1() * v1.v[0] + sp_other.exact_plane().base2() * v1.v[1];
        IkVector_2 v1v = IkVector_2(v1v3d * sp.exact_plane().base1(), v1v3d * sp.exact_plane().base2());
        IkFT f1 = v.v * v1v;

        CGAL_assertion_code(
          IkVector_3 v2v3d=sp_other.exact_plane().base1()*v2.v[0]+sp_other.exact_plane().base2()*v2.v[1];
          IkVector_2 v2v=IkVector_2(v2v3d*sp.exact_plane().base1(),v2v3d*sp.exact_plane().base2()););

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
        CGAL_assertion_code(
        IkVector_3 v1v3d = sp_other.exact_plane().base1() * v1.v[0] + sp_other.exact_plane().base2() * v1.v[1];
        IkVector_2 v1v = IkVector_2(v1v3d * sp.exact_plane().base1(), v1v3d * sp.exact_plane().base2());
        IkFT f1 = v.v * v1v;);
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

        // If the vertex hits the border edge of its face, this edge is crossed. This also holds if the hit point is an end
        // of the edge, i.e., if the vertex hits several lines at once: it is treated as hitting this line first.
        IkFT u_s = m_data.get_u_on_line(m_data.igraph().source(e), line_idx);
        IkFT u_t = m_data.get_u_on_line(m_data.igraph().target(e), line_idx);
        if (u_t < u_s)
          std::swap(u_s, u_t);
        if (u_s <= ce.u && ce.u <= u_t)
          ce.e = e;
        else
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

        if ((adj == std::size_t(-1) || vts[adj].cached_events.empty()) && v.other != std::size_t(-1))
          adj = has_adjacent_with_same_target(v.other);

        /// Does it have a constrained neighbor with the same target (= C1 event)
        if (adj != std::size_t(-1) && !vts[adj].cached_events.empty()) {
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
        CGAL_assertion(v.itarget != std::size_t(-1));
        IkPoint_2 target = sp.data().exact_plane.to_2d(m_data.igraph().point_3(v.itarget));
        IkFT t1 = intersection_time(v, target);
        // The target may be reached immediately, if the vertex was created by a vertex hitting several lines at once.
        CGAL_assertion(t1 >= time);

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
    // A support plane is counted once, even if its polygon is registered for both faces along the iedge.
    std::set<std::size_t> crossing;
    for (auto& [key, ov] : m_data.igraph().edge(edge).vertices) {
      const std::size_t sp_other = key.first;
      if (sp_other == v.sp_idx || crossing.count(sp_other))
        continue;

      Support_plane &sp = m_data.support_plane(sp_other);

      Vertex& v1 = vts[ov.first];
      Vertex& v2 = vts[ov.second];

      if (e.destination != std::size_t(-1)) {
        if (vts[ov.first].ivertex == e.destination || vts[ov.second].ivertex == e.destination)
          crossing.insert(sp_other);
      }
      else {
        IkFT u_low = l.to_vector() * (sp.to_3d(v1.p0 + (e.time - v1.t_init) * v1.v) - l.point());
        IkFT u_high = l.to_vector() * (sp.to_3d(v2.p0 + (e.time - v2.t_init) * v2.v) - l.point());
        if (u_high < u_low)
          std::swap(u_high, u_low);

        if (u_low <= e.u && e.u <= u_high)
          crossing.insert(sp_other);
      }
    }

    if (!decrement)
      return crossing.empty();

    if (signed(crossing.size()) <= v.k) {
      v.k -= static_cast<int>(crossing.size());
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
    for (auto& [key, ov] : m_data.igraph().edge(edge).vertices) {
      const std::size_t sp_other = key.first;
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
    for (auto& [key, ov] : m_data.igraph().edge(edge).vertices) {
      const std::size_t sp_other = key.first;
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
      sp.set_piece_polygon(other_face, sp.piece_polygon(vts[e.vertex].face));

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

      // The polygon in the entered face is bounded by the paired vertices on the crossed edge.
      register_on_edge(e.crossed_edge, vts[prev].other, vts[next].other);

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
    register_on_edge(e.crossed_edge, prev, next);
    //CGAL_assertion(inserted);

    vts[e.vertex].stop(e);

    for (std::size_t v : new_moving_vertices)
      calculate_events(v, e.time);

    CGAL_assertion_code(check_vertices_on_edge(e.crossed_edge);)
  }

  void case_B(const std::vector<Event>& events, std::size_t adj) {
    b_events++;
    const Event& e = events.front();
    if (events.size() > 1)
      t_events++;
    CGAL_assertion(adj != std::size_t(-1));
    Support_plane& sp = m_data.support_plane(vts[e.vertex].sp_idx);

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
    auto border_it = edge_vertices(vts[adj].constraint_edge, e.vertex);
    CGAL_assertion(border_it != m_data.igraph().edge(vts[adj].constraint_edge).vertices.end());
    border_it->second.first = *it_before_adj;
    border_it->second.second = *it;

    std::size_t other = vts[adj].other;
    std::size_t other_face;
    CGAL_assertion(other == std::size_t(-1) || vts[other].other == adj);
    if (other == std::size_t(-1)) {
      // Check if there is a vertex on the other side of vts[adj];
      other_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, vts[adj].constraint_edge, vts[e.vertex].face);
      typename std::map<IFace, std::list<std::size_t> >::iterator it;
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

        // The new vertex replaces the former one in the registration of the other face on the edge.
        replace_on_edge(vts.back().constraint_edge, *it, vts.size() - 1);

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

    bool still_moving = false;
    std::list<std::size_t>& poly = sp.data().polygons[v1.face];
    for (std::size_t i : poly)
      if (vts[i].moving) {
        still_moving = true;
        break;
      }

    if (still_moving) {
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
        auto border_v1 = edge_vertices(v1.constraint_edge, e.vertex);
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
        auto border_v2 = edge_vertices(v2.constraint_edge, events[1].vertex);
        CGAL_assertion(border_v2 != m_data.igraph().edge(v2.constraint_edge).vertices.end());
        border_v2->second.first = *it_adj;
        border_v2->second.second = events[1].vertex;
      }
    }
  }

  // Creates the piece of the polygon of the event vertex in the face `face` around the ivertex e.destination, which the
  // polygon reaches at the event. The face is bounded at the ivertex by the edges ea and eb. The input polygon the piece originates from is translated such that its centroid is at the event and
  // is clipped by the lines of ea and eb. The piece consists of a stationary vertex at the ivertex and of a vertex moving
  // along each edge. Returns the indices of the stationary vertex and of the vertices moving along ea and eb.
  std::array<std::size_t, 3> create_wedge_piece(const Event& e, IFace face, IEdge ea, IEdge eb, std::vector<std::size_t>& new_moving_vertices) {
    const std::size_t sp_idx = vts[e.vertex].sp_idx;
    Support_plane& sp = m_data.support_plane(sp_idx);
    const IVertex d = e.destination;
    const std::size_t la = m_data.igraph().edge(ea).line;
    const std::size_t lb = m_data.igraph().edge(eb).line;
    CGAL_assertion(la != lb);

    m_data.init_border(face);
    std::list<std::size_t>& poly = sp.data().polygons[face];
    CGAL_assertion(poly.empty());
    sp.set_piece_polygon(face, sp.piece_polygon(vts[e.vertex].face));

    const IVertex ta = m_data.igraph().other(ea, d);
    const IVertex tb = m_data.igraph().other(eb, d);
    // A point inside the face close to the ivertex.
    const IkPoint_2 inside = CGAL::centroid(e.p, sp.to_2d(m_data.igraph().point_3(ta)), sp.to_2d(m_data.igraph().point_3(tb)));

    const auto& input = sp.input_polygon(sp.piece_polygon(vts[e.vertex].face));
    std::vector<Vertex> vertices;
    vertices.reserve(input.vertices.size());
    const IkVector_2 shift = (e.p - input.centroid);
    for (std::size_t v = 0; v < input.vertices.size(); v++)
      vertices.emplace_back(sp_idx, input.vertices[v] + shift, input.velocities[v]);

    std::list<std::size_t> polygon(input.vertices.size()), other;
    std::iota(polygon.begin(), polygon.end(), 0);
    m_data.split_polygon(vertices, polygon, la, other);
    if (sp.data().lines[la].oriented_side(inside) == CGAL::ON_POSITIVE_SIDE)
      std::swap(polygon, other);
    other.clear();
    m_data.split_polygon(vertices, polygon, lb, other);
    if (sp.data().lines[lb].oriented_side(inside) == CGAL::ON_POSITIVE_SIDE)
      std::swap(polygon, other);

    std::array<std::size_t, 3> res = { std::size_t(-1), std::size_t(-1), std::size_t(-1) };
    for (std::size_t i : polygon) {
      if (vertices[i].constraints.empty())
        continue;
      vts.emplace_back(vertices[i]);
      const std::size_t idx = vts.size() - 1;
      Vertex& v = vts.back();
      v.t_init = e.time;
      v.p0 = e.p;
      v.face = face;
      v.other = std::size_t(-1);
      if (!v.moving) {
        CGAL_assertion(res[0] == std::size_t(-1));
        CGAL_assertion(v.constraints.size() == 2);
        v.constraint_edge = ea;
        v.other_constraint_edge = eb;
        v.ivertex = d;
        res[0] = idx;
      }
      else {
        CGAL_assertion(v.constraints.size() == 1);
        v.k = vts[e.vertex].k;
        if (*v.constraints.begin() == la) {
          CGAL_assertion(res[1] == std::size_t(-1));
          v.constraint_edge = ea;
          v.itarget = ta;
          res[1] = idx;
        }
        else {
          CGAL_assertion(*v.constraints.begin() == lb);
          CGAL_assertion(res[2] == std::size_t(-1));
          v.constraint_edge = eb;
          v.itarget = tb;
          res[2] = idx;
        }
        new_moving_vertices.push_back(idx);
      }
      poly.push_back(idx);
    }
    CGAL_assertion(res[0] != std::size_t(-1) && res[1] != std::size_t(-1) && res[2] != std::size_t(-1));

    register_on_edge(ea, res[0], res[1]);
    register_on_edge(eb, res[0], res[2]);
    return res;
  }

  // Propagates the polygon of the event vertex around the ivertex e.destination, starting with the face `face` entered via
  // entry_edge. Used if more than two lines intersect in the ivertex. The polygon enters the faces around the ivertex one
  // after the other, as long as no polygon of another support plane blocks the crossed edge.
  // The vertex `paired` moves along entry_edge in the face the polygon comes from.
  void propagate_around_ivertex(const Event& e, IEdge entry_edge, IFace face, std::size_t paired, std::vector<std::size_t>& new_moving_vertices) {
    const std::size_t sp_idx = vts[e.vertex].sp_idx;
    Support_plane& sp = m_data.support_plane(sp_idx);
    const IFace start = vts[e.vertex].face;
    bool first = true;
    while (true) {
      const IEdge exit_edge = other_face_edge(face, entry_edge, e.destination);
      CGAL_assertion(exit_edge != m_data.igraph().null_iedge());
      const auto piece = create_wedge_piece(e, face, entry_edge, exit_edge, new_moving_vertices);
      if (first) {
        // The vertex entering the face moves like its mirrored vertex.
        vts[piece[1]].v = vts[paired].v;
        first = false;
      }
      if (paired != std::size_t(-1)) {
        vts[piece[1]].other = paired;
        vts[paired].other = piece[1];
      }

      const IFace next = m_data.igraph().get_other_face(sp_idx, exit_edge, face);
      if (next == std::size_t(-1) || next == start || sp.data().polygons.find(next) != sp.data().polygons.end())
        break;
      if (!could_cross(sp_idx, face, e.destination, e.time, exit_edge))
        break;

      paired = piece[2];
      entry_edge = exit_edge;
      face = next;
    }
  }

  // A constrained vertex reaches an ivertex at the same time as its unconstrained neighbor, i.e., the edge between both
  // vertices collapses in the ivertex. The unconstrained vertex is removed and the event is handled like case_C2 for the
  // constrained vertex. The polygon edge of the removed vertex is used to compute the speed of the ricochet vertex.
  bool handle_collapse_at_ivertex(const std::vector<Event>& connected, const std::vector<std::size_t>& vertices) {
    if (vertices.size() != 2)
      return false;
    std::size_t c = vertices[0], f = vertices[1];
    if (vts[c].constraints.empty())
      std::swap(c, f);
    if (vts[c].constraints.empty() || !vts[f].constraints.empty())
      return false;

    std::vector<Event> ce;
    const Event* fe = nullptr;
    for (const Event& ev : connected)
      if (ev.vertex == c)
        ce.push_back(ev);
      else
        fe = &ev;
    if (ce.size() != 1 || fe == nullptr || ce[0].destination == m_data.igraph().null_ivertex() || fe->p != ce[0].p)
      return false;

    Support_plane& sp = m_data.support_plane(vts[c].sp_idx);
    std::list<std::size_t>& poly = sp.data().polygons[vts[c].face];
    auto it = std::find(poly.begin(), poly.end(), c);
    CGAL_assertion(it != poly.end());
    auto prev = (it == poly.begin()) ? std::prev(poly.end()) : std::prev(it);
    auto next = (std::next(it) == poly.end()) ? poly.begin() : std::next(it);
    std::list<std::size_t>::iterator fit;
    if (*prev == f)
      fit = prev;
    else if (*next == f)
      fit = next;
    else
      return false;

    const Vertex ref = vts[f];
    poly.erase(fit);
    vts[f].stop(*fe);
    case_C2(ce, &ref);
    return true;
  }

  void case_C2(const std::vector<Event>& events, const Vertex* edge_ref = nullptr) {
    c2_events++;
    const Event& e = events.front();

    if (events.size() > 2)
      t_events++;

    Support_plane &sp = m_data.support_plane(vts[e.vertex].sp_idx);
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

    dir = sp.calculate_edge_speed(edge_ref ? *edge_ref : vts[e.vertex], vts[*adjacent], e.p, sp.data().lines[hit_line_idx], e.time);
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
    register_on_edge(hit_edge, e.vertex, ricochet);

    // vertices in already covered edge property don't need adjustment as only indices are stored

    // Expanding to other side of edge is still necessary (check if there is a polygon in that face)

    std::vector<std::size_t> new_moving_vertices;
    new_moving_vertices.push_back(ricochet);

    std::size_t corner = -1;
    std::size_t prolongation = -1;

    const bool crossing = can_cross(e, hit_edge);
    IFace other_face = crossing ? m_data.igraph().get_other_face(vts[e.vertex].sp_idx, hit_edge, vts[e.vertex].face) : IFace(-1);
    const std::size_t constraint_line_idx = m_data.igraph().edge(vts[e.vertex].constraint_edge).line;
    if (crossing && m_data.igraph().edge(other_face_edge(other_face, hit_edge, e.destination)).line != constraint_line_idx) {
      // More than two lines intersect in the ivertex.
      propagate_around_ivertex(e, hit_edge, other_face, ricochet, new_moving_vertices);
    }
    else if (crossing) {
      // Initialize new face
      m_data.init_border(other_face);
      std::list<std::size_t>& new_poly = sp.data().polygons[other_face];
      CGAL_assertion(new_poly.empty());
      sp.set_piece_polygon(other_face, sp.piece_polygon(vts[e.vertex].face));

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
      // The polygon in the entered face is bounded on the hit edge by the corner and the paired ricochet vertex.
      CGAL_assertion(vts[ricochet].constraint_edge == hit_edge);
      register_on_edge(hit_edge, corner, vts.size() - 1);

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

      register_on_edge(vts.back().constraint_edge, vts.size() - 3, vts.size() - 1);

      prolongation = vts.size() - 1;

      new_poly.push_back(corner);
      new_poly.push_back(vts.size() - 2);
      new_poly.push_back(vts.size() - 1);

      new_moving_vertices.push_back(vts.size() - 2);
      new_moving_vertices.push_back(vts.size() - 1);

      // Check whether a polygon needs to be created on other side of edge
      // Edge for prolongation vertex is needed anyway
      std::size_t q3_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, vts.back().constraint_edge, other_face);
      if (q3_face != std::size_t(-1) && sp.data().polygons.find(q3_face) == sp.data().polygons.end() && could_cross(vts[e.vertex].sp_idx, other_face, e.destination, e.time, vts[prolongation].constraint_edge)) {
        // Create new polygon with prolongation vertex and ivertex
        m_data.init_border(q3_face);
        std::list<std::size_t>& poly_q3 = sp.data().polygons[q3_face];
        CGAL_assertion(poly_q3.empty());
        sp.set_piece_polygon(q3_face, sp.piece_polygon(vts[e.vertex].face));

        // Insert vertices into q3
        // The input polygon the piece originates from, translated such that its centroid is at the event.
        const auto& input = sp.input_polygon(sp.piece_polygon(vts[e.vertex].face));
        std::list<std::size_t> polygon(input.vertices.size());
        std::vector<Vertex> vertices;
        std::iota(polygon.begin(), polygon.end(), 0);

        vertices.reserve(input.vertices.size());

        IkVector_2 shift = (e.p - input.centroid);

        for (std::size_t v = 0; v < input.vertices.size(); v++)
          vertices.emplace_back(vts[e.vertex].sp_idx, input.vertices[v] + shift, input.velocities[v]);

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

        register_on_edge(vts[opposite_ricochet].constraint_edge, q3_stationary, opposite_ricochet);
        if (q3_prolongation != std::size_t(-1))
          register_on_edge(vts[q3_prolongation].constraint_edge, q3_stationary, q3_prolongation);

        std::size_t q4_face = m_data.igraph().get_other_face(vts[e.vertex].sp_idx, vts[e.vertex].constraint_edge, vts[e.vertex].face);
        if (q4_face != std::size_t(-1) && sp.data().polygons.find(q4_face) == sp.data().polygons.end() && could_cross(vts[e.vertex].sp_idx, q3_face, e.destination, e.time, vts[opposite_ricochet].constraint_edge)) {
          // Create new polygon with prolongation vertex and ivertex
          m_data.init_border(q4_face);
          std::list<std::size_t>& poly_q4 = sp.data().polygons[q4_face];
          CGAL_assertion(poly_q4.empty());
          sp.set_piece_polygon(q4_face, sp.piece_polygon(vts[e.vertex].face));

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
          if (q4_stationary != std::size_t(-1) && q4_opposite_ricochet != std::size_t(-1))
            register_on_edge(vts[q4_opposite_ricochet].constraint_edge, q4_stationary, q4_opposite_ricochet);
          if (q4_stationary != std::size_t(-1) && q4_inverse != std::size_t(-1))
            register_on_edge(vts[q4_inverse].constraint_edge, q4_stationary, q4_inverse);
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
    auto vts_it = edge_vertices(vts[e.vertex].constraint_edge, e.vertex);
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

    // Insert vertices into q3
    // The input polygon the piece originates from, translated such that its centroid is at the event.
    const auto& input = sp.input_polygon(sp.piece_polygon(vts[e.vertex].face));
    std::list<std::size_t> polygon(input.vertices.size());
    std::vector<Vertex> vertices;
    std::iota(polygon.begin(), polygon.end(), 0);

    vertices.reserve(input.vertices.size());

    IkVector_2 shift = (e.p - input.centroid);

    for (std::size_t v = 0; v < input.vertices.size(); v++)
      vertices.emplace_back(vts[e.vertex].sp_idx, input.vertices[v] + shift, input.velocities[v]);

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
    sp.set_piece_polygon(q3_face, sp.piece_polygon(vts[e.vertex].face));

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

    register_on_edge(vts[opposite_ricochet].constraint_edge, q3_stationary, opposite_ricochet);
    if (q3_prolongation != std::size_t(-1))
      register_on_edge(vts[q3_prolongation].constraint_edge, q3_stationary, q3_prolongation);


    auto it = sp.data().polygons.find(q4_face);
    if (it != sp.data().polygons.end())
      return;

    m_data.init_border(q4_face);
    std::list<std::size_t>& poly_q4 = sp.data().polygons[q4_face];
    CGAL_assertion(poly_q4.empty());
    sp.set_piece_polygon(q4_face, sp.piece_polygon(vts[e.vertex].face));

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
    if (q4_stationary != std::size_t(-1) && q4_opposite_ricochet != std::size_t(-1))
      register_on_edge(vts[q4_opposite_ricochet].constraint_edge, q4_stationary, q4_opposite_ricochet);
    if (q4_stationary != std::size_t(-1) && q4_inverse != std::size_t(-1))
      register_on_edge(vts[q4_inverse].constraint_edge, q4_stationary, q4_inverse);
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
    // If more than two lines intersect in the ivertex reached by a constrained vertex, the face beyond the hit edge is
    // not bounded by the prolongation of the constraint line at the ivertex, but by another line (edge m_e1/m_e2).
    IEdge m_e1 = m_data.igraph().null_iedge(), m_e2 = m_data.igraph().null_iedge();
    if (propagates_front) {
      const IFace front_face = m_data.igraph().get_other_face(vts[e1.vertex].sp_idx, hit_edge, vts[e1.vertex].face);
      auto fan_edge = [&](const Event& ev) {
        const IEdge prolongation = get_next_edge(vts[ev.vertex].constraint_edge, ev.destination, *vts[ev.vertex].constraints.begin());
        const IEdge bounding = other_face_edge(front_face, hit_edge, ev.destination);
        return (bounding != prolongation) ? bounding : m_data.igraph().null_iedge();
      };
      if (!vts[e1.vertex].constraints.empty())
        m_e1 = fan_edge(e1);
      if (!vts[e2.vertex].constraints.empty())
        m_e2 = fan_edge(e2);
    }
    if (propagates_front && !vts[e1.vertex].constraints.empty() && m_e1 == m_data.igraph().null_iedge())
      propagates_e1_side = expands_laterally(e1, hit_line_idx, hit_edge); // Check lateral propagation
    if (propagates_front && !vts[e2.vertex].constraints.empty() && m_e2 == m_data.igraph().null_iedge())
      propagates_e2_side = expands_laterally(e2, hit_line_idx, hit_edge); // Check lateral propagation
    std::size_t fan_slider_e1 = std::size_t(-1), fan_slider_e2 = std::size_t(-1);
    IFace front_face_idx = std::size_t(-1);

    // Local propagation
    Support_plane& sp = m_data.support_plane(vts[e1.vertex].sp_idx);
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
    CGAL_assertion(adj_e1 != std::size_t(-1));
    CGAL_assertion(adj_e2 != std::size_t(-1));

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
      if (inserted_adj_e1 == std::size_t(-1))
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
      sp.set_piece_polygon(other_face, sp.piece_polygon(vts[e1.vertex].face));

      front_face_idx = other_face;
      // Create adjacent of e1.vertex prolongation
      if (m_e1 != m_data.igraph().null_iedge()) {
        // The polygon is bounded at the ivertex by the hit line and another line through the ivertex.
        const std::size_t m_line = m_data.igraph().edge(m_e1).line;
        vts.emplace_back(vts[e1.vertex].sp_idx, e1.p, IkVector_2(0, 0), e1.time);
        vts.back().moving = false;
        vts.back().constraints.insert(m_line);
        vts.back().constraints.insert(hit_line_idx);
        vts.back().face = other_face;
        vts.back().ivertex = e1.destination;
        vts.back().constraint_edge = hit_edge;
        vts.back().other_constraint_edge = m_e1;
        new_poly.push_back(vts.size() - 1);

        vts.emplace_back(vts[e1.vertex].sp_idx, e1.p, sp.calculate_edge_speed(vts[e1.vertex], vts[e2.vertex], e1.p, sp.data().lines[m_line], e1.time), e1.time);
        vts.back().constraints.insert(m_line);
        vts.back().face = other_face;
        vts.back().k = vts[e1.vertex].k;
        vts.back().constraint_edge = m_e1;
        vts.back().itarget = m_data.igraph().other(m_e1, e1.destination);
        new_moving_vertices.push_back(vts.size() - 1);
        new_poly.push_back(vts.size() - 1);
        register_on_edge(m_e1, vts.size() - 2, vts.size() - 1);
        fan_slider_e1 = vts.size() - 1;
      }
      else if (!vts[e1.vertex].constraints.empty()) {
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
        CGAL_assertion(inserted_adj_e1 != std::size_t(-1));
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
      if (m_e1 == m_data.igraph().null_iedge()) {
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
        register_on_edge(vts.back().constraint_edge, vts.size() - 2, vts.size() - 1);
      }

      // prolongation e2.vertex
      if (m_e2 != m_data.igraph().null_iedge()) {
        const std::size_t m_line = m_data.igraph().edge(m_e2).line;
        vts.emplace_back(vts[e2.vertex].sp_idx, e2.p, sp.calculate_edge_speed(vts[e1.vertex], vts[e2.vertex], e2.p, sp.data().lines[m_line], e2.time), e2.time);
        vts.back().constraints.insert(m_line);
        vts.back().face = other_face;
        vts.back().k = vts[e2.vertex].k;
        vts.back().constraint_edge = m_e2;
        vts.back().itarget = m_data.igraph().other(m_e2, e2.destination);
        new_moving_vertices.push_back(vts.size() - 1);
        new_poly.push_back(vts.size() - 1);
        fan_slider_e2 = vts.size() - 1;

        vts.emplace_back(vts[e2.vertex].sp_idx, e2.p, IkVector_2(0, 0), e2.time);
        vts.back().moving = false;
        vts.back().constraints.insert(m_line);
        vts.back().constraints.insert(hit_line_idx);
        vts.back().face = other_face;
        vts.back().ivertex = e2.destination;
        vts.back().constraint_edge = hit_edge;
        vts.back().other_constraint_edge = m_e2;
        new_poly.push_back(vts.size() - 1);
        register_on_edge(m_e2, vts.size() - 2, vts.size() - 1);
      }
      else {
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
        register_on_edge(vts.back().other_constraint_edge, vts.size() - 2, vts.size() - 1);
      }
      else {
        CGAL_assertion(inserted_adj_e2 != std::size_t(-1));
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

      // The front polygon is bounded on the hit edge by its first and last vertex.
      register_on_edge(hit_edge, new_poly.front(), new_poly.back());
    }

    // Propagation around an ivertex in which more than two lines intersect.
    for (int side = 0; side < 2; side++) {
      const Event& ev = (side == 0) ? e1 : e2;
      const IEdge m_edge = (side == 0) ? m_e1 : m_e2;
      const std::size_t slider = (side == 0) ? fan_slider_e1 : fan_slider_e2;
      if (m_edge == m_data.igraph().null_iedge())
        continue;
      const IFace next = m_data.igraph().get_other_face(vts[ev.vertex].sp_idx, m_edge, front_face_idx);
      if (next == std::size_t(-1) || next == vts[ev.vertex].face || sp.data().polygons.find(next) != sp.data().polygons.end())
        continue;
      if (!could_cross(vts[ev.vertex].sp_idx, front_face_idx, ev.destination, ev.time, m_edge))
        continue;
      propagate_around_ivertex(ev, m_edge, next, slider, new_moving_vertices);
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
    CGAL_assertion(inserted_adj_e1 != std::size_t(-1));

    if (!vts[e2.vertex].constraints.empty())
      inserted_adj_e2 = e2.vertex;
    CGAL_assertion(inserted_adj_e2 != std::size_t(-1));

    register_on_edge(hit_edge, inserted_adj_e1, inserted_adj_e2);

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

  // Line of the support plane, the 2D line is created on first use.
  const IkLine_2& sp_line(Support_plane& sp, std::size_t line_idx) {
    auto it = sp.data().lines.find(line_idx);
    if (it == sp.data().lines.end()) {
      const auto& vertices = m_data.igraph().vertices_on_line(line_idx);
      it = sp.data().lines.insert(std::make_pair(line_idx, IkLine_2(sp.to_2d(m_data.igraph().point_3(vertices[0].second)), sp.to_2d(m_data.igraph().point_3(vertices.back().second))))).first;
    }
    return it->second;
  }

  // Edge of face f incident to ivertex x that is not e.
  IEdge other_face_edge(IFace f, IEdge e, IVertex x) const {
    for (IEdge fe : m_data.igraph().face(f).edges)
      if (fe != e && (m_data.igraph().source(fe) == x || m_data.igraph().target(fe) == x))
        return fe;
    return m_data.igraph().null_iedge();
  }

  bool ivertex_on_line(IVertex v, std::size_t line_idx) const {
    const auto& vol = m_data.igraph().vertices_on_line(line_idx);
    return std::find_if(vol.begin(), vol.end(), [&](const std::pair<IkFT, IVertex>& p) { return p.second == v; }) != vol.end();
  }

  // Does event e make its vertex reach the line line_idx?
  bool reaches_line(const Event& e, std::size_t line_idx) const {
    if (vts[e.vertex].constraints.empty())
      return e.crossed_edge != m_data.igraph().null_iedge() && m_data.igraph().edge(e.crossed_edge).line == line_idx;
    return e.destination != m_data.igraph().null_ivertex() && ivertex_on_line(e.destination, line_idx);
  }

  // Consecutive corners c_0, ..., c_{m-1} of a polygon reaching corners of their face at the same time: the polygon
  // edges between them lie on the lines of the face. Edge j connects corner j-1 and corner j and lies on line L_j.
  // In an open chain, edge 0 connects the neighbor a with c_0 and edge m connects c_{m-1} with the neighbor b,
  // a and b reach L_0 and L_m. In a closed chain, all vertices of the polygon are corners (the polygon fits its face).
  struct Corner_chain {
    bool closed = false;
    bool forward = true; // the chain follows the order of the polygon
    std::vector<std::size_t> corners;
    std::vector<std::pair<Event, Event> > corner_events; // events of corner i on L_i and on L_{i+1}
    std::vector<std::size_t> lines;                      // L_j
    std::size_t a = std::size_t(-1), b = std::size_t(-1);
    Event ea, eb;
  };

  // Extracts a chain of corners from coupled events at the same time; events that do not belong to it are returned in rest.
  bool build_corner_chain(const std::vector<Event>& events, Corner_chain& ch, std::vector<Event>& rest) {
    std::map<std::size_t, std::vector<std::size_t> > by_vertex;
    for (std::size_t i = 0; i < events.size(); i++)
      by_vertex[events[i].vertex].push_back(i);

    auto corner_lines = [&](std::size_t v, std::size_t& l1, std::size_t& l2) {
      auto it = by_vertex.find(v);
      if (it == by_vertex.end() || it->second.size() != 2 || !vts[v].constraints.empty())
        return false;
      const Event& e1 = events[it->second[0]];
      const Event& e2 = events[it->second[1]];
      if (e1.crossed_edge == m_data.igraph().null_iedge() || e2.crossed_edge == m_data.igraph().null_iedge())
        return false;
      l1 = m_data.igraph().edge(e1.crossed_edge).line;
      l2 = m_data.igraph().edge(e2.crossed_edge).line;
      return l1 != l2;
    };
    auto shared_line = [&](std::size_t p, std::size_t q, std::size_t& l) {
      std::size_t p1, p2, q1, q2;
      if (!corner_lines(p, p1, p2) || !corner_lines(q, q1, q2))
        return false;
      if (p1 == q1 || p1 == q2) { l = p1; return true; }
      if (p2 == q1 || p2 == q2) { l = p2; return true; }
      return false;
    };
    auto event_on = [&](std::size_t v, std::size_t line_idx) -> std::size_t {
      auto it = by_vertex.find(v);
      if (it == by_vertex.end())
        return std::size_t(-1);
      for (std::size_t i : it->second)
        if (reaches_line(events[i], line_idx))
          return i;
      return std::size_t(-1);
    };

    std::size_t c0 = std::size_t(-1), l1, l2, l;
    for (const auto& [v, idx] : by_vertex)
      if (corner_lines(v, l1, l2)) {
        c0 = v;
        break;
      }
    if (c0 == std::size_t(-1))
      return false;

    const std::list<std::size_t>& poly = m_data.support_plane(vts[c0].sp_idx).data().polygons[vts[c0].face];
    const std::vector<std::size_t> pv(poly.begin(), poly.end());
    const std::size_t n = pv.size();
    const std::size_t pos0 = std::find(pv.begin(), pv.end(), c0) - pv.begin();

    std::size_t s = pos0, steps = 0;
    while (steps < n && shared_line(pv[(s + n - 1) % n], pv[s], l)) {
      s = (s + n - 1) % n;
      steps++;
    }

    ch = Corner_chain();
    if (steps >= n) {
      ch.closed = true;
      for (std::size_t i = 0; i < n; i++)
        ch.corners.push_back(pv[(pos0 + i) % n]);
      for (std::size_t j = 0; j < n; j++) {
        shared_line(ch.corners[(j + n - 1) % n], ch.corners[j], l);
        ch.lines.push_back(l);
      }
    }
    else {
      ch.corners.push_back(pv[s]);
      std::size_t e = s;
      while (ch.corners.size() < n && shared_line(pv[e], pv[(e + 1) % n], l)) {
        e = (e + 1) % n;
        ch.corners.push_back(pv[e]);
      }
      ch.a = pv[(s + n - 1) % n];
      ch.b = pv[(e + 1) % n];
      if (ch.a == ch.b || std::find(ch.corners.begin(), ch.corners.end(), ch.a) != ch.corners.end()
        || std::find(ch.corners.begin(), ch.corners.end(), ch.b) != ch.corners.end())
        return false;

      const std::size_t m = ch.corners.size();
      ch.lines.resize(m + 1);
      for (std::size_t j = 1; j < m; j++)
        shared_line(ch.corners[j - 1], ch.corners[j], ch.lines[j]);
      corner_lines(ch.corners[0], l1, l2);
      if (m > 1)
        ch.lines[0] = (l1 == ch.lines[1]) ? l2 : l1;
      corner_lines(ch.corners[m - 1], l1, l2);
      if (m > 1)
        ch.lines[m] = (l1 == ch.lines[m - 1]) ? l2 : l1;
      if (m == 1) {
        ch.lines[0] = l1;
        ch.lines[1] = l2;
        if (event_on(ch.a, ch.lines[0]) == std::size_t(-1) || event_on(ch.b, ch.lines[1]) == std::size_t(-1))
          std::swap(ch.lines[0], ch.lines[1]);
      }
      const std::size_t ia = event_on(ch.a, ch.lines[0]), ib = event_on(ch.b, ch.lines[m]);
      if (ia == std::size_t(-1) || ib == std::size_t(-1))
        return false;
      ch.ea = events[ia];
      ch.eb = events[ib];

      // A single corner is processed as before (a is the neighbor with the lower index).
      if (m == 1 && ch.b < ch.a) {
        ch.forward = false;
        std::swap(ch.a, ch.b);
        std::swap(ch.ea, ch.eb);
        std::swap(ch.lines[0], ch.lines[1]);
      }
    }

    const std::size_t m = ch.corners.size();
    const std::size_t ne = ch.lines.size();
    for (std::size_t i = 0; i < m; i++) {
      const std::size_t il = event_on(ch.corners[i], ch.lines[i]), ir = event_on(ch.corners[i], ch.lines[(i + 1) % ne]);
      if (il == std::size_t(-1) || ir == std::size_t(-1) || il == ir)
        return false;
      ch.corner_events.push_back(std::make_pair(events[il], events[ir]));
    }

    // The other events are processed afterwards.
    rest.clear();
    for (const Event& e : events) {
      bool used = (!ch.closed && ((e.vertex == ch.a && e.time == ch.ea.time && e.crossed_edge == ch.ea.crossed_edge && e.destination == ch.ea.destination)
        || (e.vertex == ch.b && e.time == ch.eb.time && e.crossed_edge == ch.eb.crossed_edge && e.destination == ch.eb.destination)));
      used = used || std::find(ch.corners.begin(), ch.corners.end(), e.vertex) != ch.corners.end();
      if (!used)
        rest.push_back(e);
    }
    return true;
  }

  // Coupled events forming a chain of corners (see Corner_chain). Other events of the group are queued again.
  bool handle_corner_chain(const std::vector<Event>& events) {
    Corner_chain ch;
    std::vector<Event> rest;
    if (!build_corner_chain(events, ch, rest))
      return false;
    for (const Event& e : rest) {
      vts[e.vertex].queued_events++;
      m_event_queue.push(e);
    }
    case_corner_chain(ch);
    return true;
  }

  void case_D2(const std::vector<Event>& events) {
    if (!handle_corner_chain(events))
      KSP::internal::Log_line(m_data.log_prefix()) << "Error: Unknown event type!";
  }

  // Corners of a polygon reach corners X_i of their face at the same time (case D2 for one corner).
  // Locally, the polygon stops in its face F: the corners become fixed vertices at X_i, the end vertices a and b are
  // replaced by sliders moving away from the chain (or stop at the ends of their constraint edges if constrained).
  // The polygon may continue beyond each line L_j (face G_j, decided with can_cross(), i.e., counted against k) and
  // beyond both lines at each corner (diagonal face H_i, only through free lines like q3 in case C2). A face G_j can
  // also be reached from a diagonal face if L_j is blocked directly (like q4 in case C2).
  //
  //            G_0 L_1'| H_0               G_0: [X_0, a on L_0, a prolongation, slider on L_1']
  //                    |                   G_1: [X_0, slider on L_0', b prolongation, b on L_1]
  //     L_0 --a-------X_0-------- L_0'     H_0: [X_0, slider on L_1', c_0 prolongation, slider on L_0']
  //            F       |   G_1
  //                    b L_1
  //
  // Between two corners, G_j is [X_j, X_{j-1}, slider at X_{j-1}, slider at X_j].
  void case_corner_chain(const Corner_chain& ch) {
    d2_events++;
    const std::size_t m = ch.corners.size();
    const bool open = !ch.closed;
    const std::size_t ne = ch.lines.size();
    const std::size_t none = std::size_t(-1);
    const IEdge null_edge = m_data.igraph().null_iedge();
    const std::size_t sp_idx = vts[ch.corners[0]].sp_idx;
    const IFace F = vts[ch.corners[0]].face;
    Support_plane& sp = m_data.support_plane(sp_idx);
    const IkFT t = ch.corner_events[0].first.time;
    const std::vector<std::size_t>& L = ch.lines;

    // Edge j connects the corners left_corner(j) and right_corner(j) (none for a and b); corner i lies between edge i and edge i + 1.
    auto left_corner = [&](std::size_t j) { return open ? (j == 0 ? none : j - 1) : (j + m - 1) % m; };
    auto right_corner = [&](std::size_t j) { return (open && j == m) ? none : j; };
    auto right_edge = [&](std::size_t i) { return (i + 1) % ne; };

    // The crossed edges of the corners may have been computed in a former face of the vertex. The edges of the
    // current face meeting in the corners are used instead.
    std::vector<IVertex> X(m, m_data.igraph().null_ivertex());
    std::vector<IEdge> E(ne, null_edge);
    std::vector<IkPoint_2> pX(m);
    for (std::size_t i = 0; i < m; i++) {
      pX[i] = ch.corner_events[i].first.p;
      for (IEdge f1 : m_data.igraph().face(F).edges) {
        if (m_data.igraph().edge(f1).line != L[i])
          continue;
        for (IEdge f2 : m_data.igraph().face(F).edges) {
          if (m_data.igraph().edge(f2).line != L[right_edge(i)])
            continue;
          for (IVertex v : { m_data.igraph().source(f1), m_data.igraph().target(f1) })
            if (v == m_data.igraph().source(f2) || v == m_data.igraph().target(f2)) {
              CGAL_assertion(sp.to_2d(m_data.igraph().point_3(v)) == pX[i]);
              X[i] = v;
              CGAL_assertion(E[i] == null_edge || E[i] == f1);
              CGAL_assertion(E[right_edge(i)] == null_edge || E[right_edge(i)] == f2);
              E[i] = f1;
              E[right_edge(i)] = f2;
            }
        }
      }
      CGAL_assertion(X[i] != m_data.igraph().null_ivertex());
    }

    // End vertices of an open chain.
    const std::size_t a = ch.a, b = ch.b;
    const bool a_constrained = open && !vts[a].constraints.empty();
    const bool b_constrained = open && !vts[b].constraints.empty();
    CGAL_assertion(!a_constrained || m_data.igraph().other(E[0], X[0]) == ch.ea.destination);
    CGAL_assertion(!b_constrained || m_data.igraph().other(E[m], X[m - 1]) == ch.eb.destination);

    std::list<std::size_t>& poly = sp.data().polygons[F];
    auto cyc_next = [&](std::list<std::size_t>::iterator it) { return (std::next(it) == poly.end()) ? poly.begin() : std::next(it); };
    auto cyc_prev = [&](std::list<std::size_t>::iterator it) { return (it == poly.begin()) ? std::prev(poly.end()) : std::prev(it); };
    // Neighbor of an end vertex that is not in the chain.
    auto outer_neighbor = [&](std::size_t v, std::size_t corner) {
      auto it = std::find(poly.begin(), poly.end(), v);
      return (*cyc_prev(it) == corner) ? *cyc_next(it) : *cyc_prev(it);
    };

    // Faces and edges around the corners.
    std::vector<IFace> G(ne, IFace(-1));
    std::vector<IEdge> SL(ne, null_edge), SR(ne, null_edge); // edges of G_j at its left and right corner
    for (std::size_t j = 0; j < ne; j++) {
      if (m_data.igraph().iedge_is_on_bbox(E[j]))
        continue;
      G[j] = m_data.igraph().get_other_face(sp_idx, E[j], F);
      if (G[j] == IFace(-1))
        continue;
      if (left_corner(j) != none)
        SL[j] = other_face_edge(G[j], E[j], X[left_corner(j)]);
      if (right_corner(j) != none)
        SR[j] = other_face_edge(G[j], E[j], X[right_corner(j)]);
    }

    // Decide the propagation before changing anything.
    std::vector<bool> crossG(ne, false);
    for (std::size_t j = 0; j < ne; j++) {
      std::vector<const Event*> ends;
      if (left_corner(j) != none)
        ends.push_back(&ch.corner_events[left_corner(j)].second);
      if (right_corner(j) != none)
        ends.push_back(&ch.corner_events[right_corner(j)].first);
      if (open && j == 0)
        ends.push_back(&ch.ea);
      if (open && j == m)
        ends.push_back(&ch.eb);
      for (const Event* e : ends)
        if (can_cross(*e, E[j])) {
          crossG[j] = true;
          break;
        }
    }
    const bool lateral_a = open && crossG[0] && a_constrained && expands_laterally(ch.ea, L[0], E[0]);
    const bool lateral_b = open && crossG[m] && b_constrained && expands_laterally(ch.eb, L[m], E[m]);

    std::vector<IFace> H(m, IFace(-1));
    std::vector<bool> crossH(m, false);
    for (std::size_t i = 0; i < m; i++) {
      const std::size_t jl = i, jr = right_edge(i);
      if (SR[jl] == null_edge || SL[jr] == null_edge || m_data.igraph().edge(SR[jl]).line != L[jr] || m_data.igraph().edge(SL[jr]).line != L[jl]
        || m_data.igraph().iedge_is_on_bbox(SR[jl]) || m_data.igraph().iedge_is_on_bbox(SL[jr]))
        continue;
      H[i] = m_data.igraph().get_other_face(sp_idx, SR[jl], G[jl]);
      CGAL_assertion(H[i] == m_data.igraph().get_other_face(sp_idx, SL[jr], G[jr]));
    }
    for (bool changed = true; changed;) {
      changed = false;
      for (std::size_t i = 0; i < m; i++) {
        const std::size_t jl = i, jr = right_edge(i);
        if (H[i] == IFace(-1) || crossH[i] || sp.data().polygons.find(H[i]) != sp.data().polygons.end())
          continue;
        if ((crossG[jl] && could_cross(sp_idx, G[jl], X[i], t, SR[jl])) || (crossG[jr] && could_cross(sp_idx, G[jr], X[i], t, SL[jr]))) {
          crossH[i] = true;
          changed = true;
        }
      }
      for (std::size_t j = 0; j < ne; j++) {
        if (crossG[j] || G[j] == IFace(-1) || sp.data().polygons.find(G[j]) != sp.data().polygons.end())
          continue;
        const std::size_t lc = left_corner(j), rc = right_corner(j);
        if ((lc != none && crossH[lc] && could_cross(sp_idx, H[lc], X[lc], t, SL[j])) || (rc != none && crossH[rc] && could_cross(sp_idx, H[rc], X[rc], t, SR[j]))) {
          crossG[j] = true;
          changed = true;
        }
      }
    }

    // Speeds (computed before any vertex is stopped).
    std::vector<IkVector_2> vc(m), vL(ne, IkVector_2(0, 0)), vR(ne, IkVector_2(0, 0));
    for (std::size_t i = 0; i < m; i++)
      vc[i] = vts[ch.corners[i]].v;
    const IkVector_2 va = open ? vts[a].v : IkVector_2(0, 0), vb = open ? vts[b].v : IkVector_2(0, 0);
    IkVector_2 v_s1(0, 0), v_s2(0, 0);
    if (open && !a_constrained)
      v_s1 = sp.calculate_edge_speed(vts[a], vts[outer_neighbor(a, ch.corners[0])], ch.ea.p, sp_line(sp, L[0]), t);
    if (open && !b_constrained)
      v_s2 = sp.calculate_edge_speed(vts[b], vts[outer_neighbor(b, ch.corners[m - 1])], ch.eb.p, sp_line(sp, L[m]), t);
    for (std::size_t j = 0; j < ne; j++) {
      const std::size_t lc = left_corner(j), rc = right_corner(j);
      if (rc != none && SR[j] != null_edge && (crossG[j] || crossH[rc]))
        vR[j] = sp.calculate_edge_speed(vts[ch.corners[rc]], vts[lc != none ? ch.corners[lc] : a], pX[rc], sp_line(sp, m_data.igraph().edge(SR[j]).line), t);
      if (lc != none && SL[j] != null_edge && (crossG[j] || crossH[lc]))
        vL[j] = sp.calculate_edge_speed(vts[ch.corners[lc]], vts[rc != none ? ch.corners[rc] : b], pX[lc], sp_line(sp, m_data.igraph().edge(SL[j]).line), t);
    }

    std::vector<std::size_t> new_moving_vertices;

    auto add_fixed = [&](const IkPoint_2& p, IFace face, IVertex iv, IEdge ce, IEdge oce) {
      vts.emplace_back(sp_idx, p, IkVector_2(0, 0), t);
      vts.back().moving = false;
      vts.back().face = face;
      vts.back().ivertex = iv;
      vts.back().constraints.insert(m_data.igraph().edge(ce).line);
      vts.back().constraints.insert(m_data.igraph().edge(oce).line);
      vts.back().constraint_edge = ce;
      vts.back().other_constraint_edge = oce;
      return vts.size() - 1;
    };

    // Slider starting at ivertex 'from' (or inside edge ce) moving away from it on edge ce.
    auto add_slider = [&](const IkPoint_2& p, const IkVector_2& v, IFace face, IEdge ce, IVertex from, std::size_t twin) {
      vts.emplace_back(sp_idx, p, v, t);
      vts.back().face = face;
      vts.back().constraints.insert(m_data.igraph().edge(ce).line);
      vts.back().constraint_edge = ce;
      vts.back().itarget = m_data.igraph().other(ce, from);
      CGAL_assertion(vts.back().itarget == m_data.get_target_ivertex(vts.back(), ce));
      if (twin != none) {
        vts.back().other = twin;
        vts[twin].other = vts.size() - 1;
      }
      new_moving_vertices.push_back(vts.size() - 1);
      return vts.size() - 1;
    };

    // Continuation of an unconstrained vertex in a new face.
    auto add_free_prolongation = [&](std::size_t v_idx, const IkPoint_2& p, const IkVector_2& v, IFace face) {
      vts.emplace_back(sp_idx, p, v, t);
      vts.back().face = face;
      vts.back().k = vts[v_idx].k;
      std::swap(vts.back().known_intersections, vts[v_idx].known_intersections);
      std::swap(vts.back().cached_events, vts[v_idx].cached_events);
      new_moving_vertices.push_back(vts.size() - 1);
      return vts.size() - 1;
    };

    // Continuation of a constrained vertex beyond its destination: fixed vertex at the destination and moving vertex on the next edge.
    auto add_constrained_prolongation = [&](std::size_t v_idx, const Event& e, const IkVector_2& v, IFace face, IEdge hit_edge,
      std::size_t& fixed, std::size_t& moving) {
      std::size_t line_idx = *vts[v_idx].constraints.begin();
      IEdge prolongation_edge = get_next_edge(vts[v_idx].constraint_edge, e.destination, line_idx);
      fixed = add_fixed(e.p, face, e.destination, hit_edge, prolongation_edge);
      vts.emplace_back(sp_idx, e.p, v, t);
      vts.back().face = face;
      vts.back().k = vts[v_idx].k;
      vts.back().constraints.insert(line_idx);
      vts.back().constraint_edge = prolongation_edge;
      vts.back().itarget = m_data.igraph().other(prolongation_edge, e.destination);
      for (Cached_event& ce : vts[v_idx].cached_events)
        if (ce.d != e.destination)
          vts.back().cached_events.push_back(ce);
      moving = vts.size() - 1;
      new_moving_vertices.push_back(moving);
      register_on_edge(prolongation_edge, fixed, moving);
    };

    // Vertices are given in the order of the chain, starting with a fixed vertex.
    auto new_polygon = [&](IFace face, std::vector<std::size_t> seq) {
      m_data.init_border(face);
      std::list<std::size_t>& p = sp.data().polygons[face];
      CGAL_assertion(p.empty());
      sp.set_piece_polygon(face, sp.piece_polygon(F));
      if (!ch.forward)
        std::reverse(seq.begin(), seq.end());
      p.assign(seq.begin(), seq.end());
    };

    // Face F: a and b are replaced by sliders moving away from the chain if they are unconstrained.
    std::size_t s1 = none, s2 = none;
    if (open && !a_constrained) {
      s1 = add_slider(ch.ea.p, v_s1, F, E[0], X[0], none);
      *std::find(poly.begin(), poly.end(), a) = s1;
    }
    if (open && !b_constrained) {
      s2 = add_slider(ch.eb.p, v_s2, F, E[m], X[m - 1], none);
      *std::find(poly.begin(), poly.end(), b) = s2;
    }
    for (std::size_t j = 0; j < ne; j++) {
      const std::size_t lc = left_corner(j), rc = right_corner(j);
      const std::size_t left = (lc != none) ? ch.corners[lc] : (a_constrained ? a : s1);
      const std::size_t right = (rc != none) ? ch.corners[rc] : (b_constrained ? b : s2);
      if (rc != none && lc == none)
        register_on_edge(E[j], right, left);
      else
        register_on_edge(E[j], left, right);
    }

    // Faces G_j beyond the lines.
    std::vector<std::size_t> L_up(ne, none), R_up(ne, none);
    for (std::size_t j = 0; j < ne; j++) {
      if (!crossG[j])
        continue;
      const std::size_t lc = left_corner(j), rc = right_corner(j);
      std::size_t L_on = none, R_on = none;
      if (rc != none)
        R_on = add_fixed(pX[rc], G[j], X[rc], E[j], SR[j]);
      if (lc != none)
        L_on = add_fixed(pX[lc], G[j], X[lc], E[j], SL[j]);

      if (lc != none)
        L_up[j] = add_slider(pX[lc], vL[j], G[j], SL[j], X[lc], none);
      else if (a_constrained)
        add_constrained_prolongation(a, ch.ea, va, G[j], E[j], L_on, L_up[j]);
      else {
        L_on = add_slider(ch.ea.p, v_s1, G[j], E[j], X[rc], s1);
        L_up[j] = add_free_prolongation(a, ch.ea.p, va, G[j]);
      }

      if (rc != none)
        R_up[j] = add_slider(pX[rc], vR[j], G[j], SR[j], X[rc], none);
      else if (b_constrained)
        add_constrained_prolongation(b, ch.eb, vb, G[j], E[j], R_on, R_up[j]);
      else {
        R_on = add_slider(ch.eb.p, v_s2, G[j], E[j], X[lc], s2);
        R_up[j] = add_free_prolongation(b, ch.eb.p, vb, G[j]);
      }

      if (lc != none)
        register_on_edge(SL[j], L_on, L_up[j]);
      if (rc != none)
        register_on_edge(SR[j], R_on, R_up[j]);
      register_on_edge(E[j], L_on, R_on);

      if (rc != none)
        new_polygon(G[j], { R_on, L_on, L_up[j], R_up[j] });
      else
        new_polygon(G[j], { L_on, L_up[j], R_up[j], R_on });
    }

    // Diagonal faces H_i.
    for (std::size_t i = 0; i < m; i++) {
      if (!crossH[i])
        continue;
      const std::size_t jl = i, jr = right_edge(i);
      std::size_t x = add_fixed(pX[i], H[i], X[i], SR[jl], SL[jr]);
      std::size_t s_l = add_slider(pX[i], vR[jl], H[i], SR[jl], X[i], R_up[jl]);
      std::size_t c_prol = add_free_prolongation(ch.corners[i], pX[i], vc[i], H[i]);
      std::size_t s_r = add_slider(pX[i], vL[jr], H[i], SL[jr], X[i], L_up[jr]);
      register_on_edge(SR[jl], x, s_l);
      register_on_edge(SL[jr], x, s_r);
      new_polygon(H[i], { x, s_l, c_prol, s_r });
    }

    // Constrained a/b also continue on their constraint lines beyond the face of the corner (like in D1).
    if (lateral_a)
      propagate_laterally(ch.ea, ch.corner_events[0].first, E[0], new_moving_vertices);
    if (lateral_b)
      propagate_laterally(ch.eb, ch.corner_events[m - 1].second, E[m], new_moving_vertices);

    // Stop the event vertices.
    for (std::size_t i = 0; i < m; i++) {
      Vertex& c = vts[ch.corners[i]];
      c.stop(ch.corner_events[i].first);
      c.ivertex = X[i];
      c.constraints.insert(L[i]);
      c.constraints.insert(L[right_edge(i)]);
      c.constraint_edge = E[i];
      c.other_constraint_edge = E[right_edge(i)];
    }
    if (open) {
      vts[a].stop(ch.ea);
      vts[b].stop(ch.eb);
      if (a_constrained) {
        vts[a].constraints.insert(L[0]);
        vts[a].other_constraint_edge = E[0];
      }
      if (b_constrained) {
        vts[b].constraints.insert(L[m]);
        vts[b].other_constraint_edge = E[m];
      }
    }

    for (std::size_t v : new_moving_vertices)
      calculate_events(v, t);

    CGAL_assertion_code(for (IEdge e : E) check_vertices_on_edge(e);)
  }

  // Lines an event hits: the crossed line for an unconstrained vertex or the other lines through the destination
  // in the support plane for a constrained vertex.
  std::set<std::size_t> hit_lines(const Event& e) const {
    std::set<std::size_t> lines;
    const Vertex& v = vts[e.vertex];
    if (v.constraints.empty()) {
      if (e.crossed_edge != m_data.igraph().null_iedge())
        lines.insert(m_data.igraph().edge(e.crossed_edge).line);
    }
    else if (e.destination != m_data.igraph().null_ivertex()) {
      for (const IEdge& edge : m_data.incident_iedges(e.destination)) {
        const auto& planes = m_data.intersected_planes(edge);
        if (planes.find(v.sp_idx) == planes.end())
          continue;
        std::size_t line_idx = m_data.igraph().edge(edge).line;
        if (v.constraints.find(line_idx) == v.constraints.end())
          lines.insert(line_idx);
      }
    }
    return lines;
  }

  bool are_adjacent(std::size_t a, std::size_t b) const {
    const std::list<std::size_t>& poly = m_data.support_plane(vts[a].sp_idx).data().polygons.at(vts[a].face);
    auto it = std::find(poly.begin(), poly.end(), a);
    CGAL_assertion(it != poly.end());
    auto prev = (it == poly.begin()) ? std::prev(poly.end()) : std::prev(it);
    auto next = (std::next(it) == poly.end()) ? poly.begin() : std::next(it);
    return *prev == b || *next == b;
  }

  // Two events at the same time are coupled if they belong to the same vertex or to adjacent vertices
  // of the same polygon that reach the same destination or hit a common line.
  bool are_coupled(const Event& a, const Event& b) const {
    if (a.vertex == b.vertex)
      return true;
    if (vts[a.vertex].sp_idx != vts[b.vertex].sp_idx || vts[a.vertex].face != vts[b.vertex].face)
      return false;
    if (!are_adjacent(a.vertex, b.vertex))
      return false;
    if (!vts[a.vertex].constraints.empty() && !vts[b.vertex].constraints.empty() && a.destination == b.destination)
      return true;
    const std::set<std::size_t> la = hit_lines(a), lb = hit_lines(b);
    for (std::size_t l : la)
      if (lb.count(l))
        return true;
    return false;
  }

  // Extracts the first event and all events coupled with it (directly or indirectly) from the list of events at the same time.
  std::vector<Event> get_connected_events(std::list<Event> &events, std::vector<std::size_t> &vertices) {
    std::vector<Event> connected;
    connected.push_back(events.front());
    events.pop_front();

    bool added = true;
    while (added) {
      added = false;
      for (auto it = events.begin(); it != events.end();) {
        bool coupled = false;
        for (const Event& c : connected)
          if (are_coupled(*it, c)) {
            coupled = true;
            break;
          }
        if (coupled) {
          connected.push_back(*it);
          it = events.erase(it);
          added = true;
        }
        else
          ++it;
      }
    }

    std::set<std::size_t> unique_v;
    for (const Event& e : connected)
      unique_v.insert(e.vertex);

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
        if (handle_collapse_at_ivertex(connected, vertices))
          break;
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
        // Several corners of a polygon reach lines at the same time, e.g., a polygon fitting exactly into its face.
        if (!handle_corner_chain(connected))
          KSP::internal::Log_line(m_data.log_prefix()) << "Error: Unknown event type!";
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
