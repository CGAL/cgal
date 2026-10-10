// Copyright (c) 2023 GeometryFactory Sarl (France).
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

#ifndef CGAL_KSP_3_INITIALIZER_H
#define CGAL_KSP_3_INITIALIZER_H

#include <CGAL/license/Kinetic_space_partition.h>

// CGAL includes.
#include <CGAL/Timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/min_quadrilateral_2.h>

// Internal includes.
#include <CGAL/KSP/utils.h>
#include <CGAL/KSP/debug.h>
#include <CGAL/KSP/parameters.h>

#include <CGAL/KSP_3/Data_structure.h>

#include <CGAL/Real_timer.h>

namespace CGAL {
namespace KSP_3 {
namespace internal {

#ifdef DOXYGEN_RUNNING
#else

template<typename GeomTraits, typename IntersectionKernel>
class Initializer {

public:
  using Kernel = GeomTraits;
  using Intersection_kernel = IntersectionKernel;

private:
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using Vector_2 = typename Kernel::Vector_2;
  using Segment_2 = typename Kernel::Segment_2;
  using Segment_3 = typename Kernel::Segment_3;
  using Line_2 = typename Kernel::Line_2;
  using Transform_3 = CGAL::Aff_transformation_3<Kernel>;
  using Direction_2 = typename Kernel::Direction_2;
  using IkPoint_2 = typename Intersection_kernel::Point_2;
  using IkPoint_3 = typename Intersection_kernel::Point_3;
  using IkVector_2 = typename Intersection_kernel::Vector_2;
  using IkLine_2 = typename Intersection_kernel::Line_2;
  using IkFT = typename Intersection_kernel::FT;

  using Data_structure = KSP_3::internal::Data_structure<Kernel, Intersection_kernel>;
  using Support_plane = typename Data_structure::Support_plane;
  using Vertex = typename Support_plane::Vertex;
  using IEdge = typename Data_structure::IEdge;
  using IFace = typename Data_structure::IFace;
  using Face_property = typename Data_structure::Intersection_graph::Face_property;
  using Intersection_graph = typename Data_structure::Intersection_graph;
  using IEdge_set = typename Data_structure::IEdge_set;

  using IVertex = typename Data_structure::IVertex;

  using To_exact = CGAL::Cartesian_converter<Kernel, Intersection_kernel>;
  using From_exact = CGAL::Cartesian_converter<Intersection_kernel, Kernel>;

  using Bbox_3 = CGAL::Bbox_3;

  using Parameters = KSP::internal::Parameters_3<FT>;

  using Timer = CGAL::Real_timer;

public:
  Initializer(std::vector<std::vector<Point_3> >& input_polygons,
    std::vector<typename Intersection_kernel::Plane_3>& input_planes,
    Data_structure& data, const Parameters& parameters, bool snap)
    : m_input_polygons(input_polygons), m_input_planes(input_planes), m_data(data), m_parameters(parameters), m_snap(snap)
  {}

  Initializer(std::vector<std::vector<IkPoint_3> >& input_polygons,
    std::vector<typename Intersection_kernel::Plane_3>& input_planes, Data_structure& data, const Parameters& parameters, bool snap)
    : m_input_polygons(input_polygons), m_input_planes(input_planes), m_data(data), m_parameters(parameters), m_snap(snap)
  {}

  void initialize(const std::array<typename Intersection_kernel::Point_3, 8>& bbox, std::vector<std::size_t>& input_polygons) {
    Timer timer;
    timer.reset();
    timer.start();

    std::vector< std::vector<typename Intersection_kernel::Point_3> > bbox_faces;
    bounding_box_to_polygons(bbox, bbox_faces);
    add_polygons(bbox_faces, input_polygons);

    m_data.igraph().finished_bbox();

    if (m_parameters.verbose)
      std::cout << "* intersecting input polygons ... ";

    double t1 = timer.time();

    // Fills in the ivertices on support plane intersections inside the bbox.
    create_internal_vertices();
    double t2 = timer.time();

    if (m_parameters.debug)
      KSP_3::internal::dump(m_data, m_data.prefix() + "intersected");

    // Generation of ifaces
    create_ifaces();
    double t3 = timer.time();

    // Splitting the input polygons along intersection lines.
    initial_polygon_splitting();

    double t4 = timer.time();

    create_bbox_meshes();

    double t5 = timer.time();

    // Starting from here the intersection graph is const, it won't change anymore.
    if (m_parameters.verbose)
      std::cout << "done" << std::endl;

    CGAL_assertion(m_data.check_bbox());

    m_data.initialization_done();

    if (m_parameters.debug) {
      for (std::size_t sp = 0; sp < m_data.number_of_support_planes(); sp++)
        dump_2d_surface_mesh(m_data, sp, m_data.prefix() + "before-partition-sp" + std::to_string(sp));
    }

    if (m_parameters.verbose) {
      double t6 = timer.time();
      std::cout << "inserting polygons: " << t1 << std::endl;
      std::cout << "internal vertices:  " << (t2 - t1) << std::endl;
      std::cout << "creating ifaces:    " << (t3 - t2) << std::endl;
      std::cout << "polygon splitting:  " << (t4 - t3) << std::endl;
      std::cout << "bbox meshes:        " << (t5 - t4) << std::endl;
      std::cout << "finishing up:       " << (t6 - t5) << std::endl;
    }
  }

  void clear() {
    // to be added
  }

private:
  std::vector<std::vector<IkPoint_3> >& m_input_polygons;
  std::vector<typename Intersection_kernel::Plane_3>& m_input_planes;
  Data_structure& m_data;
  const Parameters& m_parameters;
  bool m_snap;

  void add_iface_from_iedge(std::size_t sp_idx, IEdge edge, IEdge next, bool cw) {
    IVertex s = m_data.source(edge);
    IVertex t = m_data.target(edge);

    IFace face_idx = m_data.add_iface(sp_idx);
    Face_property& face = m_data.igraph().face(face_idx);
    face.pts.push_back(m_data.support_plane(sp_idx).to_2d(m_data.igraph().point_3(s)));
    face.pts.push_back(m_data.support_plane(sp_idx).to_2d(m_data.igraph().point_3(t)));
    face.vertices.push_back(s);
    face.vertices.push_back(t);
    face.edges.push_back(edge);
    m_data.igraph().add_face(sp_idx, edge, face_idx);

    face.edges.push_back(next);
    m_data.igraph().add_face(sp_idx, next, face_idx);

    std::size_t iterations = 0;

    int dir = (cw) ? -1 : 1;
    const std::size_t uninitialized = static_cast<std::size_t>(-1);
    std::size_t inext;
    while (s != m_data.target(next) && iterations < 10000) {
      face.vertices.push_back(m_data.target(next));
      face.pts.push_back(m_data.support_plane(sp_idx).to_2d(m_data.igraph().point_3(m_data.target(next))));

      IEdge enext, eprev;
      get_prev_next(sp_idx, next, eprev, enext);

      std::vector<std::pair<IEdge, typename Intersection_kernel::Direction_2> > connected;
      m_data.get_and_sort_all_connected_iedges(sp_idx, m_data.target(next), connected);
      inext = uninitialized;
      for (std::size_t idx = 0; idx < connected.size(); idx++) {
        if (connected[idx].first == next) {
          inext = (idx + dir + connected.size()) % connected.size();
          break;
        }
      }
      CGAL_assertion(inext != uninitialized);

      next = connected[inext].first;
      face.edges.push_back(next);
      m_data.igraph().add_face(sp_idx, next, face_idx);

      iterations++;
    }

    // Loop complete, connecting face with all edges.
    for (IEdge edge : face.edges) {
      m_data.support_plane(sp_idx).add_neighbor(edge, face_idx);
      CGAL_assertion_code(IFace f1 = m_data.support_plane(sp_idx).iface(edge);)
      CGAL_assertion_code(IFace f2 = m_data.support_plane(sp_idx).other(edge, f1);)
      CGAL_assertion(f1 == face_idx || f2 == face_idx);
    }

    std::vector<typename Intersection_kernel::Point_2> pts;
    pts.reserve(face.pts.size());
    for (auto p : face.pts)
      pts.push_back(p);

    face.poly = Polygon_2<Intersection_kernel>(pts.begin(), pts.end());

    if (face.poly.orientation() != CGAL::COUNTERCLOCKWISE) {
      face.poly.reverse_orientation();
      std::reverse(face.pts.begin(), face.pts.end());
      std::reverse(face.vertices.begin(), face.vertices.end());
      std::reverse(face.edges.begin(), face.edges.end());
    }

    CGAL_assertion(face.poly.orientation() == CGAL::COUNTERCLOCKWISE);
    CGAL_assertion(face.poly.is_convex());
    CGAL_assertion(face.poly.is_simple());
  }

  void get_prev_next(std::size_t sp_idx, IEdge edge, IEdge& prev, IEdge& next) {
    CGAL_assertion(edge != Intersection_graph::null_iedge());
    CGAL_assertion(sp_idx != static_cast<std::size_t>(-1));

    std::vector<std::pair<IEdge, typename Intersection_kernel::Direction_2> > connected;
    m_data.get_and_sort_all_connected_iedges(sp_idx, m_data.target(edge), connected);

    std::size_t inext = static_cast<std::size_t>(-1), iprev = static_cast<std::size_t>(-1);
    for (std::size_t idx = 0; idx < connected.size(); idx++) {
      if (connected[idx].first == edge) {
        iprev = (idx - 1 + connected.size()) % connected.size();
        inext = (idx + 1) % connected.size();
        break;
      }
    }

    CGAL_assertion(inext != static_cast<std::size_t>(-1));
    CGAL_assertion(iprev != static_cast<std::size_t>(-1));
    prev = connected[iprev].first;
    next = connected[inext].first;
  }

  void create_ifaces() {
    for (std::size_t sp_idx = 0; sp_idx < m_data.number_of_support_planes(); sp_idx++) {
      const IEdge_set& iedges = m_data.support_plane(sp_idx).iedges();

      // Special case bbox without splits
      if (sp_idx < 6 && iedges.size() == 4) {

        IEdge next, prev;
        get_prev_next(sp_idx, *iedges.begin(), prev, next);
        add_iface_from_iedge(sp_idx, *iedges.begin(), prev, true);

        // Get first edge
        /*IEdge first = *iedges.begin();
        IEdge edge = first;
        IVertex s = m_data.source(edge);
        IVertex t = m_data.target(edge);

        // Create single IFace for unsplit bbox face
        IFace face_idx = m_data.add_iface(sp_idx);
        Face_property& face = m_data.igraph().face(face_idx);

        // Add first edge, vertices and points to face properties
        face.pts.push_back(m_data.support_plane(sp_idx).to_2d(m_data.igraph().point_3(s)));
        face.pts.push_back(m_data.support_plane(sp_idx).to_2d(m_data.igraph().point_3(t)));
        face.vertices.push_back(s);
        face.vertices.push_back(t);
        face.edges.push_back(edge);

        // Link edge and face
        m_data.igraph().add_face(sp_idx, edge, face_idx);

        // Walk around bbox face
        while (s != t) {
          auto inc_iedges = m_data.incident_iedges(t);
          for (auto next : inc_iedges) {
            // Filter edges that are not in this bbox face
            const auto iplanes = m_data.intersected_planes(next);
            if (iplanes.find(sp_idx) == iplanes.end()) {
              continue;
            }

            // Skip current edge
            if (edge == next)
              continue;

            // The only left edge is the next one.
            edge = next;
            break;
          }
          t = (m_data.target(edge) == t) ? m_data.source(edge) : m_data.target(edge);
          if (s != t) {
            face.vertices.push_back(t);
            face.pts.push_back(m_data.support_plane(sp_idx).to_2d(m_data.igraph().point_3(t)));
            face.edges.push_back(edge);
            m_data.igraph().add_face(sp_idx, edge, face_idx);
          }
        }*/

        // create polygon in proper order
      }
      else {
        bool all_on_bbox = true;
        for (auto edge : iedges) {
          bool on_edge = m_data.igraph().iedge_is_on_bbox(edge);
          // If non-bbox support plane is treated, skip all edges on bbox as they only have one face.
          if (sp_idx >= 6 && on_edge)
            continue;

          // If bbox support plane is treated, skip edges on bbox edge.
          if (sp_idx < 6 && m_data.igraph().line_is_bbox_edge(m_data.line_idx(edge)))
            continue;

          all_on_bbox = false;

          IFace n1 = m_data.support_plane(sp_idx).iface(edge);
          IFace n2 = m_data.support_plane(sp_idx).other(edge, n1);
          if (n1 != Intersection_graph::null_iface() && n2 != Intersection_graph::null_iface())
            continue;

          Face_property np1, np2;
          if (n1 != Intersection_graph::null_iface())
            np1 = m_data.igraph().face(n1);

          if (n2 != Intersection_graph::null_iface())
            np2 = m_data.igraph().face(n2);

          IEdge next, prev;
          get_prev_next(sp_idx, edge, prev, next);

          // Check if cw face already exists.
          bool skip = false;
          if (n1 != Intersection_graph::null_iface()) {
            if (np1.is_part(edge, next))
              skip = true;
          }

          if (!skip && n2 != Intersection_graph::null_iface()) {
            if (np2.is_part(edge, next))
              skip = true;
          }

          if (!skip) {
            add_iface_from_iedge(sp_idx, edge, next, false);
          }

          // Check if cw face already exists.
          skip = false;
          if (n1 != Intersection_graph::null_iface()) {
            if (np1.is_part(edge, prev))
              skip = true;
          }

          if (!skip && n2 != Intersection_graph::null_iface()) {
            if (np2.is_part(edge, prev))
              skip = true;
          }

          if (!skip) {
            add_iface_from_iedge(sp_idx, edge, prev, true);
          }
        }

        // Special case if the input polygon only intersects with the bbox.
        if (all_on_bbox) {
          IEdge next, prev;
          get_prev_next(sp_idx, *iedges.begin(), prev, next);
          add_iface_from_iedge(sp_idx, *iedges.begin(), prev, true);
        }
      }
    }
  }

  template<typename IndexRange>
  void export_poly(const IndexRange& poly, const std::string& filename) const {
    std::ofstream out(filename);
    Support_plane &sp = m_data.support_plane(m_data.kinetic_vertices()[poly.front()].sp_idx);
    out << (poly.size() + 1);
    for (const auto& i : poly) {
      const auto& v = m_data.kinetic_vertices()[i];
      IkPoint_3 p = sp.to_3d(v.p0);
      out << " " << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " " << CGAL::to_double(p.z());
    }
    IkPoint_3 p = sp.to_3d(m_data.kinetic_vertices()[poly.front()].p0);
    out << " " << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " " << CGAL::to_double(p.z()) << std::endl;

    out.close();
  }

  template<typename IndexRange>
  void export_poly_front(const IndexRange& poly, const std::string& filename) const {
    std::ofstream out(filename);

    out << (poly.size() + 1);
    for (const auto& i : poly) {
      const auto& v = m_data.kinetic_vertices()[i];
      out << " " << CGAL::to_double(v.p0.x() + (0.1 * v.v.x())) << " " << CGAL::to_double(v.p0.y() + (0.1 * v.v.y())) << " 0";
    }
    const auto& v = m_data.kinetic_vertices()[poly.front()];
    out << " " << CGAL::to_double(v.p0.x() + (0.1 * v.v.x())) << " " << CGAL::to_double(v.p0.y() + (0.1 * v.v.y())) << " 0" << std::endl;

    for (const auto& i : poly) {
      const auto& v = m_data.kinetic_vertices()[i];
      out << "2";
      out << " " << CGAL::to_double(v.p0.x()) << " " << CGAL::to_double(v.p0.y()) << " 0";
      out << " " << CGAL::to_double(v.p0.x() + (0.1 * v.v.x())) << " " << CGAL::to_double(v.p0.y() + (0.1 * v.v.y())) << " 0" << std::endl;
    }
    out.close();
  }

  template<typename IndexRange>
  void export_poly_line(const IndexRange& poly, const IkLine_2 &l, const std::string& filename) const {
    std::ofstream out(filename);
    FT fmin = FLT_MAX, fmax = -FLT_MAX;
    From_exact from_exact;
    Vector_2 dir = from_exact(l.to_vector() * (FT(1.0) / CGAL::approximate_sqrt(l.to_vector().squared_length())));

    std::vector<Vertex> &vts = m_data.kinetic_vertices();

    for (const auto& i : poly) {
      const auto& v = vts[i];
      FT t = dir * from_exact((v.p0 - l.point()));
      fmin = (std::min)(fmin, t);
      fmax = (std::max)(fmax, t);
    }
    FT fspan = fmax - fmin;
    fmin -= fspan * 0.2;
    fmax += fspan * 0.2;

    out << 2 << " " << CGAL::to_double(l.point().x() + fmin * dir.x()) << " " << CGAL::to_double(l.point().y() + fmin * dir.y()) << " 0";
    out << " " << CGAL::to_double(l.point().x() + fmax * dir.x()) << " " << CGAL::to_double(l.point().y() + fmax * dir.y()) << " 0";
    out << std::endl;

    out << (poly.size() + 1);
    for (const auto& i : poly) {
      const auto& v = vts[i];
      out << " " << CGAL::to_double(v.p0.x()) << " " << CGAL::to_double(v.p0.y()) << " 0";
    }
    out << " " << CGAL::to_double(vts[poly.front()].p0.x()) << " " << CGAL::to_double(vts[poly.front()].p0.y()) << " 0" << std::endl;

    out << (poly.size() + 1);
    for (const auto& i : poly) {
      const auto& v = vts[i];
      out << " " << CGAL::to_double(v.p0.x() + (0.1 * v.v.x())) << " " << CGAL::to_double(v.p0.y() + (0.1 * v.v.y())) << " 0";
    }
    const auto& v = vts[poly.front()];
    out << " " << CGAL::to_double(v.p0.x() + (0.1 * v.v.x())) << " " << CGAL::to_double(v.p0.y() + (0.1 * v.v.y())) << " 0" << std::endl;

    for (const auto& i : poly) {
      const auto& v = vts[i];
      out << "2";
      out << " " << CGAL::to_double(v.p0.x()) << " " << CGAL::to_double(v.p0.y()) << " 0";
      out << " " << CGAL::to_double(v.p0.x() + (0.1 * v.v.x())) << " " << CGAL::to_double(v.p0.y() + (0.1 * v.v.y())) << " 0" << std::endl;
    }

    out.close();
  }

  void initial_polygon_splitting() {
    using Vertex = typename Support_plane::Vertex;
    using Polygon = std::list<std::size_t>;
    using Polygons = std::vector<Polygon>;
    std::vector<Vertex> &vts = m_data.kinetic_vertices();

    bool plot = false;
    for (std::size_t sp_idx = 6; sp_idx < m_data.number_of_support_planes(); sp_idx++) {
      Support_plane& sp = m_data.support_plane(sp_idx);
      Polygons polygons(1);
      polygons[0].resize(sp.data().exact_vertices.size());
      std::iota(polygons[0].begin(), polygons[0].end(), vts.size());

      vts.reserve(vts.size() + sp.data().exact_vertices.size());

      // get initial polygon for support plane
      for (std::size_t v = 0; v < sp.data().exact_vertices.size(); v++)
        vts.emplace_back(sp_idx, sp.data().exact_vertices[v], sp.data().original_rays[v].to_vector());

      if (plot)
        export_poly(polygons[0], m_data.prefix() + std::to_string(sp_idx) + "_initial_polygon.polylines.txt");

      // Clip the polygon to the bbox. The sides are decided at t = 0+: a vertex on a bbox line that moves outwards
      // is replaced by vertices sliding on the bbox line where its edges cross it. Crossings on edges that lie on
      // another bbox line are fixed at the corner.
      for (auto [line_idx, line] : sp.data().bbox_lines) {
        CGAL_assertion(m_data.igraph().line_is_on_bbox(line_idx));
        const CGAL::Oriented_side inside = line.oriented_side(sp.data().ikcentroid);
        CGAL_assertion(inside != CGAL::ON_ORIENTED_BOUNDARY);

        std::list<std::size_t> pos_side;
        if (!m_data.split_polygon(polygons[0], line_idx, pos_side)) {
          CGAL_assertion_code(for (std::size_t v : polygons[0]))
            CGAL_assertion(Data_structure::side_0plus(line, vts[v].p0, vts[v].v) != CGAL::opposite(inside));
          continue;
        }

        if (inside == CGAL::ON_POSITIVE_SIDE)
          std::swap(polygons[0], pos_side);

        // The twins of the kept vertices are discarded.
        std::set<std::size_t> removed(pos_side.begin(), pos_side.end());
        for (std::size_t v : polygons[0])
          if (removed.count(vts[v].other))
            vts[v].other = std::size_t(-1);
      }

      for (auto [line_idx, line] : sp.data().lines) {
        // Skip lines on the bounding box, they won't split the polygon.
        if (m_data.igraph().line_is_on_bbox(line_idx))
          continue;

        std::size_t num_polys = polygons.size();
        for (std::size_t i = 0; i < num_polys; i++) {
          Polygon& poly = polygons[i];
          std::list<std::size_t> pos_side;
          if (m_data.split_polygon(poly, line_idx, pos_side))
            polygons.emplace_back(std::move(pos_side));
        }
      }

      // A fixed vertex may lie on more than two lines if several lines intersect in the same point.
      // Only the two lines carrying its polygon edges are kept as constraints.
      for (Polygon& poly : polygons)
        for (auto it = poly.begin(); it != poly.end(); ++it) {
          Vertex& v = vts[*it];
          if (v.moving || v.constraints.size() <= 2)
            continue;
          const Vertex& prev = vts[it == poly.begin() ? poly.back() : *std::prev(it)];
          const Vertex& next = vts[std::next(it) == poly.end() ? poly.front() : *std::next(it)];
          std::set<std::size_t> kept;
          for (std::size_t l : v.constraints)
            if (prev.constraints.count(l) || next.constraints.count(l))
              kept.insert(l);
          CGAL_assertion(kept.size() == 2);
          if (kept.size() == 2)
            v.constraints = kept;
        }

      // Identify which polygon belongs to which face and link the constrained vertices to the iedges.

      for (Polygon& poly : polygons) {
        bool has_moving = false;

        // The face is located with the centroid of the polygon at t = 0+, as the polygon may have zero area at t = 0
        // (e.g., if it touches a line in a vertex).
        IkVector_2 c0 = CGAL::NULL_VECTOR, cv = CGAL::NULL_VECTOR;
        for (std::size_t idx : poly) {
          c0 = c0 + (vts[idx].p0 - CGAL::ORIGIN);
          cv = cv + vts[idx].v;
        }
        c0 = c0 / IkFT(static_cast<unsigned int>(poly.size()));
        cv = cv / IkFT(static_cast<unsigned int>(poly.size()));
        IFace face = m_data.locate_iface_0plus(sp_idx, CGAL::ORIGIN + c0, cv);
        CGAL_assertion(face != IFace(-1));

        for (std::size_t idx : poly) {
          Vertex &v = vts[idx];
          if (v.moving)
            has_moving = true;

          if (v.constraints.size() == 1) {
            IVertex s(-1), t(-1);
            m_data.get_iedge(v, *v.constraints.begin(), s, t);
            v.constraint_edge = m_data.igraph().edge(s, t);
            v.itarget = t;

            std::pair<std::size_t, std::size_t> p1(idx, -1);

            auto pair = m_data.igraph().edge(v.constraint_edge).vertices.insert(std::make_pair(sp_idx, p1));

            // In case of paired vertices, only add one of them
            if (!pair.second && pair.first->second.first != v.other)
              pair.first->second.second = idx;
          }

          if (!v.moving) {
            std::set<std::size_t> planes;
            for (std::size_t c : v.constraints) {
              std::set<std::size_t> intersected_planes = m_data.igraph().planes_on_line(c);
              planes.insert(intersected_planes.begin(), intersected_planes.end());
            }

            v.ivertex = m_data.igraph().vertex(planes);
            if (v.ivertex == IVertex(-1))
              v.ivertex = m_data.igraph().vertex(sp.to_3d(v.p0));

            if (v.ivertex == IVertex(-1)) {
              std::cout << "vertex not found: " << sp_idx << " " << idx << std::endl;
              continue;
            }
          }
        } // end of vertex

        // Again treat fixed vertices to add them to the two iedges
        for (auto it = poly.begin();it != poly.end();++it) {
          Vertex& v = vts[*it];
          if (!v.moving) {
            Vertex& prev = vts[it == poly.begin() ? *std::prev(poly.end()) : *std::prev(it)];
            Vertex& next = vts[std::next(it) == poly.end() ? *poly.begin() : *std::next(it)];

            CGAL_assertion(v.constraints.size() == 2);
            CGAL_assertion(!prev.constraints.empty());
            CGAL_assertion(!next.constraints.empty());

            if (prev.moving)
              v.constraint_edge = prev.constraint_edge;
            else
              v.constraint_edge = m_data.igraph().edge(prev.ivertex, v.ivertex);

            if (next.moving)
              v.other_constraint_edge = next.constraint_edge;
            else
              v.other_constraint_edge = m_data.igraph().edge(v.ivertex, next.ivertex);

            CGAL_assertion(v.other_constraint_edge != m_data.igraph().null_iedge());
            CGAL_assertion(v.constraint_edge != m_data.igraph().null_iedge());

            std::pair<std::size_t, std::size_t> p1(*it, -1);

            auto pair = m_data.igraph().edge(v.constraint_edge).vertices.insert(std::make_pair(sp_idx, p1));

            // Already exists? Then fill in the second spot.
            if (!pair.second && pair.first->second.second == std::size_t(-1))
              pair.first->second.second = *it;

            auto pair2 = m_data.igraph().edge(v.other_constraint_edge).vertices.insert(std::make_pair(sp_idx, p1));

            // Already exists? Then fill in the second spot.
            if (!pair2.second && pair2.first->second.second == std::size_t(-1))
              pair2.first->second.second = *it;
          }
        }

        for (std::size_t i : poly)
          vts[i].face = face;

        // Move faces into sp data structure
        CGAL_assertion_code(auto res = ) sp.data().polygons.insert(std::make_pair(face, std::move(poly)));
        CGAL_assertion(res.second);
        sp.link_property_maps();
        m_data.init_border(face, true);
        sp.data().initial_ifaces.push_back(face);

        if(has_moving)
          sp.active_polygons++;
      } // end of poly

      for (std::size_t i = 0; i < m_data.igraph().number_of_lines(); i++) {
        const auto& line = m_data.igraph().line(i);
        const auto dir = line.to_vector();
        // Add lines lying in the support plane, e.g., when three planes intersect in the same line.
        if (CGAL::is_zero(sp.data().exact_plane.orthogonal_vector() * dir) && sp.data().exact_plane.has_on(line.point())
          && !m_data.igraph().vertices_on_line(i).empty()) {
          const std::vector<std::pair<IkFT, IVertex>>& vol = m_data.igraph().vertices_on_line(i);
          sp.data().lines.insert(std::make_pair(i, IkLine_2(sp.to_2d(m_data.igraph().point_3(vol[0].second)), sp.to_2d(m_data.igraph().point_3(vol.back().second)))));
        }
      }
    } // end of support plane

    // Check proper edge linking
    CGAL_assertion_code(
    for (std::size_t sp_idx = 6; sp_idx < m_data.number_of_support_planes(); sp_idx++) {
      Support_plane& sp = m_data.support_plane(sp_idx);
      for (auto [face, poly] : sp.data().polygons) {
        CGAL_assertion(face != IFace(-1));
        for (std::size_t i : poly) {
          CGAL_assertion(vts[i].face == face);
          if (!vts[i].moving) {
            CGAL_assertion(vts[i].ivertex != IVertex(-1));
            CGAL_assertion(vts[i].constraints.size() >= 2);
            //CGAL_assertion(vts[i].other != std::size_t(-1));
          }
          else if (!vts[i].constraints.empty()) {
            CGAL_assertion(vts[i].constraints.size() == 1);
            CGAL_assertion(vts[i].constraint_edge != m_data.igraph().null_iedge());
            //CGAL_assertion(vts[i].other != std::size_t(-1));
            auto it = m_data.igraph().edge(vts[i].constraint_edge).vertices.find(vts[i].sp_idx);
            CGAL_assertion(it != m_data.igraph().edge(vts[i].constraint_edge).vertices.end());
            Vertex& v1 = vts[it->second.first];
            Vertex& v2 = vts[it->second.second];
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
              CGAL_assertion(CGAL::is_negative(v1.v * v2.v)); // pointing in different directions
              CGAL_assertion(CGAL::is_zero(v1.v.x() * v2.v.y() - v1.v.y() * v2.v.x())); // collinear
            }
            else {
              // The moving vertex has to move away from the fixed one at t = 0+. A sliding vertex spawned at
              // a fixed vertex starts at the same position, i.e., (fixed - moving) * v is zero at t = 0 and
              // decreases with -|v|^2 afterwards.
              if (v1.moving) {
                CGAL_assertion(v1.v != CGAL::NULL_VECTOR);
                CGAL_assertion(!CGAL::is_positive((v2.p0 - v1.p0) * v1.v));
              }
              else {
                CGAL_assertion(v2.moving);
                CGAL_assertion(v2.v != CGAL::NULL_VECTOR);
                CGAL_assertion(!CGAL::is_positive((v1.p0 - v2.p0) * v2.v));
              }
            }
          }
        }
      }
    });
  }

  void bounding_box_to_polygons(const std::array<typename Intersection_kernel::Point_3, 8>& bbox, std::vector<std::vector<typename Intersection_kernel::Point_3> >& bbox_faces) const {
    bbox_faces.clear();
    bbox_faces.reserve(6);

    bbox_faces.push_back({ bbox[0], bbox[1], bbox[2], bbox[3] }); // zmin
    bbox_faces.push_back({ bbox[0], bbox[5], bbox[6], bbox[1] }); // ymin
    bbox_faces.push_back({ bbox[1], bbox[6], bbox[7], bbox[2] }); // xmax
    bbox_faces.push_back({ bbox[2], bbox[7], bbox[4], bbox[3] }); // ymax
    bbox_faces.push_back({ bbox[3], bbox[4], bbox[5], bbox[0] }); // xmin
    bbox_faces.push_back({ bbox[5], bbox[4], bbox[7], bbox[6] }); // zmax
    CGAL_assertion(bbox_faces.size() == 6);
  }

  void add_polygons(const std::vector<std::vector<typename Intersection_kernel::Point_3> >& bbox_faces, std::vector<std::size_t>& input_polygons) {
    add_bbox_faces(bbox_faces);

    // Filter input polygons
    std::vector<bool> remove(input_polygons.size(), false);
    for (std::size_t i = 0; i < 6; i++)
      for (std::size_t j = 0; j < m_input_planes.size(); j++)
        if (m_data.support_plane(i).exact_plane() == m_input_planes[j] || m_data.support_plane(i).exact_plane() == m_input_planes[j].opposite()) {
          m_data.support_plane(i).set_input_polygon(j);
          remove[j] = true;
        }

    std::size_t write = 0;
    for (std::size_t i = 0; i < input_polygons.size(); i++)
      if (!remove[i]) {
        m_input_polygons[write] = m_input_polygons[i];
        m_input_planes[write] = m_input_planes[i];
        input_polygons[write] = input_polygons[i];
        write++;
      }
    m_input_polygons.resize(write);
    m_input_planes.resize(write);
    input_polygons.resize(write);
    add_input_polygons();

    // Order vertices on lines

    std::vector<std::set<IVertex>> line2vertices(m_data.igraph().number_of_lines());
    for (IEdge e : m_data.igraph().edges()) {
      std::size_t line_idx = m_data.igraph().edge(e).line;
      CGAL_assertion(line_idx < m_data.igraph().number_of_lines());
      line2vertices[line_idx].insert(m_data.igraph().source(e));
      line2vertices[line_idx].insert(m_data.igraph().target(e));
    }

    for (std::size_t line_idx = 0; line_idx < m_data.igraph().number_of_lines(); line_idx++)
      m_data.set_vertices_on_line(line2vertices[line_idx], line_idx);
  }

  void add_bbox_faces(const std::vector< std::vector<typename Intersection_kernel::Point_3> >& bbox_faces) {
    for (const auto& bbox_face : bbox_faces)
      m_data.add_bbox_polygon(bbox_face);

    CGAL_assertion(m_data.number_of_support_planes() == 6);
    CGAL_assertion(m_data.ivertices().size() == 8);
    CGAL_assertion(m_data.iedges().size() == 12);

    if (m_parameters.verbose) {
      std::cout << "* inserted bbox faces: " << bbox_faces.size() << std::endl;
    }
  }

  void add_input_polygons() {
    using Polygon_2 = std::vector<typename Intersection_kernel::Point_2>;
    using Indices = std::vector<std::size_t>;

    std::map< std::size_t, std::pair<Polygon_2, Indices> > polygons;
    preprocess_polygons(polygons);

    for (const auto& item : polygons) {
      const std::size_t support_plane_idx = item.first;
      const auto& pair = item.second;
      const Polygon_2& polygon = pair.first;
      const Indices& input_indices = pair.second;
      m_data.add_input_polygon(support_plane_idx, input_indices, polygon);
      m_data.support_plane(support_plane_idx).set_input_polygon(static_cast<int>(item.first) - 6);
    }
    //dump_polygons(m_data, polygons, m_data.prefix() + "inserted-polygons");

    CGAL_assertion(m_data.number_of_support_planes() >= 6);
    if (m_parameters.verbose) {
      std::cout << "* provided input polygons: " << m_data.input_polygons().size() << std::endl;
      std::cout << "* inserted input polygons: " << polygons.size() << std::endl;
    }
  }

  template<typename PointRange>
  void convert_polygon(const std::size_t support_plane_idx, const PointRange& polygon_3, std::vector<typename Intersection_kernel::Point_2>& polygon_2) {
    polygon_2.clear();
    polygon_2.reserve(polygon_3.size());

    for (const auto& point : polygon_3) {
      polygon_2.push_back(m_data.support_plane(support_plane_idx).to_2d(point));
    }
    CGAL_assertion(polygon_2.size() == polygon_3.size());
  }

  void restrict_to_bbox(const std::size_t support_plane_idx, std::vector<typename Intersection_kernel::Point_2>& polygon_2) {
    // The input points are projected exactly onto their support planes before the octree clips
    // them (see Kinetic_space_partition_3::split_octree()). Hence all vertices are exactly inside
    // the bbox or on its boundary and no snapping is required.
#ifndef CGAL_NO_ASSERTIONS
    Support_plane &sp = m_data.support_plane(support_plane_idx);
    IkPoint_2 c = CGAL::centroid(polygon_2.begin(), polygon_2.end(), CGAL::Dimension_tag<0>());
    for (const auto &e : sp.data().iedges) {
      const IkLine_2 line = sp.to_2d(m_data.igraph().line(m_data.igraph().edge(e).line));
      const CGAL::Oriented_side inside = line.oriented_side(c);
      CGAL_assertion(inside != CGAL::ON_ORIENTED_BOUNDARY);
      for (const auto& p : polygon_2) {
        CGAL_assertion_code(const CGAL::Oriented_side side = line.oriented_side(p);)
        CGAL_assertion(side == inside || side == CGAL::ON_ORIENTED_BOUNDARY);
      }
    }
#else
    CGAL_USE(support_plane_idx);
#endif

    Polygon_2<Intersection_kernel> poly(polygon_2.begin(), polygon_2.end());
    polygon_2.clear();
    convex_hull_2(poly.vertices_begin(), poly.vertices_end(), std::back_inserter(polygon_2));
  }

  void preprocess_polygons(std::map< std::size_t, std::pair<std::vector<typename Intersection_kernel::Point_2>, std::vector<std::size_t> > >& polygons) {
    std::size_t input_index = 0;
    std::vector<typename Intersection_kernel::Point_2> polygon_2;
    std::vector<std::size_t> input_indices;
    for (std::size_t i = 0; i < m_input_polygons.size(); i++) {
      bool is_added = true;
      std::size_t support_plane_idx = std::size_t(-1);
      std::tie(support_plane_idx, is_added) = m_data.add_support_plane(m_input_polygons[i], false, m_input_planes[i]);
      CGAL_assertion(support_plane_idx != std::size_t(-1));
      convert_polygon(support_plane_idx, m_input_polygons[i], polygon_2);
      restrict_to_bbox(support_plane_idx, polygon_2);

      if (is_added) {
        input_indices.clear();
        input_indices.push_back(input_index);
        polygons[support_plane_idx] = std::make_pair(polygon_2, input_indices);

      }
      else {
        CGAL_assertion(polygons.find(support_plane_idx) != polygons.end());
        auto& pair = polygons.at(support_plane_idx);
        auto& other_polygon = pair.first;
        auto& other_indices = pair.second;
        other_indices.push_back(input_index);
        merge_polygons(support_plane_idx, polygon_2, other_polygon);
      }
      ++input_index;
    }
  }

  void merge_polygons(const std::size_t support_plane_idx, const std::vector<typename Intersection_kernel::Point_2>& polygon_a, std::vector<typename Intersection_kernel::Point_2>& polygon_b) {
    const bool is_debug = true;
    CGAL_assertion(support_plane_idx >= 6);
    if (is_debug) {
      std::cout << std::endl << "support plane idx: " << support_plane_idx << std::endl;
    }

    // Add points from a to b.
    auto& points = polygon_b;
    for (const auto& point : polygon_a) {
      points.push_back(point);
    }

    // Create the merged polygon.
    std::vector<typename Intersection_kernel::Point_2> merged;
    create_merged_polygon(points, merged);

    if (is_debug) {
      std::cout << "merged polygon: " << std::endl;
      From_exact from_exact;
      for (std::size_t i = 0; i < merged.size(); ++i) {
        const std::size_t ip = (i + 1) % merged.size();
        const auto& p = merged[i];
        const auto& q = merged[ip];
        std::cout << "2 " <<
          m_data.to_3d(support_plane_idx, from_exact(p)) << " " <<
          m_data.to_3d(support_plane_idx, from_exact(q)) << std::endl;
      }
    }

    // Update b with the new merged polygon.
    polygon_b = merged;
  }

  void create_merged_polygon(const std::vector<typename Intersection_kernel::Point_2>& points, std::vector<typename Intersection_kernel::Point_2>& merged) const {
    merged.clear();

    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(merged));

    CGAL_assertion(merged.size() >= 3);
  }

  void create_bbox_meshes() {
    for (std::size_t i = 0; i < 6; i++) {
      m_data.clear_pfaces(i);
      std::set<IFace> ifaces = m_data.support_plane(i).ifaces();

      for (auto iface : ifaces) {
        m_data.add_iface_to_mesh(i, iface);
      }
    }
  }

  void create_internal_vertices() {
    // First, create all transverse intersection lines.
    using Map_p2vv = std::map<std::set<std::size_t>, std::pair<IVertex, IVertex> >;
    Map_p2vv map_p2vv;

    for (const auto ivertex : m_data.ivertices()) {
      const auto key = m_data.intersected_planes(ivertex, false);
      if (key.size() < 2) {
        continue;
      }

      const auto pair = map_p2vv.insert(
        std::make_pair(key, std::make_pair(ivertex, IVertex())));
      const bool is_inserted = pair.second;
      if (!is_inserted) {
        pair.first->second.second = ivertex;
      }
    }

    std::unordered_map<std::set<std::size_t>, IVertex, boost::hash<std::set<std::size_t>>> internal_v;

    // Then, intersect these lines to find internal intersection vertices.
    using Pair_pv = std::pair< std::set<std::size_t>, std::vector<IVertex> >;

    std::vector<Pair_pv> todo;
    for (typename Map_p2vv::iterator it_a = map_p2vv.begin(); it_a != map_p2vv.end(); ++it_a) {
      const auto& set_a = it_a->first;

      todo.push_back(std::make_pair(set_a, std::vector<IVertex>()));
      auto& crossed_vertices = todo.back().second;
      crossed_vertices.push_back(it_a->second.first);

      std::unordered_set<std::set<std::size_t>, boost::hash<std::set<std::size_t>>> done;

      for (typename Map_p2vv::iterator it_b = map_p2vv.begin(); it_b != map_p2vv.end(); ++it_b) {
        if (it_a == it_b)
          continue;

        const auto& set_b = it_b->first;

        std::size_t common_plane_idx = std::size_t(-1);
        const std::function<void(const std::size_t idx)> lambda =
          [&](const std::size_t idx) {
          common_plane_idx = idx;
          };

        std::set_intersection(
          set_a.begin(), set_a.end(),
          set_b.begin(), set_b.end(),
          boost::make_function_output_iterator(lambda)
        );

        if (common_plane_idx != std::size_t(-1)) {
          auto union_set = set_a;
          union_set.insert(set_b.begin(), set_b.end());
          if (!done.insert(union_set).second)
            continue;

          IVertex v;
          auto it = internal_v.find(union_set);
          if (it != internal_v.end())
            v = it->second;
          else {
            std::size_t other_b = ((*set_b.begin()) == common_plane_idx) ? *(++set_b.begin()) : *set_b.begin();

            typename Intersection_kernel::Segment_3 seg_a(m_data.point_3(it_a->second.first), m_data.point_3(it_a->second.second));
            typename Intersection_kernel::Point_3 point;
            if (!intersection(seg_a, m_data.support_plane(other_b).exact_plane(), point))
              continue;

            v = m_data.add_ivertex(point, union_set);
            internal_v[union_set] = v;
          }

          // Several lines may intersect line a in the same vertex.
          if (std::find(crossed_vertices.begin(), crossed_vertices.end(), v) == crossed_vertices.end())
            crossed_vertices.push_back(v);
        }
      }
      crossed_vertices.push_back(it_a->second.second);
    }

    for (auto& t : todo) {
      m_data.add_iedge(t.first, t.second);
    }

    return;
  }

  template<typename Type1, typename Type2, typename ResultType>
  inline bool intersection(const Type1& t1, const Type2& t2, ResultType& result) const {

    const auto inter = CGAL::intersection(t1, t2);
    if (!inter) return false;
    if (CGAL::assign(result, inter))
      return true;

    return false;
  }
};

#endif //DOXYGEN_RUNNING

} // namespace internal
} // namespace KSP_3
} // namespace CGAL

#endif // CGAL_KSP_3_INITIALIZER_H
