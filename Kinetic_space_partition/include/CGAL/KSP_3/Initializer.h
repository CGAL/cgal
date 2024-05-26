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
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <boost/optional/optional_io.hpp>

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

  using Data_structure = KSP_3::internal::Data_structure<Kernel, Intersection_kernel>;
  using Support_plane = typename Data_structure::Support_plane;
  using IEdge = typename Data_structure::IEdge;
  using IFace = typename Data_structure::IFace;
  using Face_property = typename Data_structure::Intersection_graph::Face_property;
  using Intersection_graph = typename Data_structure::Intersection_graph;
  using IEdge_set = typename Data_structure::IEdge_set;

  using IVertex = typename Data_structure::IVertex;

  using To_exact = CGAL::Cartesian_converter<Kernel, Intersection_kernel>;
  using From_exact = CGAL::Cartesian_converter<Intersection_kernel, Kernel>;

  using Bbox_3 = CGAL::Bbox_3;
  using OBB_traits = CGAL::Oriented_bounding_box_traits_3<Kernel>;

  using Parameters = KSP::internal::Parameters_3<FT>;

  using Timer = CGAL::Real_timer;

public:
  Initializer(std::vector<std::vector<Point_3> >& input_polygons, Data_structure& data, const Parameters& parameters) :
    m_input_polygons(input_polygons), m_data(data), m_parameters(parameters)
  { }

  Initializer(std::vector<std::vector<Point_3> >& input_polygons, std::vector<typename Intersection_kernel::Plane_3>& input_planes, Data_structure& data, const Parameters& parameters) :
    m_input_polygons(input_polygons), m_input_planes(input_planes), m_data(data), m_parameters(parameters)
  { }

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

    // Fills in the ivertices on support plane intersections inside the bbox.
    make_polygons_intersection_free();

    // Generation of ifaces
    create_ifaces();

    // Splitting the input polygons along intersection lines.
    initial_polygon_iedge_intersections();

    create_bbox_meshes();

    // Starting from here the intersection graph is const, it won't change anymore.
    if (m_parameters.verbose)
      std::cout << "done" << std::endl;

    if (m_parameters.debug)
      KSP_3::internal::dump(m_data, m_data.prefix() + "intersected");

    CGAL_assertion(m_data.check_bbox());
    //m_data.set_limit_lines();
    m_data.precompute_iedge_data();

    m_data.initialization_done();

    if (m_parameters.debug) {
      for (std::size_t sp = 0; sp < m_data.number_of_support_planes(); sp++)
        dump_2d_surface_mesh(m_data, sp, m_data.prefix() + "before-partition-sp" + std::to_string(sp));
    }

    if (m_parameters.verbose) {
      std::cout << "v: " << m_data.igraph().number_of_vertices() << " f: " << m_data.igraph().number_of_faces() << std::endl;
    }
  }

  void clear() {
    // to be added
  }

private:
  std::vector<std::vector<Point_3> >& m_input_polygons;
  std::vector<typename Intersection_kernel::Plane_3>& m_input_planes;
  Data_structure& m_data;
  const Parameters& m_parameters;

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

      std::vector<std::pair<IEdge, Direction_2> > connected;
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

    std::vector<std::pair<IEdge, Direction_2> > connected;
    m_data.get_and_sort_all_connected_iedges(sp_idx, m_data.target(edge), connected);
    //if (connected.size() <= 2) ivertex is on bbox edge
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
      const IEdge_set& uiedges = m_data.support_plane(sp_idx).unique_iedges();

      // Special case bbox without splits
      if (sp_idx < 6 && uiedges.size() == 4) {
        // Get first edge
        IEdge first = *uiedges.begin();
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
          face.vertices.push_back(t);
          face.pts.push_back(m_data.support_plane(sp_idx).to_2d(m_data.igraph().point_3(t)));
          face.edges.push_back(edge);
          m_data.igraph().add_face(sp_idx, edge, face_idx);
        }

        // create polygon in proper order
      }

      bool all_on_bbox = true;
      for (auto edge : uiedges) {
        bool on_edge = m_data.igraph().iedge_is_on_bbox(edge);
        //if (m_data.igraph().iedge_is_on_bbox(edge))
        //  continue;
        //
        //Note the number of bbox lines during creation and skip all those.

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
        get_prev_next(sp_idx, *uiedges.begin(), prev, next);
        add_iface_from_iedge(sp_idx, *uiedges.begin(), prev, true);
      }
    }
  }

  void initial_polygon_iedge_intersections() {
    To_exact to_exact;
    From_exact from_exact;

    for (std::size_t sp_idx = 0; sp_idx < m_data.number_of_support_planes(); sp_idx++) {
      bool polygons_assigned = false;
      Support_plane& sp = m_data.support_plane(sp_idx);
      if (sp.is_bbox())
        continue;

      sp.mesh().clear_without_removing_property_maps();

      std::map<std::size_t, std::vector<IEdge> > line2edges;
      // Get all iedges, sort into lines and test intersection per line?
      for (const IEdge& edge : sp.unique_iedges()) {

        if (m_data.is_bbox_iedge(edge))
          continue;

        std::size_t line = m_data.igraph().line(edge);

        line2edges[line].push_back(edge);
      }

      for (auto pair : line2edges) {
        // Get line
        //Line_2 l(sp.to_2d(m_data.point_3(m_data.source(pair.second[0]))),sp.to_2d(m_data.point_3(m_data.target(pair.second[0]))));

        typename Intersection_kernel::Point_2 a(sp.to_2d(m_data.point_3(m_data.source(pair.second[0]))));
        typename Intersection_kernel::Point_2 b(sp.to_2d(m_data.point_3(m_data.target(pair.second[0]))));
        typename Intersection_kernel::Line_2 exact_line(a, b);
        Line_2 l = from_exact(exact_line);

        typename Intersection_kernel::Vector_2 ldir = exact_line.to_vector();
        ldir = (typename Intersection_kernel::FT(1.0) / CGAL::approximate_sqrt(ldir * ldir)) * ldir;
        Vector_2 dir = from_exact(ldir);

        std::vector<typename Intersection_kernel::Segment_2> crossing_polygon_segments;
        std::vector<IEdge> crossing_iedges;
        typename Intersection_kernel::FT emin = (std::numeric_limits<double>::max)();
        typename Intersection_kernel::FT emax = -(std::numeric_limits<double>::max)();
        FT min_speed = (std::numeric_limits<double>::max)(), max_speed = -(std::numeric_limits<double>::max)();

        CGAL::Oriented_side last_side = l.oriented_side(sp.data().original_vertices.back());
        Point_2 minp, maxp;
        typename Intersection_kernel::Point_2 eminp, emaxp;

        // Map polygon to line and get min&max projection
        for (std::size_t v = 0; v < sp.data().original_vertices.size(); v++) {
          const Point_2& p = sp.data().original_vertices[v];

          CGAL::Oriented_side s = l.oriented_side(p);
          if (last_side != s) {
            // Fetch former point to add segment.
            const Point_2& prev = sp.data().original_vertices[(v + sp.data().original_vertices.size() - 1) % sp.data().original_vertices.size()];
            const Vector_2 edge_dir = sp.original_edge_direction((v + sp.data().original_vertices.size() - 1) % sp.data().original_vertices.size(), v);
            typename Intersection_kernel::Segment_2 seg(to_exact(prev), to_exact(p));
            const auto result = CGAL::intersection(seg, exact_line);
            typename Intersection_kernel::Point_2 intersection;

            if (result && CGAL::assign(intersection, result)) {
              typename Intersection_kernel::FT eproj = (intersection - exact_line.point()) * ldir;
              //FT proj = to_inexact(eproj);
              if (eproj < emin) {
                eminp = intersection;
                emin = eproj;
                minp = from_exact(intersection);
                //min = proj;
                typename Intersection_kernel::FT p = dir * edge_dir;
                assert(p != 0);
                min_speed = CGAL::approximate_sqrt(edge_dir * edge_dir) / from_exact(p);
              }
              if (emax < eproj) {
                emaxp = intersection;
                emax = eproj;
                maxp = from_exact(intersection);
                //max = proj;
                typename Intersection_kernel::FT p = dir * edge_dir;
                assert(p != 0);
                max_speed = CGAL::approximate_sqrt(edge_dir * edge_dir) / from_exact(p);
              }
            }
            else std::cout << "crossing segment does not intersect line" << std::endl;
            crossing_polygon_segments.push_back(seg);
          }

          last_side = s;
        }

        // Is there any intersection?
        // As the polygon is convex there can only be one line segment on the inside of the polygon
        if (emin < emax) {
          m_data.support_plane(sp_idx).set_crossed_line(pair.first);
          // Collect crossing edges by overlapping min/max barycentric coordinates on line
          for (IEdge e : pair.second) {
            std::pair<IFace, IFace> faces;
            m_data.igraph().get_faces(sp_idx, e, faces);
            IVertex lower = m_data.source(e);
            IVertex upper = m_data.target(e);
            if (lower > upper) {
              IVertex tmp = upper;
              upper = lower;
              lower = tmp;
            }
            typename Intersection_kernel::FT s = (sp.to_2d(m_data.point_3(lower)) - exact_line.point()) * ldir;
            typename Intersection_kernel::FT t = (sp.to_2d(m_data.point_3(upper)) - exact_line.point()) * ldir;

            if (s < t) {
              if (s < emax && emin < t) {
                std::pair<IFace, IFace> faces;
                m_data.igraph().get_faces(sp_idx, e, faces);

                polygons_assigned = true;

                if (!m_data.igraph().face(faces.first).part_of_partition) {
                  auto pface = m_data.add_iface_to_mesh(sp_idx, faces.first);
                  sp.data().initial_ifaces.push_back(faces.first);
                  sp.set_initial(pface.second);
                }

                if (!m_data.igraph().face(faces.second).part_of_partition) {
                  auto pface = m_data.add_iface_to_mesh(sp_idx, faces.second);
                  sp.data().initial_ifaces.push_back(faces.second);
                  sp.set_initial(pface.second);
                }

                typename Intersection_graph::Kinetic_interval& kinetic_interval = m_data.igraph().kinetic_interval(e, sp_idx);
                crossing_iedges.push_back(e);
                if (emin > s) {
                  typename Intersection_kernel::FT bary_edge_exact = (emin - s) / (t - s);
                  FT bary_edge = from_exact((emin - s) / (t - s));
                  CGAL_assertion(bary_edge_exact >= 0);
                  FT time = CGAL::abs(from_exact(s - emin) / min_speed);
                  kinetic_interval.push_back(std::pair<FT, FT>(0, time)); // border barycentric coordinate
                  kinetic_interval.push_back(std::pair<FT, FT>(bary_edge, 0));
                }
                else {
                  kinetic_interval.push_back(std::pair<FT, FT>(0, 0));
                }

                if (t > emax) {
                  typename Intersection_kernel::FT bary_edge_exact = (emax - s) / (t - s);
                  FT bary_edge = from_exact((emax - s) / (t - s));
                  CGAL_assertion(0 <= bary_edge_exact && bary_edge_exact <= 1);
                  FT time = CGAL::abs(from_exact(emax - t) / max_speed);
                  kinetic_interval.push_back(std::pair<FT, FT>(bary_edge, 0));
                  kinetic_interval.push_back(std::pair<FT, FT>(1, time)); // border barycentric coordinate
                }
                else
                  kinetic_interval.push_back(std::pair<FT, FT>(1, 0));
              }
            }
            else if (t < emax && emin < s) {
              std::pair<IFace, IFace> faces;
              m_data.igraph().get_faces(sp_idx, e, faces);

              polygons_assigned = true;

              if (!m_data.igraph().face(faces.first).part_of_partition) {
                auto pface = m_data.add_iface_to_mesh(sp_idx, faces.first);
                sp.data().initial_ifaces.push_back(faces.first);
                sp.set_initial(pface.second);
              }

              if (!m_data.igraph().face(faces.second).part_of_partition) {
                auto pface = m_data.add_iface_to_mesh(sp_idx, faces.second);
                sp.data().initial_ifaces.push_back(faces.second);
                sp.set_initial(pface.second);
              }

              typename Intersection_graph::Kinetic_interval& kinetic_interval = m_data.igraph().kinetic_interval(e, sp_idx);
              crossing_iedges.push_back(e);
              if (s > emax) {
                typename Intersection_kernel::FT bary_edge_exact = (s - emax) / (s - t);
                FT bary_edge = from_exact((s - emax) / (s - t));
                CGAL_assertion(0 <= bary_edge_exact && bary_edge_exact <= 1);
                FT time = CGAL::abs(from_exact(emax - s) / max_speed);
                kinetic_interval.push_back(std::pair<FT, FT>(0, time)); // border barycentric coordinate
                kinetic_interval.push_back(std::pair<FT, FT>(bary_edge, 0));
              }
              else
                kinetic_interval.push_back(std::pair<FT, FT>(0, 0));

              if (emin > t) {
                typename Intersection_kernel::FT bary_edge_exact = (s - emin) / (s - t);
                FT bary_edge = from_exact(bary_edge_exact);
                CGAL_assertion(0 <= bary_edge_exact && bary_edge_exact <= 1);
                FT time = CGAL::abs(from_exact(t - emin) / min_speed);
                kinetic_interval.push_back(std::pair<FT, FT>(bary_edge, 0));
                kinetic_interval.push_back(std::pair<FT, FT>(1, time)); // border barycentric coordinate
              }
              else
                kinetic_interval.push_back(std::pair<FT, FT>(1, 0));
            }
          }
        }
      }

      // If no faces have been assigned, the input polygon lies completely inside a face.
      // Every IFace is checked whether the polygon, or just a single vertex, lies inside.
      // The implementation takes advantage of the IFace being convex.
      if (!polygons_assigned) {
        IFace face = IFace(-1);
        for (auto& f : sp.ifaces()) {
          Face_property& fp = m_data.igraph().face(f);

          typename Intersection_kernel::Point_2 p = to_exact(sp.data().centroid);
          bool outside = false;

          // poly, vertices and edges in IFace are oriented ccw
          for (std::size_t i = 0; i < fp.pts.size(); i++) {
            typename Intersection_kernel::Vector_2 ts = fp.pts[(i + fp.pts.size() - 1) % fp.pts.size()] - p;
            typename Intersection_kernel::Vector_2 tt = fp.pts[i] - p;

            bool ccw = (tt.x() * ts.y() - tt.y() * ts.x()) <= 0;
            if (!ccw) {
              outside = true;
              break;
            }
          }
          if (!outside) {
            if (face == IFace(-1))
              face = f;
            else {
              std::cout << "Two faces found for " << sp_idx << " sp, f1 " << face << " f2 " << f << std::endl;
            }
          }
        }
        if (face != IFace(-1)) {
          if (!m_data.igraph().face(face).part_of_partition) {
            auto pface = m_data.add_iface_to_mesh(sp_idx, face);
            sp.data().initial_ifaces.push_back(face);
            sp.set_initial(pface.second);
          }
        }
        else
          std::cout << "No IFace found for sp " << sp_idx << std::endl;
      }
    }
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
    using Polygon_2 = std::vector<Point_2>;
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
  void convert_polygon(const std::size_t support_plane_idx, const PointRange& polygon_3, std::vector<Point_2>& polygon_2) {
    polygon_2.clear();
    polygon_2.reserve(polygon_3.size());
    for (const auto& point : polygon_3) {
      const Point_3 converted(static_cast<FT>(point.x()), static_cast<FT>(point.y()), static_cast<FT>(point.z()));

      polygon_2.push_back(m_data.support_plane(support_plane_idx).to_2d(converted));
    }
    CGAL_assertion(polygon_2.size() == polygon_3.size());
  }

  void preprocess_polygons(std::map< std::size_t, std::pair<std::vector<Point_2>, std::vector<std::size_t> > >& polygons) {
    std::size_t input_index = 0;
    std::vector<Point_2> polygon_2;
    std::vector<std::size_t> input_indices;
    for (std::size_t i = 0; i < m_input_polygons.size(); i++) {
      bool is_added = true;
      std::size_t support_plane_idx = std::size_t(-1);
      std::tie(support_plane_idx, is_added) = m_data.add_support_plane(m_input_polygons[i], false, m_input_planes[i]);
      CGAL_assertion(support_plane_idx != std::size_t(-1));
      convert_polygon(support_plane_idx, m_input_polygons[i], polygon_2);

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

  void merge_polygons(const std::size_t support_plane_idx, const std::vector<Point_2>& polygon_a, std::vector<Point_2>& polygon_b) {
    const bool is_debug = false;
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
    std::vector<Point_2> merged;
    create_merged_polygon(points, merged);

    if (is_debug) {
      std::cout << "merged polygon: " << std::endl;
      for (std::size_t i = 0; i < merged.size(); ++i) {
        const std::size_t ip = (i + 1) % merged.size();
        const auto& p = merged[i];
        const auto& q = merged[ip];
        std::cout << "2 " <<
          m_data.to_3d(support_plane_idx, p) << " " <<
          m_data.to_3d(support_plane_idx, q) << std::endl;
      }
    }

    // Update b with the new merged polygon.
    polygon_b = merged;
  }

  void create_merged_polygon(const std::vector<Point_2>& points, std::vector<Point_2>& merged) const {
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

  void make_polygons_intersection_free() {
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

    // Then, intersect these lines to find internal intersection vertices.
    using Pair_pv = std::pair< std::set<std::size_t>, std::vector<IVertex> >;
    std::vector<Pair_pv> todo;
    for (auto it_a = map_p2vv.begin(); it_a != map_p2vv.end(); ++it_a) {
      const auto& set_a = it_a->first;

      todo.push_back(std::make_pair(set_a, std::vector<IVertex>()));
      auto& crossed_vertices = todo.back().second;
      crossed_vertices.push_back(it_a->second.first);

      std::set<std::set<std::size_t>> done;
      for (auto it_b = map_p2vv.begin(); it_b != map_p2vv.end(); ++it_b) {
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
          if (!done.insert(union_set).second) {
            continue;
          }

          typename Intersection_kernel::Point_2 point;
          typename Intersection_kernel::Segment_3 seg_a(m_data.point_3(it_a->second.first), m_data.point_3(it_a->second.second));
          typename Intersection_kernel::Segment_3 seg_b(m_data.point_3(it_b->second.first), m_data.point_3(it_b->second.second));
          if (!intersection(m_data.support_plane(common_plane_idx).to_2d(seg_a), m_data.support_plane(common_plane_idx).to_2d(seg_b), point))
            continue;

          crossed_vertices.push_back(
            m_data.add_ivertex(m_data.support_plane(common_plane_idx).to_3d(point), union_set));
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
