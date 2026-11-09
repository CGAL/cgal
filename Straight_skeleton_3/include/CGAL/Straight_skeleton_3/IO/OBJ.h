// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   db/3d/OBJFile.h
 * author Gernot Walzl
 * date   2012-05-15
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_IO_OBJ_H
#define CGAL_STRAIGHT_SKELETON_3_IO_OBJ_H

#include <CGAL/license/Straight_skeleton_3.h>

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/mark_domain_in_triangulation.h>

#include <cmath>
#include <exception>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace CGAL {

template <typename GeomTraits>
class Straight_skeleton_3;

namespace Straight_skeletons_3 {
namespace internal {
namespace HDS {

template <typename GeomTraits>
class Polyhedron;

} // namespace HDS
namespace algorithm {

template <typename GeomTraits>
class Polyhedron_transformation;

} // namespace algorithm
} // namespace internal

namespace IO {

/*
 * It is impossible for an `.obj` file to store holes inside a facet,
 * so we triangulate facets in this case.
 */
template <typename GeomTraits,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(const std::string& filename,
               std::shared_ptr<internal::HDS::Polyhedron<GeomTraits> > polyhedron,
               const CGAL_NP_CLASS& np = parameters::default_values())
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;

  using Polyhedron = internal::HDS::Polyhedron<GeomTraits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  bool do_triangulate = !choose_parameter(get_parameter(np, CGAL::internal_np::do_not_triangulate_faces), false);

  CGAL_SS3_IO_TRACE_V(8, "-- OBJ::Save(Polyhedron) to " << filename << " --");
  CGAL_SS3_IO_TRACE_V(8, "   do_triangulate: " << std::boolalpha << do_triangulate);
  // CGAL_SS3_IO_TRACE_V(8, "Saving to OBJ:\n" << polyhedron->to_string());
  CGAL_SS3_IO_TRACE_V(8, polyhedron->vertices().size() << " NV, " << polyhedron->facets().size() << " NF");

  // If we want to tolerate intersections for debugging
  // using Itag = CGAL::Exact_intersections_tag;
  // using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
  using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;

  using PK = CGAL::Projection_traits_3<GeomTraits>;
  using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
  using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
  using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
  using PTDS = CGAL::Triangulation_data_structure_2<PVb, PFb>;
  using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
  using PCDT_VH = typename PCDT::Vertex_handle;
  using PCDT_FH = typename PCDT::Face_handle;

  std::stringstream oss;
  set_stream_precision_from_NP(oss, np);

  // Map for unique vertices and their indices
  std::map<VertexSPtr, int> vertex_to_index;
  int next_index = 1; // OBJ indices start at 1

  for (const VertexSPtr& vertex : polyhedron->vertices()) {
    vertex_to_index[vertex] = next_index++;
  }

  // Write vertices
  for (const VertexSPtr& vertex : polyhedron->vertices()) {
    const Point_3& pt = vertex->point();
    oss << "v " << CGAL::to_double(pt.x()) << " "
                << CGAL::to_double(pt.y()) << " "
                << CGAL::to_double(pt.z()) << "\n";
  }

  // Write facets
  for (const FacetSPtr& facet : polyhedron->facets()) {
    bool do_triangulate_facet = do_triangulate;
    if (facet->edges().size() <= 3) {
      do_triangulate_facet = false;
    }

    if (do_triangulate_facet) {
      Vector_3 n = facet->supporting_plane().orthogonal_vector();
      CGAL_assertion(n != CGAL::NULL_VECTOR);

      PK traits(n);
      PCDT pcdt(traits);

      std::map<VertexSPtr, PCDT_VH> face_vhs;

      for (const VertexSPtr& vertex : facet->vertices()) {
        auto res = face_vhs.emplace(vertex, PCDT_VH());
        if (res.second) { // first time seeing this point
          PCDT_VH vh = pcdt.insert(vertex->point());
          vh->info() = vertex;
          res.first->second = vh;
        }
      }

      auto ne = 0;
      for (const EdgeSPtr& edge : facet->edges()) {
        VertexSPtr v0 = edge->src(facet);
        VertexSPtr v1 = edge->tgt(facet);

        if (v0->point() == v1->point()) {
          CGAL_SS3_IO_TRACE_V(64, "Warning: degenerate edge @ " << v0->point());

          CGAL_assertion(v0->degree() != 1); // @todo handle that...
          VertexSPtr vm1 = edge->prev(facet)->src(facet);

          // manually create a degenerate facet so that the resulting mesh is conforming
          oss << "f " << vertex_to_index.at(vm1) << " "
                      << vertex_to_index.at(v0) << " "
                      << vertex_to_index.at(v1) << "\n";
        }
        else
        {
          PCDT_VH vh0 = face_vhs.at(v0);
          PCDT_VH vh1 = face_vhs.at(v1);

          try
          {
            pcdt.insert_constraint(vh0, vh1);
          }
          catch(const typename PCDT::Intersection_of_constraints_exception&)
          {
            CGAL_SS3_IO_TRACE_V(1, "Warning: Intersection of constraints");
            CGAL_SS3_IO_TRACE_V(1, "While inserting " << v0->point() << " || " << v1->point());
            CGAL_SS3_IO_TRACE_V(1, facet->to_string());
            CGAL_warning_msg(false, "Intersections in CDT2 are not allowed");
            return false;
          }
          ++ne;
        }
      }

      if (pcdt.finite_vertex_handles().size() != facet->vertices().size()) {
        // CGAL_SS3_IO_TRACE("Warning: CDT #nv != facet #nv. Constraint intersection?");
        // CGAL_SS3_IO_TRACE("Facet: " << facet->to_string());
        // CGAL_SS3_IO_TRACE("CDT: " << pcdt.number_of_vertices() << " vertices, "
        //                           << pcdt.number_of_faces() << " faces");
        // CGAL_SS3_IO_TRACE("Vertices: ");
        // CGAL_SS3_IO_TRACE_CODE(for (PCDT_VH vh : pcdt.finite_vertex_handles()) {)
        // CGAL_SS3_IO_TRACE("  " << vh->point() << " -> " << vh->info()->id());
        // CGAL_SS3_IO_TRACE_CODE(})
      }

      if (ne < 3) { // degenerate facet
        CGAL_SS3_IO_TRACE_V(1, "Warning: skipping degenerate facet");
        continue;
      }

      if (do_triangulate_facet) {
        CGAL::unordered_flat_map<PCDT_FH, bool> in_domain_map;
        boost::associative_property_map<CGAL::unordered_flat_map<PCDT_FH, bool>> in_domain(in_domain_map);

        CGAL::mark_domain_in_triangulation(pcdt, in_domain);

        for (auto fh : pcdt.finite_face_handles())
        {
          if (!get(in_domain, fh)) {
            continue;
          }
          oss << "f " << vertex_to_index.at(fh->vertex(0)->info()) << " "
                      << vertex_to_index.at(fh->vertex(1)->info()) << " "
                      << vertex_to_index.at(fh->vertex(2)->info()) << "\n";
        }
      }
    }

    if (!do_triangulate_facet) {
      std::set<EdgeSPtr> visited_edges;

      for (EdgeSPtr edge : facet->edges()) {
        CGAL_SS3_DEBUG_SPTR(edge);
        if (visited_edges.find(edge) != visited_edges.end()) {
          continue; // already visited
        }

        std::vector<VertexSPtr> boundary_vertices;
        EdgeSPtr start_edge = edge;
        bool is_open = false;

        // Walk forward to collect boundary vertices
        do {
          CGAL_SS3_DEBUG_SPTR(edge);
          visited_edges.insert(edge);
          boundary_vertices.push_back(edge->src(facet));
          if (edge->tgt(facet)->degree() == 1) {
            is_open = true;
            break;
          }
          edge = edge->next(facet);
        } while (edge != start_edge);

        // If open, also walk backward to collect remaining boundary vertices
        if (is_open) {
          boundary_vertices.push_back(edge->tgt(facet));
          edge = start_edge;
          while (edge->src(facet)->degree() != 1) {
            edge = edge->prev(facet);
            visited_edges.insert(edge);
            boundary_vertices.insert(boundary_vertices.begin(), edge->src(facet));
          }
        }

        // Write the boundary as a face
        oss << "f ";
        for (const auto& vertex : boundary_vertices) {
          oss << vertex_to_index.at(vertex) << " ";
        }

        oss << "\n";
      }
    }
  }

  // result = (polyhedron->vertices().size() == vertex_id &&
  //           polyhedron->facets().size() == facet_id);

  std::ofstream ofs(filename.c_str());
  if (ofs.is_open()) {
    set_stream_precision_from_NP(oss, np);
    ofs << oss.str();
  } else {
    CGAL_SS3_IO_TRACE_V(1, "Error: failed to open file: " << filename.c_str());
    CGAL_assertion(false);
    return false;
  }

  CGAL_SS3_IO_TRACE_V(8, "-- Write OBJ end --");
  return true;
}

/*!
 * \ingroup PkgStraightSkeleton3IO
 *
 * writes a `CGAL::Straight_skeleton_3` to a Wavefront OBJ file.
 *
 * Depending on the time, the skeleton may be incomplete: skeleton edges and facets are interpolated to
 * the requested time. A skeleton facet with fewer than 3 skeleton vertices is not drawn.
 *
 * \tparam GeomTraits must be a model of `Kernel`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param filename the path to the output file
 * \param skeleton the straight skeleton to output
 * \param time the time at which to output the skeleton.
 * \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{the precision of the stream `os`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre `time <= 0`
 */
template <typename GeomTraits,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(const std::string& filename,
               const std::shared_ptr<Straight_skeleton_3<GeomTraits> >& skeleton,
               const typename GeomTraits::FT& time,
               const CGAL_NP_CLASS& np = parameters::default_values())
{
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;
  using Plane_3 = typename GeomTraits::Plane_3;

  using Polyhedron = internal::HDS::Polyhedron<GeomTraits>;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  using Skeleton = Straight_skeleton_3<GeomTraits>;
  using NodeSPtr = typename Skeleton::NodeSPtr;
  using SpokeSPtr = typename Skeleton::SpokeSPtr;
  using SheetWPtr = typename Skeleton::SheetWPtr;
  using SheetSPtr = typename Skeleton::SheetSPtr;

  using Kernel_wrapper = internal::kernel::Kernel_wrapper<GeomTraits>;
  using Transformation = internal::algorithm::Polyhedron_transformation<GeomTraits>;

  // If we want to tolerate intersections for debugging
  // using Itag = CGAL::Exact_intersections_tag;
  // using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
  using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;

  using PK = CGAL::Projection_traits_3<GeomTraits>;
  using PVb = CGAL::Triangulation_vertex_base_2<PK>;
  using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
  using PTDS = CGAL::Triangulation_data_structure_2<PVb, PFb>;
  using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
  using PCDT_FH = typename PCDT::Face_handle;

  CGAL_SS3_IO_TRACE_V(8, "-- OBJ::Save(Skeleton) to " << filename << " --");
  CGAL_SS3_IO_TRACE_V(8, skeleton->nodes().size() << " NN, " << skeleton->sheets().size() << " NS");
  CGAL_SS3_IO_TRACE_V(8, "Time = " << time);

  CGAL_precondition(time <= 0);

  std::stringstream oss;
  set_stream_precision_from_NP(oss, np);

  // Map for unique nodes and their indices
  std::map<Point_3, int> pos_to_index;
  int next_index = 1;

  // Write vertices
  for (const NodeSPtr& node : skeleton->nodes()) {
    if (node->time() < time) {
      continue;
    }

    const Point_3& point = node->point();
    pos_to_index[point] = next_index++;

    oss << "v " << CGAL::to_double(point.x()) << " "
                << CGAL::to_double(point.y()) << " "
                << CGAL::to_double(point.z()) << "\n";
  }


  // incomplete spokes get an interpolated point at the requested time
  CGAL::unordered_flat_map<SpokeSPtr, Point_3> incomplete_spoke_to_point;

  // Write sheets as faces, triangulating if requested
  for (const SheetSPtr& sheet : skeleton->sheets()) {
    CGAL_SS3_IO_TRACE_V(16, "  Sheet # " << sheet->id());
    CGAL_SS3_IO_TRACE_V(16, "  Sheet #nodes " << sheet->nodes().size());
    CGAL_SS3_IO_TRACE_V(16, "  Sheet #spokes " << sheet->spokes().size());
    CGAL_SS3_IO_TRACE_V(16, "  Sheet #contours " << sheet->contours().size());

    const auto& nodes = sheet->nodes();

    const auto first_visible_node = std::find_if(nodes.begin(), nodes.end(),
                                                 [&](const NodeSPtr& node) { return node->time() >= time; });
    if (first_visible_node == nodes.end()) {
      continue;
    }

    // Grab all the 3D positions to insert in the CDT
    std::set<Point_3> unique_points;
    for (const NodeSPtr& node : nodes) {
      if (node->time() < time) {
        continue;
      }

      CGAL_SS3_IO_TRACE_V(64, " Add point from node " << node->id() << " " << node->point() << " @ " << node->time());
      auto res = unique_points.insert(node->point());
      CGAL_warning_msg(res.second, "Node should not be found twice in the sheet's nodes");
    }

    // If the source node exists at the requested time, add incomplete spokes clipped to that time.
    Vector_3 n = sheet->supporting_plane().orthogonal_vector();
    CGAL_assertion(n != CGAL::NULL_VECTOR);

    PK traits(n);
    PCDT pcdt(traits);

    for (const SpokeSPtr& contour : sheet->contours()) {
      NodeSPtr v_src = contour->source();
      NodeSPtr v_tgt = contour->target();
      try {
        CGAL_SS3_IO_TRACE_V(64, " Try to add contour segment between " << v_src->point() << " and " << v_tgt->point());
        pcdt.insert_constraint(v_src->point(), v_tgt->point());
      } catch(const typename PCDT::Intersection_of_constraints_exception&) {
        CGAL_SS3_IO_TRACE_V(1, "Error: Intersection of constraints in sheet triangulation");
        CGAL_warning_msg(false, "Intersections in CDT2 are not allowed");
        return false;
      }
    }

    std::vector<Point_3> incomplete_points;

    for (const SpokeSPtr& spoke : sheet->spokes()) {
      Point_3 ap0, ap1;
      NodeSPtr node_src = spoke->source();
      if (node_src->time() <= time) {
        continue;
      }

      if (spoke->has_target() && spoke->target()->time() >= time) {
        NodeSPtr node_tgt = spoke->target();
        ap0 = node_src->point();
        ap1 = node_tgt->point();
        CGAL_assertion(ap0 != ap1);
      } else {
        CGAL_SS3_IO_TRACE_V(64, " Spoke from " << node_src->point() << " is partial at the requested time");
        Point_3 interpolated_point;
        if (node_src->time() == time) { // detached node: several spokes may share the same final point
          interpolated_point = node_src->point();
        } else {
          auto it = incomplete_spoke_to_point.find(spoke);
          if (it == incomplete_spoke_to_point.end()) {
            std::vector<Plane_3> shifted_planes;
            shifted_planes.reserve(3);

            // these faces are from the static polyhedron!!
            // @fixme really ought to fix this static stuff, or distinguish very clearly...
            CGAL::unordered_flat_set<FacetSPtr> ifs;
            for (const SheetWPtr& sheet_wptr : spoke->sheets()) {
              if (const SheetSPtr sheet = sheet_wptr.lock()) {
                ifs.insert(sheet->get_facet_F());
                ifs.insert(sheet->get_facet_B());
              }
            }
            CGAL_assertion(ifs.size() == 3);

            for (const FacetSPtr& f : ifs) {
              CGAL_SS3_DEBUG_SPTR(f);
              // 'f' is a face in the base polyhedron, it lives at t=0
              shifted_planes.emplace_back(Transformation::shift_plane(f, time));
            }

            CGAL_assertion(shifted_planes.size() == 3);
            const std::optional<Point_3> point = Kernel_wrapper::intersection(shifted_planes[0], shifted_planes[1], shifted_planes[2]);
            CGAL_assertion(point.has_value());
            interpolated_point = point.value();

            CGAL_assertion_code(auto res_1 =)
              incomplete_spoke_to_point.emplace(spoke, interpolated_point);
            CGAL_assertion(res_1.second);

            CGAL_assertion_code(auto res_2 =)
              pos_to_index.emplace(interpolated_point, next_index++);
            CGAL_assertion(res_2.second);

            oss << "v " << CGAL::to_double(interpolated_point.x()) << " "
                        << CGAL::to_double(interpolated_point.y()) << " "
                        << CGAL::to_double(interpolated_point.z()) << "\n";
          } else {
            interpolated_point = it->second;
          }
        }

        CGAL_SS3_IO_TRACE_V(64, " Add interpolated point from incomplete spoke "<< spoke->id() << ": " << interpolated_point);
        CGAL_SS3_IO_TRACE_V(64, " using direction: " << spoke->direction());

        incomplete_points.push_back(interpolated_point);
        ap0 = node_src->point();
        ap1 = interpolated_point;
      }

      if (ap0 == ap1) {
        CGAL_SS3_IO_TRACE_V(64, " Ignoring degenerate spoke segment at " << ap0);
        continue;
      }

      try {
        CGAL_SS3_IO_TRACE_V(64, " Try to add spoke segment between " << ap0 << " and " << ap1);
        pcdt.insert_constraint(ap0, ap1);
      } catch(const typename PCDT::Intersection_of_constraints_exception&) {
        CGAL_SS3_IO_TRACE_V(1, "Warning: Intersection of constraints in sheet triangulation (a)");
        CGAL_warning_msg(false, "Intersections in CDT2 are not allowed");
        return false;
      }
    }

    if (!incomplete_points.empty()) {
      const Plane_3 plane_B = Transformation::shift_plane(sheet->get_facet_B(), time);
      const Plane_3 plane_F = Transformation::shift_plane(sheet->get_facet_F(), time);
      const auto line = Kernel_wrapper::intersection(plane_B, plane_F);
      CGAL_assertion(line.has_value());
      const Point_3 line_origin = line->point();
      const Vector_3 line_direction = line->to_vector();
      CGAL_assertion(line_direction != CGAL::NULL_VECTOR);

      std::sort(incomplete_points.begin(), incomplete_points.end(),
                [&](const Point_3& lhs, const Point_3& rhs) {
                  const FT lhs_coord = CGAL::scalar_product(lhs - line_origin, line_direction);
                  const FT rhs_coord = CGAL::scalar_product(rhs - line_origin, line_direction);
                  return lhs_coord < rhs_coord;
                });

      for (std::size_t i = 0; i + 1 < incomplete_points.size(); i += 2) {
        const Point_3& p0 = incomplete_points[i];
        const Point_3& p1 = incomplete_points[i + 1];
        if (p0 == p1) {
          CGAL_SS3_IO_TRACE_V(64, " Ignoring degenerate closing spoke at " << p0);
          continue;
        }

        try {
          CGAL_SS3_IO_TRACE_V(64, " Try to add closing spoke between " << p0 << " and " << p1);
          pcdt.insert_constraint(p0, p1);
        } catch(const typename PCDT::Intersection_of_constraints_exception&) {
          CGAL_SS3_IO_TRACE_V(1, "Warning: Intersection of constraints in sheet triangulation (b)");
          CGAL_warning_msg(false, "Intersections in CDT2 are not allowed");
          return false;
        }
      }
    }

    if (pcdt.dimension() != 2) {
      CGAL_SS3_IO_TRACE_V(16, "Triangulation of sheet is not 2D");
      CGAL_SS3_IO_TRACE_V(16, sheet->to_string());
      continue;
    }

    // Mark domain and output triangles
    CGAL::unordered_flat_map<PCDT_FH, bool> in_domain_map;
    boost::associative_property_map<CGAL::unordered_flat_map<PCDT_FH, bool>> in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(pcdt, in_domain);

    for (auto fh : pcdt.finite_face_handles()) {
      if (!get(in_domain, fh)) {
        continue;
      }
      oss << "f " << pos_to_index.at(fh->vertex(0)->point()) << " "
                  << pos_to_index.at(fh->vertex(1)->point()) << " "
                  << pos_to_index.at(fh->vertex(2)->point()) << "\n";
    }
  }

  std::ofstream ofs(filename.c_str());
  if (ofs.is_open()) {
    set_stream_precision_from_NP(ofs, np);
    ofs << oss.str();
  } else {
    CGAL_SS3_IO_TRACE_V(1, "Error: failed to open file");
    CGAL_assertion(false);
    return false;
  }

  CGAL_SS3_IO_TRACE_V(8, "-- Write OBJ end --");
  return true;
}

} // namespace IO
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_IO_OBJ_H */
