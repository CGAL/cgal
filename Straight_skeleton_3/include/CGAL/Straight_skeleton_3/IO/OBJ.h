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

  CGAL_SS3_IO_TRACE("-- OBJ::Save(Polyhedron) to " << filename << " --");
  CGAL_SS3_IO_TRACE("   do_triangulate: " << std::boolalpha << do_triangulate);
  CGAL_SS3_IO_TRACE(polyhedron->vertices().size() << " NV, " << polyhedron->facets().size() << " NF");

  // CGAL_SS3_IO_TRACE("Saving to OBJ:\n" << polyhedron->to_string());

  // we can tolerate intersections for OBJ::save because it is only used for debugging
  // CGAL::Exact_intersections_tag
  // CGAL::No_constraint_intersection_requiring_constructions_tag
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
      Vector_3 n = facet->get_plane().orthogonal_vector();
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
          CGAL_SS3_IO_TRACE("Degenerate edge @ " << v0->point());

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
            CGAL_SS3_IO_TRACE("Warning: Intersection of constraints");
            CGAL_SS3_IO_TRACE("While inserting " << v0->point() << " || " << v1->point());
            CGAL_SS3_IO_TRACE(facet->to_string());
            CGAL_assertion_msg(false, "Intersections in CDT2 are not allowed");
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
        CGAL_SS3_IO_TRACE("Warning: skipping degenerate facet");
        continue;
      }

      if (do_triangulate_facet) {
        std::unordered_map<PCDT_FH, bool> in_domain_map;
        boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);

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
        if (visited_edges.find(edge) != visited_edges.end()) {
          continue; // already visited
        }

        std::vector<VertexSPtr> boundary_vertices;
        EdgeSPtr start_edge = edge;
        bool is_open = false;

        // Walk forward to collect boundary vertices
        do {
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
    CGAL_SS3_IO_TRACE("Error: failed to open file");
    CGAL_assertion(false);
    return false;
  }

  CGAL_SS3_IO_TRACE("-- Write OBJ end --");
  return true;
}

/*!
 * \ingroup PkgStraightSkeleton3IO
 *
 * writes a `CGAL::Straight_skeleton_3` to a Wavefront OBJ file.
 *
 * This function outputs only the complete sheets (i.e., no sheet that has partial arcs)
 * of the straight skeleton.
 *
 * \tparam GeomTraits must be a model of `Kernel`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param filename the path to the output file
 * \param skeleton the straight skeleton to output
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
 */
template <typename GeomTraits,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(const std::string& filename,
               std::shared_ptr<Straight_skeleton_3<GeomTraits> > skeleton,
               const CGAL_NP_CLASS& np = parameters::default_values())
{
  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;

  using Skeleton = Straight_skeleton_3<GeomTraits>;
  using NodeSPtr = typename Skeleton::NodeSPtr;
  using ArcSPtr = typename Skeleton::ArcSPtr;
  using SheetSPtr = typename Skeleton::SheetSPtr;

  // For triangulation
  using Itag = CGAL::Exact_intersections_tag;
  using PK = CGAL::Projection_traits_3<GeomTraits>;
  using PVbb = CGAL::Triangulation_vertex_base_with_info_2<NodeSPtr, PK>;
  using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
  using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
  using PTDS = CGAL::Triangulation_data_structure_2<PVb, PFb>;
  using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
  using PCDT_VH = typename PCDT::Vertex_handle;
  using PCDT_FH = typename PCDT::Face_handle;

  CGAL_SS3_IO_TRACE("-- OBJ::Save(Skeleton) to " << filename << " --");
  CGAL_SS3_IO_TRACE(skeleton->nodes().size() << " NN, " << skeleton->sheets().size() << " NS");

  std::stringstream oss;
  set_stream_precision_from_NP(oss, np);

  // Map for unique nodes and their indices
  std::map<NodeSPtr, int> node_to_index;
  int next_index = 1;

  // Collect all unique nodes from all sheets
  for (const NodeSPtr& node : skeleton->nodes()) {
    node_to_index[node] = next_index++;
  }

  // Write vertices
  for (const NodeSPtr& node : skeleton->nodes()) {
    const Point_3& point = node->point();
    CGAL_SS3_IO_TRACE_V(16, "Node's #arcs " << node->degree());
    CGAL_SS3_IO_TRACE_V(16, "Node's #sheets " << node->sheets().size());

    oss << "v " << CGAL::to_double(point.x()) << " "
                << CGAL::to_double(point.y()) << " "
                << CGAL::to_double(point.z()) << "\n";
  }

  // Write sheets as faces, triangulating if requested
  for (const SheetSPtr& sheet : skeleton->sheets()) {
    const auto& nodes = sheet->nodes();

    CGAL_SS3_IO_TRACE("  Sheet #nodes " << nodes.size());
    CGAL_SS3_IO_TRACE("  Sheet #arcs " << sheet->arcs().size());

    if (nodes.size() < 3) {
      CGAL_SS3_IO_TRACE("Sheet with fewer than 3 nodes found. Skipping.");
      continue;
    }

    // ignore if the sheet has incomplete arcs
    bool is_complete = true;
    for (const ArcSPtr& arc : sheet->arcs()) {
      if (!arc->has_target()) {
        is_complete = false;
        break;
      }
    }

    if (!is_complete) {
      CGAL_SS3_IO_TRACE_V(16, "  Incomplete sheet found. Skipping.");
      continue;
    }

    if (nodes.size() == 3) {
      oss << "f ";
      for (const NodeSPtr& node : nodes) {
        oss << node_to_index.at(node) << " ";
      }
      oss << "\n";
    } else {
      Vector_3 n = sheet->get_plane().orthogonal_vector();
      CGAL_assertion(n != CGAL::NULL_VECTOR);

      PK traits(n);
      PCDT pcdt(traits);

      std::map<NodeSPtr, PCDT_VH> face_vhs;
      for (const NodeSPtr& node : nodes) {
        CGAL_SS3_IO_TRACE(" CDT add node " << node->id());
        auto res = face_vhs.emplace(node, PCDT_VH());
        CGAL_warning_msg(res.second, "Node should not be found twice in the sheet's nodes");
        if (res.second) {
          PCDT_VH vh = pcdt.insert(node->point());
          vh->info() = node;
          res.first->second = vh;
        }
      }

      // Insert constraints along the boundary
      for (const ArcSPtr& arc : sheet->arcs()) {
        NodeSPtr node_src = arc->source();
        NodeSPtr node_tgt = arc->target();
        CGAL_SS3_IO_TRACE_V(32, " CDT arc between " << node_src->id() << " and " << node_tgt->id());
        CGAL_SS3_DEBUG_SPTR(node_src);
        CGAL_SS3_DEBUG_SPTR(node_tgt);
        PCDT_VH vh0 = face_vhs.at(node_src);
        PCDT_VH vh1 = face_vhs.at(node_tgt);
        try {
          pcdt.insert_constraint(vh0, vh1);
        } catch(const typename PCDT::Intersection_of_constraints_exception&) {
          CGAL_SS3_IO_TRACE("Warning: Intersection of constraints in sheet triangulation");
          CGAL_assertion_msg(false, "Intersections in CDT2 are not allowed");
          return false;
        }
      }

      for (const ArcSPtr& contour : sheet->contours()) {
        NodeSPtr v_src = contour->source();
        NodeSPtr v_tgt = contour->target();
        try {
          pcdt.insert_constraint(v_src->point(), v_tgt->point());
        } catch(const typename PCDT::Intersection_of_constraints_exception&) {
          CGAL_SS3_IO_TRACE_V(1, "Error: Intersection of constraints in sheet triangulation");
          CGAL_assertion_msg(false, "Intersections in CDT2 are not allowed");
          return false;
        }
      }

      // Mark domain and output triangles
      std::unordered_map<PCDT_FH, bool> in_domain_map;
      boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);
      CGAL::mark_domain_in_triangulation(pcdt, in_domain);

      for (auto fh : pcdt.finite_face_handles()) {
        if (!get(in_domain, fh)) {
          continue;
        }
        oss << "f " << node_to_index.at(fh->vertex(0)->info()) << " "
                    << node_to_index.at(fh->vertex(1)->info()) << " "
                    << node_to_index.at(fh->vertex(2)->info()) << "\n";
      }
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

  CGAL_SS3_IO_TRACE("-- Write OBJ end --");
  return true;
}

} // namespace IO
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_IO_OBJ_H */
