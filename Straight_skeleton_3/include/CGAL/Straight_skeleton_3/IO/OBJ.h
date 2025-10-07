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

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>

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
namespace Straight_skeletons_3 {

namespace internal {
namespace HDS {

template <typename Traits_>
class Polyhedron;

} // namespace HDS

namespace SDS {

template <typename Traits>
class StraightSkeleton;

} // namespace SDS
} // namespace internal

namespace IO {

/**
 * Wavefront obj file format
 */
class OBJFile
{
public:
  /**
    * It is impossible for an obj file to store holes inside a facet,
    * so we triangulate facet in this case. (@todo)
    */
  template <typename Traits>
  static bool save(const std::string& filename,
                   std::shared_ptr<internal::HDS::Polyhedron<Traits> > polyhedron,
                   bool do_triangulate = true,
                   bool convert_to_double = true)
  {
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Point3SPtr = std::shared_ptr<Point_3>;
    using Vector3SPtr = std::shared_ptr<Vector_3>;

    using Polyhedron = internal::HDS::Polyhedron<Traits>;
    using VertexSPtr = typename Polyhedron::VertexSPtr;
    using EdgeSPtr = typename Polyhedron::EdgeSPtr;
    using FacetSPtr = typename Polyhedron::FacetSPtr;

    using KernelFactory = internal::kernel::KernelFactory<Traits>;

    bool result = true;
    CGAL_SS3_IO_TRACE("-- OBJ::Save(Polyhedron) to " << filename << " --");
    CGAL_SS3_IO_TRACE("   do_triangulate: " << std::boolalpha << do_triangulate << "\n" <<
                      "   convert_to_double: " << convert_to_double);
    CGAL_SS3_IO_TRACE(polyhedron->vertices().size() << " NV, " << polyhedron->facets().size() << " NF");

    // CGAL_SS3_IO_TRACE("Saving to OBJ:\n" << polyhedron->toString());

    // we can tolerate intersections for OBJ::save because it is only used for debugging
    // CGAL::Exact_intersections_tag
    // CGAL::No_constraint_intersection_requiring_constructions_tag
    using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
    using PK = CGAL::Projection_traits_3<Traits>;
    using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
    using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
    using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using PTDS = CGAL::Triangulation_data_structure_2<PVb, PFb>;
    using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
    using PCDT_VH = typename PCDT::Vertex_handle;
    using PCDT_FH = typename PCDT::Face_handle;

    std::stringstream oss;
    oss.precision(17);

    // Map for unique vertices and their indices
    std::map<VertexSPtr, int> vertex_to_index;
    int next_index = 1; // OBJ indices start at 1

    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      vertex_to_index[vertex] = next_index++;
    }

    // Write vertices
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      Point3SPtr pt = vertex->getPoint();
      if (convert_to_double) {
        oss << "v " << CGAL::to_double(pt->x()) << " "
                    << CGAL::to_double(pt->y()) << " "
                    << CGAL::to_double(pt->z()) << "\n";
      } else {
        oss << "v " << pt->x().exact() << " "
                    << pt->y().exact() << " "
                    << pt->z().exact() << "\n";
      }
    }

    // Write facets
    for (const FacetSPtr& facet : polyhedron->facets()) {
      bool do_triangulate_facet = do_triangulate;
      if (facet->edges().size() <= 3) {
        do_triangulate_facet = false;
      }

      if (do_triangulate_facet) {
        Vector3SPtr n = KernelFactory::createVector3(facet->plane());
        CGAL_assertion(*n != CGAL::NULL_VECTOR);

        PK traits(*n);
        PCDT pcdt(traits);

        std::map<VertexSPtr, PCDT_VH> face_vhs;

        for (const VertexSPtr& vertex : facet->vertices()) {
          auto res = face_vhs.emplace(vertex, PCDT_VH());
          if (res.second) { // first time seeing this point
            PCDT_VH vh = pcdt.insert(*(vertex->getPoint()));
            vh->info() = vertex;
            res.first->second = vh;
          }
        }

        auto ne = 0;
        for (const EdgeSPtr& edge : facet->edges()) {
          VertexSPtr v0 = edge->src(facet);
          VertexSPtr v1 = edge->dst(facet);

          if (*(v0->getPoint()) == *(v1->getPoint())) {
            CGAL_SS3_IO_TRACE("Warning: encountered degenerate edge @ " << *(v0->getPoint()));

            CGAL_assertion(v0->degree() != 1); // @todo handle that...
            VertexSPtr vm1 = edge->prev(facet)->src(facet);

            // manually create a degenerate facet so that the resulting mesh is conforming
            oss << "f " << vertex_to_index[vm1] << " "
                        << vertex_to_index[v0] << " "
                        << vertex_to_index[v1] << "\n";
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
              CGAL_SS3_IO_TRACE("While inserting " << *(v0->getPoint()) << " || " << *(v1->getPoint()));
              CGAL_SS3_IO_TRACE(facet->toString());
              CGAL_assertion_msg(false, "Intersections in CDT2 are not allowed");
              return false;
            }
            ++ne;
          }
        }

        if (pcdt.finite_vertex_handles().size() != facet->vertices().size()) {
          // CGAL_SS3_IO_TRACE("Warning: CDT #nv != facet #nv. Constraint intersection?");
          // CGAL_SS3_IO_TRACE("Facet: " << facet->toString());
          // CGAL_SS3_IO_TRACE("CDT: " << pcdt.number_of_vertices() << " vertices, "
          //                           << pcdt.number_of_faces() << " faces");
          // CGAL_SS3_IO_TRACE("Vertices: ");
          // CGAL_SS3_IO_TRACE_CODE(for (PCDT_VH vh : pcdt.finite_vertex_handles()) {)
          // CGAL_SS3_IO_TRACE("  " << vh->point() << " -> " << vh->info()->getID());
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
            oss << "f " << vertex_to_index[fh->vertex(0)->info()] << " "
                        << vertex_to_index[fh->vertex(1)->info()] << " "
                        << vertex_to_index[fh->vertex(2)->info()] << "\n";
          }
        }
      }

      if (!do_triangulate_facet) {
        std::set<EdgeSPtr> visited_edges;

        for (const EdgeSPtr& edge : facet->edges()) {
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
            if (edge->dst(facet)->degree() == 1) {
              is_open = true;
              break;
            }
              edge = edge->next(facet);
          } while (edge != start_edge);

          // If open, also walk backward to collect remaining boundary vertices
          if (is_open) {
            boundary_vertices.push_back(edge->dst(facet));
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
            oss << vertex_to_index[vertex] << " ";
          }

          oss << "\n";
        }
      }
    }

    // result = (polyhedron->vertices().size() == vertex_id &&
    //           polyhedron->facets().size() == facet_id);

    std::ofstream ofs(filename.c_str());
    if (ofs.is_open()) {
        ofs.precision(17);
        ofs << oss.str();
    } else {
        CGAL_SS3_IO_TRACE("Error: failed to open file");
        CGAL_assertion(false);
    }

    CGAL_SS3_IO_TRACE("-- Write OBJ end --");
    return result;
  }
};

} // namespace IO
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_IO_OBJ_H */
