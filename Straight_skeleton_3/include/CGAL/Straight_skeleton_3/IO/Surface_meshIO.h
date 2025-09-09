// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_IO_SURFACE_MESH_IO_H
#define CGAL_STRAIGHT_SKELETON_3_IO_SURFACE_MESH_IO_H

#include <CGAL/Straight_skeleton_3/Configuration.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>
#if 0
# include <CGAL/Polygon_mesh_processing/region_growing.h>
#endif

#include <cmath>
#include <exception>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <unordered_map>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace IO {

template <typename Traits>
class Surface_meshIO
{
  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;
  using Vector_3 = typename Traits::Vector_3;
  using Plane_3 = typename Traits::Plane_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Vector3SPtr = std::shared_ptr<Vector_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;

private:
  using Mesh = CGAL::Surface_mesh<Point_3>;

  using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  using edge_descriptor = typename boost::graph_traits<Mesh>::edge_descriptor;
  using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

private:
  using Polyhedron = internal::HDS::Polyhedron<Traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Vertex = typename Polyhedron::template Vertex<Traits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using Facet = typename Polyhedron::template Facet<Traits>;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  using SkelFacetData = typename Polyhedron::SkelFacetData;
  using SkelFacetDataSPtr = typename Polyhedron::SkelFacetDataSPtr;

private:
  using KernelFactory = internal::kernel::KernelFactory<Traits>;
  using PolyhedronTransformation = internal::algorithm::PolyhedronTransformation<Traits>;

public:
  template <typename NamedParameters>
  static PolyhedronSPtr load(const Mesh& sm,
                             std::map<edge_descriptor, EdgeWPtr>& e2e,
                             const NamedParameters& np)
  {
    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::get_parameter;

    CGAL_warning(CGAL::is_triangle_mesh(sm));

    PolyhedronSPtr result = Polyhedron::create();

    unsigned int vertex_id_new = 0;
    std::vector<VertexSPtr> vertices;

    for (vertex_descriptor vi : CGAL::vertices(sm)) {
      ++vertex_id_new;
      Point3SPtr point = KernelFactory::createPoint3(sm.point(vi));
      VertexSPtr vertex = Vertex::create(point);
      vertex->setID(vertex_id_new);
      result->addVertex(vertex);
    }

    vertices = std::vector<VertexSPtr>(result->vertices().begin(), result->vertices().end());

    auto weight_pmap = choose_parameter(get_parameter(np, CGAL::internal_np::face_weight),
                                        CGAL::Constant_property_map<std::size_t, double>(1.0));

    int facet_id_new = -1;
    for (face_descriptor fi : faces(sm)) {
      ++facet_id_new;

      unsigned int num_vertices = sm.degree(fi);
      CGAL_SS3_IO_TRACE("new: F" << facet_id_new << " with " << num_vertices << " vertices");
      CGAL_assertion(num_vertices > 2);

      std::vector<VertexSPtr> poly_vertices(num_vertices);
      Vector3SPtr normal_sum;
      for (unsigned int i = 0; i < num_vertices; ++i) {
        poly_vertices[i] = VertexSPtr();
      }

      unsigned int pos = 0;
      for (halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(fi, sm), sm)) {
        unsigned int vertex_id = source(h, sm);
        if (vertex_id < vertices.size()) {
          poly_vertices[pos++] = vertices[vertex_id];
          CGAL_SS3_IO_TRACE("  V" << vertices[vertex_id]->getID() << "; "
                                  << *(vertices[vertex_id]->getPoint()));
        } else {
          std::stringstream whatstream;
          whatstream << "Vertex with id=" << vertex_id << " does not exist.";
          throw std::runtime_error(whatstream.str());
        }
      }

      FacetSPtr facet = Facet::create(poly_vertices);
      facet->setID(facet_id_new);

      // Correspondence between the edges of the original mesh and the new edges
      // in the polyhedron
      // poly_vertices is filled starting at source() of the first edge
      // Facet::create() creates the i-th edge between vertices[i] and vertices[i+1]
      halfedge_descriptor h = halfedge(fi, sm);
      typename std::list<EdgeSPtr>::iterator eit = facet->edges().begin();
      for (unsigned int i = 0; i < num_vertices; ++i) {
        EdgeSPtr e = *eit++;
        e2e[edge(h, sm)] = e;
        h = next(h, sm);
      }

      if (num_vertices == 3) {
        Plane3SPtr plane = KernelFactory::createPlane3(poly_vertices[0]->getPoint(),
                                                       poly_vertices[1]->getPoint(),
                                                       poly_vertices[2]->getPoint());
        facet->setPlane(plane);
      } else {
        // @todo is there a point handling non triangulated inputs here and everywhere...?
        return { };
      }
      result->addFacet(facet);

      FT weight = 1;
      if (weight_pmap) {
        weight = get(weight_pmap, fi);
      } else {
        std::cerr << "Warning: no weights in Surface_mesh being read?" << std::endl;
      }

      CGAL_assertion(weight != 0);
      SkelFacetDataSPtr data = SkelFacetData::create(facet);
      data->setSpeed(weight);
    }

    for (EdgeSPtr edge : result->edges()) {
      if (!(edge->getFacetL() && edge->getFacetR())) {
        CGAL_SS3_IO_TRACE("Warning: Polyhedron has no closed boundary.");
        CGAL_SS3_IO_TRACE(edge->toString());
      }
    }

    CGAL_postcondition(bool(result));
    CGAL_postcondition(result->isConsistent());

    ConfigurationSPtr config = Configuration::getInstance();
    std::string section("db_3d_PLYFile");
    if (config->isLoaded() &&
      config->contains(section, "merge_coplanar_faces") &&
      config->getBool(section, "merge_coplanar_faces")) {
      double epsilon = 0.;
      std::string key("epsilon_coplanarity");
      if (config->contains(section, key)) {
        epsilon = config->getDouble(section, key);
      }
      PolyhedronTransformation::mergeCoplanarFacets(result, epsilon);
    }

    PolyhedronTransformation::removeVerticesDegLt3(result);

    CGAL_postcondition(bool(result));
    CGAL_postcondition(result->isConsistent());

    return result;
  }

  template <typename NamedParameters = CGAL::parameters::Default_named_parameters>
  static PolyhedronSPtr load(const Mesh& sm,
                             const NamedParameters& np = CGAL::parameters::default_values())
  {
    std::map<edge_descriptor, EdgeWPtr> unused_e2e;
    return load(sm, unused_e2e, np);
  }

  template <typename NamedParameters = CGAL::parameters::Default_named_parameters>
  static bool save(const PolyhedronSPtr& polyhedron,
                   Mesh& sm,
                   const NamedParameters& np = CGAL::parameters::default_values())
  {
    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::get_parameter;

    using Itag = CGAL::Exact_intersections_tag;
    using PK = CGAL::Projection_traits_3<Traits>;
    using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
    using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
    using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using PTDS = CGAL::Triangulation_data_structure_2<PVb, PFb>;
    using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
    using PCDT_VH = typename PCDT::Vertex_handle;
    using PCDT_FH = typename PCDT::Face_handle;

    CGAL_SS3_IO_TRACE_V(8, "Save polyhedron with " << polyhedron->vertices().size() << " vertices and "
                                                   << polyhedron->facets().size() << " facets");

    // @todo do not systematically triangulate, but use this NP and if it is false,
    // only triangulate what is not representable otherwise (see code in PMP::remesh_planar_faces)
    // bool do_triangulate = !choose_parameter(get_parameter(np, CGAL::internal_np::do_not_triangulate_faces), false);

    // Vertices
    std::unordered_map<VertexSPtr, vertex_descriptor> v_map;
    for (VertexSPtr v : polyhedron->vertices()) {
      vertex_descriptor idx = sm.add_vertex(*(v->getPoint()));
      v_map[v] = idx;
    }

    // Write facets
    auto weight_pmap = choose_parameter(get_parameter(np, CGAL::internal_np::face_weight),
                                        CGAL::Constant_property_map<std::size_t, double>(1.0));

    for (FacetSPtr facet : polyhedron->facets()) {
      CGAL_assertion(facet->getID() != -1);
      CGAL_assertion(facet->edges().size() >= 3);

      double speed = 1.;
      if (facet->hasData()) {
        SkelFacetDataSPtr skel_data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
        if (skel_data) {
          speed = CGAL::to_double(skel_data->getSpeed());
        }
      }

      Vector3SPtr n = KernelFactory::createVector3(facet->plane());
      CGAL_assertion(*n != CGAL::NULL_VECTOR);

      PK traits(*n);
      PCDT pcdt(traits);

      std::map<VertexSPtr, PCDT_VH> face_vhs;

      for (VertexSPtr vertex : facet->vertices()) {
        auto res = face_vhs.emplace(vertex, PCDT_VH());
        if (res.second) { // first time seeing this point
          PCDT_VH vh = pcdt.insert(*(vertex->getPoint()));
          res.first->second = vh;
          vh->info() = vertex;
        }
      }

      unsigned int ne = 0;
      for (EdgeSPtr edge : facet->edges()) {
        VertexSPtr v0 = edge->src(facet);
        VertexSPtr v1 = edge->dst(facet);

        if(*(v0->getPoint()) == *(v1->getPoint()))
        {
          CGAL_SS3_IO_TRACE("Warning: degenerate edge at " << *(v0->getPoint()));
          break;
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
            CGAL_SS3_IO_TRACE("Error: Intersection of constraints");
            CGAL_SS3_IO_TRACE("While inserting " << *(v0->getPoint()) << " || " << *(v1->getPoint()));
            CGAL_SS3_IO_TRACE(facet->toString());
            CGAL_assertion_msg(false, "Intersections in CDT2 are not allowed");
            return false;
          }
          ++ne;
        }
      }

      if(ne < 3) // degenerate face
      {
        CGAL_SS3_IO_TRACE("Warning: skipping degenerate face");
        continue;
      }

      std::unordered_map<PCDT_FH, bool> in_domain_map;
      boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);
      CGAL::mark_domain_in_triangulation(pcdt, in_domain);

      for (auto fh : pcdt.finite_face_handles()) {
        if(!get(in_domain, fh))
          continue;

        auto vr = CGAL::make_array(v_map[fh->vertex(0)->info()],
                                   v_map[fh->vertex(1)->info()],
                                   v_map[fh->vertex(2)->info()]);
        face_descriptor sm_f = CGAL::Euler::add_face(vr, sm);
        if (sm_f == boost::graph_traits<Mesh>::null_face()) {
          std::cerr << "Error: failed to add face to surface mesh (2)" << std::endl;
          std::cerr << "Face:\n" << fh->vertex(0)->point() << " "
                                 << fh->vertex(1)->point() << " "
                                 << fh->vertex(2)->point() << std::endl;
          CGAL::IO::write_polygon_mesh("results/failed.off", sm, CGAL::parameters::stream_precision(17));

          CGAL_assertion_code(bool cannot_add = CGAL::Euler::can_add_face(vr, sm, true /*verbose*/);)
          CGAL_assertion(!cannot_add);

          return false;
        }

        put(weight_pmap, sm_f, speed);
      }
    }

    return true;
  }

  template <typename NamedParameters = CGAL::parameters::Default_named_parameters>
  static PolyhedronSPtr convert(const CGAL::Surface_mesh<Point_3>& sm,
                                const NamedParameters& np = CGAL::parameters::default_values())
  {
    CGAL_SS3_TRANSF_TRACE("Converting mesh...");

    bool merge_faces = false;

    ConfigurationSPtr config = Configuration::getInstance();
    std::string section("main");
    if (config->isLoaded() &&
      config->contains(section, "merge_coplanar_faces") &&
      config->getBool(section, "merge_coplanar_faces")) {
      merge_faces = true;
    }

    if (!merge_faces) {
      return IO::Surface_meshIO<Traits>::load(sm, np);
    }

#if 1
    PolyhedronSPtr polyhedron = IO::Surface_meshIO<Traits>::load(sm, np);
    PolyhedronTransformation::mergeCoplanarFacets(polyhedron);
#else
    namespace PMP = CGAL::Polygon_mesh_processing;

    CGAL::Bbox_3 bbox = PMP::bbox(sm);
    const FT diag_length = CGAL::approximate_sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                                  CGAL::square(bbox.ymax() - bbox.ymin()) +
                                                  CGAL::square(bbox.zmax() - bbox.zmin()));

    // Use shape detection to analyze the mesh
    std::vector<std::size_t> region_ids(num_faces(sm));
    boost::vector_property_map<Plane_3> plane_map; // supporting planes of the regions detected

    const FT cos_of_max_angle = 0.98;
    const FT max_distance = 0.0001 * diag_length;

    // detect planar regions in the mesh
    // @todo growing should:
    // - use the .ini value of 'epsilon_coplanarity'
    // - stop if it merges faces with different weights
    // - give an error for adjacent coplanar faces that have different weights
    std::size_t nb_regions =
        PMP::region_growing_of_planes_on_faces(sm,
                                               CGAL::make_random_access_property_map(region_ids),
                                               CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle)
                                                               .region_primitive_map(plane_map)
                                                               .maximum_distance(max_distance));

    static int region_dump_id = -1;
    utils::save_colored_mesh(sm, region_ids, "results/regions_" + std::to_string(++region_dump_id) + ".ply");

    // detect corner vertices on the boundary of planar regions
    std::vector<std::size_t> corner_ids(num_vertices(sm), -1); // corner status of vertices
    std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions

    std::size_t nb_corners =
        PMP::detect_corners_of_regions(sm,
                                      CGAL::make_random_access_property_map(region_ids),
                                      nb_regions,
                                      CGAL::make_random_access_property_map(corner_ids),
                                      CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle).
                                                        maximum_distance(max_distance).
                                                        edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

    CGAL_SS3_TRANSF_TRACE_CODE(for (face_descriptor f : faces(sm)))
    CGAL_SS3_TRANSF_TRACE("facet " << f << " is in region " << region_ids[f]);

    // the almost-coplanar merge is performed after the conversion to the Polyhedron
    // data structure because we want to be able to create faces that have holes,
    // which the CGAL::Surface_mesh class does not support
    std::map<edge_descriptor, EdgeWPtr> e2e;
    PolyhedronSPtr polyhedron = db::_3d::Surface_meshIO::load(sm, e2e, np);

    // merge the facets incident to an unconstrained edge (i.e., the edge is interior to a region)
    for (edge_descriptor e: edges(sm)) {
        if (ecm[e]) {
          continue;
        }

        EdgeSPtr edge = e2e[e].lock();
        if (!edge) {
          continue;
        }

        CGAL_SS3_TRANSF_TRACE("Merging facets " << edge->getFacetL()->getID() << " and " << edge->getFacetR()->getID());
        CGAL_assertion(sm.point(source(e, sm)) == *(edge->getVertexSrc()->getPoint()));
        CGAL_assertion(sm.point(target(e, sm)) == *(edge->getVertexDst()->getPoint()));

        // @fixme it seems like intermediate states are somewhat unsound during edge merging
        mergeFacets(edge, polyhedron);
    }

    polyhedron->initializeAllIDs();

    sanitize(polyhedron);
#endif

    CGAL_SS3_TRANSF_TRACE("Converted, " << polyhedron->facets().size() << " facets");

    return polyhedron;
  }
};

} // namespace IO
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_IO_SURFACE_MESH_IO_H */
