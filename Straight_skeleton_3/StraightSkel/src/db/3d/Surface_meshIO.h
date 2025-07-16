// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef DB_3D_SURFACE_MESH_H
#define DB_3D_SURFACE_MESH_H

#include "debug.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Triangle.h"
#include "data/3d/KernelFactory.h"
#include "db/3d/AbstractFile.h"
#include "util/StringFuncs.h"
#include "util/Configuration.h"

#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SkelFacetData.h"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/property_map.h>

#include <cmath>
#include <exception>
#include <fstream>
#include <vector>
#include <unordered_map>

namespace db { namespace _3d {

using namespace data::_3d;

/**
 * PLY file format
 */
class Surface_meshIO : public AbstractFile
{
    using Mesh = CGAL::Surface_mesh<Point3>;

    using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using edge_descriptor = typename boost::graph_traits<Mesh>::edge_descriptor;
    using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

public:
    virtual ~Surface_meshIO();

    template <typename NamedParameters>
    static PolyhedronSPtr load(const Mesh& sm,
                               std::map<edge_descriptor, EdgeWPtr>& e2e,
                               const NamedParameters& np);

    template <typename NamedParameters = CGAL::parameters::Default_named_parameters>
    static PolyhedronSPtr load(const Mesh& sm,
                               const NamedParameters& np = CGAL::parameters::default_values());

    template <typename NamedParameters = CGAL::parameters::Default_named_parameters>
    static bool save(const PolyhedronSPtr& polyhedron,
                     CGAL::Surface_mesh<Point3>& sm,
                     const NamedParameters& np = CGAL::parameters::default_values());

protected:
    Surface_meshIO();
};

template <typename NamedParameters>
PolyhedronSPtr Surface_meshIO::load(const CGAL::Surface_mesh<Point3>& sm,
                                    std::map<edge_descriptor, EdgeWPtr>& e2e,
                                    const NamedParameters& np)
{
    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::get_parameter;

    using Face_index = CGAL::Surface_mesh<Point3>::Face_index;

    CGAL_warning(CGAL::is_triangle_mesh(sm));

    PolyhedronSPtr result = Polyhedron::create();

    unsigned int vertex_id_new = 0;
    std::vector<VertexSPtr> vertices;

    for(vertex_descriptor vi : CGAL::vertices(sm)) {
        vertex_id_new++;
        Point3SPtr point = KernelFactory::createPoint3(sm.point(vi));
        VertexSPtr vertex = Vertex::create(point);
        vertex->setID(vertex_id_new);
        result->addVertex(vertex);
    }

    vertices = std::vector<VertexSPtr>(result->vertices().begin(), result->vertices().end());

    auto weight_pmap = choose_parameter(get_parameter(np, CGAL::internal_np::face_weight),
                                        CGAL::Constant_property_map<std::size_t, double>(1.0));

    unsigned int facet_id_new = 0;
    for(face_descriptor fi : faces(sm)) {
        facet_id_new++;
        unsigned int num_vertices = sm.degree(fi);
        // std::cout << "new face of size " << num_vertices << std::endl;
        CGAL_assertion(num_vertices > 2);

        VertexSPtr poly_vertices[num_vertices];
        Vector3SPtr normal_sum;
        for (unsigned int i = 0; i < num_vertices; i++) {
            poly_vertices[i] = VertexSPtr();
        }

        unsigned int pos = 0;
        for (halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(fi, sm), sm)) {
            unsigned int vertex_id = source(h, sm);
            if (vertex_id < vertices.size()) {
                poly_vertices[pos++] = vertices[vertex_id];
            } else {
                std::stringstream whatstream;
                whatstream << "Vertex with id="
                        << vertex_id
                        << " does not exist.";
                throw std::runtime_error(whatstream.str());
            }
        }

        FacetSPtr facet = Facet::create(num_vertices, poly_vertices);
        facet->setID(facet_id_new);

        // Correspondence between the edges of the original mesh and the new edges
        // in the polyhedron
        // poly_vertices is filled starting at source() of the first edge
        // Facet::create() creates the i-th edge between vertices[i] and vertices[i+1]
        halfedge_descriptor h = halfedge(fi, sm);
        std::list<EdgeSPtr>::iterator eit = facet->edges().begin();
        for (unsigned int i = 0; i < num_vertices; ++i) {
            EdgeSPtr e = *eit++;
            e2e[edge(h, sm)] = e;
            h = next(h, sm);
        }

        if (num_vertices == 3) {
            Plane3SPtr plane = KernelFactory::createPlane3(
                    poly_vertices[0]->getPoint(),
                    poly_vertices[1]->getPoint(),
                    poly_vertices[2]->getPoint());
            facet->setPlane(plane);
        } else {
            // @todo is there a point handling non triangulated inputs here and everywhere...?
            return { };
        }
        result->addFacet(facet);

        CGAL::FT weight = 1;
        if (weight_pmap) {
            weight = get(weight_pmap, fi);
        } else {
            std::cerr << "Warning: no weights in Surface_mesh being read?" << std::endl;
        }

        CGAL_assertion(weight != 0);
        data::_3d::skel::SkelFacetDataSPtr data = data::_3d::skel::SkelFacetData::create(facet);
        data->setSpeed(weight);
    }

    std::list<EdgeSPtr>::iterator it_e = result->edges().begin();
    while (it_e != result->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (!(edge->getFacetL() && edge->getFacetR())) {
            DEBUG_PRINT("Warning: Polyhedron has no closed boundary.");
            DEBUG_PRINT(edge->toString());
        }
    }

    util::ConfigurationSPtr config = util::Configuration::getInstance();
    std::string section("db_3d_PLYFile");
    if (config->isLoaded() &&
        config->contains(section, "merge_coplanar_faces") &&
        config->getBool(section, "merge_coplanar_faces")) {
        double epsilon = 0.;
        std::string key("epsilon_coplanarity");
        if (config->contains(section, key)) {
            epsilon = config->getDouble(section, key);
        }
        mergeCoplanarFacets(result, epsilon);
    }

    removeVerticesDegLt3(result);

    CGAL_postcondition(result->isConsistent());

    return result;
}

template <typename NamedParameters>
PolyhedronSPtr Surface_meshIO::load(const CGAL::Surface_mesh<Point3>& sm,
                                    const NamedParameters& np)
{
  std::map<edge_descriptor, EdgeWPtr> unused_e2e;
  return load(sm, unused_e2e, np);
}

template <typename NamedParameters>
bool Surface_meshIO::save(const PolyhedronSPtr& polyhedron,
                          CGAL::Surface_mesh<Point3>& sm,
                          const NamedParameters& np)
{
    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::get_parameter;

   // @fixme factorize with OBJ/PLY save()
    using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
    using PK = CGAL::Projection_traits_3<CGAL::K>;
    using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
    using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
    using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using PTDS = CGAL::Triangulation_data_structure_2<PVb, PFb>;
    using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
    using PCDT_VH = PCDT::Vertex_handle;
    using PCDT_FH = PCDT::Face_handle;

    std::cout << "--> SM of polyhedron with " << polyhedron->vertices().size() << " vertices and "
              << polyhedron->facets().size() << " facets" << std::endl;

    bool do_triangulate = !choose_parameter(get_parameter(np, CGAL::internal_np::do_not_triangulate_faces), false);

    // Vertices
    std::unordered_map<VertexSPtr, vertex_descriptor> v_map;
    for (VertexSPtr v : polyhedron->vertices()) {
        vertex_descriptor idx = sm.add_vertex(*(v->getPoint()));
        v_map[v] = idx;
    }

    // Write facets
    auto weight_pmap = choose_parameter(get_parameter(np, CGAL::internal_np::face_weight),
                                        CGAL::Constant_property_map<std::size_t, double>(1.0));

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        CGAL_assertion(facet->getID() != -1);
        // std::cout << "handle facet " << facet->getID() << std::endl;

        double speed = 1.;
        if (facet->hasData()) {
            auto skel_data = std::dynamic_pointer_cast<data::_3d::skel::SkelFacetData>(facet->getData());
            if (skel_data) {
                speed = CGAL::to_double(skel_data->getSpeed());
            }
        }

        bool do_triangulate_face = do_triangulate;
        if (facet->edges().size() < 3)
            do_triangulate_face = false;

        if (do_triangulate_face)
        {
            Vector3SPtr n = KernelFactory::createVector3(facet->plane());

            // @todo might have to do something fancier than staring from a single vertex
            // for degenerate faces with zigzagging edges...
            CGAL_assertion(*n != CGAL::NULL_VECTOR);

            PK traits(*n);
            PCDT pcdt(traits);

            std::map<VertexSPtr, PCDT_VH> face_vhs;

            std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
            while (it_v != facet->vertices().end()) {
                VertexSPtr vertex = *it_v++;
                auto res = face_vhs.emplace(vertex, PCDT_VH());
                if(res.second) // first time seeing this point
                {
                    PCDT_VH vh = pcdt.insert(*(vertex->getPoint()));
                    vh->info() = vertex;
                    res.first->second = vh;
                }
            }

            auto ne = 0;
            std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
            while (it_e != facet->edges().end()) {
                EdgeSPtr edge = *it_e++;
                VertexSPtr v0 = edge->src(facet);
                VertexSPtr v1 = edge->dst(facet);

                if(*(v0->getPoint()) == *(v1->getPoint()))
                {
                    // std::cerr << "W: encountered degenerate edge @ " << *(v0->getPoint()) << std::endl;

                    CGAL_assertion(v0->degree() != 1);
                    VertexSPtr vm1 = edge->prev(facet)->src(facet);

                    face_descriptor sm_f = sm.add_face(v_map[vm1], v_map[v0], v_map[v1]);
                    if (sm_f == boost::graph_traits<CGAL::Surface_mesh<Point3>>::null_face()) {
                        std::cerr << "Error: failed to add face to surface mesh." << std::endl;
                        return false;
                    }

                    put(weight_pmap, sm_f, speed);
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
                        DEBUG_PRINT("Error: Intersection of constraints");
                        DEBUG_PRINT("While inserting " << *(v0->getPoint()) << " || " << *(v1->getPoint()));
                        DEBUG_PRINT(facet->toString());
                        CGAL_warning_msg(false, "Intersections in CDT2 not allowed");
                        do_triangulate_face = false;
                        break;
                    }
                    ++ne;
                }
            }

            if(ne < 3) // degenerate face
            {
                std::cerr << "Warning: skipping degenerate face" << std::endl;
                continue;
            }

            if(do_triangulate_face)
            {
                std::unordered_map<PCDT_FH, bool> in_domain_map;
                boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);

                CGAL::mark_domain_in_triangulation(pcdt, in_domain);

                for(auto fh : pcdt.finite_face_handles())
                {
                    if(!get(in_domain, fh))
                      continue;

                    face_descriptor sm_f = sm.add_face(v_map[fh->vertex(0)->info()],
                                                       v_map[fh->vertex(1)->info()],
                                                       v_map[fh->vertex(2)->info()]);
                    if (sm_f == boost::graph_traits<CGAL::Surface_mesh<Point3>>::null_face()) {
                        std::cerr << "Error: failed to add face to surface mesh." << std::endl;
                        return false;
                    }

                    put(weight_pmap, sm_f, speed);
                }
            }
        }

        if(!do_triangulate_face)
        {
            std::set<EdgeSPtr> visited_edges;

            std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
            while (it_e != facet->edges().end()) {
                EdgeSPtr edge = *it_e++;
                if (visited_edges.find(edge) != visited_edges.end()) {
                    continue; // already visited
                }

                std::vector<vertex_descriptor> boundary_vertices;
                EdgeSPtr start_edge = edge;
                // std::cout << "start a walk @ " << start_edge->src(facet)->getID() << " " << start_edge->dst(facet)->getID() << std::endl;
                bool is_open = false;

                // Walk forward to collect boundary vertices
                do {
                    visited_edges.insert(edge);
                    boundary_vertices.push_back(v_map[edge->src(facet)]);
                    if (edge->dst(facet)->degree() == 1) {
                        is_open = true;
                        break;
                    }
                    edge = edge->next(facet);
                } while (edge != start_edge);

                // If open, also walk backward to collect remaining boundary vertices
                if (is_open) {
                    boundary_vertices.push_back(v_map[edge->dst(facet)]);
                    edge = start_edge;
                    while (edge->src(facet)->degree() != 1) {
                        edge = edge->prev(facet);
                        visited_edges.insert(edge);
                        boundary_vertices.insert(boundary_vertices.begin(), v_map[edge->src(facet)]);
                    }
                }

                // Write the boundary as a face
                face_descriptor sm_f = sm.add_face(boundary_vertices);
                if (sm_f == boost::graph_traits<CGAL::Surface_mesh<Point3>>::null_face()) {
                    std::cerr << "Error: failed to add face to surface mesh." << std::endl;
                    return false;
                }
                put(weight_pmap, sm_f, speed);
            }
        }
    }

    return true;
}

} }

#endif /* DB_3D_SURFACE_MESH_H */
