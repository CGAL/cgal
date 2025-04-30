// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#include "db/3d/PLYFile.h"

#include "debug.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Triangle.h"
#include "data/3d/KernelFactory.h"
#include "util/StringFuncs.h"
#include "util/Configuration.h"

#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SkelFacetData.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO/PLY.h>

#include <cmath>
#include <exception>
#include <fstream>
#include <sstream>
#include <vector>

namespace db { namespace _3d {

PLYFile::PLYFile() {
    // intentionally does nothing
}

PLYFile::~PLYFile() {
    // intentionally does nothing
}

PolyhedronSPtr PLYFile::load(const CGAL::Surface_mesh<Point3>& sm,
                             const std::string& description,
                             std::map<edge_descriptor, EdgeWPtr>& e2e)
{
    typedef CGAL::Surface_mesh<Point3>::Face_index Face_index;

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

    auto weight_pmap_res = sm.property_map<Face_index, double>("f:weight");

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

        // std::cout << "Final edges:" << std::endl;
        // for(edge_descriptor e : facet->edges())
        //     std::cout << e->toString() << std::endl;

        if (num_vertices == 3) {
            Plane3SPtr plane = KernelFactory::createPlane3(
                    poly_vertices[0]->getPoint(),
                    poly_vertices[1]->getPoint(),
                    poly_vertices[2]->getPoint());
            facet->setPlane(plane);
        } else if (normal_sum && num_vertices > 3) {
            // vertex normals are used as a hint for facet normal
            // @fixme no reason for these 3 points not to be collinear?
            Plane3SPtr plane = KernelFactory::createPlane3(
                    poly_vertices[0]->getPoint(),
                    poly_vertices[1]->getPoint(),
                    poly_vertices[2]->getPoint());
            Vector3SPtr normal_plane = KernelFactory::createVector3(plane);
            double angle = 0.0;
            CGAL::FT arg = 0.0;
#ifdef USE_CGAL
            // @fixme tolerating this sqrt for now, but this should find an extremum vertex
            // in the face plane, and then do an orientation test with prev/next
            arg = ((*normal_plane)*(*normal_sum)) /
                    CGAL::sqrt_with_warning(normal_plane->squared_length() * normal_sum->squared_length());
#else
            arg = ((*normal_plane)*(*normal_sum)) /
                    sqrt(normal_plane->squared_length() * normal_sum->squared_length());
#endif
            // fixes issues with floating point precision
            // @fixme rewrite this without trigonometry
            if (arg <= -1.0) {
                angle = M_PI;
            } else if (arg >= 1.0) {
                angle = 0.0;
            } else {
                angle = acos(CGAL::to_double(arg));
            }
            if (angle > M_PI/2.0) {
                plane = KernelFactory::createPlane3(
                        poly_vertices[2]->getPoint(),
                        poly_vertices[1]->getPoint(),
                        poly_vertices[0]->getPoint());
            }
            facet->setPlane(plane);
            facet->makeFirstConvex();
        }
        else {
            facet->initPlane(); // num_vertices > 3 but no normals provided
            facet->makeFirstConvex();
        }
        result->addFacet(facet);

        CGAL::FT weight = 1;
        if (weight_pmap_res) {
            weight = get(*weight_pmap_res, fi);
        } else {
            std::cerr << "Error: face was not associated a weight?" << std::endl;
            return {};
        }

        CGAL_assertion(weight != 0);
        data::_3d::skel::SkelFacetDataSPtr data = data::_3d::skel::SkelFacetData::create(facet);
        data->setSpeed(weight);
    }

    result->setDescription(description);
    std::list<EdgeSPtr>::iterator it_e = result->edges().begin();
    while (it_e != result->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (!(edge->getFacetL() && edge->getFacetR())) {
            DEBUG_VAL("Warning: Polyhedron has no closed boundary.");
            DEBUG_VAR(edge->toString());
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

PolyhedronSPtr PLYFile::load(const CGAL::Surface_mesh<Point3>& sm,
                             const std::string& description)
{
  std::map<edge_descriptor, EdgeWPtr> unused_e2e;
  return load(sm, description, unused_e2e);
}

// This particular reader can read weights if they are stored as a property in the PLY file
PolyhedronSPtr PLYFile::load(const std::string& filename) {
    PolyhedronSPtr result = PolyhedronSPtr();

    std::ifstream ifs(filename.c_str());
    if (ifs.is_open()) {
        // @todo avoid building the intermediate mesh (but still use CGAL readers)
        CGAL::Surface_mesh<Point3> sm;
        bool success = CGAL::IO::read_PLY(ifs, sm);
        if (!success) {
          return {};
        }

        result = load(sm, "filename='"+filename+"'; ");
    }

    return result;
}

} }
