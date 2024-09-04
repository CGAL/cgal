// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   db/3d/OBJFile.cpp
 * @author Gernot Walzl
 * @date   2012-05-15
 */

#include "db/3d/OBJFile.h"

#include "debug.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Triangle.h"
#include "data/3d/KernelFactory.h"
#include "util/StringFuncs.h"
#include "util/Configuration.h"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/mark_domain_in_triangulation.h>

#include <cmath>
#include <exception>
#include <fstream>
#include <sstream>
#include <vector>

namespace db { namespace _3d {

OBJFile::OBJFile() {
    // intentionally does nothing
}

OBJFile::~OBJFile() {
    // intentionally does nothing
}

PolyhedronSPtr OBJFile::load(const std::string& filename) {
    PolyhedronSPtr result = PolyhedronSPtr();
    std::ifstream ifs(filename.c_str());
    if (ifs.is_open()) {
        result = Polyhedron::create();
        unsigned int vertex_id_new = 0;
        unsigned int facet_id_new = 0;
        std::vector<VertexSPtr> vertices;
        std::vector<Vector3SPtr> normals;
        while (ifs.good()) {
            std::string line;
            getline(ifs, line);
            if (line.find("v ") == 0) {  // vertex
                std::vector<std::string> strs = util::StringFuncs::split(line, " \t", false);
                if (strs.size() >= 4) {
                    vertex_id_new++;
                    Point3SPtr point = KernelFactory::createPoint3(
                            atof(strs[1].c_str()), atof(strs[2].c_str()), atof(strs[3].c_str()));
                    VertexSPtr vertex = Vertex::create(point);
                    vertex->setID(vertex_id_new);
                    result->addVertex(vertex);
                }
            } else if (line.find("vn ") == 0) {  // vertex normal
                std::vector<std::string> strs = util::StringFuncs::split(line, " \t", false);
                if (strs.size() >= 4) {
                    Vector3SPtr normal = KernelFactory::createVector3(
                            atof(strs[1].c_str()), atof(strs[2].c_str()), atof(strs[3].c_str()));
                    normals.push_back(normal);
                }
            } else if (line.find("f ") == 0) {  // face
                if (vertices.size() != result->vertices().size()) {  // initialize vector
                    vertices = std::vector<VertexSPtr>(
                            result->vertices().begin(), result->vertices().end());
                }
                std::vector<std::string> strs = util::StringFuncs::split(line, " \t", false);
                if (strs.size() >= 4) {
                    facet_id_new++;
                    unsigned int num_vertices = strs.size() - 1;
                    std::cout << "new face of size " << num_vertices << std::endl;
                    CGAL_assertion(num_vertices > 2);

                    VertexSPtr poly_vertices[num_vertices];
                    Vector3SPtr normal_sum;
                    for (unsigned int i = 0; i < num_vertices; i++) {
                        poly_vertices[i] = VertexSPtr();
                    }
                    for (unsigned int i = 0; i < num_vertices; i++) {
                        std::vector<std::string> strsvertex =
                                util::StringFuncs::split(strs[i+1], "/", true);
                        if (strsvertex.size() > 0) {
                            unsigned int vertex_id = atoi(strsvertex[0].c_str());
                            if (0 < vertex_id && vertex_id <= vertices.size()) {
                                poly_vertices[i] = vertices[vertex_id - 1];
                            } else {
                                std::stringstream whatstream;
                                whatstream << "Vertex with id="
                                        << vertex_id
                                        << " does not exist.";
                                throw std::runtime_error(whatstream.str());
                            }
                            if (strsvertex.size() > 2) {
                                unsigned int normal_id = atoi(strsvertex[2].c_str());
                                if (0 < normal_id && normal_id <= normals.size()) {
                                    Vector3SPtr normal = normals[normal_id - 1];
                                    if (normal_sum) {
                                        normal_sum = KernelFactory::createVector3(*normal_sum + *normal);
                                    } else {
                                        normal_sum = normal;
                                    }
                                }
                            }
                        }
                    }

                    FacetSPtr facet = Facet::create(num_vertices, poly_vertices);
                    facet->setID(facet_id_new);

                    std::cout << "Final edges:" << std::endl;
                    for(auto e : facet->edges())
                      std::cout << e->toString() << std::endl;

                    if (num_vertices == 3) {
                        Triangle::create(facet, poly_vertices);
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
                }
            }
        }
        ifs.close();
        result->setDescription("filename='"+filename+"'; ");
        std::list<EdgeSPtr>::iterator it_e = result->edges().begin();
        while (it_e != result->edges().end()) {
            EdgeSPtr edge = *it_e++;
            if (!(edge->getFacetL() && edge->getFacetR())) {
                DEBUG_VAL("Warning: Polyhedron has no closed boundary.");
                DEBUG_VAR(edge->toString());
            }
        }

        std::cout << "Loaded OBJ: " << result->toString() << std::endl;

        util::ConfigurationSPtr config = util::Configuration::getInstance();
        if (config->contains("main", "merge_coplanar_faces") &&
            config->getBool("main", "merge_coplanar_faces")) {
            double epsilon = 0.;
            std::string section("db_3d_OBJFile");
            std::string key("epsilon_coplanarity");
            if (config->contains(section, key)) {
                epsilon = config->getDouble(section, key);
            }
            mergeCoplanarFacets(result, epsilon);
        }

        removeVerticesDegLt3(result);

        // @tmp with the degree checking, this is false if the input is not triangulated
        // assert(result->isConsistent());
    }

    std::cout << "Final load: " << result->toString() << std::endl;

    return result;
}

bool OBJFile::save(const std::string& filename, PolyhedronSPtr polyhedron) {
    bool result = false;
    std::ofstream ofs(filename.c_str());
    if (ofs.is_open()) {
        WriteLock l(polyhedron->mutex());
        unsigned int vertex_id = 0;
        std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
        while (it_v != polyhedron->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            vertex_id++;
            vertex->setID(vertex_id);
            ofs << "v " << CGAL::to_double(vertex->getX()) << " "
                        << CGAL::to_double(vertex->getY()) << " "
                        << CGAL::to_double(vertex->getZ()) << "\n";
        }
        unsigned int facet_id = 0;
        std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
        while (it_f != polyhedron->facets().end()) {
            FacetSPtr facet = *it_f++;
            facet->makeFirstConvex();
            facet_id++;
            facet->setID(facet_id);

            using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
            using PK = CGAL::Projection_traits_3<CGAL::K>;
            using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
            using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
            using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
            using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
            using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
            using PCDT_VH = PCDT::Vertex_handle;
            using PCDT_FH = PCDT::Face_handle;

            auto eit = facet->edges().begin();
            const VertexSPtr p0 = (*eit++)->src(facet);
            const VertexSPtr p1 = (*eit++)->src(facet);
            const VertexSPtr p2 = (*eit)->src(facet);

            const auto n = CGAL::cross_product(*(p1->getPoint()) - *(p0->getPoint()),
                                               *(p2->getPoint()) - *(p0->getPoint()));
            PK traits(n);
            PCDT pcdt(traits);

            std::map<VertexSPtr, PCDT_VH> face_vhs;

            auto ins = [&pcdt, &face_vhs] (const VertexSPtr v) -> PCDT_VH
            {
              auto res = face_vhs.emplace(v, PCDT_VH());
              if(res.second) // first time seeing this point
              {
                PCDT_VH vh = pcdt.insert(*(v->getPoint()));
                vh->info() = v;
                res.first->second = vh;
              }
              return res.first->second;
            };

            bool use_default = false;

            for(auto e : facet->edges_)
            {
              PCDT_VH vh0 = ins(e->src(facet));
              PCDT_VH vh1 = ins(e->dst(facet));

              try
              {
                pcdt.insert_constraint(vh0, vh1);
              }
              catch(const typename PCDT::Intersection_of_constraints_exception&)
              {
                std::cerr << "Warning: Failed to triangulate face" << std::endl;
                DEBUG_VAR(facet->toString());
                use_default = true;
              }
            }

            if(!use_default)
            {
              std::unordered_map<PCDT_FH, bool> in_domain_map;
              boost::associative_property_map< std::unordered_map<PCDT_FH, bool> > in_domain(in_domain_map);

              CGAL::mark_domain_in_triangulation(pcdt, in_domain);

              unsigned int num_f = 0;
              for(auto fh : pcdt.finite_face_handles())
              {
                if(!get(in_domain, fh))
                  continue;

                ofs << "f " << fh->vertex(0)->info()->getID() << " "
                            << fh->vertex(1)->info()->getID() << " "
                            << fh->vertex(2)->info()->getID() << "\n";

                ++num_f;
              }
              // std::cout << "face split into " << num_f << " faces" << std::endl;
            }
            else
            {
              ofs << "f ";
              unsigned int num_edges = 0;
              EdgeSPtr first = EdgeSPtr();
              EdgeSPtr edge = facet->edges().front();
              while (first != edge) {
                  if (!first) {
                      first = edge;
                  } else {
                      ofs << " ";
                  }
                  VertexSPtr vertex = edge->src(facet);
                  ofs << vertex->getID();
                  edge = edge->next(facet);
                  num_edges++;
              }
              if (num_edges != facet->edges().size()) {
                  DEBUG_VAL("Warning: Facet does not consist of connected edges only (" << num_edges << " VS " << facet->edges().size() << ")");
                  DEBUG_VAL("Warning: It is impossible for an obj file to store holes inside a facet.");
                  DEBUG_VAR(facet->toString());
              }
              ofs << "\n";
            }
        }
        result = (polyhedron->vertices().size() == vertex_id &&
                  polyhedron->facets().size() == facet_id);
        ofs.close();
    }
    std::cout << "-- Write OBJ end --" << std::endl;
    return result;
}

} }
