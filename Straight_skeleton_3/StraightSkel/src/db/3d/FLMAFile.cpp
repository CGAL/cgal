/**
 * @file   db/3d/FLMAFile.cpp
 * @author Gernot Walzl
 * @date   2013-12-16
 */

#include "db/3d/FLMAFile.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/Triangle.h"
#include "data/3d/KernelFactory.h"
#include "util/Configuration.h"
#include "util/StringFuncs.h"
#include <fstream>
#include <vector>

namespace db { namespace _3d {

FLMAFile::FLMAFile() {
    // intentionally does nothing
}

FLMAFile::~FLMAFile() {
    // intentionally does nothing
}

PolyhedronSPtr FLMAFile::load(const std::string& filename) {
    PolyhedronSPtr result = PolyhedronSPtr();
    std::ifstream ifs(filename.c_str());
    if (ifs.is_open()) {
        unsigned int l = 0;
        result = Polyhedron::create();
        unsigned int num_vertices = 0;
        unsigned int num_facets = 0;
        unsigned int poly_id = 0;
        std::vector<VertexSPtr> vertices;
        while (ifs.good()) {
            l++;
            std::string line;
            getline(ifs, line);
            if (l == 1) {
                num_vertices = atoi(line.c_str());
            }
            if (l == 2) {
                std::vector<std::string> str_coords =
                        util::StringFuncs::split(line, " \t", false);
                if (str_coords.size()%3 != 0) {
                    DEBUG_VAL("Error: " << filename << ":" << l);
                }
                for (unsigned int i = 0; i < num_vertices; i++) {
                    double x = atof(str_coords[i*3].c_str());
                    double y = atof(str_coords[i*3+1].c_str());
                    double z = atof(str_coords[i*3+2].c_str());
                    Point3SPtr point = KernelFactory::createPoint3(x,y,z);
                    VertexSPtr vertex = Vertex::create(point);
                    vertex->setID(i);
                    result->addVertex(vertex);
                }
                vertices = std::vector<VertexSPtr>(
                        result->vertices().begin(), result->vertices().end());
            }
            if (l == 3) {
                num_facets = atoi(line.c_str());
            }
            if (l > 3 && l%2 == 1 && poly_id < num_facets) {
                std::vector<std::string> str_vs = util::StringFuncs::split(line, " \t", false);
                unsigned int num_poly_vertices = str_vs.size();
                VertexSPtr poly_vertices[num_poly_vertices];
                for (unsigned int i = 0; i < num_poly_vertices; i++) {
                    poly_vertices[i] = VertexSPtr();
                }
                for (unsigned int i = 0; i < num_poly_vertices; i++) {
                    unsigned int vertex_id = atoi(str_vs[i].c_str());
                    if (vertex_id < num_vertices) {
                        poly_vertices[num_poly_vertices-i-1] = vertices[vertex_id];
                    } else {
                        DEBUG_VAL("Error: " << filename << ":" << l);
                    }
                }
                FacetSPtr facet = Facet::create(num_poly_vertices, poly_vertices);
                facet->setID(poly_id);
                poly_id++;
                if (num_poly_vertices == 3) {
                    Triangle::create(facet, poly_vertices);
                }
                if (num_poly_vertices >= 3) {
                    facet->initPlane();
                }
                result->addFacet(facet);
            }
        }
        ifs.close();
        result->setDescription("filename='"+filename+"'; ");
        double epsilon = 0.0001;
        util::ConfigurationSPtr config = util::Configuration::getInstance();
        std::string section("db_3d_FLMAFile");
        std::string key("epsilon_coplanarity");
        if (config->contains(section, key)) {
            epsilon = config->getDouble(section, key);
        }
        mergeCoplanarFacets(result, epsilon);
        removeVerticesDegLt3(result);
        assert(result->isConsistent());
    }
    return result;
}

bool FLMAFile::save(const std::string& filename, PolyhedronSPtr polyhedron) {
    bool result = false;
    // TODO
    return result;
}

} }
