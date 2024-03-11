/**
 * @file   db/2d/FLMAFile.cpp
 * @author Gernot Walzl
 * @date   2013-12-16
 */

#include "db/2d/FLMAFile.h"

#include "data/2d/KernelFactory.h"
#include "data/2d/Polygon.h"
#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "util/StringFuncs.h"
#include "util/Configuration.h"
#include <fstream>
#include <vector>

namespace db { namespace _2d {

FLMAFile::FLMAFile() {
    // intentionally does nothing
}

FLMAFile::~FLMAFile() {
    // intentionally does nothing
}

PolygonSPtr FLMAFile::load(const std::string& filename) {
    PolygonSPtr result = PolygonSPtr();
    std::ifstream ifs(filename.c_str());
    if (ifs.is_open()) {
        unsigned int l = 0;
        result = Polygon::create();
        unsigned int num_vertices = 0;
        unsigned int num_edges = 0;
        unsigned int edge_id = 0;
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
                    //double z = atof(str_coords[i*3+2].c_str());
                    Point2SPtr point = KernelFactory::createPoint2(x,y);
                    VertexSPtr vertex = Vertex::create(point);
                    vertex->setID(i);
                    result->addVertex(vertex);
                }
                vertices = std::vector<VertexSPtr>(
                        result->vertices().begin(), result->vertices().end());
            }
            if (l == 3) {
                num_edges = atoi(line.c_str());
            }
            if (l > 3 && l%2 == 1 && edge_id < num_edges) {
                std::vector<std::string> str_vs = util::StringFuncs::split(line, " \t", false);
                if (str_vs.size() >= 2) {
                    unsigned int vertex_src_id = atoi(str_vs[0].c_str());
                    unsigned int vertex_dst_id = atoi(str_vs[1].c_str());
                    if (vertex_src_id < num_vertices &&
                            vertex_dst_id < num_vertices) {
                        VertexSPtr vertex_src = vertices[vertex_src_id];
                        VertexSPtr vertex_dst = vertices[vertex_dst_id];
                        if (vertex_src->getEdgeOut() || vertex_dst->getEdgeIn()) {
                            DEBUG_VAL("Error: " << filename << ":" << l);
                        }
                        EdgeSPtr edge = Edge::create(vertex_src, vertex_dst);
                        edge->setID(edge_id);
                        edge_id++;
                        result->addEdge(edge);
                    } else {
                        DEBUG_VAL("Error: " << filename << ":" << l);
                    }
                } else {
                    DEBUG_VAL("Error: " << filename << ":" << l);
                }
            }
        }
        ifs.close();
        result->setDescription("filename='"+filename+"'; ");
        double epsilon = 0.0001;
        util::ConfigurationSPtr config = util::Configuration::getInstance();
        std::string section("db_2d_FLMAFile");
        std::string key("epsilon_collinearity");
        if (config->contains(section, key)) {
            epsilon = config->getDouble(section, key);
        }
        mergeCollinearEdges(result, epsilon);
        assert(result->isConsistent());
    }
    return result;
}

bool FLMAFile::save(const std::string& filename, PolygonSPtr polygon) {
    bool result = false;
    // TODO
    return result;
}

} }
