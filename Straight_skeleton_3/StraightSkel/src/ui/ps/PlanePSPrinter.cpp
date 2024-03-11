/**
 * @file   ui/ps/PlanePSPrinter.cpp
 * @author Gernot Walzl
 * @date   2012-11-14
 */

#include "ui/ps/PlanePSPrinter.h"

#include "debug.h"
#include "typedefs_thread.h"
#include "data/2d/Polygon.h"
#include "data/2d/Edge.h"
#include "data/2d/Vertex.h"
#include "data/2d/skel/Node.h"
#include "data/2d/skel/Arc.h"
#include "data/2d/skel/StraightSkeleton.h"
#include "data/2d/mesh/Mesh.h"
#include "data/2d/mesh/MeshCell.h"
#include "data/2d/mesh/MeshVertex.h"
#include <cmath>
#include <limits>
#include <list>

namespace ui { namespace ps {

PlanePSPrinter::PlanePSPrinter() {
    for (unsigned int i = 0; i < 2; i++) {
        translate_[i] = 0.0f;
    }
    scale_ = 20.0f;
}

PlanePSPrinter::~PlanePSPrinter() {
    // intentionally does nothing
}

PlanePSPrinterSPtr PlanePSPrinter::create() {
    return PlanePSPrinterSPtr(new PlanePSPrinter());
}

void PlanePSPrinter::boundingBoxMin(PolygonSPtr polygon, vec2f& out) {
    for (unsigned int i = 0; i < 2; i++) {
        out[i] = std::numeric_limits<float>::max();
    }
    std::list<data::_2d::VertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        data::_2d::VertexSPtr vertex = *it_v++;
        if ((float)vertex->getX() < out[0]) {
            out[0] = (float)vertex->getX();
        }
        if ((float)vertex->getY() < out[1]) {
            out[1] = (float)vertex->getY();
        }
    }
}

void PlanePSPrinter::boundingBoxMax(PolygonSPtr polygon, vec2f& out) {
    for (unsigned int i = 0; i < 2; i++) {
        out[i] = -std::numeric_limits<float>::max();
    }
    std::list<data::_2d::VertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        data::_2d::VertexSPtr vertex = *it_v++;
        if ((float)vertex->getX() > out[0]) {
            out[0] = (float)vertex->getX();
        }
        if ((float)vertex->getY() > out[1]) {
            out[1] = (float)vertex->getY();
        }
    }
}

float PlanePSPrinter::getScale() const {
    return scale_;
}

void PlanePSPrinter::setScale(float scale) {
    scale_ = scale;
}

void PlanePSPrinter::initBoundingBox(PolygonSPtr polygon) {
    ReadLock l(polygon->mutex());
    vec2f polygon_min;
    vec2f polygon_max;
    boundingBoxMin(polygon, polygon_min);
    boundingBoxMax(polygon, polygon_max);
    for (unsigned int i = 0; i < 2; i++) {
        translate_[i] = -polygon_min[i];
    }
    vec2f paper_max;
    toPaper(polygon_max, paper_max);
    bounding_box_[0] = 0;
    bounding_box_[1] = 0;
    bounding_box_[2] = (int)ceilf(paper_max[0]);
    bounding_box_[3] = (int)ceilf(paper_max[1]);
}

void PlanePSPrinter::toPaper(const vec2f in, vec2f& out) {
    for (unsigned int i = 0; i < 2; i++) {
        out[i] = (in[i] + translate_[i]) * scale_;
    }
}

void PlanePSPrinter::printPolygon(PolygonSPtr polygon, std::ostream& out) {
    if (!polygon) {
        DEBUG_PRINT("Warning: polygon is null.");
        return;
    }
    ReadLock l(polygon->mutex());
    std::list<data::_2d::EdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        data::_2d::EdgeSPtr edge = *it_e++;
        data::_2d::VertexSPtr vertex_src = edge->getVertexSrc();
        data::_2d::VertexSPtr vertex_dst = edge->getVertexDst();
        vec2f src = {(float)vertex_src->getX(),
                     (float)vertex_src->getY()};
        vec2f dst = {(float)vertex_dst->getX(),
                     (float)vertex_dst->getY()};
        vec2f paper_src;
        vec2f paper_dst;
        toPaper(src, paper_src);
        toPaper(dst, paper_dst);
        printLine(paper_src, paper_dst, out);
    }
}

void PlanePSPrinter::printSkel(data::_2d::skel::StraightSkeletonSPtr skel, std::ostream& out) {
    ReadLock l(skel->mutex());
    std::list<data::_2d::skel::ArcSPtr>::iterator it_a = skel->arcs().begin();
    while (it_a != skel->arcs().end()) {
        data::_2d::skel::ArcSPtr arc = *it_a++;
        if (!arc->hasNodeDst()) {
            continue;
        }
        data::_2d::skel::NodeSPtr node_src = arc->getNodeSrc();
        data::_2d::skel::NodeSPtr node_dst = arc->getNodeDst();
        vec2f src = {(float)node_src->getX(),
                     (float)node_src->getY()};
        vec2f dst = {(float)node_dst->getX(),
                     (float)node_dst->getY()};
        vec2f paper_src;
        vec2f paper_dst;
        toPaper(src, paper_src);
        toPaper(dst, paper_dst);
        printLine(paper_src, paper_dst, out);
    }
}

void PlanePSPrinter::printMesh(data::_2d::mesh::MeshSPtr mesh, std::ostream& out) {
    ReadLock l(mesh->mutex());
    std::list<data::_2d::mesh::MeshCellSPtr>::iterator it_c = mesh->cells().begin();
    while (it_c != mesh->cells().end()) {
        data::_2d::mesh::MeshCellSPtr cell = *it_c++;
        unsigned int num_points = cell->vertices().size();
        vec2f points[num_points];
        unsigned int i = 0;
        std::list<data::_2d::mesh::MeshVertexSPtr>::iterator it_v = cell->vertices().begin();
        while (it_v != cell->vertices().end()) {
            data::_2d::mesh::MeshVertexSPtr vertex = *it_v++;
            vec2f point = {(float)vertex->getX(),
                           (float)vertex->getY()};
            toPaper(point, points[i]);
            i++;
        }
        printPath(num_points, points, true, false, out);
    }

}

} }
