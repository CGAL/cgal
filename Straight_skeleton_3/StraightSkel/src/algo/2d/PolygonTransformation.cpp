/**
 * @file   algo/2d/PolygonTransformation.cpp
 * @author Gernot Walzl
 * @date   2012-09-18
 */

#include "algo/2d/PolygonTransformation.h"

#include "algo/2d/KernelWrapper.h"
#include "data/2d/Edge.h"
#include "data/2d/Vertex.h"
#include "data/2d/Polygon.h"
#include "util/StringFactory.h"
#include <cstdlib>
#include <list>

namespace algo { namespace _2d {

PolygonTransformation::PolygonTransformation() {
    // intentionally does nothing.
}

PolygonTransformation::~PolygonTransformation() {
    // intentionally does nothing.
}

bool PolygonTransformation::hasParallelLines(PolygonSPtr polygon) {
    bool result = false;
    std::list<EdgeSPtr>::iterator it_e1 = polygon->edges().begin();
    while (it_e1 != polygon->edges().end()) {
        EdgeSPtr edge1 = *it_e1++;
        std::list<EdgeSPtr>::iterator it_e2 = it_e1;
        while (it_e2 != polygon->edges().end()) {
            EdgeSPtr edge2 = *it_e2++;
            if (!KernelWrapper::intersection(
                    edge1->line(), edge2->line())) {
                result = true;
                break;
            }
        }
        if (result) {
            break;
        }
    }
    return result;
}

Vector2SPtr PolygonTransformation::randVec(double min, double max) {
    double rval[2];
    for (unsigned int i = 0; i < 2; i++) {
        rval[i] = (max-min)*((double)rand()/(double)RAND_MAX) + min;
    }
    Vector2SPtr result = KernelFactory::createVector2(rval[0], rval[1]);
    return result;
}

void PolygonTransformation::randMovePoints(PolygonSPtr polygon, double range) {
    // srand(time(NULL));
    srand(0);   // set seed to a const value to reproduce errors
    std::list<VertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point2SPtr p = vertex->getPoint();
        Vector2SPtr v_r = randVec(-range/2.0, range/2.0);
        Point2SPtr p_t = KernelFactory::createPoint2(*p + *v_r);
        vertex->setPoint(p_t);
    }

    polygon->appendDescription("rand_move_points_range=" +
            util::StringFactory::fromDouble(range) + "; ");
}


Point2SPtr PolygonTransformation::boundingBoxMin(PolygonSPtr polygon) {
    double p_min[2];
    for (unsigned int i = 0; i < 2; i++) {
        p_min[i] = std::numeric_limits<double>::max();
    }
    std::list<VertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point2SPtr p = vertex->getPoint();
        for (unsigned int i = 0; i < 2; i++) {
            if ((*p)[i] < p_min[i]) {
                p_min[i] = (*p)[i];
            }
        }
    }
    Point2SPtr result = KernelFactory::createPoint2(p_min[0],p_min[1]);
    return result;
}

Point2SPtr PolygonTransformation::boundingBoxMax(PolygonSPtr polygon) {
    double p_max[2];
    for (unsigned int i = 0; i < 2; i++) {
        p_max[i] = -std::numeric_limits<double>::max();
    }
    std::list<VertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point2SPtr p = vertex->getPoint();
        for (unsigned int i = 0; i < 2; i++) {
            if ((*p)[i] > p_max[i]) {
                p_max[i] = (*p)[i];
            }
        }
    }
    Point2SPtr result = KernelFactory::createPoint2(p_max[0],p_max[1]);
    return result;
}

} }
