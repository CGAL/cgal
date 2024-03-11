/**
 * @file   algo/3d/PolyhedronTransformation.cpp
 * @author Gernot Walzl
 * @date   2012-09-01
 */

#include "algo/3d/PolyhedronTransformation.h"

#include "algo/3d/KernelWrapper.h"
#include "data/3d/Vertex.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "util/StringFactory.h"
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <list>
#include <limits>

namespace algo { namespace _3d {

PolyhedronTransformation::PolyhedronTransformation() {
    // intentionally does nothing.
}

PolyhedronTransformation::~PolyhedronTransformation() {
    // intentionally does nothing.
}

void PolyhedronTransformation::translate(PolyhedronSPtr polyhedron, Vector3SPtr v_t) {
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        Point3SPtr p_t = KernelFactory::createPoint3(*p + *v_t);
        vertex->setPoint(p_t);
    }
    polyhedron->initPlanes();

    polyhedron->appendDescription("translate=<" +
            util::StringFactory::fromDouble((*v_t)[0]) + ", " +
            util::StringFactory::fromDouble((*v_t)[1]) + ", " +
            util::StringFactory::fromDouble((*v_t)[2]) + ">; ");
}

void PolyhedronTransformation::scale(PolyhedronSPtr polyhedron, Vector3SPtr v_s) {
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        Point3SPtr p_s = KernelFactory::createPoint3(
                (*p)[0] * (*v_s)[0], (*p)[1] * (*v_s)[1], (*p)[2] * (*v_s)[2]);
        vertex->setPoint(p_s);
    }
    polyhedron->initPlanes();

    polyhedron->appendDescription("scale=<" +
            util::StringFactory::fromDouble((*v_s)[0]) + ", " +
            util::StringFactory::fromDouble((*v_s)[1]) + ", " +
            util::StringFactory::fromDouble((*v_s)[2]) + ">; ");
}

bool PolyhedronTransformation::hasParallelPlanes(PolyhedronSPtr polyhedron) {
    bool result = false;
    std::list<FacetSPtr>::iterator it_f1 = polyhedron->facets().begin();
    while (it_f1 != polyhedron->facets().end()) {
        FacetSPtr facet1 = *it_f1++;
        std::list<FacetSPtr>::iterator it_f2 = it_f1;
        while (it_f2 != polyhedron->facets().end()) {
            FacetSPtr facet2 = *it_f2++;
            if (!KernelWrapper::intersection(
                    facet1->plane(), facet2->plane())) {
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

bool PolyhedronTransformation::doAll3PlanesIntersect(PolyhedronSPtr polyhedron) {
    bool result = true;
    std::list<FacetSPtr>::iterator it_f1 = polyhedron->facets().begin();
    while (it_f1 != polyhedron->facets().end()) {
        FacetSPtr facet1 = *it_f1++;
        std::list<FacetSPtr>::iterator it_f2 = it_f1;
        while (it_f2 != polyhedron->facets().end()) {
            FacetSPtr facet2 = *it_f2++;
            std::list<FacetSPtr>::iterator it_f3 = it_f2;
            while (it_f3 != polyhedron->facets().end()) {
                FacetSPtr facet3 = *it_f3++;
                if (!KernelWrapper::intersection(
                        facet1->plane(), facet2->plane(), facet3->plane())) {
                    result = false;
                    break;
                }
            }
            if (!result) {
                break;
            }
        }
        if (!result) {
            break;
        }
    }
    return result;
}

Vector3SPtr PolyhedronTransformation::randVec(double min, double max) {
    double rval[3];
    for (unsigned int i = 0; i < 3; i++) {
        rval[i] = (max-min)*((double)rand()/(double)RAND_MAX) + min;
    }
    Vector3SPtr result = KernelFactory::createVector3(rval[0],rval[1],rval[2]);
    return result;
}

void PolyhedronTransformation::randMovePoints(PolyhedronSPtr polyhedron, double range) {
    // srand(time(NULL));
    srand(0);   // set seed to a const value to reproduce errors
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        Vector3SPtr v_r = randVec(-range/2.0, range/2.0);
        Point3SPtr p_t = KernelFactory::createPoint3(*p + *v_r);
        vertex->setPoint(p_t);
    }
    polyhedron->initPlanes();

    polyhedron->appendDescription("rand_move_points_range=" +
            util::StringFactory::fromDouble(range) + "; ");
}

Point3SPtr PolyhedronTransformation::boundingBoxMin(PolyhedronSPtr polyhedron) {
    double p_min[3];
    for (unsigned int i = 0; i < 3; i++) {
        p_min[i] = std::numeric_limits<double>::max();
    }
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        for (unsigned int i = 0; i < 3; i++) {
            if ((*p)[i] < p_min[i]) {
                p_min[i] = (*p)[i];
            }
        }
    }
    Point3SPtr result = KernelFactory::createPoint3(p_min[0],p_min[1],p_min[2]);
    return result;
}

Point3SPtr PolyhedronTransformation::boundingBoxMax(PolyhedronSPtr polyhedron) {
    double p_max[3];
    for (unsigned int i = 0; i < 3; i++) {
        p_max[i] = -std::numeric_limits<double>::max();
    }
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        for (unsigned int i = 0; i < 3; i++) {
            if ((*p)[i] > p_max[i]) {
                p_max[i] = (*p)[i];
            }
        }
    }
    Point3SPtr result = KernelFactory::createPoint3(p_max[0],p_max[1],p_max[2]);
    return result;
}

void PolyhedronTransformation::translateNscale(PolyhedronSPtr polyhedron,
        Point3SPtr p_box_min, Point3SPtr p_box_max) {
    Vector3SPtr v_box_min = KernelFactory::createVector3(p_box_min);
    Vector3SPtr v_box_max = KernelFactory::createVector3(p_box_max);
    Vector3SPtr v_size = KernelFactory::createVector3(*v_box_max - *v_box_min);
    Vector3SPtr v_center = KernelFactory::createVector3(
            (*v_box_min + *v_box_max) / 2.0);

    Point3SPtr p_box_min_curr = boundingBoxMin(polyhedron);
    Point3SPtr p_box_max_curr = boundingBoxMax(polyhedron);
    Vector3SPtr v_box_min_curr = KernelFactory::createVector3(p_box_min_curr);
    Vector3SPtr v_box_max_curr = KernelFactory::createVector3(p_box_max_curr);
    Vector3SPtr v_size_curr = KernelFactory::createVector3(
            *v_box_max_curr - *v_box_min_curr);
    Vector3SPtr v_center_curr = KernelFactory::createVector3(
            (*v_box_min_curr + *v_box_max_curr) / 2.0);

    double scale_factor = std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < 3; i++) {
        double s = (*v_size)[i]/(*v_size_curr)[i];
        if (scale_factor > s) {
            scale_factor = s;
        }
    }
    scale_factor = floor(scale_factor*1000.0)/1000.0;
    DEBUG_VAR(scale_factor);
    Vector3SPtr v_s = KernelFactory::createVector3(
            scale_factor, scale_factor, scale_factor);

    Vector3SPtr v_t = KernelFactory::createVector3((*v_center_curr) * -1.0);
    if (v_t->squared_length() > 0.0) {
        translate(polyhedron, v_t);
    }
    if (scale_factor != 1.0) {
        scale(polyhedron, v_s);
    }
    if (v_center->squared_length() > 0.0) {
        translate(polyhedron, v_center);
    }
}

bool PolyhedronTransformation::isInsideBox(PolyhedronSPtr polyhedron,
        Point3SPtr p_box_min, Point3SPtr p_box_max) {
    bool result = true;
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        for (unsigned int i = 0; i < 3; i++) {
            if (!((*p_box_min)[i] <= (*p)[i] &&
                    (*p)[i] <= (*p_box_max)[i])) {
                result = false;
                break;
            }
        }
        if (!result) {
            break;
        }
    }
    return result;
}

} }
