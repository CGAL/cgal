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
            util::StringFactory::fromDouble(CGAL::to_double((*v_t)[0])) + ", " +
            util::StringFactory::fromDouble(CGAL::to_double((*v_t)[1])) + ", " +
            util::StringFactory::fromDouble(CGAL::to_double((*v_t)[2])) + ">; ");
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
            util::StringFactory::fromDouble(CGAL::to_double((*v_s)[0])) + ", " +
            util::StringFactory::fromDouble(CGAL::to_double((*v_s)[1])) + ", " +
            util::StringFactory::fromDouble(CGAL::to_double((*v_s)[2])) + ">; ");
}

// @todo this function cannot deal with degree 1 vertices
Point3SPtr PolyhedronTransformation::shiftPoint(VertexSPtr vertex,
                                                CGAL::FT offset)
{
    // std::cout << "shift point, " << vertex->facets().size() << " incident facets" << std::endl;

    Plane3SPtr planes[3];
    unsigned int i = 0;
    std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (i < 3 && it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (!facet_wptr.expired()) {
            FacetSPtr facet = FacetSPtr(facet_wptr);
            Plane3SPtr plane = facet->plane();
            // std::cout << "Facet #" << facet->getID() << std::endl;
            // std::cout << "  Plane: " << *plane << std::endl;

            // @fixme do a proper independence by checking
            // at i = 1, that cross(u, v) != null_vector
            // at i = 2, that det(u, v, w) != 0
            bool independent = true;
            for(unsigned int j=0; j<i; ++j) {
                if (CGAL::parallel(*plane, *(planes[j]))) {
                    independent = false;
                    break;
                }
            }

            if (!independent) {
                continue;
            }

            CGAL::FT speed = 1.0;
            if (facet->hasData()) {
                speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
            }

            planes[i] = KernelWrapper::offsetPlane(plane, offset*speed);
            // std::cout << "  Offset Plane[" << i << "] = " << *(planes[i]) << std::endl;

            ++i;
        }
    }

    if (i < 3) {
        std::cerr << "Warning: Couldn't find three independent planes" << std::endl;
        return { };
    }

    Point3SPtr point = KernelWrapper::intersection(planes[0], planes[1], planes[2]);

    if (!point) {
#if 1 // Temporarily allowing points to not exist for simultaneous events:
        std::cerr << "Warning: triplet of planes doesn't define a point!" << std::endl;
#else
        Point3SPtr result = PolyhedronSPtr();
        DEBUG_SPTR(result);
        return result;
#endif
    }

    return point;
}

// @todo plenty of needless recomputations
PolyhedronSPtr PolyhedronTransformation::shiftFacets(PolyhedronSPtr polyhedron,
                                                     CGAL::FT offset,
                                                     const bool recompute_positions,
                                                     const bool with_sanity_checks)
{
    std::cout << "~~~~ Shift by " << offset << std::endl;

    std::ofstream shift_out("results/last_shift.polylines.txt");
    shift_out.precision(17);

    PolyhedronSPtr result = Polyhedron::create();

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr point = vertex->getPoint();

        // those are treated in the next loop
        if (vertex->degree() == 1) {
            continue;
        }

        if (offset != 0 || recompute_positions) {
            point = shiftPoint(vertex, offset);
            if (!point) {
                std::cerr << "Warning: Failed to create shifted polyhedron" << std::endl;
                return { };
            }
        }

        VertexSPtr offset_vertex = Vertex::create(point);
        // SkelVertexData for each vertex should be created by init
        SkelVertexDataSPtr data;
        if (vertex->hasData()) {
            data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
            SkelVertexDataSPtr offset_data = SkelVertexData::create(offset_vertex);
            offset_data->setArc(data->getArc());
        } else {
            data = SkelVertexData::create(vertex);
        }
        data->setOffsetVertex(offset_vertex);
        result->addVertex(offset_vertex);

        // std::cout << *(vertex->getPoint()) << " to " << *point << std::endl;
        shift_out << "2 " << *(vertex->getPoint()) << " " << *point << std::endl;
    }

    // Now, deal with degree 1 vertices.
    //
    // Do NOT merge the two vertices loops together: this function assumes that degree 1 vertices
    // are adjacent to degree 3+ vertices that have been offset in the first loop, and offsets
    // degree 1 vertices using the offset of the adjacent degree 3+ vertex.
    it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr point = vertex->getPoint();

        EdgeSPtr edge;
        unsigned int i = 0;
        std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
        while (it_e != vertex->edges().end()) {
            EdgeWPtr edge_wptr = *it_e++;
            if (!edge_wptr.expired()) {
                edge = EdgeSPtr(edge_wptr);
                i++;
            }
        }
        if (i == 1) {
            VertexSPtr vertex_other = (edge->getVertexSrc() == vertex) ? edge->getVertexDst()
                                                                       : edge->getVertexSrc();

            // otherwise vertex_other has not been offset in the previous vertex loop
            CGAL_assertion(vertex_other->degree() >= 3);

            if (offset != 0 || recompute_positions) {
                SkelVertexDataSPtr data_other = std::dynamic_pointer_cast<SkelVertexData>(vertex_other->getData());
                VertexSPtr offset_vertex_other = data_other->getOffsetVertex();
                Vector3 direction = *(offset_vertex_other->getPoint()) - *(vertex_other->getPoint());
                point = KernelFactory::createPoint3(*(vertex->getPoint()) + direction);
            }

            VertexSPtr offset_vertex = Vertex::create(point);
            SkelVertexDataSPtr data;
            if (vertex->hasData()) {
                data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
            } else {
                data = SkelVertexData::create(vertex);
            }
            data->setOffsetVertex(offset_vertex);
            result->addVertex(offset_vertex);

            std::cout << *(vertex->getPoint()) << " to " << *point << std::endl;
            shift_out << "2 " << *(vertex->getPoint()) << " " << *point << std::endl;
        } else {
            CGAL_assertion_msg(i != 0, "no edge on degree 1 vertex?");
        }
    }

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        SkelVertexDataSPtr vertex_src_data = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexSrc()->getData());
        SkelVertexDataSPtr vertex_dst_data = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexDst()->getData());
        if (vertex_src_data && vertex_dst_data) {
            VertexSPtr offset_vertex_src = vertex_src_data->getOffsetVertex();
            VertexSPtr offset_vertex_dst = vertex_dst_data->getOffsetVertex();
            EdgeSPtr offset_edge = Edge::create(offset_vertex_src, offset_vertex_dst);
            SkelEdgeDataSPtr data;
            if (edge->hasData()) {
                data = std::dynamic_pointer_cast<SkelEdgeData>(edge->getData());
                SkelEdgeDataSPtr offset_data = SkelEdgeData::create(offset_edge);
                offset_data->setSheet(data->getSheet());
            } else {
                data = SkelEdgeData::create(edge);
            }
            data->setOffsetEdge(offset_edge);
            result->addEdge(offset_edge);
        }
    }

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        FacetSPtr offset_facet = Facet::create();
        SkelFacetDataSPtr data;
        CGAL::FT speed = 1.0;
        if (facet->hasData()) {
            data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
            speed = data->getSpeed();
            SkelFacetDataSPtr data_offset = SkelFacetData::create(offset_facet);
            data_offset->setFacetOrigin(data->getFacetOrigin());
            data_offset->setSpeed(speed);
        } else {
            data = SkelFacetData::create(facet);
            data->setSpeed(speed);
        }
        Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->plane(), offset*speed);
        offset_facet->setPlane(offset_plane);
        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            if (vertex->hasData()) {
                SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex->getData());
                VertexSPtr offset_vertex = vertex_data->getOffsetVertex();
                if (offset_vertex) {
                    offset_facet->addVertex(offset_vertex);

                    // @fixme only a warning because sometimes I still run perturbed, non-triangulated cases
                    CGAL_warning(offset_plane->has_on(*(offset_vertex->getPoint())));

                    if (with_sanity_checks) {
                        std::cout << KernelWrapper::squared_distance(facet->getPlane(),
                                                                     offset_vertex->getPoint())
                                  << " VS1 " << CGAL::square(offset*speed) << std::endl;
                        CGAL_assertion(CGAL::abs(KernelWrapper::squared_distance(
                                           facet->getPlane(), offset_vertex->getPoint()) - CGAL::square(offset*speed)) < 1e-5);
                        CGAL_assertion(KernelWrapper::squared_distance(
                                           facet->getPlane(), offset_vertex->getPoint()) == CGAL::square(offset*speed));
                    }
                }
            }
        }
        std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                edge->getData());
            EdgeSPtr offset_edge = edge_data->getOffsetEdge();
            if (facet == edge->getFacetL()) {
                offset_edge->setFacetL(offset_facet);
                offset_facet->addEdge(offset_edge);
            }
            if (facet == edge->getFacetR()) {
                offset_edge->setFacetR(offset_facet);
                offset_facet->addEdge(offset_edge);
            }
        }
        data->setOffsetFacet(offset_facet);
        result->addFacet(offset_facet);
    }

    assert(polyhedron->edges().size() == result->edges().size());
    assert(polyhedron->facets().size() == result->facets().size());

    return result;
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
    CGAL::FT p_min[3];
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
    CGAL::FT p_max[3];
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

    CGAL::FT scale_factor = std::numeric_limits<double>::max(); // do not put FT
    for (unsigned int i = 0; i < 3; i++) {
        CGAL::FT s = (*v_size)[i]/(*v_size_curr)[i];
        if (scale_factor > s) {
            scale_factor = s;
        }
    }
    scale_factor = floor(CGAL::to_double(scale_factor)*1000.0)/1000.0; // @fixme interval
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
