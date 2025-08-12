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
 * @file   algo/3d/AbstractVertexSplitter.cpp
 * @author Gernot Walzl
 * @date   2012-10-17
 */

#include "algo/3d/AbstractVertexSplitter.h"

#include "algo/3d/KernelWrapper.h"
#include "algo/3d/PolyhedronTransformation.h"
#include "algo/3d/SelfIntersection.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/CircularArc.h"
#include "data/3d/skel/SkelVertexData.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "data/3d/skel/SkelFacetData.h"

#include <limits>

namespace algo { namespace _3d {

AbstractVertexSplitter::AbstractVertexSplitter() {
    type_ = 0;
}

AbstractVertexSplitter::~AbstractVertexSplitter() {
    // intentionally does nothing
}


int AbstractVertexSplitter::getType() const {
    return type_;
}


PolyhedronSPtr AbstractVertexSplitter::splitConvexVertex(VertexSPtr vertex) {
    assert(vertex->isConvex());
    PolyhedronSPtr result = vertex->getPolyhedron();
    while (vertex->degree() > 3) {
        CGAL::FT sq_speed_max = 0.0;
        FacetSPtr facet_1;
        FacetSPtr facet_2;
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (FacetSPtr facet = facet_wptr.lock()) {
                FacetSPtr facet_prev = facet->prev(vertex);
                FacetSPtr facet_next = facet->next(vertex);
                CGAL::FT speed = 1.0;
                if (facet->hasData()) {
                    speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
                }
                CGAL::FT speed_prev = 1.0;
                if (facet_prev->hasData()) {
                    speed_prev = std::dynamic_pointer_cast<SkelFacetData>(facet_prev->getData())->getSpeed();
                }
                CGAL::FT speed_next = 1.0;
                if (facet_next->hasData()) {
                    speed_next = std::dynamic_pointer_cast<SkelFacetData>(facet_next->getData())->getSpeed();
                }

                // @todo I don't think the quick split of this function works for weighted cases
                CGAL_assertion(speed == 1.0 && speed_prev == 1.0 && speed_next == 1.0);

                Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->plane(), - speed);
                Plane3SPtr offset_plane_prev = KernelWrapper::offsetPlane(facet_prev->plane(), -speed_prev);
                Plane3SPtr offset_plane_next = KernelWrapper::offsetPlane(facet_next->plane(), -speed_next);

                Point3SPtr offset_point = KernelWrapper::intersection(offset_plane_prev, offset_plane, offset_plane_next);

                CGAL::FT sq_speed = KernelWrapper::squared_distance(vertex->getPoint(), offset_point);
                if (sq_speed > sq_speed_max) {
                    facet_1 = facet_next;
                    facet_2 = facet_prev;
                    sq_speed_max = sq_speed;
                }
            }
        }
        VertexSPtr vertex_splitted = vertex->split(facet_1, facet_2);
        if (vertex->hasData()) {
            SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex->getData());
            SkelVertexDataSPtr data_splitted = SkelVertexData::create(
                    vertex_splitted);
            data_splitted->setNode(data->getNode());
        }
    }
    return result;
}

PolyhedronSPtr AbstractVertexSplitter::splitReflexVertex(VertexSPtr vertex) {
    assert(vertex->isReflex());
    PolyhedronSPtr result = vertex->getPolyhedron();
    while (vertex->degree() > 3) {
        CGAL::FT speed_min = std::numeric_limits<double>::max(); // do not put FT
        FacetSPtr facet_1;
        FacetSPtr facet_2;
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (FacetSPtr facet = facet_wptr.lock()) {
                FacetSPtr facet_prev = facet->prev(vertex);
                FacetSPtr facet_next = facet->next(vertex);
                CGAL::FT speed = 1.0;
                if (facet->hasData()) {
                    speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
                }
                CGAL::FT speed_prev = 1.0;
                if (facet_prev->hasData()) {
                    speed_prev = std::dynamic_pointer_cast<SkelFacetData>(facet_prev->getData())->getSpeed();
                }
                CGAL::FT speed_next = 1.0;
                if (facet_next->hasData()) {
                    speed_next = std::dynamic_pointer_cast<SkelFacetData>(facet_next->getData())->getSpeed();
                }

                // @todo I don't think the quick split of this function works for weighted cases
                CGAL_assertion(speed == 1.0 && speed_prev == 1.0 && speed_next == 1.0);

                Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->plane(), -speed);
                Plane3SPtr offset_plane_prev = KernelWrapper::offsetPlane(facet_prev->plane(), -speed_prev);
                Plane3SPtr offset_plane_next = KernelWrapper::offsetPlane(facet_next->plane(), -speed_next);

                Point3SPtr offset_point = KernelWrapper::intersection(
                        offset_plane_prev, offset_plane, offset_plane_next);

                CGAL::FT v_speed = KernelWrapper::distance(vertex->getPoint(), offset_point);
                if (v_speed < speed_min) {
                    facet_1 = facet_next;
                    facet_2 = facet_prev;
                    speed_min = v_speed;
                }
            }
        }
        VertexSPtr vertex_splitted = vertex->split(facet_1, facet_2);
        if (vertex->hasData()) {
            SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex->getData());
            SkelVertexDataSPtr data_splitted = SkelVertexData::create(
                    vertex_splitted);
            data_splitted->setNode(data->getNode());
        }
    }
    return result;
}

bool AbstractVertexSplitter::checkSplitted(PolyhedronSPtr polyhedron) {
    bool result = false;
    PolyhedronSPtr polyhedron_offset = PolyhedronTransformation::shiftFacets(polyhedron, -1.0);
    if (polyhedron_offset) {
        result = !SelfIntersection::hasSelfIntersectingSurface(polyhedron_offset);
    }
    return result;
}


std::string AbstractVertexSplitter::toString() const {
    std::string result;
    switch (getType()) {
        case ANGLE_VERTEX_SPLITTER:
            result = "AngleVertexSplitter";
            break;
        case COMBI_VERTEX_SPLITTER:
            result = "CombiVertexSplitter";
            break;
        case CONVEX_VERTEX_SPLITTER:
            result = "ConvexVertexSplitter";
            break;
        case VOLUME_VERTEX_SPLITTER:
            result = "VolumeVertexSplitter";
            break;
        case WEIGHT_VERTEX_SPLITTER:
            result = "WeightVertexSplitter";
            break;
        case SPHERE_VERTEX_SPLITTER:
            result = "SphereVertexSplitter";
            break;
        default:
            result = "AbstractVertexSplitter";
    }
    return result;
}

} }

