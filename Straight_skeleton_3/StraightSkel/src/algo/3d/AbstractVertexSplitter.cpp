/**
 * @file   algo/3d/AbstractVertexSplitter.cpp
 * @author Gernot Walzl
 * @date   2012-10-17
 */

#include "algo/3d/AbstractVertexSplitter.h"

#include "algo/3d/KernelWrapper.h"
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
        double speed_max = 0.0;
        FacetSPtr facet_1;
        FacetSPtr facet_2;
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (facet_wptr.expired()) {
                continue;
            }
            FacetSPtr facet(facet_wptr);
            FacetSPtr facet_prev = facet->prev(vertex);
            FacetSPtr facet_next = facet->next(vertex);
            Plane3SPtr offset_plane =
                    KernelWrapper::offsetPlane(facet->plane(), -1.0);
            Plane3SPtr offset_plane_prev =
                    KernelWrapper::offsetPlane(facet_prev->plane(), -1.0);
            Plane3SPtr offset_plane_next =
                    KernelWrapper::offsetPlane(facet_next->plane(), -1.0);
            Point3SPtr offset_point = KernelWrapper::intersection(
                    offset_plane_prev, offset_plane, offset_plane_next);
            double speed = KernelWrapper::distance(
                    vertex->getPoint(), offset_point);
            if (speed > speed_max) {
                facet_1 = facet_next;
                facet_2 = facet_prev;
                speed_max = speed;
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
        double speed_min = std::numeric_limits<double>::max();
        FacetSPtr facet_1;
        FacetSPtr facet_2;
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (facet_wptr.expired()) {
                continue;
            }
            FacetSPtr facet(facet_wptr);
            FacetSPtr facet_prev = facet->prev(vertex);
            FacetSPtr facet_next = facet->next(vertex);
            Plane3SPtr offset_plane =
                    KernelWrapper::offsetPlane(facet->plane(), -1.0);
            Plane3SPtr offset_plane_prev =
                    KernelWrapper::offsetPlane(facet_prev->plane(), -1.0);
            Plane3SPtr offset_plane_next =
                    KernelWrapper::offsetPlane(facet_next->plane(), -1.0);
            Point3SPtr offset_point = KernelWrapper::intersection(
                    offset_plane_prev, offset_plane, offset_plane_next);
            double speed = KernelWrapper::distance(
                    vertex->getPoint(), offset_point);
            if (speed < speed_min) {
                facet_1 = facet_next;
                facet_2 = facet_prev;
                speed_min = speed;
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


PolyhedronSPtr AbstractVertexSplitter::shiftFacets(PolyhedronSPtr polyhedron, double offset) {
    PolyhedronSPtr result = Polyhedron::create();

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Plane3SPtr planes[3];
        unsigned int i = 0;
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (i < 3 && it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet = FacetSPtr(facet_wptr);
                double speed = 1.0;
                if (facet->hasData()) {
                    speed = std::dynamic_pointer_cast<SkelFacetData>(
                            facet->getData())->getSpeed();
                }
                planes[i] = KernelWrapper::offsetPlane(facet->plane(), offset*speed);
                i++;
            }
        }
        if (i >= 3) {
            Point3SPtr point = KernelWrapper::intersection(planes[0], planes[1], planes[2]);
            if (!point) {
                result = PolyhedronSPtr();
                DEBUG_SPTR(result);
                return result;
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
        }
    }

    // the following code is used by VertexSplitter::checkSplitted(...)
    it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
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
            VertexSPtr vertex_other;
            if (edge->getVertexSrc() == vertex) {
                vertex_other = edge->getVertexDst();
            } else if (edge->getVertexDst() == vertex) {
                vertex_other = edge->getVertexSrc();
            }
            SkelVertexDataSPtr data_other = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_other->getData());
            VertexSPtr offset_vertex_other = data_other->getOffsetVertex();
            Vector3 direction =
                    *(offset_vertex_other->getPoint()) - *(vertex_other->getPoint());
            Point3SPtr point = KernelFactory::createPoint3(
                    *(vertex->getPoint()) + direction);
            VertexSPtr offset_vertex = Vertex::create(point);
            SkelVertexDataSPtr data;
            if (vertex->hasData()) {
                data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
            } else {
                data = SkelVertexData::create(vertex);
            }
            data->setOffsetVertex(offset_vertex);
            result->addVertex(offset_vertex);
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
        double speed = 1.0;
        if (facet->hasData()) {
            data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
            speed = data->getSpeed();
            SkelFacetDataSPtr data_offset = SkelFacetData::create(offset_facet);
            data_offset->setFacetOrigin(data->getFacetOrigin());
            data_offset->setSpeed(data->getSpeed());
        } else {
            data = SkelFacetData::create(facet);
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

bool AbstractVertexSplitter::checkSplitted(PolyhedronSPtr polyhedron) {
    bool result = false;
    PolyhedronSPtr polyhedron_offset = shiftFacets(polyhedron, -1.0);
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

