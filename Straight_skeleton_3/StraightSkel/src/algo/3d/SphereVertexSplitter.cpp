/**
 * @file   algo/3d/SphereVertexSplitter.cpp
 * @author Gernot Walzl
 * @date   2012-12-20
 */

#include "algo/3d/SphereVertexSplitter.h"

#include "algo/Controller.h"
#include "algo/3d/KernelWrapper.h"
#include "algo/3d/AbstractSimpleSphericalSkel.h"
#include "algo/3d/ProjSimpleSphericalSkel.h"
#include "algo/3d/RotSimpleSphericalSkel.h"
#include "algo/3d/TransSimpleSphericalSkel.h"
#include "algo/3d/SpeedSimpleSphericalSkel.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/SphericalPolygon.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/CircularArc.h"
#include "data/3d/skel/SkelVertexData.h"
#include "data/3d/skel/SphericalSkelVertexData.h"
#include "data/3d/skel/SphericalSkelEdgeData.h"
#include "util/Configuration.h"
#include <list>
#include <map>

namespace algo { namespace _3d {

SphereVertexSplitter::SphereVertexSplitter() {
    type_ = AbstractVertexSplitter::SPHERE_VERTEX_SPLITTER;
    controller_ = ControllerSPtr();
}

SphereVertexSplitter::SphereVertexSplitter(ControllerSPtr controller) {
    type_ = AbstractVertexSplitter::SPHERE_VERTEX_SPLITTER;
    controller_ = controller;
}

SphereVertexSplitter::~SphereVertexSplitter() {
    controller_.reset();
}

SphereVertexSplitterSPtr SphereVertexSplitter::create() {
    SphereVertexSplitterSPtr result =
            SphereVertexSplitterSPtr(new SphereVertexSplitter());
    return result;
}

SphereVertexSplitterSPtr SphereVertexSplitter::create(ControllerSPtr controller) {
    SphereVertexSplitterSPtr result =
            SphereVertexSplitterSPtr(new SphereVertexSplitter(controller));
    return result;
}


FacetSPtr SphereVertexSplitter::findFacet(EdgeSPtr edge_1, EdgeSPtr edge_2) {
    FacetSPtr result;
    FacetSPtr facet_1l = edge_1->getFacetL();
    FacetSPtr facet_1r = edge_1->getFacetR();
    FacetSPtr facet_2l = edge_2->getFacetL();
    FacetSPtr facet_2r = edge_2->getFacetR();
    if (edge_1->getVertexSrc() == edge_2->getVertexSrc() ||
            edge_1->getVertexSrc() == edge_2->getVertexDst()) {
        FacetSPtr facet_tmp = facet_1r;
        facet_1r = facet_1l;
        facet_1l = facet_tmp;
    }
    if (edge_2->getVertexDst() == edge_1->getVertexSrc() ||
            edge_2->getVertexDst() == edge_1->getVertexDst()) {
        FacetSPtr facet_tmp = facet_2r;
        facet_2r = facet_2l;
        facet_2l = facet_tmp;
    }
    if (facet_1r == facet_2l || facet_1r == facet_2r) {
        result = facet_1r;
    } else if (facet_1l == facet_2l || facet_1l == facet_2r) {
        result = facet_1l;
    }
    DEBUG_SPTR(result);
    return result;
}


SphericalPolygonSPtr SphereVertexSplitter::intersectPolyhedron(VertexSPtr vertex) {
    Point3SPtr p_center = vertex->getPoint();
    double radius = 1.0;
    Sphere3SPtr sphere = KernelFactory::createSphere3(p_center, radius);
    SphericalPolygonSPtr result = SphericalPolygon::create(sphere);

    std::list<EdgeSPtr> edges;
    std::list<EdgeWPtr>::iterator it_ew = vertex->edges().begin();
    while (it_ew != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_ew++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
            edges.push_back(edge);
        }
    }

    while (edges.size() > 0) {
        CircularVertexSPtr cvertex_begin;
        CircularVertexSPtr cvertex_prev;
        CircularVertexSPtr cvertex_current;
        EdgeSPtr edge_begin;
        EdgeSPtr edge_prev;
        EdgeSPtr edge = edges.front();
        while (edge != edge_begin) {
            if (!edge_begin) {
                edge_begin = edge;
            }
            Point3SPtr p_dst;
            if (edge->getVertexSrc() == vertex) {
                p_dst = edge->getVertexDst()->getPoint();
            } else {
                p_dst = edge->getVertexSrc()->getPoint();
            }
            Vector3SPtr dir = KernelWrapper::normalize(
                    KernelFactory::createVector3(*p_dst - *p_center));
            Point3SPtr p_sphere = KernelFactory::createPoint3(
                    *p_center + ((*dir) * radius));
            cvertex_current = CircularVertex::create(p_sphere);
            SphericalSkelVertexDataSPtr cvertex_data =
                    SphericalSkelVertexData::create(cvertex_current);
            cvertex_data->setEdgeOrigin(edge);
            result->addVertex(cvertex_current);
            if (!cvertex_begin) {
                cvertex_begin = cvertex_current;
            }
            if (cvertex_prev) {
                CircularEdgeSPtr cedge = CircularEdge::create(cvertex_prev, cvertex_current);
                FacetSPtr facet = findFacet(edge_prev, edge);
                cedge->setSupportingPlane(facet->plane());
                SphericalSkelEdgeDataSPtr data = SphericalSkelEdgeData::create(cedge);
                data->setFacetOrigin(facet);
                result->addEdge(cedge);
            }
            cvertex_prev = cvertex_current;
            edge_prev = edge;
            edge = edge->next(vertex);
            edges.remove(edge_prev);
        }
        if (cvertex_current && cvertex_begin) {
            CircularEdgeSPtr cedge = CircularEdge::create(cvertex_current, cvertex_begin);
            FacetSPtr facet = findFacet(edge_prev, edge);
            cedge->setSupportingPlane(facet->plane());
            SphericalSkelEdgeDataSPtr data = SphericalSkelEdgeData::create(cedge);
            data->setFacetOrigin(facet);
            result->addEdge(cedge);
        }
    }
    return result;
}


PolyhedronSPtr SphereVertexSplitter::splitVertex(VertexSPtr vertex) {
    PolyhedronSPtr polyhedron = vertex->getPolyhedron();
    if (vertex->degree() <= 3) {
        return polyhedron;
    }
    if (vertex->hasData()) {
        vertex->getData()->setHighlight(true);
    }
    if (controller_) {
        controller_->wait();
    }
    SphericalPolygonSPtr polygon = intersectPolyhedron(vertex);
    DEBUG_VAR(polygon->toString());

    AbstractSimpleSphericalSkelSPtr algo_sphericalskel;
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        std::string s_algo_sphericalskel = config->getString(
                "algo_3d_SphereVertexSplitter", "algo_sphericalskel");
        if (s_algo_sphericalskel.compare("ProjSimpleSphericalSkel") == 0) {
            algo_sphericalskel = ProjSimpleSphericalSkel::create(polygon, controller_);
        } else if (s_algo_sphericalskel.compare("RotSimpleSphericalSkel") == 0) {
            algo_sphericalskel = RotSimpleSphericalSkel::create(polygon, controller_);
        } else if (s_algo_sphericalskel.compare("TransSimpleSphericalSkel") == 0) {
            algo_sphericalskel = TransSimpleSphericalSkel::create(polygon, controller_);
        } else if (s_algo_sphericalskel.compare("SpeedSimpleSphericalSkel") == 0) {
            algo_sphericalskel = SpeedSimpleSphericalSkel::create(polygon, controller_);
        } else {
            DEBUG_VAL("Warning: " << s_algo_sphericalskel << " not found.");
            DEBUG_VAL("Using ProjSimpleSphericalSkel.");
            algo_sphericalskel = ProjSimpleSphericalSkel::create(polygon, controller_);
        }
    } else {
        algo_sphericalskel = ProjSimpleSphericalSkel::create(polygon, controller_);
    }
    SphericalSkeletonSPtr sphericalskel = algo_sphericalskel->getResult();
    if (controller_) {
        controller_->setDispSphericalPolygon(polygon);
        controller_->setDispSphericalSkel(sphericalskel);
    }
    algo_sphericalskel->run();
    if (controller_) {
        controller_->haltSkip();
        controller_->wait();
    }

    WriteLock l(polyhedron->mutex());

    NodeSPtr node_vertex;
    if (vertex->hasData()) {
        skel::SkelVertexDataSPtr data =
                std::dynamic_pointer_cast<skel::SkelVertexData>(vertex->getData());
        node_vertex = data->getNode();
    }
    std::map<CircularNodeSPtr, VertexSPtr> vertices;
    std::list<CircularNodeSPtr>::iterator it_n = sphericalskel->nodes().begin();
    while (it_n != sphericalskel->nodes().end()) {
        CircularNodeSPtr node = *it_n++;
        if (node->degree() > 1) {
            VertexSPtr vertex_split = Vertex::create(vertex->getPoint());
            vertices[node] = vertex_split;
            polyhedron->addVertex(vertex_split);
            skel::SkelVertexDataSPtr data_split =
                    skel::SkelVertexData::create(vertex_split);
            data_split->setNode(node_vertex);
        }
    }

    std::list<CircularArcSPtr>::iterator it_a = sphericalskel->arcs().begin();
    while (it_a != sphericalskel->arcs().end()) {
        CircularArcSPtr arc = *it_a++;
        CircularNodeSPtr node_src = arc->getNodeSrc();
        CircularNodeSPtr node_dst = arc->getNodeDst();
        CircularEdgeSPtr edge_l = arc->getEdgeLeft();
        CircularEdgeSPtr edge_r = arc->getEdgeRight();
        SphericalSkelEdgeDataSPtr data_l =
                std::dynamic_pointer_cast<SphericalSkelEdgeData>(edge_l->getData());
        SphericalSkelEdgeDataSPtr data_r =
                std::dynamic_pointer_cast<SphericalSkelEdgeData>(edge_r->getData());
        FacetSPtr facet_l = data_l->getFacetOrigin();
        FacetSPtr facet_r = data_r->getFacetOrigin();
        if (node_src->degree() == 1) {
            if (node_dst->degree() == 1) {
                // 2 edges with the same adjacent facets have to be merged
                EdgeSPtr edge_1;
                EdgeSPtr edge_2;
                std::list<EdgeSPtr>::iterator it_e = facet_l->edges().begin();
                while (it_e != facet_l->edges().end()) {
                    EdgeSPtr edge = *it_e++;
                    if (edge->getFacetL() == facet_r ||
                            edge->getFacetR() == facet_r) {
                        if (!edge_1) {
                            edge_1 = edge;
                        } else {
                            edge_2 = edge;
                            break;
                        }
                    }
                }

                if (edge_1->getVertexSrc() == vertex) {
                    if (edge_2->getVertexSrc() == vertex) {
                        edge_1->replaceVertexSrc(edge_2->getVertexDst());
                    } else if (edge_2->getVertexDst() == vertex) {
                        edge_1->replaceVertexSrc(edge_2->getVertexSrc());
                    }
                } else if (edge_1->getVertexDst() == vertex) {
                    if (edge_2->getVertexSrc() == vertex) {
                        edge_1->replaceVertexDst(edge_2->getVertexDst());
                    } else if (edge_2->getVertexDst() == vertex) {
                        edge_1->replaceVertexDst(edge_2->getVertexSrc());
                    }
                }
                facet_l->removeEdge(edge_2);
                facet_r->removeEdge(edge_2);
                polyhedron->removeEdge(edge_2);
            } else {
                EdgeSPtr edge;
                std::list<CircularVertexSPtr>::iterator it_cv = polygon->vertices().begin();
                while (it_cv != polygon->vertices().end()) {
                    CircularVertexSPtr cvertex = *it_cv++;
                    SphericalSkelVertexDataSPtr cvertex_data =
                            std::dynamic_pointer_cast<SphericalSkelVertexData>(
                            cvertex->getData());
                    if (cvertex_data->getNode() == node_src) {
                        edge = cvertex_data->getEdgeOrigin();
                        break;
                    }
                }
                VertexSPtr vertex_dst = vertices[node_dst];
                if (edge->getVertexSrc() == vertex) {
                    edge->replaceVertexSrc(vertex_dst);
                }
                if (edge->getVertexDst() == vertex) {
                    edge->replaceVertexDst(vertex_dst);
                }
                if (!facet_l->containsVertex(vertex_dst)) {
                    facet_l->addVertex(vertex_dst);
                }
                if (!facet_r->containsVertex(vertex_dst)) {
                    facet_r->addVertex(vertex_dst);
                }
            }
        } else {
            EdgeSPtr edge = Edge::create(vertices[node_src], vertices[node_dst]);
            edge->setFacetL(facet_l);
            edge->setFacetR(facet_r);
            facet_l->addEdge(edge);
            facet_r->addEdge(edge);
            polyhedron->addEdge(edge);
        }
    }

    std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (!facet_wptr.expired()) {
            FacetSPtr facet(facet_wptr);
            facet->removeVertex(vertex);
        }
    }
    polyhedron->removeVertex(vertex);
    return polyhedron;
}

} }
