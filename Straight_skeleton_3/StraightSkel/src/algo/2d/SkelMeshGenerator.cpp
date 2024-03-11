/**
 * @file   algo/2d/SkelMeshGenerator.cpp
 * @author Gernot Walzl
 * @date   2014-01-28
 */

#include "algo/2d/SkelMeshGenerator.h"

#include "debug.h"
#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "data/2d/Polygon.h"
#include "data/2d/skel/Node.h"
#include "data/2d/skel/Arc.h"
#include "data/2d/skel/StraightSkeleton.h"
#include "data/2d/skel/SkelVertexData.h"
#include "data/2d/skel/AbstractEvent.h"
#include "data/2d/mesh/Mesh.h"
#include "data/2d/mesh/MeshEdgeData.h"
#include "data/2d/mesh/MeshVertex.h"
#include "data/2d/mesh/MeshCell.h"
#include "data/2d/mesh/MeshRay.h"
#include "algo/Controller.h"
#include "algo/2d/KernelWrapper.h"
#include "algo/2d/MeshModifier.h"
#include <list>

namespace algo { namespace _2d {

SkelMeshGenerator::SkelMeshGenerator(StraightSkeletonSPtr skel) {
    skel_ = skel;
    mesh_result_ = Mesh::create();
}

SkelMeshGenerator::SkelMeshGenerator(StraightSkeletonSPtr skel, ControllerSPtr controller) {
    skel_ = skel;
    controller_ = controller;
    mesh_result_ = Mesh::create();
}

SkelMeshGenerator::~SkelMeshGenerator() {
    skel_.reset();
    controller_.reset();
    mesh_result_.reset();
}

SkelMeshGeneratorSPtr SkelMeshGenerator::create(StraightSkeletonSPtr skel) {
    return SkelMeshGeneratorSPtr(new SkelMeshGenerator(skel));
}

SkelMeshGeneratorSPtr SkelMeshGenerator::create(StraightSkeletonSPtr skel, ControllerSPtr controller) {
    return SkelMeshGeneratorSPtr(new SkelMeshGenerator(skel, controller));
}


void SkelMeshGenerator::initEdgeDatas() {
    PolygonSPtr polygon = skel_->getPolygon();
    WriteLock l(polygon->mutex());
    std::list<EdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        EdgeSPtr edge = *it_e++;
        MeshEdgeData::create(edge);
    }
    std::list<ArcSPtr>::iterator it_a = skel_->arcs().begin();
    while (it_a != skel_->arcs().end()) {
        ArcSPtr arc = *it_a++;
        MeshEdgeDataSPtr data_l = std::dynamic_pointer_cast<MeshEdgeData>(
                arc->getEdgeLeft()->getData());
        MeshEdgeDataSPtr data_r = std::dynamic_pointer_cast<MeshEdgeData>(
                arc->getEdgeRight()->getData());
        data_l->addArc(arc);
        data_r->addArc(arc);
    }
    it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        EdgeSPtr edge = *it_e++;
        sortArcs(edge);
    }
}

void SkelMeshGenerator::sortArcs(EdgeSPtr edge) {
    MeshEdgeDataSPtr data = std::dynamic_pointer_cast<MeshEdgeData>(edge->getData());
    std::list<ArcSPtr> arcs;
    std::list<ArcWPtr>::iterator it_a_w = data->arcs().begin();
    while (it_a_w != data->arcs().end()) {
        ArcWPtr arc_wptr = *it_a_w++;
        if (!arc_wptr.expired()) {
            ArcSPtr arc(arc_wptr);
            arcs.push_back(arc);
        }
    }

    std::list<ArcSPtr> arcs_sorted;
    SkelVertexDataSPtr data_vertex = std::dynamic_pointer_cast<SkelVertexData>(
            edge->getVertexDst()->getData());
    NodeSPtr node_first = data_vertex->getNode();
    NodeSPtr node_current = node_first;
    while (!arcs.empty()) {
        std::list<ArcSPtr>::iterator it_a = arcs.begin();
        while (it_a != arcs.end()) {
            std::list<ArcSPtr>::iterator it_current = it_a;
            ArcSPtr arc = *it_a++;
            if (arc->getNodeSrc() == node_current) {
                node_current = arc->getNodeDst();
                arcs_sorted.push_back(arc);
                arcs.erase(it_current);
                break;
            } else if (arc->getNodeDst() == node_current) {
                node_current = arc->getNodeSrc();
                arcs_sorted.push_back(arc);
                arcs.erase(it_current);
                break;
            }
        }
    }

    data->arcs().clear();
    std::list<ArcSPtr>::iterator it_a = arcs_sorted.begin();
    while (it_a != arcs_sorted.end()) {
        ArcSPtr arc = *it_a++;
        data->addArc(arc);
    }
}

void SkelMeshGenerator::initCells() {
    WriteLock l(mesh_result_->mutex());
    PolygonSPtr polygon = skel_->getPolygon();
    std::list<EdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        EdgeSPtr edge = *it_e++;
        MeshEdgeDataSPtr data = std::dynamic_pointer_cast<MeshEdgeData>(
                edge->getData());
        std::list<Point2SPtr> points;
        Point2SPtr p_last;
        std::list<ArcWPtr>::iterator it_a = data->arcs().begin();
        while (it_a != data->arcs().end()) {
            ArcWPtr arc_wptr = *it_a++;
            if (!arc_wptr.expired()) {
                ArcSPtr arc(arc_wptr);
                Point2SPtr p_src = arc->getNodeSrc()->getPoint();
                Point2SPtr p_dst = arc->getNodeDst()->getPoint();
                if (!p_last) {
                    points.push_back(p_src);
                    points.push_back(p_dst);
                    p_last = p_dst;
                } else {
                    if (p_src == p_last) {
                        points.push_back(p_dst);
                        p_last = p_dst;
                    } else if (p_dst == p_last) {
                        points.push_back(p_src);
                        p_last = p_src;
                    } else {
                        DEBUG_PRINT("Error: Arcs are not sorted.");
                    }
                }
            }
        }
        std::list<Point2SPtr>::iterator it_p = points.begin();
        unsigned int num_vertices = points.size();
        MeshVertexSPtr vertices[num_vertices];
        for (unsigned int i = 0; i < num_vertices; i++) {
            Point2SPtr point = *it_p++;
            MeshVertexSPtr vertex = mesh_result_->getVertex(point);
            if (!vertex) {
                vertex = MeshVertex::create(point);
                mesh_result_->addVertex(vertex);
            }
            vertices[i] = vertex;
        }
        MeshCellSPtr cell = MeshCell::create(num_vertices, vertices);
        mesh_result_->addCell(cell);
        data->addCell(cell);
    }
    assert(mesh_result_->isConsistent());
}

void SkelMeshGenerator::createRays(EdgeSPtr edge) {
    WriteLock l(mesh_result_->mutex());
    MeshEdgeDataSPtr data = std::dynamic_pointer_cast<MeshEdgeData>(edge->getData());

    std::list<NodeSPtr> nodes;
    NodeSPtr node_prev;
    std::list<ArcWPtr>::iterator it_a = data->arcs().begin();
    while (it_a != data->arcs().end()) {
        ArcWPtr arc_wptr = *it_a++;
        if (!arc_wptr.expired()) {
            ArcSPtr arc(arc_wptr);
            if (nodes.empty()) {
                nodes.push_back(arc->getNodeSrc());
                nodes.push_back(arc->getNodeDst());
                node_prev = arc->getNodeDst();
            } else {
                if (arc->getNodeSrc() == node_prev) {
                    nodes.push_back(arc->getNodeDst());
                    node_prev = arc->getNodeDst();
                } else if (arc->getNodeDst() == node_prev) {
                    nodes.push_back(arc->getNodeSrc());
                    node_prev = arc->getNodeSrc();
                }
            }
        }
    }

    Point2SPtr p_src = edge->getVertexSrc()->getPoint();
    Point2SPtr p_dst = edge->getVertexDst()->getPoint();
    Vector2SPtr v_dir = KernelFactory::createVector2(*p_dst - *p_src);
    std::list<NodeSPtr> nodes_sorted;
    while (!nodes.empty()) {
        NodeSPtr node_min = nodes.front();
        std::list<NodeSPtr>::iterator it_min = nodes.begin();
        std::list<NodeSPtr>::iterator it_n = nodes.begin();
        while (it_n != nodes.end()) {
            std::list<NodeSPtr>::iterator it_current = it_n;
            NodeSPtr node = *it_n++;
            if (KernelWrapper::compatePoints(v_dir,
                    node_min->getPoint(), node->getPoint()) > 0) {
                node_min = node;
                it_min = it_current;
            }
        }
        nodes_sorted.push_back(node_min);
        nodes.erase(it_min);
    }

    nodes_sorted.pop_front();
    nodes_sorted.pop_back();

    std::list<NodeSPtr>::iterator it_n = nodes_sorted.begin();
    while (it_n != nodes_sorted.end()) {
        NodeSPtr node = *it_n++;
        MeshVertexSPtr vertex = mesh_result_->getVertex(node->getPoint());
        MeshRaySPtr ray = MeshRay::create(edge, vertex);
        data->addRay(ray);
        mesh_result_->addRay(ray);
    }
}

void SkelMeshGenerator::findRayDsts() {
    WriteLock l(mesh_result_->mutex());
    std::list<MeshRaySPtr> rays;
    std::list<MeshRaySPtr>::iterator it_r = mesh_result_->rays().begin();
    while (it_r != mesh_result_->rays().end()) {
        MeshRaySPtr ray = *it_r++;
        if (!ray->getDst()) {
            rays.push_back(ray);
        }
    }

    while (!rays.empty()) {
        MeshRaySPtr ray = rays.front();
        rays.pop_front();
        EdgeSPtr edge = ray->getEdge();
        Line2SPtr l_edge = edge->line();
        Vector2SPtr dir_ray = KernelWrapper::perpendicular(
                KernelFactory::createVector2(
                        KernelWrapper::opposite(l_edge)));
        Point2SPtr p_ray_src = ray->getSrc()->getPoint();
        Line2SPtr l_ray = KernelFactory::createLine2(p_ray_src, dir_ray);
        Point2SPtr p_intersect = KernelWrapper::intersection(l_ray, l_edge);
        Point2SPtr p_src = edge->getVertexSrc()->getPoint();
        Point2SPtr p_dst = edge->getVertexDst()->getPoint();
        if (p_ray_src != p_src && p_ray_src != p_dst &&
                KernelWrapper::isInside(p_intersect, p_src, p_dst)) {
            MeshVertexSPtr vertex_dst = MeshVertex::create(p_intersect);
            mesh_result_->addVertex(vertex_dst);
            ray->setDst(vertex_dst);
        } else {
            MeshEdgeDataSPtr data = std::dynamic_pointer_cast<MeshEdgeData>(
                    edge->getData());
            double dist_max = 0.0;
            ArcSPtr arc_max;
            std::list<ArcWPtr>::iterator it_a = data->arcs().begin();
            while (it_a != data->arcs().end()) {
                ArcWPtr arc_wptr = *it_a++;
                if (!arc_wptr.expired()) {
                    ArcSPtr arc(arc_wptr);
                    Line2SPtr l_arc = arc->line();
                    Point2SPtr p_arc_inter = KernelWrapper::intersection(l_ray, l_arc);
                    p_src = arc->getNodeSrc()->getPoint();
                    p_dst = arc->getNodeDst()->getPoint();
                    if (p_ray_src != p_src && p_ray_src != p_dst &&
                            KernelWrapper::isInside(p_arc_inter, p_src, p_dst)) {
                        double dist = KernelWrapper::distance(p_ray_src, p_arc_inter);
                        if (dist >= dist_max) {
                            dist_max = dist;
                            arc_max = arc;
                            p_intersect = p_arc_inter;
                        }
                    }
                }
            }
            if (arc_max) {
                MeshVertexSPtr vertex = MeshVertex::create(p_intersect);
                mesh_result_->addVertex(vertex);
                ray->setDst(vertex);
                EdgeSPtr edge_next = arc_max->getEdgeLeft();
                if (edge_next == edge) {
                    edge_next = arc_max->getEdgeRight();
                }
                MeshRaySPtr ray_next = MeshRay::create(edge_next, vertex);
                MeshEdgeDataSPtr data_next = std::dynamic_pointer_cast<MeshEdgeData>(
                        edge_next->getData());
                data_next->addRay(ray_next);
                mesh_result_->addRay(ray_next);
                rays.push_back(ray_next);
            } else {
                DEBUG_PRINT("Error: Unable to find intersection");
            }
        }
    }
}

void SkelMeshGenerator::sortRays(EdgeSPtr edge) {
    MeshEdgeDataSPtr data = std::dynamic_pointer_cast<MeshEdgeData>(edge->getData());
    Point2SPtr p_src = edge->getVertexSrc()->getPoint();
    Point2SPtr p_dst = edge->getVertexDst()->getPoint();
    Vector2SPtr v_dir = KernelFactory::createVector2(*p_dst - *p_src);
    std::list<MeshRaySPtr> rays;
    std::list<MeshRayWPtr>::iterator it_r_wptr = data->rays().begin();
    while (it_r_wptr != data->rays().end()) {
        MeshRayWPtr ray_wptr = *it_r_wptr++;
        if (!ray_wptr.expired()) {
            MeshRaySPtr ray(ray_wptr);
            rays.push_back(ray);
        }
    }
    data->rays().clear();
    while (!rays.empty()) {
        MeshRaySPtr ray_min = rays.front();
        std::list<MeshRaySPtr>::iterator it_min = rays.begin();
        std::list<MeshRaySPtr>::iterator it_r = rays.begin();
        while (it_r != rays.end()) {
            std::list<MeshRaySPtr>::iterator it_current = it_r;
            MeshRaySPtr ray = *it_r++;
            if (KernelWrapper::compatePoints(v_dir,
                    ray_min->getSrc()->getPoint(), ray->getSrc()->getPoint()) > 0) {
                ray_min = ray;
                it_min = it_current;
            }
        }
        MeshRayWPtr ray_min_wptr(ray_min);
        data->rays().push_back(ray_min_wptr);
        rays.erase(it_min);
    }
}

void SkelMeshGenerator::splitCellR(EdgeSPtr edge) {
    WriteLock l(mesh_result_->mutex());
    MeshEdgeDataSPtr data = std::dynamic_pointer_cast<MeshEdgeData>(edge->getData());
    MeshCellSPtr cell = MeshCellSPtr(data->cells().front());
    Vector2SPtr dir_line = KernelFactory::createVector2(edge->line());
    Vector2SPtr dir_inside = KernelWrapper::perpendicular(dir_line);
    std::list<MeshRayWPtr>::iterator it_r = data->rays().begin();
    while (it_r != data->rays().end()) {
        MeshRayWPtr ray_wptr = *it_r++;
        if (!ray_wptr.expired()) {
            MeshRaySPtr ray(ray_wptr);
            Vector2SPtr dir_ray = KernelFactory::createVector2(
                    *(ray->getDst()->getPoint()) - *(ray->getSrc()->getPoint()));
            MeshCellSPtr cell_left;
            if ((*dir_ray * *dir_inside) > 0) {  // angle < M_PI/2.0
                cell_left = splitCell(cell, ray->getSrc(), ray->getDst());
            } else {
                cell_left = splitCell(cell, ray->getDst(), ray->getSrc());
            }
            data->addCell(cell_left);
        }
    }
}


MeshVertexSPtr SkelMeshGenerator::findNearestIntersection(MeshCellSPtr cell,
        MeshVertexSPtr vertex, Vector2SPtr direction) {
    MeshVertexSPtr result;
    Point2SPtr point = vertex->getPoint();
    Line2SPtr line = KernelFactory::createLine2(point, direction);
    double dist_min = std::numeric_limits<double>::max();
    MeshVertexSPtr vertex_current;
    MeshVertexSPtr vertex_prev = cell->vertices().back();
    std::list<MeshVertexSPtr>::iterator it_v = cell->vertices().begin();
    while (it_v != cell->vertices().end()) {
        vertex_current = *it_v++;
        if (vertex_current == vertex) {
            return vertex_current;
        }
        Point2SPtr p_src = vertex_prev->getPoint();
        Point2SPtr p_dst = vertex_current->getPoint();
        int side_src = KernelWrapper::side(line, p_src);
        int side_dst = KernelWrapper::side(line, p_dst);
        if ((side_src + side_dst) == 0) {
            Line2SPtr l_current = KernelFactory::createLine2(p_src, p_dst);
            Point2SPtr p_intersect = KernelWrapper::intersection(line, l_current);
            double dist = KernelWrapper::distance(point, p_intersect);
            if (dist < dist_min) {
                dist_min = dist;
                result = vertex_current;
            }
        }
        vertex_prev = vertex_current;
    }
    return result;
}

MeshCellSPtr SkelMeshGenerator::splitCell(MeshCellSPtr cell,
        MeshVertexSPtr v_src, MeshVertexSPtr v_dst) {
    Vector2SPtr dir = KernelFactory::createVector2(
            *(v_dst->getPoint()) - *(v_src->getPoint()));
    if (!cell->containsVertex(v_src)) {
        MeshVertexSPtr v_edge_dst = findNearestIntersection(cell, v_src, dir);
        MeshVertexSPtr v_edge_src = v_edge_dst->prev(cell);
        MeshModifier::splitEdge(v_edge_src, v_edge_dst, v_src);
    }
    if (!cell->containsVertex(v_dst)) {
        MeshVertexSPtr v_edge_dst = findNearestIntersection(cell, v_dst, dir);
        MeshVertexSPtr v_edge_src = v_edge_dst->prev(cell);
        MeshModifier::splitEdge(v_edge_src, v_edge_dst, v_dst);
    }
    MeshCellSPtr result = MeshModifier::splitCell(cell, v_src, v_dst);
    return result;
}


void SkelMeshGenerator::offsetEdges() {
    std::list<AbstractEventSPtr>::iterator it_ev = skel_->events().begin();
    while (it_ev != skel_->events().end()) {
        AbstractEventSPtr event = *it_ev++;
        if (event->getType() == AbstractEvent::CONST_OFFSET_EVENT) {
            PolygonSPtr polygon_offset = event->getPolygonResult();
            std::list<EdgeSPtr>::iterator it_e = polygon_offset->edges().begin();
            while (it_e != polygon_offset->edges().end()) {
                EdgeSPtr edge_offset = *it_e++;
                SkelVertexDataSPtr data_src =
                        std::dynamic_pointer_cast<SkelVertexData>(
                        edge_offset->getVertexSrc()->getData());
                EdgeSPtr edge_origin = data_src->getArc()->getEdgeRight();
                MeshEdgeDataSPtr data_origin =
                        std::dynamic_pointer_cast<MeshEdgeData>(
                        edge_origin->getData());
                data_origin->addEdge(edge_offset);
            }
        }
    }
}

void SkelMeshGenerator::splitCellsE(EdgeSPtr edge, double radius_snap) {
    WriteLock l(mesh_result_->mutex());
    MeshEdgeDataSPtr data = std::dynamic_pointer_cast<MeshEdgeData>(edge->getData());
    std::list<MeshCellSPtr> cells;
    std::list<MeshCellWPtr>::iterator it_c = data->cells().begin();
    while (it_c != data->cells().end()) {
        MeshCellWPtr cell_wptr = *it_c++;
        if (!cell_wptr.expired()) {
            MeshCellSPtr cell(cell_wptr);
            cells.push_back(cell);
        }
    }

    while (!cells.empty()) {
        MeshCellSPtr cell = cells.front();
        cells.pop_front();
        std::list<EdgeWPtr>::iterator it_e = data->edges().begin();
        while (it_e != data->edges().end()) {
            EdgeWPtr edge_wptr = *it_e++;
            if (edge_wptr.expired()) continue;
            EdgeSPtr edge(edge_wptr);
            cell = splitCell(cell, edge->line(), radius_snap);
            if (!cell) break;
            if (cell->vertices().size() == 0) break;
        }
    }
}

MeshCellSPtr SkelMeshGenerator::splitCell(MeshCellSPtr cell, Line2SPtr line,
        double radius_snap) {
    int side = 0;
    MeshCellSPtr result;
    MeshVertexSPtr v_intersect_src;
    MeshVertexSPtr v_intersect_dst;
    MeshVertexSPtr pos_insert_src;
    MeshVertexSPtr pos_insert_dst;
    MeshVertexSPtr vertex_current;
    MeshVertexSPtr vertex_prev = cell->vertices().back();
    std::list<MeshVertexSPtr>::iterator it_v = cell->vertices().begin();
    while (it_v != cell->vertices().end()) {
        vertex_current = *it_v++;
        Point2SPtr p_src = vertex_prev->getPoint();
        Point2SPtr p_dst = vertex_current->getPoint();
        // TODO: improve floating point kung fu
        Vector2SPtr dir = KernelWrapper::normalize(
                KernelFactory::createVector2(*p_dst - *p_src));
        int side_src = KernelWrapper::side(line,
                KernelFactory::createPoint2(*p_src - *dir * (radius_snap/1000.0)));
        int side_dst = KernelWrapper::side(line,
                KernelFactory::createPoint2(*p_dst + *dir * (radius_snap/1000.0)));
        if ((side_src + side_dst) == 0) {
            Line2SPtr l_current = KernelFactory::createLine2(p_src, p_dst);
            Point2SPtr p_intersect = KernelWrapper::intersection(line, l_current);
            if (side_src > 0) {
                if (KernelWrapper::distance(p_src, p_intersect) < radius_snap) {
                    v_intersect_src = vertex_prev;
                } else if (KernelWrapper::distance(p_dst, p_intersect) < radius_snap) {
                    v_intersect_src = vertex_current;
                } else {
                    v_intersect_src = MeshVertex::create(p_intersect);
                    pos_insert_src = vertex_current;
                }
            } else {
                if (KernelWrapper::distance(p_src, p_intersect) < radius_snap) {
                    v_intersect_dst = vertex_prev;
                } else if (KernelWrapper::distance(p_dst, p_intersect) < radius_snap) {
                    v_intersect_dst = vertex_current;
                } else {
                    v_intersect_dst = MeshVertex::create(p_intersect);
                    pos_insert_dst = vertex_current;
                }
            }
        }
        side += side_dst;
        vertex_prev = vertex_current;
    }

    if (v_intersect_src && v_intersect_dst &&
            v_intersect_src != v_intersect_dst) {
        if (pos_insert_src && pos_insert_dst) {
            mesh_result_->addVertex(v_intersect_src);
            mesh_result_->addVertex(v_intersect_dst);
            MeshModifier::splitEdge(pos_insert_src->prev(cell), pos_insert_src, v_intersect_src);
            MeshModifier::splitEdge(pos_insert_dst->prev(cell), pos_insert_dst, v_intersect_dst);
            result = MeshModifier::splitCell(cell, v_intersect_src, v_intersect_dst);
        } else if (pos_insert_src) {
            if (pos_insert_src != v_intersect_dst->next(cell)) {
                mesh_result_->addVertex(v_intersect_src);
                MeshModifier::splitEdge(pos_insert_src->prev(cell), pos_insert_src, v_intersect_src);
                result = MeshModifier::splitCell(cell, v_intersect_src, v_intersect_dst);
            } else if (side > 0) {
                result = cell;
            }
        } else if (pos_insert_dst) {
            if (pos_insert_dst != v_intersect_src) {
                mesh_result_->addVertex(v_intersect_dst);
                MeshModifier::splitEdge(pos_insert_dst->prev(cell), pos_insert_dst, v_intersect_dst);
                result = MeshModifier::splitCell(cell, v_intersect_src, v_intersect_dst);
            } else if (side > 0) {
                result = cell;
            }
        } else if (v_intersect_src->next(cell) != v_intersect_dst &&
                v_intersect_src->prev(cell) != v_intersect_dst) {
            result = MeshModifier::splitCell(cell, v_intersect_src, v_intersect_dst);
        } else if (side > 0) {
            result = cell;
        }
    } else if (side > 0) {
        result = cell;
    }
    return result;
}


void SkelMeshGenerator::mergeAdjTriangles() {
    WriteLock l(mesh_result_->mutex());
    std::list<NodeSPtr>::iterator it_n = skel_->nodes().begin();
    while (it_n != skel_->nodes().end()) {
        NodeSPtr node = *it_n++;
        if (node->degree() > 1) {
            MeshVertexSPtr vertex = mesh_result_->getVertex(node->getPoint());
            MeshCellSPtr cell_current = vertex->firstCell();
            MeshCellSPtr cell_prev = cell_current->prev(vertex);
            MeshCellSPtr cell_next;
            MeshCellSPtr cell_begin;
            while (cell_current != cell_begin) {
                if (!cell_begin) {
                    cell_begin = cell_current;
                }
                cell_next = cell_current->next(vertex);
                if (cell_prev->vertices().size() == 3 &&
                        cell_current->vertices().size() == 3) {
                    if (cell_current == cell_begin) {
                        cell_begin = cell_prev;
                    }
                    MeshModifier::mergeCells(cell_prev, cell_current);
                }
                cell_prev = cell_current;
                cell_current = cell_next;
            }
        }
    }
    std::list<MeshCellSPtr>::iterator it_c = mesh_result_->cells().begin();
    while (it_c != mesh_result_->cells().end()) {
        MeshCellSPtr cell = *it_c;
        if (cell->vertices().size() == 3) {
            std::list<MeshVertexSPtr>::iterator it_v = cell->vertices().begin();
            while (it_v != cell->vertices().end()) {
                MeshVertexSPtr vertex = *it_v++;
                MeshCellSPtr cell_next = cell->next(vertex);
                if (cell_next->vertices().size() == 3) {
                    MeshModifier::mergeCells(cell, cell_next);
                    break;
                }
            }
        }
        it_c++;
    }
}

void SkelMeshGenerator::mergeVertices() {
    WriteLock l(mesh_result_->mutex());
    std::list<MeshVertexSPtr> vertices;
    std::map<Point2SPtr, MeshVertexSPtr>::iterator it_mv = mesh_result_->vertices().begin();
    while (it_mv != mesh_result_->vertices().end()) {
        MeshVertexSPtr vertex = it_mv->second;
        if (vertex->countCells() == 3) {
            vertices.push_back(vertex);
        }
        it_mv++;
    }
    std::list<MeshVertexSPtr>::iterator it_v = vertices.begin();
    while (it_v != vertices.end()) {
        MeshVertexSPtr vertex = *it_v;
        std::list<MeshCellWPtr>::iterator it_c = vertex->cells().begin();
        while (it_c != vertex->cells().end()) {
            MeshCellWPtr cell_wptr = *it_c++;
            if (cell_wptr.expired()) continue;
            MeshCellSPtr cell(cell_wptr);
            MeshVertexSPtr vertex_next = vertex->next(cell);
            if (vertex_next->countCells() == 3) {
                MeshModifier::mergeVertices(vertex, vertex_next);
                vertices.remove(vertex_next);
                break;
            }
        }
        it_v++;
    }
}


void SkelMeshGenerator::run() {
    DEBUG_PRINT("== Skeleton Mesh Generator started ==");
    initEdgeDatas();
    initCells();
    if (controller_) controller_->wait();
    PolygonSPtr polygon = skel_->getPolygon();
    std::list<EdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        EdgeSPtr edge = *it_e++;
        createRays(edge);
    }
    findRayDsts();
    it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        EdgeSPtr edge = *it_e++;
        sortRays(edge);
        splitCellR(edge);
    }
    assert(mesh_result_->isConsistent());
    if (controller_) controller_->wait();
    offsetEdges();
    it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        EdgeSPtr edge = *it_e++;
        splitCellsE(edge, 0.2);
    }
    assert(mesh_result_->isConsistent());
    if (controller_) controller_->wait();
    mergeAdjTriangles();
    assert(mesh_result_->isConsistent());
    if (controller_) controller_->wait();
    mergeVertices();
    assert(mesh_result_->isConsistent());
    DEBUG_PRINT("== Skeleton Mesh Generator finished ==");
    DEBUG_VAR(mesh_result_->toString());
}

ThreadSPtr SkelMeshGenerator::startThread() {
    return ThreadSPtr(new std::thread(
            std::bind(&SkelMeshGenerator::run, this)));
}

MeshSPtr SkelMeshGenerator::getResult() const {
    return this->mesh_result_;
}

} }
