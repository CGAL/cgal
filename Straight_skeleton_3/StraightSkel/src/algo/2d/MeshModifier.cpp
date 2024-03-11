/**
 * @file   algo/2d/MeshModifier.cpp
 * @author Gernot Walzl
 * @date   2014-04-27
 */

#include "algo/2d/MeshModifier.h"

#include "debug.h"
#include "algo/2d/KernelWrapper.h"
#include "data/2d/mesh/Mesh.h"
#include "data/2d/mesh/MeshCell.h"
#include "data/2d/mesh/MeshVertex.h"
#include "data/2d/KernelFactory.h"

namespace algo { namespace _2d {

MeshModifier::MeshModifier() {
    // intentionally does nothing
}

MeshModifier::~MeshModifier() {
    // intentionally does nothing
}

void MeshModifier::splitEdge(MeshVertexSPtr v_src, MeshVertexSPtr v_dst,
        MeshVertexSPtr v_insert) {
    std::list<MeshCellSPtr> cells_common;
    std::list<MeshCellWPtr>::iterator it_c1 = v_src->cells().begin();
    while (it_c1 != v_src->cells().end()) {
        MeshCellWPtr cell_src_wptr = *it_c1++;
        if (cell_src_wptr.expired()) continue;
        MeshCellSPtr cell_src(cell_src_wptr);
        std::list<MeshCellWPtr>::iterator it_c2 = v_dst->cells().begin();
        while (it_c2 != v_dst->cells().end()) {
            MeshCellWPtr cell_dst_wptr = *it_c2++;
            if (cell_dst_wptr.expired()) continue;
            MeshCellSPtr cell_dst(cell_dst_wptr);
            if (cell_src == cell_dst) {
                cells_common.push_back(cell_src);
            }
        }
    }
    std::list<MeshCellSPtr>::iterator it_c = cells_common.begin();
    while (it_c != cells_common.end()) {
        MeshCellSPtr cell = *it_c++;
        if (v_src->next(cell) == v_dst) {
            cell->addVertexBefore(v_dst, v_insert);
        } else {
            cell->addVertexBefore(v_src, v_insert);
        }
    }
}

MeshCellSPtr MeshModifier::splitCell(MeshCellSPtr cell,
        MeshVertexSPtr v_src, MeshVertexSPtr v_dst) {
    MeshCellSPtr result;
    if (cell->containsVertex(v_src) && cell->containsVertex(v_dst)) {
        result = MeshCell::create();
        Line2SPtr line = KernelFactory::createLine2(v_src->getPoint(), v_dst->getPoint());
        std::list<MeshVertexSPtr>::iterator it_v = cell->vertices().begin();
        while (it_v != cell->vertices().end()) {
            MeshVertexSPtr vertex = *it_v++;
            if (vertex == v_src || vertex == v_dst) {
                result->addVertex(vertex);
                continue;
            }
            if (KernelWrapper::side(line, vertex->getPoint()) > 0) {
                cell->removeVertex(vertex);
                result->addVertex(vertex);
            }
        }
        MeshSPtr mesh = cell->getMesh();
        if (mesh) {
            mesh->addCell(result);
        }
    }
    DEBUG_SPTR(result);
    return result;
}

void MeshModifier::mergeCells(MeshCellSPtr cell_1, MeshCellSPtr cell_2) {
    MeshVertexSPtr vertex_common_1;
    MeshVertexSPtr vertex_common_2;
    MeshVertexSPtr vertex_1_prev;
    std::list<MeshVertexSPtr>::iterator it_v1 = cell_1->vertices().begin();
    while (it_v1 != cell_1->vertices().end()) {
        MeshVertexSPtr vertex_1 = *it_v1++;
        std::list<MeshVertexSPtr>::iterator it_v2 = cell_2->vertices().begin();
        while (it_v2 != cell_2->vertices().end()) {
            MeshVertexSPtr vertex_2 = *it_v2++;
            if (vertex_1 == vertex_2) {
                if (vertex_common_1) {
                    vertex_common_2 = vertex_1;
                } else {
                    vertex_common_1 = vertex_1;
                }
                break;
            }
        }
        if (vertex_common_1 && vertex_common_2) {
            if (!(vertex_common_1 == vertex_1_prev &&
                    vertex_common_2 == vertex_1)) {
                MeshVertexSPtr vertex_tmp = vertex_common_1;
                vertex_common_1 = vertex_common_2;
                vertex_common_2 = vertex_tmp;
            }
            break;
        }
        vertex_1_prev = vertex_1;
    }
    if (vertex_common_1 && vertex_common_2) {
        std::list<MeshVertexSPtr> vertices_tomove;
        MeshVertexSPtr vertex_current = vertex_common_1->next(cell_2);
        while (vertex_current != vertex_common_2) {
            vertices_tomove.push_back(vertex_current);
            vertex_current = vertex_current->next(cell_2);
        }
        std::list<MeshVertexSPtr>::iterator it_v = vertices_tomove.begin();
        while (it_v != vertices_tomove.end()) {
            MeshVertexSPtr vertex = *it_v++;
            cell_2->removeVertex(vertex);
            cell_1->addVertexBefore(vertex_common_2, vertex);
        }
        cell_2->removeVertex(vertex_common_1);
        cell_2->removeVertex(vertex_common_1);
        cell_2->getMesh()->removeCell(cell_2);
    } else {
        DEBUG_PRINT("Error: MeshCells are not adjacent.");
    }
}

void MeshModifier::mergeVertices(MeshVertexSPtr vertex_1,
        MeshVertexSPtr vertex_2) {
    MeshSPtr mesh = vertex_1->getMesh();
    Point2SPtr p_1 = vertex_1->getPoint();
    Point2SPtr p_2 = vertex_2->getPoint();
    Vector2SPtr dir = KernelFactory::createVector2(*p_2 - *p_1);
    Point2SPtr p_1_moved = KernelFactory::createPoint2(*p_1 + (*dir * 0.5));
    mesh->vertices().erase(p_1);
    vertex_1->setPoint(p_1_moved);
    mesh->vertices()[p_1_moved] = vertex_1;

    MeshCellSPtr cell_common_1;
    MeshCellSPtr cell_common_2;
    std::list<MeshCellWPtr>::iterator it_c1 = vertex_1->cells().begin();
    while (it_c1 != vertex_1->cells().end()) {
        MeshCellWPtr cell_1_wptr = *it_c1++;
        if (cell_1_wptr.expired()) continue;
        MeshCellSPtr cell_1(cell_1_wptr);
        std::list<MeshCellWPtr>::iterator it_c2 = vertex_2->cells().begin();
        while (it_c2 != vertex_2->cells().end()) {
            MeshCellWPtr cell_2_wptr = *it_c2++;
            if (cell_2_wptr.expired()) continue;
            MeshCellSPtr cell_2(cell_2_wptr);
            if (cell_1 == cell_2) {
                if (cell_common_1) {
                    cell_common_2 = cell_1;
                } else {
                    cell_common_1 = cell_1;
                }
                break;
            }
        }
        if (cell_common_1 && cell_common_2) {
            break;
        }
    }
    if (cell_common_1 && cell_common_2) {
        if (cell_common_1->next(vertex_1) != cell_common_2) {
            MeshCellSPtr cell_tmp = cell_common_1;
            cell_common_1 = cell_common_2;
            cell_common_2 = cell_tmp;
        }
        std::list<MeshCellSPtr> cells_tomove;
        MeshCellSPtr cell_current = cell_common_1->next(vertex_2);
        while (cell_current != cell_common_2) {
            cells_tomove.push_back(cell_current);
            cell_current = cell_current->next(vertex_2);
        }
        std::list<MeshCellSPtr>::iterator it_c = cells_tomove.begin();
        while (it_c != cells_tomove.end()) {
            MeshCellSPtr cell = *it_c++;
            cell->addVertexBefore(vertex_2, vertex_1);
            cell->removeVertex(vertex_2);
        }
        cell_common_1->removeVertex(vertex_2);
        cell_common_2->removeVertex(vertex_2);
        mesh->removeVertex(vertex_2);
    } else {
        DEBUG_PRINT("Error: MeshVertices are not adjacent.");
    }
}

} }
