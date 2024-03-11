/**
 * @file   data/2d/mesh/MeshCell.cpp
 * @author Gernot Walzl
 * @date   2014-01-24
 */

#include "data/2d/mesh/MeshCell.h"

#include "debug.h"
#include "data/2d/mesh/MeshVertex.h"
#include "util/StringFactory.h"
#include <algorithm>

namespace data { namespace _2d { namespace mesh {

MeshCell::MeshCell() {
}

MeshCell::~MeshCell() {
}

MeshCellSPtr MeshCell::create() {
    return MeshCellSPtr(new MeshCell());
}

MeshCellSPtr MeshCell::create(unsigned int num_vertices, MeshVertexSPtr vertices[]) {
    MeshCellSPtr result = MeshCellSPtr(new MeshCell());;
    for (unsigned int i = 0; i < num_vertices; i++) {
        result->addVertex(vertices[i]);
    }
    return result;
}

MeshSPtr MeshCell::getMesh() const {
    DEBUG_WPTR(mesh_);
    if (this->mesh_.expired())
        return MeshSPtr();
    else
        return MeshSPtr(this->mesh_);
}

void MeshCell::setMesh(MeshSPtr mesh) {
    this->mesh_ = mesh;
}

std::list<MeshCellSPtr>::iterator MeshCell::getListIt() const {
    return this->list_it_;
}

void MeshCell::setListIt(std::list<MeshCellSPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

void MeshCell::addVertex(MeshVertexSPtr vertex) {
    vertices_.insert(vertices_.end(), vertex);
    vertex->addCell(shared_from_this());
}

bool MeshCell::removeVertex(MeshVertexSPtr vertex) {
    bool result = false;
    std::list<MeshVertexSPtr>::iterator it_v =
            std::find(vertices_.begin(), vertices_.end(), vertex);
    if (it_v != vertices_.end()) {
        vertices_.erase(it_v);
        vertex->removeCell(shared_from_this());
        result = true;
    }
    return result;
}

bool MeshCell::addVertexBefore(MeshVertexSPtr position, MeshVertexSPtr vertex) {
    bool result = false;
    std::list<MeshVertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        std::list<MeshVertexSPtr>::iterator it_current = it_v;
        MeshVertexSPtr v_current = *it_v++;
        if (v_current == position) {
            vertices_.insert(it_current, vertex);
            vertex->addCell(shared_from_this());
            result = true;
            break;
        }
    }
    return result;
}

bool MeshCell::containsVertex(MeshVertexSPtr vertex) const {
    bool result = false;
    std::list<MeshVertexSPtr>::const_iterator it_v =
            std::find(vertices_.begin(), vertices_.end(), vertex);
    if (it_v != vertices_.end()) {
        result = true;
    }
    return result;
}

MeshCellSPtr MeshCell::next(MeshVertexSPtr vertex) {
    MeshCellSPtr result;
    MeshVertexSPtr v_prev = vertex->prev(shared_from_this());
    std::list<MeshCellWPtr>::const_iterator it_c1 = vertex->cells().begin();
    while (it_c1 != vertex->cells().end()) {
        MeshCellWPtr cell_1_wptr = *it_c1++;
        if (!cell_1_wptr.expired()) {
            MeshCellSPtr cell_1(cell_1_wptr);
            if (cell_1.get() == this) {
                continue;
            }
            std::list<MeshCellWPtr>::const_iterator it_c2 = v_prev->cells().begin();
            while (it_c2 != v_prev->cells().end()) {
                MeshCellWPtr cell_2_wptr = *it_c2++;
                if (!cell_2_wptr.expired()) {
                    MeshCellSPtr cell_2(cell_2_wptr);
                    if (cell_1 == cell_2) {
                        result = cell_1;
                    }
                }
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

MeshCellSPtr MeshCell::prev(MeshVertexSPtr vertex) {
    MeshCellSPtr result;
    MeshVertexSPtr v_next = vertex->next(shared_from_this());
    std::list<MeshCellWPtr>::const_iterator it_c1 = vertex->cells().begin();
    while (it_c1 != vertex->cells().end()) {
        MeshCellWPtr cell_1_wptr = *it_c1++;
        if (!cell_1_wptr.expired()) {
            MeshCellSPtr cell_1(cell_1_wptr);
            if (cell_1.get() == this) {
                continue;
            }
            std::list<MeshCellWPtr>::const_iterator it_c2 = v_next->cells().begin();
            while (it_c2 != v_next->cells().end()) {
                MeshCellWPtr cell_2_wptr = *it_c2++;
                if (!cell_2_wptr.expired()) {
                    MeshCellSPtr cell_2(cell_2_wptr);
                    if (cell_1 == cell_2) {
                        result = cell_1;
                    }
                }
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

std::list<MeshVertexSPtr>& MeshCell::vertices() {
    return this->vertices_;
}

std::string MeshCell::toString() const {
    std::string result("MeshCell(");
    result += util::StringFactory::fromPointer(this) + ",\n";
    std::list<MeshVertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        MeshVertexSPtr vertex = *it_v++;
        result += vertex->toString() + "\n";
    }
    result += ")";
    return result;
}

} } }
