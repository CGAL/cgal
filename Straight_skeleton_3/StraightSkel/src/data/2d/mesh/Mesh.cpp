/**
 * @file   data/2d/mesh/Mesh.cpp
 * @author Gernot Walzl
 * @date   2014-01-24
 */

#include "data/2d/mesh/Mesh.h"

#include "data/2d/mesh/MeshVertex.h"
#include "data/2d/mesh/MeshCell.h"
#include "data/2d/mesh/MeshRay.h"
#include "debug.h"
#include "util/StringFactory.h"
#include <sstream>

namespace data { namespace _2d { namespace mesh {

Mesh::Mesh() {
    // intentionally does nothing
}

Mesh::~Mesh() {
    vertices_.clear();
    cells_.clear();
}

MeshSPtr Mesh::create() {
    return MeshSPtr(new Mesh());
}

bool Mesh::addVertex(MeshVertexSPtr vertex) {
    bool result = false;
    Point2SPtr point = vertex->getPoint();
    if (vertices_.find(point) == vertices_.end()) {
        vertices_[point] = vertex;
        vertex->setMesh(shared_from_this());
        result = true;
    }
    return result;
}

bool Mesh::removeVertex(MeshVertexSPtr vertex) {
    bool result = false;
    if (vertex->getMesh() == shared_from_this()) {
        Point2SPtr point = vertex->getPoint();
        std::map<Point2SPtr, MeshVertexSPtr>::iterator it = vertices_.find(point);
        vertices_.erase(it);
        vertex->setMesh(MeshSPtr());
        std::list<MeshCellWPtr>::iterator it_c = vertex->cells().begin();
        while (it_c != vertex->cells().end()) {
            MeshCellWPtr cell_wptr = *it_c++;
            if (!cell_wptr.expired()) {
                MeshCellSPtr cell(cell_wptr);
                this->removeCell(cell);
            }
        }
        result = true;
    }
    return result;
}

MeshVertexSPtr Mesh::getVertex(Point2SPtr point) const {
    MeshVertexSPtr result;
    std::map<Point2SPtr, MeshVertexSPtr>::const_iterator it = vertices_.find(point);
    if (it != vertices_.end()) {
        result = it->second;
    }
    //DEBUG_SPTR(result);
    return result;
}

void Mesh::addCell(MeshCellSPtr cell) {
    std::list<MeshCellSPtr>::iterator it = cells_.insert(cells_.end(), cell);
    cell->setMesh(shared_from_this());
    cell->setListIt(it);
    std::list<MeshVertexSPtr>::iterator it_v = cell->vertices().begin();
    while (it_v != cell->vertices().end()) {
        MeshVertexSPtr vertex = *it_v++;
        if (vertex->getMesh() != shared_from_this()) {
            this->addVertex(vertex);
        }
    }
}

bool Mesh::removeCell(MeshCellSPtr cell) {
    bool result = false;
    if (cell->getMesh() == shared_from_this()) {
        cells_.erase(cell->getListIt());
        cell->setMesh(MeshSPtr());
        cell->setListIt(std::list<MeshCellSPtr>::iterator());
        result = true;
    }
    return result;
}

void Mesh::addRay(MeshRaySPtr ray) {
    std::list<MeshRaySPtr>::iterator it = rays_.insert(rays_.end(), ray);
    ray->setMesh(shared_from_this());
    ray->setListIt(it);
    MeshVertexSPtr vertex = ray->getSrc();
    if (vertex->getMesh() != shared_from_this()) {
        this->addVertex(vertex);
    }
    vertex = ray->getDst();
    if (vertex) {
        if (vertex->getMesh() != shared_from_this()) {
            this->addVertex(vertex);
        }
    }
}

bool Mesh::removeRay(MeshRaySPtr ray) {
    bool result = false;
    if (ray->getMesh() == shared_from_this()) {
        rays_.erase(ray->getListIt());
        ray->setMesh(MeshSPtr());
        ray->setListIt(std::list<MeshRaySPtr>::iterator());
        result = true;
    }
    return result;
}

SharedMutex& Mesh::mutex() {
    return this->mutex_;
}

std::map<Point2SPtr, MeshVertexSPtr>& Mesh::vertices() {
    return this->vertices_;
}

std::list<MeshCellSPtr>& Mesh::cells() {
    return this->cells_;
}

std::list<MeshRaySPtr>& Mesh::rays() {
    return this->rays_;
}

bool Mesh::isConsistent() const {
    bool result = true;

    std::map<Point2SPtr, MeshVertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        Point2SPtr point = it_v->first;
        MeshVertexSPtr vertex = it_v->second;
        it_v++;
        if (vertex->getPoint() != point) {
            result = false;
            break;
        }
        if (vertex->getMesh() != shared_from_this()) {
            result = false;
            break;
        }
        std::list<MeshCellWPtr>::const_iterator it_c = vertex->cells().begin();
        while (it_c != vertex->cells().end()) {
            MeshCellWPtr cell_wptr = *it_c++;
            if (!cell_wptr.expired()) {
                MeshCellSPtr cell(cell_wptr);
                if (!cell->containsVertex(vertex)) {
                    result = false;
                    break;
                }
            }
        }
    }

    std::list<MeshCellSPtr>::const_iterator it_c = cells_.begin();
    while (it_c != cells_.end()) {
        MeshCellSPtr cell = *it_c++;
        if (cell->getMesh() != shared_from_this()) {
            result = false;
            break;
        }
        std::list<MeshVertexSPtr>::const_iterator it_v = cell->vertices().begin();
        while (it_v != cell->vertices().end()) {
            MeshVertexSPtr vertex = *it_v++;
            if (!vertex->containsCell(cell)) {
                result = false;
                break;
            }
        }
    }

    return result;
}

std::string Mesh::toString() const {
    std::stringstream sstr;
    sstr << "Mesh(";
    sstr << util::StringFactory::fromPointer(this) << ", \n";
    std::map<unsigned int, unsigned int> cell_sizes;
    std::list<MeshCellSPtr>::const_iterator it_m = cells_.begin();
    while (it_m != cells_.end()) {
        MeshCellSPtr cell = *it_m++;
        unsigned int size = cell->vertices().size();
        if (cell_sizes.find(size) == cell_sizes.end()) {
            cell_sizes[size] = 0;
        }
        cell_sizes[size] += 1;
    }
    std::map<unsigned int, unsigned int>::iterator it = cell_sizes.begin();
    while (it != cell_sizes.end()) {
        sstr << it->first << ": " << it->second << "\n";
        it++;
    }
    sstr << ")\n";
    return sstr.str();
}

} } }
