/**
 * @file   data/2d/mesh/MeshVertex.cpp
 * @author Gernot Walzl
 * @date   2014-01-24
 */

#include "data/2d/mesh/MeshVertex.h"

#include "data/2d/mesh/MeshCell.h"
#include "debug.h"
#include "util/StringFactory.h"

namespace data { namespace _2d { namespace mesh {

MeshVertex::MeshVertex(Point2SPtr point) {
    point_ = point;
}

MeshVertex::~MeshVertex() {
    point_.reset();
}

MeshVertexSPtr MeshVertex::create(Point2SPtr point) {
    MeshVertexSPtr result = MeshVertexSPtr(new MeshVertex(point));
    return result;
}

Point2SPtr MeshVertex::getPoint() const {
    return point_;
}

void MeshVertex::setPoint(Point2SPtr point) {
    this->point_ = point;
}

MeshSPtr MeshVertex::getMesh() const {
    DEBUG_WPTR(mesh_);
    if (this->mesh_.expired())
        return MeshSPtr();
    else
        return MeshSPtr(this->mesh_);
}

void MeshVertex::setMesh(MeshSPtr mesh) {
    this->mesh_ = mesh;
}

void MeshVertex::addCell(MeshCellSPtr cell) {
    cells_.insert(cells_.end(), MeshCellWPtr(cell));
}

bool MeshVertex::removeCell(MeshCellSPtr cell) {
    bool result = false;
    std::list<MeshCellWPtr>::iterator it = cells_.begin();
    while (it != cells_.end()) {
        std::list<MeshCellWPtr>::iterator it_current = it;
        MeshCellWPtr cell_wptr = *it++;
        if (!cell_wptr.expired()) {
            if (cell_wptr.lock() == cell) {
                cells_.erase(it_current);
                result = true;
                break;
            }
        }
    }
    return result;
}

bool MeshVertex::containsCell(MeshCellSPtr cell) const {
    bool result = false;
    std::list<MeshCellWPtr>::const_iterator it = cells_.begin();
    while (it != cells_.end()) {
        MeshCellWPtr cell_wptr = *it++;
        if (!cell_wptr.expired()) {
            if (cell_wptr.lock() == cell) {
                result = true;
                break;
            }
        }
    }
    return result;
}


MeshCellSPtr MeshVertex::firstCell() const {
    MeshCellSPtr result;
    std::list<MeshCellWPtr>::const_iterator it_c = cells_.begin();
    while (it_c != cells_.end()) {
        MeshCellWPtr cell_wptr = *it_c++;
        if (!cell_wptr.expired()) {
            result = MeshCellSPtr(cell_wptr);
            break;
        }
    }
    DEBUG_SPTR(result);
    return result;
}


MeshVertexSPtr MeshVertex::next(MeshCellSPtr cell) const {
    MeshVertexSPtr result;
    std::list<MeshVertexSPtr>::const_iterator it_v = cell->vertices().begin();
    while (it_v != cell->vertices().end()) {
        MeshVertexSPtr vertex = *it_v++;
        if (vertex.get() == this) {
            if (it_v == cell->vertices().end()) {
                result = cell->vertices().front();
            } else {
                result = *it_v;
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

MeshVertexSPtr MeshVertex::prev(MeshCellSPtr cell) const {
    MeshVertexSPtr result;
    std::list<MeshVertexSPtr>::const_reverse_iterator it_v = cell->vertices().rbegin();
    while (it_v != cell->vertices().rend()) {
        MeshVertexSPtr vertex = *it_v++;
        if (vertex.get() == this) {
            if (it_v == cell->vertices().rend()) {
                result = cell->vertices().back();
            } else {
                result = *it_v;
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

std::list<MeshCellWPtr>& MeshVertex::cells() {
    return this->cells_;
}

unsigned int MeshVertex::countCells() const {
    unsigned int result = 0;
    std::list<MeshCellWPtr>::const_iterator it_c = cells_.begin();
    while (it_c != cells_.end()) {
        MeshCellWPtr cell_wptr = *it_c++;
        if (!cell_wptr.expired()) {
            result++;
        }
    }
    return result;
}

double MeshVertex::getX() const {
#ifdef USE_CGAL
    return this->point_->x();
#else
    return this->point_->getX();
#endif
}

double MeshVertex::getY() const {
#ifdef USE_CGAL
    return this->point_->y();
#else
    return this->point_->getY();
#endif
}

std::string MeshVertex::toString() const {
    std::string result("MeshVertex(");
    result += util::StringFactory::fromPointer(this) + ", ";
    result += "<" + util::StringFactory::fromDouble(getX()) + ", ";
    result += util::StringFactory::fromDouble(getY()) + ">)";
    return result;
}

} } }
