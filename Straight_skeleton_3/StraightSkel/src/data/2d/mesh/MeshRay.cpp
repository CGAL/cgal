/**
 * @file   data/2d/mesh/MeshRay.cpp
 * @author Gernot Walzl
 * @date   2014-02-03
 */

#include "data/2d/mesh/MeshRay.h"

#include "debug.h"

namespace data { namespace _2d { namespace mesh {

MeshRay::MeshRay(EdgeSPtr edge, MeshVertexSPtr src) {
    edge_ = edge;
    src_ = src;
}

MeshRay::~MeshRay() {
    // intentionally does nothing
}

MeshRaySPtr MeshRay::create(EdgeSPtr edge, MeshVertexSPtr src) {
    MeshRaySPtr result = MeshRaySPtr(new MeshRay(edge, src));
    return result;
}

MeshSPtr MeshRay::getMesh() const {
    DEBUG_WPTR(mesh_);
    if (this->mesh_.expired())
        return MeshSPtr();
    else
        return MeshSPtr(this->mesh_);
}

void MeshRay::setMesh(MeshSPtr mesh) {
    this->mesh_ = mesh;
}

std::list<MeshRaySPtr>::iterator MeshRay::getListIt() const {
    return this->list_it_;
}

void MeshRay::setListIt(std::list<MeshRaySPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

EdgeSPtr MeshRay::getEdge() const {
    DEBUG_SPTR(edge_);
    return edge_;
}

void MeshRay::setEdge(EdgeSPtr edge) {
    edge_ = edge;
}

MeshVertexSPtr MeshRay::getSrc() const {
    DEBUG_SPTR(src_);
    return src_;
}

void MeshRay::setSrc(MeshVertexSPtr src) {
    src_ = src;
}

MeshVertexSPtr MeshRay::getDst() const {
    //DEBUG_SPTR(dst_);
    return dst_;
}

void MeshRay::setDst(MeshVertexSPtr dst) {
    dst_ = dst;
}

} } }
