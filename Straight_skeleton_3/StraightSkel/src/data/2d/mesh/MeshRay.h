/**
 * @file   data/2d/mesh/MeshRay.h
 * @author Gernot Walzl
 * @date   2014-02-03
 */

#ifndef DATA_2D_MESH_MESHRAY_H
#define DATA_2D_MESH_MESHRAY_H

#include "data/2d/ptrs.h"
#include "data/2d/mesh/ptrs.h"
#include <list>

namespace data { namespace _2d { namespace mesh {

class MeshRay {
public:
    virtual ~MeshRay();
    static MeshRaySPtr create(EdgeSPtr edge, MeshVertexSPtr src);

    MeshSPtr getMesh() const;
    void setMesh(MeshSPtr mesh);
    std::list<MeshRaySPtr>::iterator getListIt() const;
    void setListIt(std::list<MeshRaySPtr>::iterator list_it);

    EdgeSPtr getEdge() const;
    void setEdge(EdgeSPtr edge);
    MeshVertexSPtr getSrc() const;
    void setSrc(MeshVertexSPtr src);
    MeshVertexSPtr getDst() const;
    void setDst(MeshVertexSPtr dst);

protected:
    MeshRay(EdgeSPtr edge, MeshVertexSPtr src);
    MeshWPtr mesh_;
    std::list<MeshRaySPtr>::iterator list_it_;
    EdgeSPtr edge_;
    MeshVertexSPtr src_;
    MeshVertexSPtr dst_;
};

} } }

#endif /* DATA_2D_MESH_MESHRAY_H */

