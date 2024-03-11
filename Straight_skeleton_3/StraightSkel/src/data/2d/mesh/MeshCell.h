/**
 * @file   data/2d/mesh/MeshCell.h
 * @author Gernot Walzl
 * @date   2014-01-24
 */

#ifndef DATA_2D_MESH_MESHCELL_H
#define DATA_2D_MESH_MESHCELL_H

#include "data/2d/ptrs.h"
#include "data/2d/mesh/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _2d { namespace mesh {

class MeshCell : public std::enable_shared_from_this<MeshCell> {
public:
    virtual ~MeshCell();
    static MeshCellSPtr create();

    static MeshCellSPtr create(unsigned int num_vertices, MeshVertexSPtr vertices[]);

    MeshSPtr getMesh() const;
    void setMesh(MeshSPtr mesh);
    std::list<MeshCellSPtr>::iterator getListIt() const;
    void setListIt(std::list<MeshCellSPtr>::iterator list_it);

    void addVertex(MeshVertexSPtr vertex);
    bool removeVertex(MeshVertexSPtr vertex);
    bool addVertexBefore(MeshVertexSPtr position, MeshVertexSPtr vertex);
    bool containsVertex(MeshVertexSPtr vertex) const;

    MeshCellSPtr next(MeshVertexSPtr vertex);
    MeshCellSPtr prev(MeshVertexSPtr vertex);

    std::list<MeshVertexSPtr>& vertices();

    std::string toString() const;

protected:
    MeshCell();
    MeshWPtr mesh_;
    std::list<MeshCellSPtr>::iterator list_it_;
    std::list<MeshVertexSPtr> vertices_;
};

} } }

#endif /* DATA_2D_MESH_MESHCELL_H */
