/**
 * @file   data/2d/mesh/MeshVertex.h
 * @author Gernot Walzl
 * @date   2014-01-24
 */

#ifndef DATA_2D_MESH_MESHVERTEX_H
#define DATA_2D_MESH_MESHVERTEX_H

#include "data/2d/ptrs.h"
#include "data/2d/mesh/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _2d { namespace mesh {

class MeshVertex {
public:
    virtual ~MeshVertex();
    static MeshVertexSPtr create(Point2SPtr point);

    Point2SPtr getPoint() const;
    void setPoint(Point2SPtr point);
    MeshSPtr getMesh() const;
    void setMesh(MeshSPtr mesh);

    void addCell(MeshCellSPtr cell);
    bool removeCell(MeshCellSPtr cell);
    bool containsCell(MeshCellSPtr cell) const;

    MeshCellSPtr firstCell() const;

    MeshVertexSPtr next(MeshCellSPtr cell) const;
    MeshVertexSPtr prev(MeshCellSPtr cell) const;

    std::list<MeshCellWPtr>& cells();

    unsigned int countCells() const;

    double getX() const;
    double getY() const;

    std::string toString() const;

protected:
    MeshVertex(Point2SPtr point);
    Point2SPtr point_;
    MeshWPtr mesh_;
    std::list<MeshCellWPtr> cells_;
};

} } }

#endif /* DATA_2D_MESH_MESHVERTEX_H */
