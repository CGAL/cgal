/**
 * @file   data/2d/mesh/Mesh.h
 * @author Gernot Walzl
 * @date   2014-01-24
 */

#ifndef DATA_2D_MESH_MESH_H
#define DATA_2D_MESH_MESH_H

#include "typedefs_thread.h"
#include "data/2d/ptrs.h"
#include "data/2d/mesh/ptrs.h"
#include <list>
#include <map>
#include <string>

namespace data { namespace _2d { namespace mesh {

class Mesh : public std::enable_shared_from_this<Mesh> {
public:
    virtual ~Mesh();
    static MeshSPtr create();

    bool addVertex(MeshVertexSPtr vertex);
    bool removeVertex(MeshVertexSPtr vertex);
    MeshVertexSPtr getVertex(Point2SPtr point) const;
    void addCell(MeshCellSPtr cell);
    bool removeCell(MeshCellSPtr cell);
    void addRay(MeshRaySPtr ray);
    bool removeRay(MeshRaySPtr ray);

    SharedMutex& mutex();

    std::map<Point2SPtr, MeshVertexSPtr>& vertices();
    std::list<MeshCellSPtr>& cells();
    std::list<MeshRaySPtr>& rays();

    bool isConsistent() const;

    std::string toString() const;

protected:
    Mesh();
    mutable SharedMutex mutex_;
    std::map<Point2SPtr, MeshVertexSPtr> vertices_;
    std::list<MeshCellSPtr> cells_;
    std::list<MeshRaySPtr> rays_;
};

} } }

#endif /* DATA_2D_MESH_MESH_H */
