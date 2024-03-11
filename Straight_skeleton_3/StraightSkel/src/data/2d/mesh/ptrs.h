/**
 * @file   data/2d/mesh/ptrs.h
 * @author Gernot Walzl
 * @date   2014-01-24
 */

#ifndef DATA_2D_MESH_PTRS_H
#define DATA_2D_MESH_PTRS_H

#include "smarter_ptr.h"

namespace data { namespace _2d { namespace mesh {

class Mesh;
class MeshVertex;
class MeshCell;
class MeshRay;

class MeshEdgeData;

typedef SHARED_PTR<Mesh> MeshSPtr;
typedef WEAK_PTR<Mesh> MeshWPtr;
typedef SHARED_PTR<MeshVertex> MeshVertexSPtr;
typedef WEAK_PTR<MeshVertex> MeshVertexWPtr;
typedef SHARED_PTR<MeshCell> MeshCellSPtr;
typedef WEAK_PTR<MeshCell> MeshCellWPtr;
typedef SHARED_PTR<MeshRay> MeshRaySPtr;
typedef WEAK_PTR<MeshRay> MeshRayWPtr;

typedef SHARED_PTR<MeshEdgeData> MeshEdgeDataSPtr;
typedef WEAK_PTR<MeshEdgeData> MeshEdgeDataWPtr;

} } }

#endif /* DATA_2D_MESH_PTRS_H */
