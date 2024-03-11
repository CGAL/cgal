/**
 * @file   data/2d/mesh/MeshEdgeData.h
 * @author Gernot Walzl
 * @date   2014-03-28
 */

#ifndef DATA_2D_MESH_MESHEDGEDATA_H
#define DATA_2D_MESH_MESHEDGEDATA_H

#include "data/2d/ptrs.h"
#include "data/2d/EdgeData.h"
#include "data/2d/skel/ptrs.h"
#include "data/2d/mesh/ptrs.h"
#include <list>

namespace data { namespace _2d { namespace mesh {

using namespace data::_2d;
using namespace data::_2d::skel;

class MeshEdgeData : public EdgeData {
public:
    virtual ~MeshEdgeData();

    static MeshEdgeDataSPtr create(EdgeSPtr edge);

    void addArc(ArcSPtr arc);
    bool removeArc(ArcSPtr arc);
    void addRay(MeshRaySPtr ray);
    bool removeRay(MeshRaySPtr ray);
    void addCell(MeshCellSPtr cell);
    bool removeCell(MeshCellSPtr cell);

    void addEdge(EdgeSPtr edge);
    bool removeEdge(EdgeSPtr edge);

    std::list<ArcWPtr>& arcs();
    std::list<MeshRayWPtr>& rays();
    std::list<MeshCellWPtr>& cells();

    std::list<EdgeWPtr>& edges();

protected:
    MeshEdgeData();
    std::list<ArcWPtr> arcs_;
    std::list<MeshRayWPtr> rays_;
    std::list<MeshCellWPtr> cells_;

    std::list<EdgeWPtr> edges_;
};

} } }

#endif /* DATA_3D_MESH_MESHEDGEDATA_H */
