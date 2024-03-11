/**
 * @file   data/3d/skel/SkelVertexData.h
 * @author Gernot Walzl
 * @date   2012-04-05
 */

#ifndef DATA_3D_SKEL_SKELVERTEXDATA_H
#define DATA_3D_SKEL_SKELVERTEXDATA_H

#include "data/3d/ptrs.h"
#include "data/3d/VertexData.h"
#include "data/3d/skel/ptrs.h"

namespace data { namespace _3d { namespace skel {

class SkelVertexData : public VertexData {
public:
    virtual ~SkelVertexData();

    static SkelVertexDataSPtr create(VertexSPtr vertex);

    ArcSPtr getArc() const;
    void setArc(ArcSPtr arc);
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    VertexSPtr getOffsetVertex() const;
    void setOffsetVertex(VertexSPtr offset_vertex);

protected:
    SkelVertexData();
    ArcWPtr arc_;
    NodeWPtr node_;
    VertexWPtr offset_vertex_;
};

} } }

#endif /* DATA_3D_SKEL_SKELVERTEXDATA_H */

