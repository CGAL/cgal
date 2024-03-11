/**
 * @file   data/2d/skel/SkelVertexData.h
 * @author Gernot Walzl
 * @date   2012-02-09
 */

#ifndef DATA_2D_SKEL_SKELVERTEXDATA_H
#define DATA_2D_SKEL_SKELVERTEXDATA_H

#include "debug.h"
#include "data/2d/ptrs.h"
#include "data/2d/VertexData.h"
#include "data/2d/skel/ptrs.h"

namespace data { namespace _2d { namespace skel {

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

#endif /* DATA_2D_SKEL_SKELVERTEXDATA_H */

