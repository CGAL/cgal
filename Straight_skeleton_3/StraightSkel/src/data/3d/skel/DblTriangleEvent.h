/**
 * @file   data/3d/skel/DblTriangleEvent.h
 * @author Gernot Walzl
 * @date   2012-09-11
 */

#ifndef DATA_3D_SKEL_DBLTRIANGLEEVENT_H
#define DATA_3D_SKEL_DBLTRIANGLEEVENT_H


#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class DblTriangleEvent : public AbstractEvent {
public:
    virtual ~DblTriangleEvent();
    static DblTriangleEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    double getOffset() const;
    EdgeSPtr getEdge() const;
    void setEdge(EdgeSPtr edge);
    void getVertices(VertexSPtr out[4]) const;
    void getEdges(EdgeSPtr out[5]) const;
    void setHighlight(bool highlight);
protected:
    DblTriangleEvent();
    NodeSPtr node_;
    EdgeSPtr edge_;
};

} } }

#endif /* DATA_3D_SKEL_DBLTRIANGLEEVENT_H */

