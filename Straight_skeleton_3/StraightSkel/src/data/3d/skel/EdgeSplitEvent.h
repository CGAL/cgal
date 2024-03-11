/**
 * @file   data/3d/skel/EdgeSplitEvent.h
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#ifndef DATA_3D_SKEL_EDGESPLITEVENT_H
#define DATA_3D_SKEL_EDGESPLITEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class EdgeSplitEvent : public AbstractEvent {
public:
    virtual ~EdgeSplitEvent();
    static EdgeSplitEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    double getOffset() const;
    EdgeSPtr getEdge1() const;
    void setEdge1(EdgeSPtr edge1);
    EdgeSPtr getEdge2() const;
    void setEdge2(EdgeSPtr edge2);
    void setHighlight(bool highlight);
protected:
    EdgeSplitEvent();
    NodeSPtr node_;
    EdgeSPtr edge1_;
    EdgeSPtr edge2_;
};

} } }

#endif /* DATA_3D_SKEL_EDGESPLITEVENT_H */

