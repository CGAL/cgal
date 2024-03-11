/**
 * @file   data/3d/skel/SurfaceEvent.h
 * @author Gernot Walzl
 * @date   2012-09-10
 */

#ifndef DATA_3D_SKEL_SURFACEEVENT_H
#define DATA_3D_SKEL_SURFACEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SurfaceEvent : public AbstractEvent {
public:
    virtual ~SurfaceEvent();
    static SurfaceEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    double getOffset() const;
    EdgeSPtr getEdge1() const;
    void setEdge1(EdgeSPtr edge1);
    EdgeSPtr getEdge2() const;
    void setEdge2(EdgeSPtr edge2);
    void setHighlight(bool highlight);
protected:
    SurfaceEvent();
    NodeSPtr node_;
    EdgeSPtr edge1_;
    EdgeSPtr edge2_;
};

} } }

#endif /* DATA_3D_SKEL_SURFACEEVENT_H */

