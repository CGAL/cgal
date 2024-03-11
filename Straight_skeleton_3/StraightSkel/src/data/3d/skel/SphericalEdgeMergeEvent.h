/**
 * @file   data/3d/skel/SphericalEdgeMergeEvent.h
 * @author Gernot Walzl
 * @date   2013-02-08
 */

#ifndef DATA_3D_SKEL_SPHERICALEDGEMERGEEVENT_H
#define DATA_3D_SKEL_SPHERICALEDGEMERGEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalAbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SphericalEdgeMergeEvent : public SphericalAbstractEvent {
public:
    virtual ~SphericalEdgeMergeEvent();
    static SphericalEdgeMergeEventSPtr create();
    CircularNodeSPtr getNode() const;
    void setNode(CircularNodeSPtr node);
    double getOffset() const;
    CircularEdgeSPtr getEdge1() const;
    void setEdge1(CircularEdgeSPtr edge_1);
    CircularEdgeSPtr getEdge2() const;
    void setEdge2(CircularEdgeSPtr edge_2);
    void setHighlight(bool highlight);
protected:
    SphericalEdgeMergeEvent();
    CircularNodeSPtr node_;
    CircularEdgeSPtr edge_1_;
    CircularEdgeSPtr edge_2_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALEDGEMERGEEVENT_H */

