/**
 * @file   data/3d/skel/SphericalSplitEvent.h
 * @author Gernot Walzl
 * @date   2012-11-29
 */

#ifndef DATA_3D_SKEL_SPHERICALSPLITEVENT_H
#define DATA_3D_SKEL_SPHERICALSPLITEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalAbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SphericalSplitEvent : public SphericalAbstractEvent {
public:
    virtual ~SphericalSplitEvent();
    static SphericalSplitEventSPtr create();
    CircularNodeSPtr getNode() const;
    void setNode(CircularNodeSPtr node);
    double getOffset() const;
    CircularVertexSPtr getVertex() const;
    void setVertex(CircularVertexSPtr vertex);
    CircularEdgeSPtr getEdge() const;
    void setEdge(CircularEdgeSPtr edge);
    void setHighlight(bool highlight);
protected:
    SphericalSplitEvent();
    CircularNodeSPtr node_;
    CircularVertexSPtr vertex_;
    CircularEdgeSPtr edge_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALSPLITEVENT_H */
