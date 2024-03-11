/**
 * @file   data/3d/skel/SphericalTriangleEvent.h
 * @author Gernot Walzl
 * @date   2012-11-29
 */

#ifndef DATA_3D_SKEL_SPHERICALTRIANGLEEVENT_H
#define DATA_3D_SKEL_SPHERICALTRIANGLEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalAbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SphericalTriangleEvent : public SphericalAbstractEvent {
public:
    virtual ~SphericalTriangleEvent();
    static SphericalTriangleEventSPtr create();
    CircularNodeSPtr getNode() const;
    void setNode(CircularNodeSPtr node);
    double getOffset() const;
    CircularEdgeSPtr getEdgeBegin() const;
    void setEdgeBegin(CircularEdgeSPtr edge_begin);
    void getVertices(CircularVertexSPtr out[3]) const;
    void getEdges(CircularEdgeSPtr out[3]) const;
    void setHighlight(bool highlight);
protected:
    SphericalTriangleEvent();
    CircularNodeSPtr node_;
    CircularEdgeSPtr edge_begin_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALTRIANGLEEVENT_H */

