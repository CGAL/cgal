/**
 * @file   data/3d/skel/SphericalDblEdgeEvent.h
 * @author Gernot Walzl
 * @date   2013-02-22
 */

#ifndef DATA_3D_SKEL_SPHERICALDBLEDGEEVENT_H
#define DATA_3D_SKEL_SPHERICALDBLEDGEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalAbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SphericalDblEdgeEvent : public SphericalAbstractEvent {
public:
    virtual ~SphericalDblEdgeEvent();
    static SphericalDblEdgeEventSPtr create();
    double getOffset() const;
    void setOffset(double offset);
    CircularEdgeSPtr getEdge1() const;
    void setEdge1(CircularEdgeSPtr edge_1);
    CircularEdgeSPtr getEdge2() const;
    void setEdge2(CircularEdgeSPtr edge_2);
    void setHighlight(bool highlight);
protected:
    SphericalDblEdgeEvent();
    double offset_;
    CircularEdgeSPtr edge_1_;
    CircularEdgeSPtr edge_2_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALDBLEDGEEVENT_H */
