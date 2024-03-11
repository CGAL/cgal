/**
 * @file   data/2d/skel/EdgeEvent.h
 * @author Gernot Walzl
 * @date   2012-02-03
 */

#ifndef DATA_2D_SKEL_EDGEEVENT_H
#define DATA_2D_SKEL_EDGEEVENT_H

#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"
#include "data/2d/skel/AbstractEvent.h"

namespace data { namespace _2d { namespace skel {

class EdgeEvent : public AbstractEvent {
public:
    virtual ~EdgeEvent();
    static EdgeEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    double getOffset() const;
    EdgeSPtr getEdge() const;
    void setEdge(EdgeSPtr edge);
    void setHighlight(bool highlight);
protected:
    EdgeEvent();
    NodeSPtr node_;
    EdgeSPtr edge_;
};

} } }

#endif /* DATA_2D_SKEL_EDGEEVENT_H */

