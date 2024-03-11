/**
 * @file   data/3d/skel/TriangleEvent.h
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#ifndef DATA_3D_TRIANGLEEVENT_H
#define DATA_3D_TRIANGLEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class TriangleEvent : public AbstractEvent {
public:
    virtual ~TriangleEvent();
    static TriangleEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    double getOffset() const;
    FacetSPtr getFacet() const;
    void setFacet(FacetSPtr facet);
    EdgeSPtr getEdgeBegin() const;
    void setEdgeBegin(EdgeSPtr edge_begin);
    void getVertices(VertexSPtr out[3]) const;
    void getEdges(EdgeSPtr out[3]) const;
    void setHighlight(bool highlight);
protected:
    TriangleEvent();
    NodeSPtr node_;
    FacetSPtr facet_;
    EdgeSPtr edge_begin_;
};

} } }

#endif /* DATA_3D_SKEL_TRIANGLEEVENT_H */

