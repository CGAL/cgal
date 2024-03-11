/**
 * @file   data/3d/skel/DblEdgeMergeEvent.h
 * @author Gernot Walzl
 * @date   2012-10-30
 */

#ifndef DATA_3D_SKEL_DBLEDGEMERGEEVENT_H
#define DATA_3D_SKEL_DBLEDGEMERGEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class DblEdgeMergeEvent : public AbstractEvent {
public:
    virtual ~DblEdgeMergeEvent();
    static DblEdgeMergeEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    double getOffset() const;
    FacetSPtr getFacet1() const;
    void setFacet1(FacetSPtr facet_1);
    EdgeSPtr getEdge11() const;
    void setEdge11(EdgeSPtr edge_11);
    EdgeSPtr getEdge12() const;
    void setEdge12(EdgeSPtr edge_12);
    FacetSPtr getFacet2() const;
    void setFacet2(FacetSPtr facet_2);
    EdgeSPtr getEdge21() const;
    void setEdge21(EdgeSPtr edge_21);
    EdgeSPtr getEdge22() const;
    void setEdge22(EdgeSPtr edge_22);
    void getVertices(VertexSPtr out[4]) const;
    void getEdges(EdgeSPtr out[4]) const;
    void setHighlight(bool highlight);
protected:
    DblEdgeMergeEvent();
    NodeSPtr node_;
    FacetSPtr facet_1_;
    EdgeSPtr edge_11_;
    EdgeSPtr edge_12_;
    FacetSPtr facet_2_;
    EdgeSPtr edge_21_;
    EdgeSPtr edge_22_;
};

} } }

#endif /* DATA_3D_SKEL_DBLEDGEMERGEEVENT_H */

