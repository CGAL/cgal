/**
 * @file   data/3d/skel/TetrahedronEvent.h
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#ifndef DATA_3D_SKEL_TETRAHEDRONEVENT_H
#define DATA_3D_SKEL_TETRAHEDRONEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class TetrahedronEvent : public AbstractEvent {
public:
    virtual ~TetrahedronEvent();
    static TetrahedronEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    double getOffset() const;
    EdgeSPtr getEdgeBegin() const;
    void setEdgeBegin(EdgeSPtr edge_begin);
    void getVertices(VertexSPtr out[4]) const;
    void getEdges(EdgeSPtr out[6]) const;
    void getFacets(FacetSPtr out[4]) const;
    void setHighlight(bool highlight);
protected:
    TetrahedronEvent();
    NodeSPtr node_;
    EdgeSPtr edge_begin_;
};

} } }

#endif /* DATA_3D_SKEL_TETRAHEDRONEVENT_H */

