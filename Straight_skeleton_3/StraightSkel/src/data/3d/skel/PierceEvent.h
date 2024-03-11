/**
 * @file   data/3d/skel/PierceEvent.h
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#ifndef DATA_3D_SKEL_PIERCEEVENT_H
#define DATA_3D_SKEL_PIERCEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class PierceEvent : public AbstractEvent {
public:
    virtual ~PierceEvent();
    static PierceEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    double getOffset() const;
    FacetSPtr getFacet() const;
    void setFacet(FacetSPtr facet);
    VertexSPtr getVertex() const;
    void setVertex(VertexSPtr vertex);
    void setHighlight(bool highlight);
protected:
    PierceEvent();
    NodeSPtr node_;
    FacetSPtr facet_;
    VertexSPtr vertex_;
};

} } }

#endif /* DATA_3D_SKEL_PIERCEEVENT_H */
