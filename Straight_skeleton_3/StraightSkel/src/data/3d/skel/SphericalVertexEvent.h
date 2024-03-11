/**
 * @file   data/3d/skel/SphericalVertexEvent.h
 * @author Gernot Walzl
 * @date   2013-02-20
 */

#ifndef DATA_3D_SKEL_SPHERICALVERTEXEVENT_H
#define DATA_3D_SKEL_SPHERICALVERTEXEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalAbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SphericalVertexEvent : public SphericalAbstractEvent {
public:
    virtual ~SphericalVertexEvent();
    static SphericalVertexEventSPtr create();
    CircularNodeSPtr getNode() const;
    void setNode(CircularNodeSPtr node);
    double getOffset() const;
    CircularVertexSPtr getVertex1() const;
    void setVertex1(CircularVertexSPtr vertex_1);
    CircularVertexSPtr getVertex2() const;
    void setVertex2(CircularVertexSPtr vertex_2);
    void setHighlight(bool highlight);
protected:
    SphericalVertexEvent();
    CircularNodeSPtr node_;
    CircularVertexSPtr vertex_1_;
    CircularVertexSPtr vertex_2_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALVERTEXEVENT_H */

