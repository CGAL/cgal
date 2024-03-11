/**
 * @file   data/3d/skel/SphericalReturnEvent.h
 * @author Gernot Walzl
 * @date   2013-01-28
 */

#ifndef DATA_3D_SKEL_SPHERICALRETURNEVENT_H
#define DATA_3D_SKEL_SPHERICALRETURNEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalAbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SphericalReturnEvent : public SphericalAbstractEvent {
public:
    virtual ~SphericalReturnEvent();

    static SphericalReturnEventSPtr create();

    double getOffset() const;
    void setOffset(double offset);
    CircularVertexSPtr getVertex() const;
    void setVertex(CircularVertexSPtr vertex);

    void setHighlight(bool highlight);

protected:
    SphericalReturnEvent();
    double offset_;
    CircularVertexSPtr vertex_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALRETURNEVENT_H */

