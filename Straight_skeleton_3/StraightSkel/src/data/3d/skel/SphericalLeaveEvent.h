/**
 * @file   data/3d/skel/SphericalLeaveEvent.h
 * @author Gernot Walzl
 * @date   2013-01-28
 */

#ifndef DATA_3D_SKEL_SPHERICALLEAVEEVENT_H
#define DATA_3D_SKEL_SPHERICALLEAVEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalAbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SphericalLeaveEvent : public SphericalAbstractEvent {
public:
    virtual ~SphericalLeaveEvent();

    static SphericalLeaveEventSPtr create();

    CGAL::FT getOffset() const;
    void setOffset(CGAL::FT offset);
    CircularVertexSPtr getVertex() const;
    void setVertex(CircularVertexSPtr vertex);

    void setHighlight(bool highlight);

protected:
    SphericalLeaveEvent();
    CGAL::FT offset_;
    CircularVertexSPtr vertex_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALLEAVEEVENT_H */
