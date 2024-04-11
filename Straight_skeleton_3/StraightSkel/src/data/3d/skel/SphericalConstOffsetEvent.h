/**
 * @file   data/3d/skel/SphericalConstOffsetEvent.h
 * @author Gernot Walzl
 * @date   2012-11-29
 */

#ifndef DATA_3D_SKEL_SPHERICALCONSTOFFSETEVENT_H
#define DATA_3D_SKEL_SPHERICALCONSTOFFSETEVENT_H

#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalAbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SphericalConstOffsetEvent : public SphericalAbstractEvent {
public:
    virtual ~SphericalConstOffsetEvent();
    static SphericalConstOffsetEventSPtr create(CGAL::FT offset);

    CGAL::FT getOffset() const;
    void setOffset(CGAL::FT offset);

protected:
    SphericalConstOffsetEvent(CGAL::FT offset);
    CGAL::FT offset_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALCONSTOFFSETEVENT_H */
