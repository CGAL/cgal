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
    static SphericalConstOffsetEventSPtr create(double offset);

    double getOffset() const;
    void setOffset(double offset);

protected:
    SphericalConstOffsetEvent(double offset);
    double offset_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALCONSTOFFSETEVENT_H */
