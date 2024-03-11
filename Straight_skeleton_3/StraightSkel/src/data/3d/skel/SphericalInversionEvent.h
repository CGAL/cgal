/**
 * @file   data/3d/skel/SphericalInversionEvent.h
 * @author Gernot Walzl
 * @date   2013-03-13
 */

#ifndef DATA_3D_SKEL_SPHERICALINVERSIONEVENT_H
#define DATA_3D_SKEL_SPHERICALINVERSIONEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalAbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SphericalInversionEvent : public SphericalAbstractEvent {
public:
    virtual ~SphericalInversionEvent();
    static SphericalInversionEventSPtr create();
    double getOffset() const;
    void setOffset(double offset);
    SphericalPolygonSPtr getPolygon() const;
    void setPolygon(SphericalPolygonSPtr polygon);
protected:
    SphericalInversionEvent();
    double offset_;
    SphericalPolygonSPtr polygon_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALINVERSIONEVENT_H */

