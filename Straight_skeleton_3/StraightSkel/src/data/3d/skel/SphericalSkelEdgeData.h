/**
 * @file   data/3d/skel/SphericalSkelEdgeData.h
 * @author Gernot Walzl
 * @date   2012-11-30
 */

#ifndef DATA_3D_SKEL_SPHERICALSKELEDGEDATA_H
#define DATA_3D_SKEL_SPHERICALSKELEDGEDATA_H

#include "data/3d/ptrs.h"
#include "data/3d/CircularEdgeData.h"
#include "data/3d/skel/ptrs.h"

namespace data { namespace _3d { namespace skel {

class SphericalSkelEdgeData : public CircularEdgeData {
public:
    virtual ~SphericalSkelEdgeData();

    static SphericalSkelEdgeDataSPtr create(CircularEdgeSPtr edge);

    CircularEdgeSPtr getOffsetEdge() const;
    void setOffsetEdge(CircularEdgeSPtr offset_edge);

    double getSpeed() const;
    void setSpeed(double speed);

    FacetSPtr getFacetOrigin() const;
    void setFacetOrigin(FacetSPtr facet_origin);

    Line3SPtr getRotationAxis() const;
    void setRotationAxis(Line3SPtr rotation_axis);

protected:
    SphericalSkelEdgeData();
    CircularEdgeWPtr offset_edge_;
    double speed_;
    FacetWPtr facet_origin_;
    Line3SPtr rotation_axis_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALSKELEDGEDATA_H */

