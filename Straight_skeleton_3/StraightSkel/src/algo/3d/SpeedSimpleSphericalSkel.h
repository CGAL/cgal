/**
 * @file   algo/3d/SpeedSimpleSphericalSkel.h
 * @author Gernot Walzl
 * @date   2013-09-04
 */

#ifndef ALGO_3D_SPEEDSIMPLESPHERICALSKEL_H
#define ALGO_3D_SPEEDSIMPLESPHERICALSKEL_H

#include "typedefs_thread.h"
#include "algo/ptrs.h"
#include "algo/3d/ptrs.h"
#include "algo/3d/AbstractSimpleSphericalSkel.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

class SpeedSimpleSphericalSkel : public AbstractSimpleSphericalSkel {
public:
    virtual ~SpeedSimpleSphericalSkel();

    static SpeedSimpleSphericalSkelSPtr create(SphericalPolygonSPtr polygon);
    static SpeedSimpleSphericalSkelSPtr create(SphericalPolygonSPtr polygon, ControllerSPtr controller);

    virtual void run();

    virtual bool isReflex(CircularVertexSPtr vertex);

    virtual bool init(SphericalPolygonSPtr polygon);

    void initSpeeds(SphericalPolygonSPtr polygon);

    double speed(CircularVertexSPtr vertex);

    static Point3SPtr vanishesAt(CircularEdgeSPtr edge);
    static Point3SPtr crashAt(CircularVertexSPtr vertex, CircularEdgeSPtr edge);

    double offsetTo(CircularVertexSPtr vertex, Point3SPtr point);

    SphericalEdgeEventSPtr nextEdgeEvent(SphericalPolygonSPtr polygon);
    SphericalSplitEventSPtr nextSplitEvent(SphericalPolygonSPtr polygon);
    SphericalTriangleEventSPtr nextTriangleEvent(SphericalPolygonSPtr polygon);

    SphericalAbstractEventSPtr nextEvent(SphericalPolygonSPtr polygon);

    SphericalPolygonSPtr copyPolygon(SphericalPolygonSPtr polygon);

    void handleEdgeEvent(SphericalEdgeEventSPtr event, SphericalPolygonSPtr polygon);
    void handleSplitEvent(SphericalSplitEventSPtr event, SphericalPolygonSPtr polygon);
    void handleTriangleEvent(SphericalTriangleEventSPtr event, SphericalPolygonSPtr polygon);

protected:
    SpeedSimpleSphericalSkel(SphericalPolygonSPtr polygon);
    SpeedSimpleSphericalSkel(SphericalPolygonSPtr polygon, ControllerSPtr controller);
};

} }

#endif /* ALGO_3D_SPEEDSIMPLESPHERICALSKEL_H */
