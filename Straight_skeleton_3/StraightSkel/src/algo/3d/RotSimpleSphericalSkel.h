/**
 * @file   algo/3d/RotSimpleSphericalSkel.h
 * @author Gernot Walzl
 * @date   2012-12-28
 */

#ifndef ALGO_3D_ROTSIMPLESPHERICALSKEL_H
#define ALGO_3D_ROTSIMPLESPHERICALSKEL_H

#include "algo/ptrs.h"
#include "algo/3d/ptrs.h"
#include "algo/3d/AbstractSimpleSphericalSkel.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <list>

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

class RotSimpleSphericalSkel : public AbstractSimpleSphericalSkel {
public:
    virtual ~RotSimpleSphericalSkel();

    static RotSimpleSphericalSkelSPtr create(SphericalPolygonSPtr polygon);
    static RotSimpleSphericalSkelSPtr create(SphericalPolygonSPtr polygon, ControllerSPtr controller);

    virtual void run();

    virtual bool isReflex(CircularVertexSPtr vertex);

    virtual bool init(SphericalPolygonSPtr polygon);

    void initSpeeds(SphericalPolygonSPtr polygon);
    void updateSpeeds(SphericalPolygonSPtr polygon_prev);

    static double approxOffsetTo(CircularVertexSPtr vertex, Point3SPtr point);

    SphericalEdgeEventSPtr nextEdgeEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalSplitEventSPtr nextSplitEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalTriangleEventSPtr nextTriangleEvent(SphericalPolygonSPtr polygon, double offset);

    SphericalAbstractEventSPtr nextEvent(SphericalPolygonSPtr polygon, double offset);

    CircularEdgeSPtr findLongestEdge(std::list<CircularEdgeSPtr> edges);

    double distanceOffset(CircularEdgeSPtr edge_begin, double offset, double speed_dst);
    double findMinDistance(CircularEdgeSPtr edge_begin, double offset);

    SphericalPolygonSPtr shiftEdges(SphericalPolygonSPtr polygon, double offset);
    SphericalPolygonSPtr shiftEdges2(SphericalPolygonSPtr polygon, double offset);

    void checkAngles(SphericalPolygonSPtr polygon);

    void handleEdgeEvent(SphericalEdgeEventSPtr event, SphericalPolygonSPtr polygon);
    void handleSplitEvent(SphericalSplitEventSPtr event, SphericalPolygonSPtr polygon);
    void handleTriangleEvent(SphericalTriangleEventSPtr event, SphericalPolygonSPtr polygon);

protected:
    RotSimpleSphericalSkel(SphericalPolygonSPtr polygon);
    RotSimpleSphericalSkel(SphericalPolygonSPtr polygon, ControllerSPtr controller);
};

} }

#endif /* ALGO_3D_ROTSIMPLESPHERICALSKEL_H */
