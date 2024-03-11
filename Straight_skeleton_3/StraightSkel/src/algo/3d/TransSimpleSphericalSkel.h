/**
 * @file   algo/3d/TransSimpleSphericalSkel.h
 * @author Gernot Walzl
 * @date   2013-01-11
 */

#ifndef ALGO_3D_TRANSSIMPLESPHERICALSKEL_H
#define ALGO_3D_TRANSSIMPLESPHERICALSKEL_H

#include "algo/ptrs.h"
#include "algo/3d/ptrs.h"
#include "algo/3d/AbstractSimpleSphericalSkel.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

class TransSimpleSphericalSkel : public AbstractSimpleSphericalSkel {
public:
    virtual ~TransSimpleSphericalSkel();

    static TransSimpleSphericalSkelSPtr create(SphericalPolygonSPtr polygon);
    static TransSimpleSphericalSkelSPtr create(SphericalPolygonSPtr polygon, ControllerSPtr controller);

    virtual void run();

    virtual bool isReflex(CircularVertexSPtr vertex);

    virtual CircularArcSPtr createArc(CircularVertexSPtr vertex);

    virtual bool init(SphericalPolygonSPtr polygon);

    Point3SPtr offsetPoint(CircularVertexSPtr vertex, double offset);

    double offsetTo(CircularEdgeSPtr edge, Point3SPtr point);

    double angleOf(CircularVertexSPtr vertex);

    bool isEdgeEvent(CircularEdgeSPtr edge, double offset);

    SphericalEdgeEventSPtr nextEdgeEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalSplitEventSPtr nextSplitEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalTriangleEventSPtr nextTriangleEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalDblEdgeEventSPtr nextDblEdgeEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalLeaveEventSPtr nextLeaveEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalReturnEventSPtr nextReturnEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalDblLeaveEventSPtr nextDblLeaveEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalDblReturnEventSPtr nextDblReturnEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalVertexEventSPtr nextVertexEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalEdgeMergeEventSPtr nextEdgeMergeEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalInversionEventSPtr nextInversionEvent(SphericalPolygonSPtr polygon, double offset);

    SphericalAbstractEventSPtr nextEvent(SphericalPolygonSPtr polygon, double offset);

    SphericalPolygonSPtr shiftEdges(SphericalPolygonSPtr polygon, double offset);

    void handleEdgeEvent(SphericalEdgeEventSPtr event, SphericalPolygonSPtr polygon);
    void handleSplitEvent(SphericalSplitEventSPtr event, SphericalPolygonSPtr polygon);
    void handleTriangleEvent(SphericalTriangleEventSPtr event, SphericalPolygonSPtr polygon);
    void handleDblEdgeEvent(SphericalDblEdgeEventSPtr event, SphericalPolygonSPtr polygon);
    void handleLeaveEvent(SphericalLeaveEventSPtr event, SphericalPolygonSPtr polygon);
    void handleReturnEvent(SphericalReturnEventSPtr event, SphericalPolygonSPtr polygon);
    void handleDblLeaveEvent(SphericalDblLeaveEventSPtr event, SphericalPolygonSPtr polygon);
    void handleDblReturnEvent(SphericalDblReturnEventSPtr event, SphericalPolygonSPtr polygon);
    void handleVertexEvent(SphericalVertexEventSPtr event, SphericalPolygonSPtr polygon);
    void handleEdgeMergeEvent(SphericalEdgeMergeEventSPtr event, SphericalPolygonSPtr polygon);
    void handleInversionEvent(SphericalInversionEventSPtr event, SphericalPolygonSPtr polygon);

protected:
    TransSimpleSphericalSkel(SphericalPolygonSPtr polygon);
    TransSimpleSphericalSkel(SphericalPolygonSPtr polygon, ControllerSPtr controller);
    double epsilon_;   // TODO: epsilon environment is not good
};

} }

#endif /* ALGO_3D_TRANSSIMPLESPHERICALSKEL_H */
