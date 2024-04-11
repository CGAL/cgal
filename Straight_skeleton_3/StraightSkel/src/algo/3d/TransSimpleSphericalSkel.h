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

    Point3SPtr offsetPoint(CircularVertexSPtr vertex, CGAL::FT offset);

    CGAL::FT offsetTo(CircularEdgeSPtr edge, Point3SPtr point);

    double angleOf(CircularVertexSPtr vertex);

    bool isEdgeEvent(CircularEdgeSPtr edge, CGAL::FT offset);

    SphericalEdgeEventSPtr nextEdgeEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);
    SphericalSplitEventSPtr nextSplitEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);
    SphericalTriangleEventSPtr nextTriangleEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);
    SphericalDblEdgeEventSPtr nextDblEdgeEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);
    SphericalLeaveEventSPtr nextLeaveEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);
    SphericalReturnEventSPtr nextReturnEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);
    SphericalDblLeaveEventSPtr nextDblLeaveEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);
    SphericalDblReturnEventSPtr nextDblReturnEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);
    SphericalVertexEventSPtr nextVertexEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);
    SphericalEdgeMergeEventSPtr nextEdgeMergeEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);
    SphericalInversionEventSPtr nextInversionEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);

    SphericalAbstractEventSPtr nextEvent(SphericalPolygonSPtr polygon, CGAL::FT offset);

    SphericalPolygonSPtr shiftEdges(SphericalPolygonSPtr polygon, CGAL::FT offset);

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
