/**
 * @file   algo/3d/ProjSimpleSphericalSkel.h
 * @author Gernot Walzl
 * @date   2012-12-03
 */

#ifndef ALGO_3D_PROJSIMPLESPHERICALSKEL_H
#define ALGO_3D_PROJSIMPLESPHERICALSKEL_H

#include "algo/ptrs.h"
#include "algo/3d/ptrs.h"
#include "data/3d/ptrs.h"
#include "algo/3d/AbstractSimpleSphericalSkel.h"
#include "data/3d/skel/ptrs.h"

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

class ProjSimpleSphericalSkel : public AbstractSimpleSphericalSkel {
public:
    virtual ~ProjSimpleSphericalSkel();

    static ProjSimpleSphericalSkelSPtr create(SphericalPolygonSPtr polygon);
    static ProjSimpleSphericalSkelSPtr create(SphericalPolygonSPtr polygon, ControllerSPtr controller);

    virtual void run();

    virtual bool isReflex(CircularVertexSPtr vertex);

    /**
     * Initializes a rotation axis for each CircularEdge.
     * All rotation axes need to be coplanar. (plane of rotation axes)
     */
    void initRotationAxes(SphericalPolygonSPtr polygon);

    /**
     * Constant speed factors are determined by the angles between
     * the plane of rotation axes and adjacent faces of the polyhedron.
     */
    void initConstSpeeds(SphericalPolygonSPtr polygon);

    /**
     * Creates CircularNodes and CircularArcs for each CircularVertex.
     */
    virtual bool init(SphericalPolygonSPtr polygon);

    /**
     * CircularEdges that cross the plane of rotation axes have no valid rotation axis.
     */
    bool hasValidRotationAxis(CircularEdgeSPtr edge);

    /**
     * Determines the angle between the plane of rotation axes and
     * the supporting plane of the given CircularEdge.
     */
    double angleTo(CircularEdgeSPtr edge);

    /**
     * Determines the offset until the given CircularEdge reaches the point.
     */
    double offsetTo(CircularEdgeSPtr edge, Point3SPtr point);

    SphericalEdgeEventSPtr nextEdgeEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalSplitEventSPtr nextSplitEvent(SphericalPolygonSPtr polygon, double offset);
    SphericalTriangleEventSPtr nextTriangleEvent(SphericalPolygonSPtr polygon, double offset);

    SphericalAbstractEventSPtr nextEvent(SphericalPolygonSPtr polygon, double offset);

    /**
     * Negative offset shrinks the given polygon.
     */
    SphericalPolygonSPtr shiftEdges(SphericalPolygonSPtr polygon, double offset);

    void handleEdgeEvent(SphericalEdgeEventSPtr event, SphericalPolygonSPtr polygon);
    void handleSplitEvent(SphericalSplitEventSPtr event, SphericalPolygonSPtr polygon);
    void handleTriangleEvent(SphericalTriangleEventSPtr event, SphericalPolygonSPtr polygon);

protected:
    ProjSimpleSphericalSkel(SphericalPolygonSPtr polygon);
    ProjSimpleSphericalSkel(SphericalPolygonSPtr polygon, ControllerSPtr controller);
};

} }

#endif /* ALGO_3D_PROJSIMPLESPHERICALSKEL_H */
