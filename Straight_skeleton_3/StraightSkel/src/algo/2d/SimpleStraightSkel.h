/**
 * @file   algo/2d/SimpleStraightSkel.h
 * @author Gernot Walzl
 * @date   2012-02-06
 */

#ifndef ALGO_2D_SIMPLESTRAIGHTSKEL_H
#define ALGO_2D_SIMPLESTRAIGHTSKEL_H

#include "typedefs_thread.h"
#include "algo/ptrs.h"
#include "algo/2d/ptrs.h"
#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"
#include <list>

namespace algo { namespace _2d {

using namespace data::_2d;
using namespace data::_2d::skel;

class SimpleStraightSkel {
public:
    virtual ~SimpleStraightSkel();

    static SimpleStraightSkelSPtr create(PolygonSPtr polygon);
    static SimpleStraightSkelSPtr create(PolygonSPtr polygon, ControllerSPtr controller);

    void run();
    ThreadSPtr startThread();

    static NodeSPtr createNode(VertexSPtr vertex);
    static ArcSPtr createArc(VertexSPtr vertex);
    bool init(PolygonSPtr polygon);
    static std::list<EdgeEventSPtr> nextEdgeEvent(PolygonSPtr polygon, double offset);
    static Point2SPtr crashAt(VertexSPtr vertex, EdgeSPtr edge);
    static std::list<SplitEventSPtr> nextSplitEvent(PolygonSPtr polygon, double offset);
    static std::list<AbstractEventSPtr> nextEvent(PolygonSPtr polygon, double offset);

    static PolygonSPtr shiftEdges(PolygonSPtr polygon, double offset);

    void appendEventNode(NodeSPtr node);

    void handleEdgeEvent(EdgeEventSPtr event, PolygonSPtr polygon);
    void handleSplitEvent(SplitEventSPtr event, PolygonSPtr polygon);
    void handleTriangleEvent(TriangleEventSPtr event, PolygonSPtr polygon);

    StraightSkeletonSPtr getResult() const;

protected:
    SimpleStraightSkel(PolygonSPtr polygon);
    SimpleStraightSkel(PolygonSPtr polygon, ControllerSPtr controller);

    PolygonSPtr polygon_;
    ControllerSPtr controller_;
    StraightSkeletonSPtr skel_result_;
};

} }

#endif /* ALGO_2D_SIMPLESTRAIGHTSKEL_H */
