/**
 * @file   algo/2d/FastStraightSkel.h
 * @author Gernot Walzl
 * @date   2015-09-28
 */

#ifndef ALGO_2D_FASTSTRAIGHTSKEL_H
#define ALGO_2D_FASTSTRAIGHTSKEL_H

#include "typedefs_thread.h"
#include "algo/ptrs.h"
#include "algo/2d/ptrs.h"
#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"

namespace algo { namespace _2d {

using namespace data::_2d;
using namespace data::_2d::skel;

/**
 * This algorithm is just an idea.
 */
class FastStraightSkel {
public:
    virtual ~FastStraightSkel();

    static FastStraightSkelSPtr create(PolygonSPtr polygon);
    static FastStraightSkelSPtr create(PolygonSPtr polygon, ControllerSPtr controller);

    void run();
    ThreadSPtr startThread();

    static NodeSPtr createNode(VertexSPtr vertex);
    static ArcSPtr createArc(VertexSPtr vertex);
    bool init(PolygonSPtr polygon);

    static EdgeEventSPtr nextEvent(PolygonSPtr polygon, double offset, int direction);

    static PolygonSPtr shiftEdges(PolygonSPtr polygon, double offset);

    void appendEventNode(NodeSPtr node);
    void handleEvent(EdgeEventSPtr event, PolygonSPtr polygon);

    StraightSkeletonSPtr getResult() const;

protected:
    FastStraightSkel(PolygonSPtr polygon);
    FastStraightSkel(PolygonSPtr polygon, ControllerSPtr controller);

    PolygonSPtr polygon_;
    ControllerSPtr controller_;
    StraightSkeletonSPtr skel_result_;
};

} }

#endif /* ALGO_2D_FASTSTRAIGHTSKEL_H */
