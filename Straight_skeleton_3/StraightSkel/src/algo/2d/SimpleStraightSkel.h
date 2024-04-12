// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

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
    static std::list<EdgeEventSPtr> nextEdgeEvent(PolygonSPtr polygon, CGAL::FT offset);
    static Point2SPtr crashAt(VertexSPtr vertex, EdgeSPtr edge);
    static std::list<SplitEventSPtr> nextSplitEvent(PolygonSPtr polygon, CGAL::FT offset);
    static std::list<AbstractEventSPtr> nextEvent(PolygonSPtr polygon, CGAL::FT offset);

    static PolygonSPtr shiftEdges(PolygonSPtr polygon, CGAL::FT offset);

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
