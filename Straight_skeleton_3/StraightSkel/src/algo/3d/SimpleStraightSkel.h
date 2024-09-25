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
 * @file   algo/3d/SimpleStraightSkel.h
 * @author Gernot Walzl
 * @date   2012-03-08
 */

#ifndef ALGO_3D_SIMPLESTRAIGHTSKEL_H
#define ALGO_3D_SIMPLESTRAIGHTSKEL_H

#include "typedefs_thread.h"

#include "algo/ptrs.h"
#include "algo/3d/ptrs.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"
#include "data/3d/skel/ptrs.h"

#include <list>
#include <optional>
#include <queue>
#include <utility>

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

class SimpleStraightSkel {
public:
    using PQ = std::priority_queue<AbstractEventSPtr,
                                   std::vector<AbstractEventSPtr>,
                                   AbstractEventSPtrCompare>;

    virtual ~SimpleStraightSkel();

    static SimpleStraightSkelSPtr create(PolyhedronSPtr polyhedron);
    static SimpleStraightSkelSPtr create(PolyhedronSPtr polyhedron, ControllerSPtr controller);
    static SimpleStraightSkelSPtr create(PolyhedronSPtr polyhedron, ControllerSPtr controller, const std::list<CGAL::FT>& save_offsets);

    void initVertexSplitter();
    void initEdgeEvent();

    static bool isReflex(EdgeSPtr edge);
    static bool isReflex(VertexSPtr vertex);
    static bool isConvex(VertexSPtr vertex);

    static Line3SPtr line(EdgeSPtr edge);

    bool run();
    ThreadSPtr startThread();

    /**
     * All vertices of the input polyhedron have to have a maximum degree of 3.
     * The splitting is done by using the speed during the shrinking process.
     * This works for convex vertices only.
     */
    static void splitConvexVertex(VertexSPtr vertex);

    static void splitReflexVertex(VertexSPtr vertex);

    /**
     * Creates a new node for the vertex data.
     * Used by init(...) only.
     */
    static NodeSPtr createNode(VertexSPtr vertex);

    /**
     * Creates a new arc for the vertex data.
     * The node of the vertex data has to be set before.
     */
    static ArcSPtr createArc(VertexSPtr vertex);

    /**
     * Creates a new sheet for the edge data.
     * The arcs of the vertices have to be set before.
     */
    static SheetSPtr createSheet(EdgeSPtr edge);

    /**
     * Split all vertices with degree > 3 and
     * initializes the data variables of all edges and vertices.
     */
    bool init(PolyhedronSPtr polyhedron);

    /**
     * Checks if the given edge is part of a triangle
     * on the surface of the polyhedron.
     * Used by: nextEdgeEvent, nextTriangleEvent
     */
    static bool isTriangle(FacetSPtr facet, EdgeSPtr edge_begin);

    /**
     * Checks if the given edge is part of a tetrahedron.
     * Used by: nextTriangleEvent, nextTetrahedronEvent
     */
    static bool isTetrahedron(EdgeSPtr edge_begin);

    /**
     * Checks if the point `point` is on the correct side of the bisector arc
     * emaneting from the vertex common to `edge` and the edge
     * shared by `edge->next(f)` and `f_third`
     */
    static bool check_bisector(EdgeSPtr edge,
                               FacetSPtr f,
                               CGAL::FT t,
                               FacetSPtr f_third,
                               Point3SPtr point);

    static std::pair<Point3SPtr, CGAL::FT> vanishesAtGeneric(FacetSPtr facet_0,
                                                             FacetSPtr facet_1,
                                                             FacetSPtr facet_2,
                                                             FacetSPtr facet_3);
    static std::pair<Point3SPtr, CGAL::FT> vanishesAtOnePairOpposite(FacetSPtr facet_0,
                                                                     FacetSPtr facet_1,
                                                                     FacetSPtr facet_2,
                                                                     FacetSPtr facet_3);
    static std::pair<Point3SPtr, CGAL::FT> vanishesAtOnePairContiguous(FacetSPtr facet_0,
                                                                       FacetSPtr facet_1,
                                                                       FacetSPtr facet_2,
                                                                       FacetSPtr facet_3);
    static std::pair<Point3SPtr, CGAL::FT> vanishesAtTwoPairs(FacetSPtr facet_0,
                                                              FacetSPtr facet_1,
                                                              FacetSPtr facet_2,
                                                              FacetSPtr facet_3);

    /**
     * Returns the point where the edge will vanish.
     */
    // static Point3SPtr vanishesAt(EdgeSPtr edge);
    static std::pair<Point3SPtr, CGAL::FT> vanishesAt(EdgeSPtr edge);

    /**
     * Returns the point where 2 edges will crash into each other.
     *
     * If `offset_max` is passed, ignore the crash if it happens in the future.
     */
    static std::pair<Point3SPtr, CGAL::FT> crashAt(EdgeSPtr edge_1, EdgeSPtr edge_2,
                                                   const std::optional<CGAL::FT> offset_max = std::nullopt);

    /**
     * Returns the offset (time) when the facet will reach the given point.
     */
    static CGAL::FT offsetDist(FacetSPtr facet, Point3SPtr point);

    /**
     * Edge flip event.
     */
    static void collectEdgeEvents(PolyhedronSPtr polyhedron,
                                  const CGAL::FT current_offset,
                                  CGAL::FT& curr_time_to_next_event,
                                  PQ& queue);

    static void collectEdgeMergeEvents(PolyhedronSPtr polyhedron,
                                       const CGAL::FT current_offset,
                                       CGAL::FT& curr_time_to_next_event,
                                       PQ& queue);

    /**
     * The triangle on the surface vanishes.
     */
    static void collectTriangleEvents(PolyhedronSPtr polyhedron,
                                      const CGAL::FT current_offset,
                                      CGAL::FT& curr_time_to_next_event,
                                      PQ& queue);

    static void collectDblEdgeMergeEvents(PolyhedronSPtr polyhedron,
                                          const CGAL::FT current_offset,
                                          CGAL::FT& curr_time_to_next_event,
                                          PQ& queue);

    static void collectDblTriangleEvents(PolyhedronSPtr polyhedron,
                                         const CGAL::FT current_offset,
                                         CGAL::FT& curr_time_to_next_event,
                                         PQ& queue);

    /**
     * A tetrahedron causes one final event only.
     */
    static void collectTetrahedronEvents(PolyhedronSPtr polyhedron,
                                         const CGAL::FT current_offset,
                                         CGAL::FT& curr_time_to_next_event,
                                         PQ& queue);

    /**
     * Two vertices crash into each other.
     */
    static void collectVertexEvents(PolyhedronSPtr polyhedron,
                                    const CGAL::FT current_offset,
                                    CGAL::FT& curr_time_to_next_event,
                                    PQ& queue);

    static void collectFlipVertexEvents(PolyhedronSPtr polyhedron,
                                        const CGAL::FT current_offset,
                                        CGAL::FT& curr_time_to_next_event,
                                        PQ& queue);

    /**
     * Split event on the surface.
     * Edges do not need to be reflex.
     */
    static void collectSurfaceEvents(PolyhedronSPtr polyhedron,
                                     const CGAL::FT current_offset,
                                     CGAL::FT& curr_time_to_next_event,
                                     PQ& queue);

    /**
     * This event occurs when two edges collide.
     * The first edge is always reflex.
     */
    static void collectPolyhedronSplitEvents(PolyhedronSPtr polyhedron,
                                             const CGAL::FT current_offset,
                                             CGAL::FT& curr_time_to_next_event,
                                             PQ& queue);

    static void collectSplitMergeEvents(PolyhedronSPtr polyhedron,
                                        const CGAL::FT current_offset,
                                        CGAL::FT& curr_time_to_next_event,
                                        PQ& queue);

    /**
     * This event occurs when two edges collide.
     * The first edge is always reflex.
     */
    static void collectEdgeSplitEvents(PolyhedronSPtr polyhedron,
                                       const CGAL::FT current_offset,
                                       CGAL::FT& curr_time_to_next_event,
                                       PQ& queue);

    /**
     * A reflex vertex reaches a facet.
     */
    static void collectPierceEvents(PolyhedronSPtr polyhedron,
                                    const CGAL::FT current_offset,
                                    CGAL::FT& curr_time_to_next_event,
                                    PQ& queue);

    /**
     * Collect all events for the polyhedron
     */
    void collectEvents(PolyhedronSPtr polyhedron,
                       const CGAL::FT current_offset,
                       PQ& queue);

    /**
     * Determines the next event.
     */
    std::pair<AbstractEventSPtr, bool> nextEvent(PQ& queue);

    /**
     * Normalize facet planes
    */
    static void normalizeFacetPlanes(PolyhedronSPtr polyhedron);

    /**
     * Normalize facet planes, ensuring parallel facets receive the same plane coefficients.
    */
    static void harmonizeFacetPlanes(PolyhedronSPtr polyhedron);

    PolyhedronSPtr enablePerturbedMode(PolyhedronSPtr polyhedron, const CGAL::FT simultaneousTime);

    PolyhedronSPtr disablePerturbedMode(PolyhedronSPtr polyhedron,
                                        CGAL::FT currentOffset,
                                        CGAL::FT nextEventOffset);

    /**
     * Appends a node of an event to the skeleton.
     * It links all adjacent arcs and sheets to this node.
     */
    void appendEventNode(NodeSPtr node);

    void handleEdgeEvent(EdgeEventSPtr event, PolyhedronSPtr polyhedron);
    void handleEdgeMergeEvent(EdgeMergeEventSPtr event, PolyhedronSPtr polyhedron);
    void handleTriangleEvent(TriangleEventSPtr event, PolyhedronSPtr polyhedron);
    void handleDblEdgeMergeEvent(DblEdgeMergeEventSPtr event, PolyhedronSPtr polyhedron);
    void handleDblTriangleEvent(DblTriangleEventSPtr event, PolyhedronSPtr polyhedron);
    void handleTetrahedronEvent(TetrahedronEventSPtr event, PolyhedronSPtr polyhedron);
    void handleVertexEvent(VertexEventSPtr event, PolyhedronSPtr polyhedron);
    void handleFlipVertexEvent(FlipVertexEventSPtr event, PolyhedronSPtr polyhedron);
    void handleSurfaceEvent(SurfaceEventSPtr event, PolyhedronSPtr polyhedron);
    void handlePolyhedronSplitEvent(PolyhedronSplitEventSPtr event, PolyhedronSPtr polyhedron);
    void handleSplitMergeEvent(SplitMergeEventSPtr event, PolyhedronSPtr polyhedron);
    void handleEdgeSplitEvent(EdgeSplitEventSPtr event, PolyhedronSPtr polyhedron);
    void handlePierceEvent(PierceEventSPtr event, PolyhedronSPtr polyhedron);

    StraightSkeletonSPtr getResult() const;

protected:
    SimpleStraightSkel(PolyhedronSPtr polyhedron);
    SimpleStraightSkel(PolyhedronSPtr polyhedron, ControllerSPtr controller);
    SimpleStraightSkel(PolyhedronSPtr polyhedron, ControllerSPtr controller, const std::list<CGAL::FT>& save_offsets);

    static FacetSPtr getFacetSrc(EdgeSPtr edge);
    static FacetSPtr getFacetDst(EdgeSPtr edge);

    PolyhedronSPtr polyhedron_;
    ControllerSPtr controller_;
    std::list<CGAL::FT> save_offsets_;
    bool use_fast_vertex_splitter_;
    AbstractVertexSplitterSPtr vertex_splitter_;
    int edge_event_;
    StraightSkeletonSPtr skel_result_;

    bool usingPerturbedMode_ = false;
    CGAL::FT simultaneousOffset_ = 0;
};

} }

#endif /* ALGO_3D_SIMPLESTRAIGHTSKEL_H */
