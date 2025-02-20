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

#include <boost/functional/hash.hpp>

#include <filesystem>
#include <list>
#include <optional>
#include <queue>
#include <utility>

namespace std {
template <>
struct hash<std::array<std::size_t, 4>> {
    size_t operator()(const std::array<std::size_t, 4>& arr) const {
        size_t seed = 0;
        for (std::size_t value : arr) {
            boost::hash_combine(seed, value); // Use boost::hash_combine
        }
        return seed;
    }
};
}

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
    static SimpleStraightSkelSPtr create(PolyhedronSPtr polyhedron, ControllerSPtr controller,
                                         const std::list<CGAL::FT>& save_offsets,
                                         const std::filesystem::path& save_path);

    void initVertexSplitter();
    void initEdgeEvent();

    static bool isReflex(EdgeSPtr edge);
    static bool isReflex(VertexSPtr vertex);
    static bool isConvex(VertexSPtr vertex);

    static Line3SPtr line(EdgeSPtr edge);

    /**
     * Save offset, attempt to un-tilt if there is a tilt
     */
    bool savePolyhedron(PolyhedronSPtr polyhedron,
                        CGAL::FT current_offset,
                        const bool do_triangulate = true,
                        const bool dump_exact = true,
                        const bool attempt_untilting = false);

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

    /**
     * @todo hamornize time and offset nomenclature
     * Return the intersection point and time of the 4 shifting planes.
     * Query the cache first, and fill it if it was not.
     */
    std::pair<Point3SPtr, CGAL::FT> intersectionPointAndTimeOffsetPlanes(FacetSPtr facet_0,
                                                                         FacetSPtr facet_1,
                                                                         FacetSPtr facet_2,
                                                                         FacetSPtr facet_3,
                                                                         const CGAL::FT& past_bound,
                                                                         const CGAL::FT& future_bound);

    /**
     * Returns the point where the edge will vanish.
     */
    // static Point3SPtr vanishesAt(EdgeSPtr edge);
    std::pair<Point3SPtr, CGAL::FT> vanishesAt(EdgeSPtr edge,
                                               const CGAL::FT offset_past_bound,
                                               const CGAL::FT offset_future_bound);

    /**
     * Returns the point where 2 edges will crash into each other.
     *
     * If `offset_max` is passed, ignore the crash if it happens in the future.
     */
    std::pair<Point3SPtr, CGAL::FT> crashAt(EdgeSPtr edge_1, EdgeSPtr edge_2,
                                            const CGAL::FT offset_past_bound,
                                            const CGAL::FT offset_future_bound);
    /**
     * Returns the offset (time) when the facet will reach the given point.
     */
    static CGAL::FT offsetDist(FacetSPtr facet, Point3SPtr point);

    /**
     * Collects the facets that are affected by the event.
     */
    static std::list<FacetSPtr> eventFacets(AbstractEventSPtr event);

    /**
     * Returns `true` if the event has facets that have a "step_id" that is
     * greater than the event's "step_id".
     */
    static bool isEventPotentiallyObsolete(AbstractEventSPtr event);

    /**
     * Vanish events.
     */
    void collectVanishEvents(PolyhedronSPtr polyhedron,
                             const CGAL::FT current_offset,
                             CGAL::FT& offset_of_nearest_event,
                             PQ& queue);

    /**
     * Edge flip event.
     */
    void collectEdgeEvents(const std::list<EdgeSPtr>& edges,
                           PolyhedronSPtr polyhedron,
                           const CGAL::FT current_offset,
                           CGAL::FT& offset_of_nearest_event,
                           PQ& queue);
     void collectEdgeEvents(PolyhedronSPtr polyhedron,
                           const CGAL::FT current_offset,
                           CGAL::FT& offset_of_nearest_event,
                           PQ& queue);

    /**
     * Edge Merge event.
     */
    void collectEdgeMergeEvents(const std::list<EdgeSPtr>& edges,
                                PolyhedronSPtr polyhedron,
                                const CGAL::FT current_offset,
                                CGAL::FT& offset_of_nearest_event,
                                PQ& queue);
    void collectEdgeMergeEvents(PolyhedronSPtr polyhedron,
                                const CGAL::FT current_offset,
                                CGAL::FT& offset_of_nearest_event,
                                PQ& queue);

    /**
     * The triangle on the surface vanishes.
     */
    void collectTriangleEvents(const std::list<EdgeSPtr>& edges,
                               PolyhedronSPtr polyhedron,
                               const CGAL::FT current_offset,
                               CGAL::FT& offset_of_nearest_event,
                               PQ& queue);
     void collectTriangleEvents(PolyhedronSPtr polyhedron,
                               const CGAL::FT current_offset,
                               CGAL::FT& offset_of_nearest_event,
                               PQ& queue);

    /**
     * Dbl Edge Merge.
     */
     void collectDblEdgeMergeEvents(const std::list<EdgeSPtr>& edges,
                                    PolyhedronSPtr polyhedron,
                                    const CGAL::FT current_offset,
                                    CGAL::FT& offset_of_nearest_event,
                                    PQ& queue);
     void collectDblEdgeMergeEvents(PolyhedronSPtr polyhedron,
                                   const CGAL::FT current_offset,
                                   CGAL::FT& offset_of_nearest_event,
                                   PQ& queue);

    /**
     * Dbl Edge Merge.
     */
    void collectDblTriangleEvents(const std::list<EdgeSPtr>& edges,
                                  PolyhedronSPtr polyhedron,
                                  const CGAL::FT current_offset,
                                  CGAL::FT& offset_of_nearest_event,
                                  PQ& queue);
     void collectDblTriangleEvents(PolyhedronSPtr polyhedron,
                                  const CGAL::FT current_offset,
                                  CGAL::FT& offset_of_nearest_event,
                                  PQ& queue);

    /**
     * A tetrahedron causes one final event only.
     */
    void collectTetrahedronEvents(const std::list<EdgeSPtr>& edges,
                                  PolyhedronSPtr polyhedron,
                                  const CGAL::FT current_offset,
                                  CGAL::FT& offset_of_nearest_event,
                                  PQ& queue);
     void collectTetrahedronEvents(PolyhedronSPtr polyhedron,
                                  const CGAL::FT current_offset,
                                  CGAL::FT& offset_of_nearest_event,
                                  PQ& queue);

    /**
     * Two vertices crash into each other.
     */
    void collectVertexEvents(const std::list<VertexSPtr>& vertices,
                             PolyhedronSPtr polyhedron,
                             const CGAL::FT current_offset,
                             CGAL::FT& offset_of_nearest_event,
                             PQ& queue);
    void collectVertexEvents(PolyhedronSPtr polyhedron,
                             const CGAL::FT current_offset,
                             CGAL::FT& offset_of_nearest_event,
                             PQ& queue);

    /**
     * Flip vertex event
     */
    void collectFlipVertexEvents(const std::list<VertexSPtr>& vertices,
                                 PolyhedronSPtr polyhedron,
                                 const CGAL::FT current_offset,
                                 CGAL::FT& offset_of_nearest_event,
                                 PQ& queue);
    void collectFlipVertexEvents(PolyhedronSPtr polyhedron,
                                 const CGAL::FT current_offset,
                                 CGAL::FT& offset_of_nearest_event,
                                 PQ& queue);

    /**
     * Split event on the surface.
     * Edges do not need to be reflex.
     */
    void collectSurfaceEvents(const std::list<EdgeSPtr>& edges,
                              PolyhedronSPtr polyhedron,
                              const CGAL::FT current_offset,
                              CGAL::FT& offset_of_nearest_event,
                              PQ& queue);
     void collectSurfaceEvents(PolyhedronSPtr polyhedron,
                              const CGAL::FT current_offset,
                              CGAL::FT& offset_of_nearest_event,
                              PQ& queue);

    /**
     * This event occurs when two edges collide.
     * The first edge is always reflex.
     */
    void collectPolyhedronSplitEvents(const std::list<EdgeSPtr>& edges,
                                      PolyhedronSPtr polyhedron,
                                      const CGAL::FT current_offset,
                                      CGAL::FT& offset_of_nearest_event,
                                      PQ& queue);
     void collectPolyhedronSplitEvents(PolyhedronSPtr polyhedron,
                                      const CGAL::FT current_offset,
                                      CGAL::FT& offset_of_nearest_event,
                                      PQ& queue);

    /**
     * Split Merge event
     */
     void collectSplitMergeEvents(const std::list<VertexSPtr>& vertices,
                                  PolyhedronSPtr polyhedron,
                                  const CGAL::FT current_offset,
                                  CGAL::FT& offset_of_nearest_event,
                                  PQ& queue);
     void collectSplitMergeEvents(PolyhedronSPtr polyhedron,
                                 const CGAL::FT current_offset,
                                 CGAL::FT& offset_of_nearest_event,
                                 PQ& queue);

    /**
     * This event occurs when two edges collide.
     * The first edge is always reflex.
     */
    void collectEdgeSplitEvents(const std::list<EdgeSPtr>& edges_1,
                                const std::list<EdgeSPtr>& edges_2,
                                PolyhedronSPtr polyhedron,
                                const CGAL::FT current_offset,
                                CGAL::FT& offset_of_nearest_event,
                                PQ& queue);
    void collectEdgeSplitEvents(PolyhedronSPtr polyhedron,
                                const CGAL::FT current_offset,
                                CGAL::FT& offset_of_nearest_event,
                                PQ& queue);

    /**
     * A reflex vertex reaches a facet.
     */
    void collectPierceEvents(const std::list<VertexSPtr>& vertices,
                             const std::list<FacetSPtr>& facets,
                             PolyhedronSPtr polyhedron,
                             const CGAL::FT current_offset,
                             CGAL::FT& offset_of_nearest_event,
                             PQ& queue);
    void collectPierceEvents(PolyhedronSPtr polyhedron,
                             const CGAL::FT current_offset,
                             CGAL::FT& offset_of_nearest_event,
                             PQ& queue);

    void printQueue(const PQ& queue);

    /**
     * Checks if the queue updated with local events is displaying the same events
     * as if it had been built from scratch.
     */
    bool checkQueueCorrectness(const PQ& queue,
                               PolyhedronSPtr polyhedron,
                               const CGAL::FT current_offset);

    /**
     * Collect all events for the polyhedron
     */
    void collectLocalEvents(const std::set<FacetSPtr>& facets,
                            PolyhedronSPtr polyhedron,
                            const CGAL::FT current_offset,
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
    AbstractEventSPtr nextEvent(PQ& queue, CGAL::FT offset_current);

    /**
     * Appends a node of an event to the skeleton.
     * It links all adjacent arcs and sheets to this node.
     */
    void appendEventNode(NodeSPtr node);

    static PolyhedronSPtr soup_to_polyhedron(const std::vector<Point3>& points,
                                             const std::vector<std::vector<std::size_t> >& triangles,
                                             const std::vector<Plane3SPtr>& planes,
                                             const std::vector<CGAL::FT>& speeds);

    void tagFacetsWithStepID(const std::set<FacetSPtr>& facets);

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
    SimpleStraightSkel(PolyhedronSPtr polyhedron, ControllerSPtr controller,
                       const std::list<CGAL::FT>& save_offsets,
                       const std::filesystem::path& save_path);

    PolyhedronSPtr polyhedron_;

    ControllerSPtr controller_;

    std::list<CGAL::FT> save_offsets_;
    std::filesystem::path save_path_;
    bool use_fast_vertex_splitter_;
    AbstractVertexSplitterSPtr vertex_splitter_;
    int edge_event_;

    StraightSkeletonSPtr skel_result_;


    std::vector<Plane3SPtr> basePlanes_;

    int step_id_;
    std::set<FacetSPtr> post_op_facets_; // @tmp here just to simplify my life

#ifndef CGAL_SS3_NO_CACHING
    std::unordered_map<std::array<std::size_t, 4>,
                       std::pair<Point3SPtr, CGAL::FT> > intersectionCache_;
#endif
};

} }

#endif /* ALGO_3D_SIMPLESTRAIGHTSKEL_H */
