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
#include <vector>

namespace std {

template <>
struct hash<std::array<int, 4>> {
    size_t operator()(const std::array<int, 4>& arr) const {
        size_t seed = 0;
        for (std::size_t value : arr) {
            boost::hash_combine(seed, value);
        }
        return seed;
    }
};

} // namespace std

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

namespace utils {

struct Base_mesh_offset_visitor
{
    virtual bool go_further(int, PolyhedronSPtr, CGAL::FT) = 0;
    virtual void before_offset_event(PolyhedronSPtr, CGAL::FT, AbstractEventSPtr) = 0;
    virtual void on_save_offset_event(PolyhedronSPtr, CGAL::FT) = 0;
    virtual void after_offset_event(PolyhedronSPtr, CGAL::FT) = 0;
};

struct Default_mesh_offset_visitor
    : public Base_mesh_offset_visitor
{
    bool go_further(int, PolyhedronSPtr, CGAL::FT) override { return true; }
    void before_offset_event(PolyhedronSPtr, CGAL::FT, AbstractEventSPtr) override { }
    void on_save_offset_event(PolyhedronSPtr, CGAL::FT) override { }
    void after_offset_event(PolyhedronSPtr, CGAL::FT) override { }
};

} // namespace utils

class SimpleStraightSkel {
public:
    using PQ = std::priority_queue<AbstractEventSPtr,
                                   std::vector<AbstractEventSPtr>,
                                   AbstractEventSPtrCompare>;

    virtual ~SimpleStraightSkel();

    static SimpleStraightSkelSPtr create(PolyhedronSPtr polyhedron);
    static SimpleStraightSkelSPtr create(PolyhedronSPtr polyhedron, ControllerSPtr controller);
    static SimpleStraightSkelSPtr create(PolyhedronSPtr polyhedron, ControllerSPtr controller,
                                         const std::vector<CGAL::FT>& save_offsets,
                                         const std::filesystem::path& save_path);

    void setVisitor(utils::Base_mesh_offset_visitor* visitor) {
        visitor_ = visitor;
    }

    void initVertexSplitter();
    void initEdgeEvent();

    static bool isReflex(EdgeSPtr edge, const bool future_facing = true);
    static bool isReflex(VertexSPtr vertex);
    static bool isConvex(VertexSPtr vertex);

    static Line3SPtr line(EdgeSPtr edge);

    static Point3SPtr getFinalPoint(VertexSPtr vertex, const CGAL::FT& offset_future_bound);
    static Plane3SPtr getFinalPlane(FacetSPtr facet, const CGAL::FT& offset_future_bound);

    /**
     * Save offset, attempt to un-tilt if there is a tilt
     */
    bool savePolyhedron(PolyhedronSPtr polyhedron,
                        const CGAL::FT& current_offset,
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
     * Store within each facet the coefficients of the plane at t=0
     */
    void cacheBasePlanes(PolyhedronSPtr polyhedron);

    /**
     * Split all vertices with degree > 3 and
     * initializes the data variables of all edges and vertices.
     */
    static bool init(PolyhedronSPtr polyhedron,
                      AbstractVertexSplitterSPtr vertex_splitter,
                      const bool use_fast_vertex_splitter = false,
                      ControllerSPtr controller = nullptr);

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
                               const CGAL::FT& t,
                               FacetSPtr f_third,
                               Point3SPtr point);

    static bool check_bisectors(EdgeSPtr edge1,
                                EdgeSPtr edge2,
                                Point3SPtr point,
                                const CGAL::FT& t);

    /**
     * @todo hamornize time and offset nomenclature
     * Return the intersection point and time of the 4 shifting planes.
     * Query the cache first, and fill it if it was not.
     */
    std::pair<Point3SPtr, CGAL::FT> intersectionPointAndTimeOffsetPlanes(FacetSPtr facet_0,
                                                                         FacetSPtr facet_1,
                                                                         FacetSPtr facet_2,
                                                                         FacetSPtr facet_3,
                                                                         const std::optional<CGAL::FT>& past_bound = std::nullopt,
                                                                         const std::optional<CGAL::FT>& future_bound = std::nullopt);

    /**
     * Returns the point where the edge will vanish.
     */
    // static Point3SPtr vanishesAt(EdgeSPtr edge);
    std::pair<Point3SPtr, CGAL::FT> vanishesAt(EdgeSPtr edge,
                                               const std::optional<CGAL::FT>& offset_past_bound = std::nullopt,
                                               const std::optional<CGAL::FT>& offset_future_bound = std::nullopt);

    /**
     * Returns the point where 2 edges will crash into each other.
     *
     * If `offset_max` is passed, ignore the crash if it happens in the future.
     */
    std::pair<Point3SPtr, CGAL::FT> crashAt(EdgeSPtr edge_1, EdgeSPtr edge_2,
                                            const std::optional<CGAL::FT>& offset_past_bound = std::nullopt,
                                            const std::optional<CGAL::FT>& offset_future_bound = std::nullopt);
    /**
     * Returns the offset (time) when the facet will reach the given point.
     */
    static CGAL::FT offsetDist(FacetSPtr facet, Point3SPtr point);

    /**
     * Returns `true` if the event is in the past
     */
    static bool isEventInThePast(const CGAL::FT& current_offset,
                                 AbstractEventSPtr event);

    /**
     * Returns `true` if the neighborhood of an event has changed.
     */
    static bool isEventObsolete(AbstractEventSPtr event);

    /**
     * Some combinatorial or geometric checks are very expensive to perform,
     * so delay them until the event is the best in the queue.
     * The gain is that the event can be invalidated combinatorially by the time
     * it gets popped.
     */
    static bool isActualEvent(const CGAL::FT& current_offset,
                              AbstractEventSPtr event,
                              PolyhedronSPtr polyhedron);

    static bool isActualVertexEvent(VertexEventSPtr event,
                                    PolyhedronSPtr polyhedron);
    static bool isActualFlipVertexEvent(FlipVertexEventSPtr event,
                                        PolyhedronSPtr polyhedron);
    static bool isActualSurfaceEvent(SurfaceEventSPtr event,
                                     PolyhedronSPtr polyhedron);
    static bool isActualPolyhedronSplitEvent(PolyhedronSplitEventSPtr event,
                                             const CGAL::FT& current_offset,
                                             PolyhedronSPtr polyhedron);
    static bool isActualSplitMergeEvent(SplitMergeEventSPtr event,
                                        PolyhedronSPtr polyhedron);
    static bool isActualPierceEvent(PierceEventSPtr event,
                                    const CGAL::FT& current_offset,
                                    PolyhedronSPtr polyhedron);

    /**
     * Vanish events.
     */
    void collectVanishEvents(const std::list<EdgeSPtr>& edges,
                             PolyhedronSPtr polyhedron,
                             const CGAL::FT& current_offset,
                             const std::optional<CGAL::FT>& offset_future_bound,
                             PQ& queue);
     void collectVanishEvents(PolyhedronSPtr polyhedron,
                             const CGAL::FT& current_offset,
                             const std::optional<CGAL::FT>& offset_future_bound,
                             PQ& queue);

    /**
     * Edge flip event.
     */
    void collectEdgeEvents(const std::list<EdgeSPtr>& edges,
                           PolyhedronSPtr polyhedron,
                           const CGAL::FT& current_offset,
                           const std::optional<CGAL::FT>& offset_future_bound,
                           PQ& queue);
     void collectEdgeEvents(PolyhedronSPtr polyhedron,
                           const CGAL::FT& current_offset,
                           const std::optional<CGAL::FT>& offset_future_bound,
                           PQ& queue);

    /**
     * Edge Merge event.
     */
    void collectEdgeMergeEvents(const std::list<EdgeSPtr>& edges,
                                PolyhedronSPtr polyhedron,
                                const bool use_canonical_event_reps,
                                const CGAL::FT& current_offset,
                                const std::optional<CGAL::FT>& offset_future_bound,
                                PQ& queue);
    void collectEdgeMergeEvents(PolyhedronSPtr polyhedron,
                                const CGAL::FT& current_offset,
                                const std::optional<CGAL::FT>& offset_future_bound,
                                PQ& queue);

    /**
     * The triangle on the surface vanishes.
     */
    void collectTriangleEvents(const std::list<EdgeSPtr>& edges,
                               PolyhedronSPtr polyhedron,
                               const bool use_canonical_event_reps,
                               const CGAL::FT& current_offset,
                               const std::optional<CGAL::FT>& offset_future_bound,
                               PQ& queue);
     void collectTriangleEvents(PolyhedronSPtr polyhedron,
                               const CGAL::FT& current_offset,
                               const std::optional<CGAL::FT>& offset_future_bound,
                               PQ& queue);

    /**
     * Dbl Edge Merge.
     */
     void collectDblEdgeMergeEvents(const std::list<EdgeSPtr>& edges,
                                    PolyhedronSPtr polyhedron,
                                    const bool use_canonical_event_reps,
                                    const CGAL::FT& current_offset,
                                    const std::optional<CGAL::FT>& offset_future_bound,
                                    PQ& queue);
     void collectDblEdgeMergeEvents(PolyhedronSPtr polyhedron,
                                   const CGAL::FT& current_offset,
                                   const std::optional<CGAL::FT>& offset_future_bound,
                                   PQ& queue);

    /**
     * Dbl Edge Merge.
     */
    void collectDblTriangleEvents(const std::list<EdgeSPtr>& edges,
                                  PolyhedronSPtr polyhedron,
                                  const bool use_canonical_event_reps,
                                  const CGAL::FT& current_offset,
                                  const std::optional<CGAL::FT>& offset_future_bound,
                                  PQ& queue);
     void collectDblTriangleEvents(PolyhedronSPtr polyhedron,
                                  const CGAL::FT& current_offset,
                                  const std::optional<CGAL::FT>& offset_future_bound,
                                  PQ& queue);

    /**
     * A tetrahedron causes one final event only.
     */
    void collectTetrahedronEvents(const std::list<EdgeSPtr>& edges,
                                  PolyhedronSPtr polyhedron,
                                  const bool use_canonical_event_reps,
                                  const CGAL::FT& current_offset,
                                  const std::optional<CGAL::FT>& offset_future_bound,
                                  PQ& queue);
     void collectTetrahedronEvents(PolyhedronSPtr polyhedron,
                                  const CGAL::FT& current_offset,
                                  const std::optional<CGAL::FT>& offset_future_bound,
                                  PQ& queue);

    /**
     * Two vertices crash into each other.
     */
    void collectVertexEvents(const std::list<VertexSPtr>& vertices,
                             PolyhedronSPtr polyhedron,
                             const bool use_canonical_event_reps,
                             const CGAL::FT& current_offset,
                             const std::optional<CGAL::FT>& offset_future_bound,
                             PQ& queue);
    void collectVertexEvents(PolyhedronSPtr polyhedron,
                             const CGAL::FT& current_offset,
                             const std::optional<CGAL::FT>& offset_future_bound,
                             PQ& queue);

    /**
     * Flip vertex event
     */
    void collectFlipVertexEvents(const std::list<VertexSPtr>& vertices,
                                 PolyhedronSPtr polyhedron,
                                 const bool use_canonical_event_reps,
                                 const CGAL::FT& current_offset,
                                 const std::optional<CGAL::FT>& offset_future_bound,
                                 PQ& queue);
    void collectFlipVertexEvents(PolyhedronSPtr polyhedron,
                                 const CGAL::FT& current_offset,
                                 const std::optional<CGAL::FT>& offset_future_bound,
                                 PQ& queue);

    /**
     * Split Merge event
     */
    void collectSplitMergeEvents(const std::list<VertexSPtr>& vertices,
                                 PolyhedronSPtr polyhedron,
                                 const bool use_canonical_event_reps,
                                 const CGAL::FT& current_offset,
                                 const std::optional<CGAL::FT>& offset_future_bound,
                                 PQ& queue);
    void collectSplitMergeEvents(PolyhedronSPtr polyhedron,
                                 const CGAL::FT& current_offset,
                                 const std::optional<CGAL::FT>& offset_future_bound,
                                 PQ& queue);

    /**
     * Split event on the surface.
     * Edges do not need to be reflex.
     */
    void collectSurfaceEvent(EdgeSPtr edge_1,
                             EdgeSPtr edge_2,
                             PolyhedronSPtr polyhedron,
                             const CGAL::FT& current_offset,
                             const std::optional<CGAL::FT>& offset_future_bound,
                             PQ& queue);
     void collectSurfaceEvents(const std::list<EdgeSPtr>& edges,
                              PolyhedronSPtr polyhedron,
                              const CGAL::FT& current_offset,
                              const std::optional<CGAL::FT>& offset_future_bound,
                              PQ& queue);
     void collectSurfaceEvents(PolyhedronSPtr polyhedron,
                              const CGAL::FT& current_offset,
                              const std::optional<CGAL::FT>& offset_future_bound,
                              PQ& queue);

    /**
     * This event occurs when two edges collide.
     * The first edge is always reflex.
     */
    void collectPolyhedronSplitEvent(EdgeSPtr edge_1,
                                     EdgeSPtr edge_2,
                                     PolyhedronSPtr polyhedron,
                                     const CGAL::FT& current_offset,
                                     const std::optional<CGAL::FT>& offset_future_bound,
                                     PQ& queue);
    void collectPolyhedronSplitEvents(const std::list<EdgeSPtr>& edges,
                                      PolyhedronSPtr polyhedron,
                                      const CGAL::FT& current_offset,
                                      const std::optional<CGAL::FT>& offset_future_bound,
                                      PQ& queue);
     void collectPolyhedronSplitEvents(PolyhedronSPtr polyhedron,
                                      const CGAL::FT& current_offset,
                                      const std::optional<CGAL::FT>& offset_future_bound,
                                      PQ& queue);

    /**
     * This event occurs when two edges collide.
     * The first edge is always reflex.
     */
    void collectEdgeSplitEvents(const std::list<EdgeSPtr>& edges_1,
                                const std::list<EdgeSPtr>& edges_2,
                                PolyhedronSPtr polyhedron,
                                const bool use_canonical_event_reps,
                                const CGAL::FT& current_offset,
                                const std::optional<CGAL::FT>& offset_future_bound,
                                PQ& queue);
    void collectEdgeSplitEvents(PolyhedronSPtr polyhedron,
                                const CGAL::FT& current_offset,
                                const std::optional<CGAL::FT>& offset_future_bound,
                                PQ& queue);

    /**
     * A reflex vertex reaches a facet.
     */
    void collectPierceEvents(const std::list<VertexSPtr>& vertices,
                             const std::list<FacetSPtr>& facets,
                             PolyhedronSPtr polyhedron,
                             const CGAL::FT& current_offset,
                             const std::optional<CGAL::FT>& offset_future_bound,
                             PQ& queue);
    void collectPierceEvents(PolyhedronSPtr polyhedron,
                             const CGAL::FT& current_offset,
                             const std::optional<CGAL::FT>& offset_future_bound,
                             PQ& queue);

    /**
     * Collect all events for a subset of the polyhedron
     */
    void collectLocalEvents(PolyhedronSPtr polyhedron,
                            const CGAL::FT& current_offset,
                            const std::optional<CGAL::FT>& offset_future_bound,
                            PQ& queue);

    /**
     * Collect all events for the polyhedron
     */
    void collectEvents(PolyhedronSPtr polyhedron,
                       const CGAL::FT& current_offset,
                       const std::optional<CGAL::FT>& offset_future_bound,
                       PQ& queue);

    void printQueue(const PQ& queue);

    /**
     * Checks if the queue updated with local events is displaying the same events
     * as if it had been built from scratch.
     */
    bool checkQueueCorrectness(const PQ& queue,
                               PolyhedronSPtr polyhedron,
                               const CGAL::FT& current_offset,
                               const std::optional<CGAL::FT>& offset_future_bound);

    /**
     * Determines the next event.
     */
    AbstractEventSPtr nextEvent(PQ& queue, const CGAL::FT& offset_current);

    /**
     * Appends a node of an event to the skeleton.
     * It links all adjacent arcs and sheets to this node.
     */
    void appendEventNode(NodeSPtr node);

    PolyhedronSPtr shiftToEventOffset(PolyhedronSPtr polyhedron,
                                      const CGAL::FT& start_offset,
                                      const CGAL::FT& target_offset);

    enum class EventStatus {
      NON_EVENT = 0,
      EVENT_HANDLED,
      EVENT_NOT_HANDLED
    };

    EventStatus handleSaveOffsetEvent(SaveOffsetEventSPtr event,
                                      const CGAL::FT& current_offset,
                                      PolyhedronSPtr polyhedron);
    EventStatus handleConstOffsetEvent(ConstOffsetEventSPtr event,
                                       const CGAL::FT& current_offset,
                                       PolyhedronSPtr polyhedron);
    EventStatus handleVanishEvent(VanishEventSPtr event,
                                  const CGAL::FT& current_offset,
                                  const std::optional<CGAL::FT>& offset_future_bound,
                                  PolyhedronSPtr polyhedron);
    EventStatus handleEdgeEvent(EdgeEventSPtr event,
                                const CGAL::FT& current_offset,
                                const std::optional<CGAL::FT>& offset_future_bound,
                                PolyhedronSPtr polyhedron);
    EventStatus handleEdgeMergeEvent(EdgeMergeEventSPtr event,
                                     const CGAL::FT& current_offset,
                                     const std::optional<CGAL::FT>& offset_future_bound,
                                     PolyhedronSPtr polyhedron);
    EventStatus handleTriangleEvent(TriangleEventSPtr event,
                                    const CGAL::FT& current_offset,
                                    const std::optional<CGAL::FT>& offset_future_bound,
                                    PolyhedronSPtr polyhedron);
    EventStatus handleDblEdgeMergeEvent(DblEdgeMergeEventSPtr event,
                                        const CGAL::FT& current_offset,
                                        const std::optional<CGAL::FT>& offset_future_bound,
                                        PolyhedronSPtr polyhedron);
    EventStatus handleDblTriangleEvent(DblTriangleEventSPtr event,
                                       const CGAL::FT& current_offset,
                                       const std::optional<CGAL::FT>& offset_future_bound,
                                       PolyhedronSPtr polyhedron);
    EventStatus handleTetrahedronEvent(TetrahedronEventSPtr event,
                                       const CGAL::FT& current_offset,
                                       const std::optional<CGAL::FT>& offset_future_bound,
                                       PolyhedronSPtr polyhedron);
    EventStatus handleVertexEvent(VertexEventSPtr event,
                                  const CGAL::FT& current_offset,
                                  const std::optional<CGAL::FT>& offset_future_bound,
                                  PolyhedronSPtr polyhedron);
    EventStatus handleFlipVertexEvent(FlipVertexEventSPtr event,
                                      const CGAL::FT& current_offset,
                                      const std::optional<CGAL::FT>& offset_future_bound,
                                      PolyhedronSPtr polyhedron);
    EventStatus handleSurfaceEvent(SurfaceEventSPtr event,
                                   const CGAL::FT& current_offset,const std::optional<CGAL::FT>& offset_future_bound,
                                   PolyhedronSPtr polyhedron);
    EventStatus handlePolyhedronSplitEvent(PolyhedronSplitEventSPtr event,
                                           const CGAL::FT& current_offset,
                                           const std::optional<CGAL::FT>& offset_future_bound,
                                           PolyhedronSPtr polyhedron);
    EventStatus handleSplitMergeEvent(SplitMergeEventSPtr event,
                                      const CGAL::FT& current_offset,
                                      const std::optional<CGAL::FT>& offset_future_bound,
                                      PolyhedronSPtr polyhedron);
    EventStatus handleEdgeSplitEvent(EdgeSplitEventSPtr event,
                                     const CGAL::FT& current_offset,
                                     const std::optional<CGAL::FT>& offset_future_bound,
                                     PolyhedronSPtr polyhedron);
    EventStatus handlePierceEvent(PierceEventSPtr event,
                                  const CGAL::FT& current_offset,
                                  const std::optional<CGAL::FT>& offset_future_bound,
                                  PolyhedronSPtr polyhedron);
    EventStatus handleEvent(AbstractEventSPtr event,
                            const CGAL::FT& current_offset,
                            const std::optional<CGAL::FT>& offset_future_bound,
                            PolyhedronSPtr polyhedron);

    StraightSkeletonSPtr getResult() const;

protected:
    SimpleStraightSkel(PolyhedronSPtr polyhedron);
    SimpleStraightSkel(PolyhedronSPtr polyhedron, ControllerSPtr controller);
    SimpleStraightSkel(PolyhedronSPtr polyhedron, ControllerSPtr controller,
                       const std::vector<CGAL::FT>& save_offsets,
                       const std::filesystem::path& save_path);

    PolyhedronSPtr polyhedron_;

    ControllerSPtr controller_;

    utils::Base_mesh_offset_visitor* visitor_ = nullptr;

    std::vector<CGAL::FT> save_offsets_;
    std::filesystem::path save_path_;
    bool use_fast_vertex_splitter_;
    AbstractVertexSplitterSPtr vertex_splitter_;
    int edge_event_;

    StraightSkeletonSPtr skel_result_;

    int step_id_;

    std::set<VertexSPtr> post_op_vertices_;
    std::set<EdgeSPtr> post_op_edges_;
    std::set<FacetSPtr> post_op_facets_;

    std::set<VertexSPtr> post_op_vertices_VV_;
    std::set<VertexSPtr> post_op_vertices_pierce_;
    std::set<EdgeSPtr> post_op_edges_edgesplit_;
};

} }

#endif /* ALGO_3D_SIMPLESTRAIGHTSKEL_H */
