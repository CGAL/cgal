// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   algo/3d/SimpleStraightSkel.h
 * author Gernot Walzl
 * date   2012-03-08
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_STRAIGHT_SKELETON_BUILDER_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_STRAIGHT_SKELETON_BUILDER_H

// @fixme yesterday:
//   What about using Shape detection then...
// - 362913 queue correctness assertion
// - 2119792 with perf weights => segfaults
// data/wolfgang/724_star-Butterfly.obj
// data/wolfgang/726_star-PNSplit.obj
// data/wolfgang/chinese_lion_174.obj
// data/wolfgang/venus_267.obj
// -save-time -inf works?

// @fixme:
// - more combinatorial checks should happen at pop time (?)
// - re-using items in event handlers goes against detecting obsolete events with expired pointers
// - perturbation depth-limit without hardcoded digit size values
// - merge vertices with equal position in output (save offset at an event time)

// @fixme later:

// @fixme latest:
// - Fix simultaneous events still happening sometimes (likely the same event multiple times
//   since we don't check the queue before pushing)
// - EPECK -> EPICK could create self-intersections

// @speed
// - don't actually shift at all intermediate steps: we can do everything with base planes,
//   including the computations that use point positions
// - avoid duplicate computations in checkBisectors()
// - delay bisector checks till after cheap combinatorial checks
// - edge split events: add an event ptr within each edge for the next best edge split (might not
//   be valid anymore with pop time checks)
// - store vanish event ptr in edges as to:
//   - avoid recomputing it if the edge is already associated a vanish event
//   - get a tigther future bound on events: e.g., if the edge vanishes at T0, it's pointless
//     to investigate contact events of that edge at after T0
// - Fix and use CGAL_SS3_VV_VERTEX_2_WALK_FACES_FOR_DETECTION
// - If an edge is growing, there is no point computing its vanish event
// - For contact events: exit early if the 4 planes are clearly not intersecting (diametral spheres
//   around the edges of size something?)

// @todo: cleaning
// - HdsUtils / Polyhedron transformation cyclic dependency
// - Do not accept non manifold inputs (clarify doc)
// - clean up all the code related to local queue (horrible variable names, duplicates, etc.)
// - clean up macros
// - check for overly shared objects, redundant function calls (plane normalization, for example)
// - CamelCase in header name, classes, variables (?)
// - tone down the useless shared ptr like the SDS builder

// @todo
// - add tests; doc figures

// @todo later:
// - check if it costly to recompute the time and position at event pop time as make events lighter
// - Re-introduce graph checker
// - add navigation to the skeleton data structure (?)
// - if checking perturbation fails because of self-intersections, use a smaller epsilon
// - perform facet merging using CGAL's region growing and remesh_planar_faces()
//   which only ensures that points are on supporting planes, but does NOT perturb the planes.
// - re-enable the option to translate and scale (?)
// - tolerate non triangulated inputs

// @todo latest:
// - use traits' functors
// - get rid of the exact construction requirement? At least if we do not have to split high degree
//   vertices, it should be possible, but that would require writing filtered predicates like SLS2's.
// - factorize the three VV events but be very careful with the tiny differences
// - lighter polyhedron data structures
// - write a sanitization algorithm without perturbaton, something akin to: randTiltPlanesv3(p, eps=0),
// - implement lazy perturbations where we only perturb if we encounter an issue?
//   * Would cost more (detection of simultaneous events, etc.)
//   * How to detect "hidden" simultaneous events?
// - splitting high degree vertices in reasonable time

// ----

/*
  As to not waste energy building the skeleton if we do not care about it.
*/
#define CGAL_SS3_NO_SKELETON_DS

/*
  Some events can be detected from multime elements. Reduce that to a single element.
  Not used when adding local elements to the queue because we might encounter
  only *some* representatives in the subset, excludin the canonical one.
*/
#define CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS

/*
  Use CGAL::box_d to filter edge pairs
*/
#define CGAL_SS3_DETECT_EDGE_SPLIT_EVENTS_WITH_BOX_D

// use v2 (faster) is inside
#define CGAL_SLS3_NEW_IS_INSIDE

// ----

// #define CGAL_SS3_DEBUG_PRINT_QUEUE

// ----

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/SDS/Straight_skeleton.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Geom_utils.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/HDS_utils.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/vertex_splitters.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/events.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_perturbation.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_self_intersection.h>
#include <CGAL/Straight_skeleton_3/IO/OBJ.h>

#include <CGAL/assertions.h>
#ifdef CGAL_SS3_DETECT_EDGE_SPLIT_EVENTS_WITH_BOX_D
# include <CGAL/box_intersection_d.h>
#endif
#ifdef CGAL_SS3_RUN_TIMERS
# include <CGAL/Real_timer.h>
#endif

#include <array>
#include <filesystem>
#include <limits>
#include <list>
#include <optional>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
struct Base_mesh_offset_visitor
{
  using FT = typename Traits::FT;
  using PolyhedronSPtr = std::shared_ptr<HDS::Polyhedron<Traits> >;
  using AbstractEventSPtr = std::shared_ptr<AbstractEvent<Traits> >;

  virtual bool go_further(int, PolyhedronSPtr, FT) = 0;
  virtual void before_event(PolyhedronSPtr, FT, AbstractEventSPtr) = 0;
  virtual void on_save_event(PolyhedronSPtr, FT) = 0;
  virtual void after_event(PolyhedronSPtr, FT) = 0;
};

template <typename Traits>
struct Default_mesh_offset_visitor
  : public Base_mesh_offset_visitor<Traits>
{
  using FT = typename Traits::FT;
  using PolyhedronSPtr = std::shared_ptr<HDS::Polyhedron<Traits> >;
  using AbstractEventSPtr = std::shared_ptr<AbstractEvent<Traits> >;

  bool go_further(int, PolyhedronSPtr, FT) override { return true; }
  void before_event(PolyhedronSPtr, FT, AbstractEventSPtr) override { }
  void on_save_event(PolyhedronSPtr, FT) override { }
  void after_event(PolyhedronSPtr, FT) override { }
};

template <typename Traits>
class SimpleStraightSkel
{
  using SimpleStraightSkelSPtr = std::shared_ptr<SimpleStraightSkel<Traits> >;

private:
  // Geometry
  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;
  using Segment_3 = typename Traits::Segment_3;
  using Vector_3 = typename Traits::Vector_3;
  using Line_3 = typename Traits::Line_3;
  using Plane_3 = typename Traits::Plane_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Segment3SPtr = std::shared_ptr<Segment_3>;
  using Vector3SPtr = std::shared_ptr<Vector_3>;
  using Line3SPtr = std::shared_ptr<Line_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;

private:
  // Polyhedron Data Structure
  using Polyhedron = HDS::Polyhedron<Traits>;
  using PolyhedronWPtr = typename Polyhedron::PolyhedronWPtr;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Vertex = typename Polyhedron::template Vertex<Traits>;
  using VertexWPtr = typename Polyhedron::VertexWPtr;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using VertexData = typename Polyhedron::VertexData;
  using VertexDataSPtr = typename Polyhedron::VertexDataSPtr;
  using SkelVertexData = typename Polyhedron::SkelVertexData;
  using SkelVertexDataSPtr = typename Polyhedron::SkelVertexDataSPtr;
  using Edge = typename Polyhedron::template Edge<Traits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using EdgeData = typename Polyhedron::EdgeData;
  using EdgeDataSPtr = typename Polyhedron::EdgeDataSPtr;
  using SkelEdgeData = typename Polyhedron::SkelEdgeData;
  using SkelEdgeDataSPtr = typename Polyhedron::SkelEdgeDataSPtr;
  using Facet = typename Polyhedron::template Facet<Traits>;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;
  using FacetData = typename Polyhedron::FacetData;
  using FacetDataSPtr = typename Polyhedron::FacetDataSPtr;
  using SkelFacetData = typename Polyhedron::SkelFacetData;
  using SkelFacetDataSPtr = typename Polyhedron::SkelFacetDataSPtr;

private:
  // Straight Skeleton Data Structure
  using StraightSkeleton = SDS::StraightSkeleton<Traits>;
  using StraightSkeletonWPtr = typename StraightSkeleton::StraightSkeletonWPtr;
  using StraightSkeletonSPtr = typename StraightSkeleton::StraightSkeletonSPtr;

  using Node = typename StraightSkeleton::Node;
  using NodeSPtr = typename StraightSkeleton::NodeSPtr;
  using Arc = typename StraightSkeleton::Arc;
  using ArcWPtr = typename StraightSkeleton::ArcWPtr;
  using ArcSPtr = typename StraightSkeleton::ArcSPtr;
  using Sheet = typename StraightSkeleton::Sheet;
  using SheetWPtr = typename StraightSkeleton::SheetWPtr;
  using SheetSPtr = typename StraightSkeleton::SheetSPtr;

private:
  // Vertex Splitters
  using AbstractVertexSplitter = algorithm::AbstractVertexSplitter<Traits>;
  using AbstractVertexSplitterSPtr = std::shared_ptr<AbstractVertexSplitter>;
  using CombiVertexSplitter = algorithm::CombiVertexSplitter<Traits>;
  using ConvexVertexSplitter = algorithm::ConvexVertexSplitter<Traits>;

private:
  // Events
  using AbstractEvent = algorithm::AbstractEvent<Traits>;
  using AbstractEventSPtr = std::shared_ptr<AbstractEvent>;

  using ConstTimeEvent = algorithm::ConstTimeEvent<Traits>;
  using ConstTimeEventSPtr = std::shared_ptr<ConstTimeEvent>;
  using SaveEvent = algorithm::SaveEvent<Traits>;
  using SaveEventSPtr = std::shared_ptr<SaveEvent>;

  using VanishEvent = algorithm::VanishEvent<Traits>;
  using VanishEventSPtr = std::shared_ptr<VanishEvent>;
  using EdgeEvent = algorithm::EdgeEvent<Traits>;
  using EdgeEventSPtr = std::shared_ptr<EdgeEvent>;
  using EdgeMergeEvent = algorithm::EdgeMergeEvent<Traits>;
  using EdgeMergeEventSPtr = std::shared_ptr<EdgeMergeEvent>;
  using TriangleEvent = algorithm::TriangleEvent<Traits>;
  using TriangleEventSPtr = std::shared_ptr<TriangleEvent>;
  using DblEdgeMergeEvent = algorithm::DblEdgeMergeEvent<Traits>;
  using DblEdgeMergeEventSPtr = std::shared_ptr<DblEdgeMergeEvent>;
  using DblTriangleEvent = algorithm::DblTriangleEvent<Traits>;
  using DblTriangleEventSPtr = std::shared_ptr<DblTriangleEvent>;
  using TetrahedronEvent = algorithm::TetrahedronEvent<Traits>;
  using TetrahedronEventSPtr = std::shared_ptr<TetrahedronEvent>;

  using VertexEvent = algorithm::VertexEvent<Traits>;
  using VertexEventSPtr = std::shared_ptr<VertexEvent>;
  using FlipVertexEvent = algorithm::FlipVertexEvent<Traits>;
  using FlipVertexEventSPtr = std::shared_ptr<FlipVertexEvent>;
  using SurfaceEvent = algorithm::SurfaceEvent<Traits>;
  using SurfaceEventSPtr = std::shared_ptr<SurfaceEvent>;
  using PolyhedronSplitEvent = algorithm::PolyhedronSplitEvent<Traits>;
  using PolyhedronSplitEventSPtr = std::shared_ptr<PolyhedronSplitEvent>;
  using SplitMergeEvent = algorithm::SplitMergeEvent<Traits>;
  using SplitMergeEventSPtr = std::shared_ptr<SplitMergeEvent>;
  using EdgeSplitEvent = algorithm::EdgeSplitEvent<Traits>;
  using EdgeSplitEventSPtr = std::shared_ptr<EdgeSplitEvent>;
  using PierceEvent = algorithm::PierceEvent<Traits>;
  using PierceEventSPtr = std::shared_ptr<PierceEvent>;

  enum class EventStatus {
    NON_EVENT = 0,
    EVENT_HANDLED,
    EVENT_NOT_HANDLED
  };

private:
  using KernelFactory = kernel::KernelFactory<Traits>;
  using KernelWrapper = kernel::KernelWrapper<Traits>;
  using GeomUtils = algorithm::GeomUtils<Traits>;
  using HdsUtils = algorithm::HdsUtils<Traits>;
  using Transformation = algorithm::PolyhedronTransformation<Traits>;
  using Perturbation = algorithm::PolyhedronPerturbation<Traits>;
  using SelfIntersection = algorithm::SelfIntersection<Traits>;

private:
  using PQ = std::priority_queue<AbstractEventSPtr,
                                 std::vector<AbstractEventSPtr>,
                                 AbstractEventSPtrCompare<Traits> >;

public:
  SimpleStraightSkel(PolyhedronSPtr polyhedron)
    : polyhedron_(polyhedron),
      save_path_(std::filesystem::current_path()),
      skel_result_(StraightSkeleton::create())
  {
    initVertexSplitter();
    initEdgeEvent();
  }

  SimpleStraightSkel(PolyhedronSPtr polyhedron,
                     const std::vector<FT>& save_times,
                     const std::filesystem::path& save_path)
    : polyhedron_(polyhedron),
      save_times_(save_times), // intentional copy
      save_path_(save_path),
      skel_result_(StraightSkeleton::create())
  {
    std::sort(save_times_.begin(), save_times_.end(),
              [](const FT& a, const FT& b) { return CGAL::abs(a) < CGAL::abs(b); });

    initVertexSplitter();
    initEdgeEvent();
  }

  ~SimpleStraightSkel()
  {
    polyhedron_.reset();
    vertex_splitter_.reset();
    events_.clear();
    skel_result_.reset();
  }

  static SimpleStraightSkelSPtr create(PolyhedronSPtr polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    return std::make_shared<SimpleStraightSkel>(polyhedron);
  }

  static SimpleStraightSkelSPtr create(PolyhedronSPtr polyhedron,
                                       const std::vector<FT>& save_times)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    return std::make_shared<SimpleStraightSkel>(polyhedron, save_times);
  }

  static SimpleStraightSkelSPtr create(PolyhedronSPtr polyhedron,
                                       const std::vector<FT>& save_times,
                                       const std::filesystem::path& save_path)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    return std::make_shared<SimpleStraightSkel>(polyhedron, save_times, save_path);
  }

  void setVisitor(Base_mesh_offset_visitor<Traits>* visitor)
  {
    visitor_ = visitor;
  }

  void initVertexSplitter()
  {
    ConfigurationSPtr config = Configuration::getInstance();
    std::string s_vertex_splitter;
    if (config->isLoaded()) {
      s_vertex_splitter = config->getString("Algorithm", "vertex_splitter");
      if (s_vertex_splitter.compare("CombiVertexSplitter") == 0) {
        vertex_splitter_ = CombiVertexSplitter::create();
      } else if (s_vertex_splitter.compare("ConvexVertexSplitter") == 0) {
        vertex_splitter_ = ConvexVertexSplitter::create();
      } else {
        CGAL_SS3_SPLITTER_TRACE("Warning: option '" << s_vertex_splitter << "' not found.");
        CGAL_SS3_SPLITTER_TRACE("Using 'CombiVertexSplitter'.");
        vertex_splitter_ = CombiVertexSplitter::create();
      }
    } else {
      vertex_splitter_ = CombiVertexSplitter::create();
    }
  }

  void initEdgeEvent()
  {
    ConfigurationSPtr config = Configuration::getInstance();
    std::string s_edge_event;
    if (config->isLoaded()) {
      s_edge_event = Configuration::getInstance()->getString("Algorithm", "edge_event");
      if (s_edge_event.compare("convex") == 0) {
        edge_event_ = 0;
      } else if (s_edge_event.compare("reflex") == 0) {
        edge_event_ = 1;
      } else if (s_edge_event.compare("flip") == 0) {
        edge_event_ = 2;
      } else {
        CGAL_SS3_CORE_TRACE("Warning: option '" << s_edge_event << "' not found.");
        CGAL_SS3_CORE_TRACE("Using 'convex'.");
        edge_event_ = 0;
        s_edge_event = "convex";
      }
    } else {
      edge_event_ = 0;
      s_edge_event = "convex";
    }
  }

  bool run()
  {
    CGAL_SS3_CORE_TRACE_V(1, "== Straight Skeleton 3D started ==");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    CGAL_assertion(bool(polyhedron_));
    CGAL_assertion(polyhedron_->isConsistent());

#ifdef CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("results/input.obj", polyhedron_, false /*do not triangulate*/);
#endif

    CGAL_SS3_CORE_TRACE_V(1, polyhedron_->vertices().size() << " NV " << polyhedron_->facets().size() << " NF");

    CGAL_assertion(Perturbation::doAll2PlanesIntersect(polyhedron_));
    CGAL_assertion(Perturbation::doAll3PlanesIntersect(polyhedron_));
    CGAL_assertion(!SelfIntersection::hasSelfIntersectingSurface(polyhedron_));

    skel_result_->setPolyhedron(polyhedron_); // skeleton's polyhedron is fixed
    PolyhedronSPtr polyhedron = polyhedron_->clone();

    // store base plane coefficients
    cacheBasePlanes(polyhedron);

// @tmp some hardcoded weights for specific inputs
// #define CGAL_SS3_ACUTE_WEIGHTS
// #define CGAL_SS3_MERGING_WEIGHTS
// #define CGAL_SS3_PERFORMANCE_WEIGHTS
#if defined(CGAL_SS3_ACUTE_WEIGHTS) || defined(CGAL_SS3_MERGING_WEIGHTS) || defined(CGAL_SS3_PERFORMANCE_WEIGHTS)
# ifdef CGAL_SS3_ACUTE_WEIGHTS
    const FT x_speed = 20;
    const FT y_speed = 20;
    const FT z_speed = 20;
    const FT other_speed = 18.7939;
# elif defined(CGAL_SS3_MERGING_WEIGHTS)
    const FT x_speed = 20;
    const FT y_speed = 20;
    const FT z_speed = 20;
    const FT other_speed = 19.8777;
# elif defined(CGAL_SS3_PERFORMANCE_WEIGHTS)
    const FT x_speed = 5;
    const FT y_speed = 5;
    const FT z_speed = 2;
    const FT other_speed = 5;
# else
#  error
# endif

    for (const FacetSPtr& facet : polyhedron->facets()) {
      FT speed = other_speed;
      const auto pl = facet->getPlane();
      const auto normal = KernelFactory::createVector3(pl);
      // DEBUG_PRINT("SP X " << CGAL::scalar_product(*normal, Vector_3(1,0,0)));
      // DEBUG_PRINT("SP Y " << CGAL::scalar_product(*normal, Vector_3(0,1,0)));
      // DEBUG_PRINT("SP Z " << CGAL::scalar_product(*normal, Vector_3(0,0,1)));
      if (CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector_3(1,0,0))) - 1) < 1e-3)
        speed = x_speed;
      if (CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector_3(0,1,0))) - 1) < 1e-3)
        speed = y_speed;
      if (CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector_3(0,0,1))) - 1) < 1e-3)
        speed = z_speed;

      HdsUtils::setSpeed(facet, speed);
      // DEBUG_PRINT("speed to " << speed);
    }
#endif

    if (!init(polyhedron, vertex_splitter_)) {
      CGAL_SS3_CORE_TRACE_V(8, "Error: failed to initialize polyhedron");
      return false;
    }

    // If we stop immediately after the last save event, there is no point investigating events
    // that are farther away than the last save event's time
    std::optional<FT> time_future_bound;
    if (!save_times_.empty()) {
      ConfigurationSPtr config = Configuration::getInstance();
      if (config->isLoaded()) {
        if ((config->contains("Algorithm", "stop_after_last_save_event") &&
              config->getBool("Algorithm", "stop_after_last_save_event"))) {
          time_future_bound = save_times_.back();
        }
      }
    }

    step_id_ = -1;
    FT current_time = 0;
    FT upcoming_event_time;

    CGAL_assertion_code(const bool is_emptiness_expected = save_times_.empty();)

    PQ queue;
    collectEvents(polyhedron, current_time, time_future_bound, queue);

    for (;;) {
      ++step_id_;

      CGAL_SS3_CORE_TRACE_V(2, "\n=========== ITERATION #" << step_id_ << " AT TIME " << current_time);
      CGAL_SS3_CORE_TRACE_V(2, polyhedron->vertices().size() << " NV " << polyhedron->facets().size() << " NF");

      if (visitor_) {
        if (!visitor_->go_further(step_id_, polyhedron, current_time)) {
          CGAL_SS3_CORE_TRACE_V(2, "Stopping on visitor request");
          break;
        }
      }

      CGAL_assertion_code(for (const FacetSPtr& facet : polyhedron->facets()) {)
      CGAL_assertion(facet->getPlane()->a() == HdsUtils::getBasePlane(facet)->a());
      CGAL_assertion(facet->getPlane()->b() == HdsUtils::getBasePlane(facet)->b());
      CGAL_assertion(facet->getPlane()->c() == HdsUtils::getBasePlane(facet)->c());
      CGAL_assertion_code(FT speed = HdsUtils::getSpeed(facet);)
      CGAL_assertion(facet->getPlane()->d() == HdsUtils::getBasePlane(facet)->d() - speed * current_time);
      CGAL_assertion_code(})

      AbstractEventSPtr event = nextEvent(queue, polyhedron, current_time);
      if (!event) {
        CGAL_SS3_CORE_TRACE_V(2, "No more events to treat");
        break;
      }

      CGAL_SS3_CORE_TRACE_V(2, "popped E" << event->getID() << " Type [" << event->getType() << "]");

      static int event_id = -1;
      CGAL_SS3_CORE_TRACE_V(2, "--> Accepted event #" << ++event_id << " " << event->toString() << " --");

      upcoming_event_time = event->getTime();

      // the next event should be at a time that is further away than the current one
      CGAL_assertion(upcoming_event_time < current_time);

#ifdef CGAL_SS3_RUN_TIMERS
      CGAL_SS3_CORE_TRACE_V(2, "current elapsed time: " << timer.time());
#endif

      if (visitor_) {
        visitor_->before_event(polyhedron, current_time, event);
      }

      // Event treatment
      EventStatus es = handleEvent(event, current_time, time_future_bound, polyhedron);
      CGAL_assertion(es != EventStatus::EVENT_NOT_HANDLED);
      if (es == EventStatus::NON_EVENT) {
        continue;
      }

      current_time = upcoming_event_time;

#ifdef CGAL_SS3_DUMP_FILES
      IO::OBJFile::save("results/event_" + std::to_string(event_id) + ".obj", polyhedron, false /*do_triangulate*/);
      IO::OBJFile::save("results/event_" + std::to_string(event_id) + "_triangulated.obj", polyhedron);
#endif

      if (visitor_ && event->getType() == AbstractEvent::SAVE_EVENT) {
        visitor_->on_save_event(polyhedron, current_time);
      }

      CGAL_SS3_CORE_TRACE_V(2, skel_result_->toString());

#ifdef CGAL_SS3_DUMP_FILES
      // Dump skeleton nodes in an .xyz file
      std::ofstream nodes_out("final_nodes.xyz");
      nodes_out.precision(17);
      for (NodeSPtr node : skel_result_->nodes()) {
        nodes_out << *(node->getPoint()) << "\n";
      }
      nodes_out.close();

      // Dump skeleton arcs as CGAL polylines
      std::ofstream arcs_out("final_arcs.polylines.txt");
      arcs_out.precision(17);
      for (ArcSPtr arc : skel_result_->arcs()) {
        arcs_out << "2 ";
        arcs_out << *(arc->getNodeSrc()->getPoint()) << " ";
        if (arc->hasNodeDst()) {
          arcs_out << *(arc->getNodeDst()->getPoint()) << "\n";
        } else {
          Point3SPtr src_pt = arc->getNodeSrc()->getPoint();
          Vector3SPtr dir = arc->getDirection();
          constexpr double ray_length = 0.1; // @todo relative value
          Point_3 ray_pt = *src_pt + ray_length * (*dir);
          arcs_out << ray_pt << "\n";
        }
      }
      arcs_out.close();
#endif

      CGAL_postcondition(polyhedron->isConsistent());
      CGAL_postcondition(skel_result_->isConsistent());

      if (visitor_) {
        visitor_->after_event(polyhedron, current_time);
      }

      // If we are only interested in specific times, there is no point going further
      if (event->getType() == AbstractEvent::SAVE_EVENT && save_times_.empty()) {
        ConfigurationSPtr config = Configuration::getInstance();
        if (config->isLoaded() &&
            config->contains("Algorithm", "stop_after_last_save_event") &&
            config->getBool("Algorithm", "stop_after_last_save_event")) {
          break;
        }
      }

      // Update the event priority queue
      collectLocalEvents(polyhedron, current_time, time_future_bound, queue);

      post_op_vertices_.clear();
      post_op_edges_.clear();
      post_op_facets_.clear();

      post_op_vertices_VV_.clear();
      post_op_vertices_pierce_.clear();
      post_op_edges_edgesplit_.clear();
    }

    CGAL_SS3_CORE_TRACE_V(8, "== Straight Skeleton 3D finished ==");

    CGAL_warning(!is_emptiness_expected || polyhedron->empty());

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
#endif

    CGAL_SS3_CORE_TRACE_V(2, eventSummary());


    CGAL_assertion(skel_result_->isConsistent(false /*is_partial*/));
#ifdef CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("final_skeleton.obj", skel_result_, true /*convert_to_double*/);
#endif

    return true;
  }

  /**
    * Creates a new node for the vertex data.
    * Used by init(...) only.
    */
  static NodeSPtr createNode(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    NodeSPtr result = Node::create();
    result->setTime(0);
    result->setPoint(vertex->getPoint());
    HdsUtils::setNode(vertex, result);
    return result;
  }

  /**
    * Creates a new arc for the vertex data.
    * The node of the vertex data has to be set before.
    */
  static ArcSPtr createArc(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->degree() == 3);
    ArcSPtr result = ArcSPtr();
    CGAL_precondition(vertex->hasData());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());

    std::array<FacetSPtr, 3> facets;
    unsigned int i = 0;
    typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (i < 3 && it_f != vertex->facets().end()) {
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        facets[i] = facet;
        ++i;
      }
    }

    CGAL_assertion(i == 3);

    Vector3SPtr direction;
    Plane3SPtr plane_1 = facets[0]->getPlane();
    Plane3SPtr plane_2 = facets[1]->getPlane();
    Plane3SPtr plane_3 = facets[2]->getPlane();

    const FT& speed_1 = HdsUtils::getSpeed(facets[0]);
    const FT& speed_2 = HdsUtils::getSpeed(facets[1]);
    const FT& speed_3 = HdsUtils::getSpeed(facets[2]);

    Vector3SPtr n_1 = KernelFactory::createVector3(plane_1);
    Vector3SPtr n_2 = KernelFactory::createVector3(plane_2);
    Vector3SPtr n_3 = KernelFactory::createVector3(plane_3);

    // @fixme likely wrong
    direction = KernelFactory::createVector3(speed_1 * (*n_1) + speed_2 * (*n_2) + speed_3 * (*n_3));

    if (direction) {
      result = Arc::create(data->getNode(), direction);
      data->setArc(result);
    }

    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  /**
    * Creates a new sheet for the edge data.
    * The arcs of the vertices have to be set before.
    */
  static SheetSPtr createSheet(const EdgeSPtr& edge)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    CGAL_precondition(edge->getVertexSrc() != edge->getVertexDst());

    SheetSPtr result = SheetSPtr();

    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    CGAL_SS3_DEBUG_SPTR(facet_l);
    CGAL_SS3_DEBUG_SPTR(facet_r);

    Plane3SPtr plane_l = facet_l->getPlane();
    Plane3SPtr plane_r = facet_r->getPlane();
    FacetSPtr facet_b = HdsUtils::getFacetOrigin(facet_l);
    FT speed_l = HdsUtils::getSpeed(facet_l);
    FacetSPtr facet_f = HdsUtils::getFacetOrigin(facet_r);
    FT speed_r = HdsUtils::getSpeed(facet_r);

    Plane3SPtr plane_sheet;
    if (speed_l == speed_r) {
      if (HdsUtils::isReflex(edge)) {
        plane_sheet = KernelWrapper::bisector(KernelWrapper::opposite(plane_l), plane_r);
      } else {
        plane_sheet = KernelWrapper::bisector(plane_l, KernelWrapper::opposite(plane_r));
      }
    } else {
      Line3SPtr line = KernelWrapper::intersection(plane_l, plane_r);
      Plane3SPtr offset_l = GeomUtils::offsetPlane(plane_l, -speed_l);
      Plane3SPtr offset_r = GeomUtils::offsetPlane(plane_r, -speed_r);
      Line3SPtr line_offset = KernelWrapper::intersection(offset_l, offset_r);
      Point3SPtr point_1 = KernelFactory::createPoint3(line->point());
      Vector3SPtr direction = KernelFactory::createVector3(line);
      Point3SPtr point_2 = KernelFactory::createPoint3(*point_1 + *direction);
      Point3SPtr point_3 = KernelFactory::createPoint3(line_offset->point());
      if (HdsUtils::isReflex(edge)) {
        plane_sheet = KernelFactory::createPlane3(point_3, point_2, point_1);
      } else {
        plane_sheet = KernelFactory::createPlane3(point_1, point_2, point_3);
      }
    }

    result = Sheet::create();
    result->setPlane(plane_sheet);
    result->setFacetB(facet_b);
    result->setFacetF(facet_f);

    HdsUtils::setSheet(edge, result);

    NodeSPtr node_src = HdsUtils::getNode(edge->getVertexSrc());
    NodeSPtr node_dst = HdsUtils::getNode(edge->getVertexDst());
    result->addNode(node_src);
    if (node_src != node_dst) {
      result->addNode(node_dst);
    }
    result->addArc(HdsUtils::getArc(edge->getVertexSrc()));
    result->addArc(HdsUtils::getArc(edge->getVertexDst()));

    return result;
  }

  void mergeSheets(EdgeSPtr edge_into,
                   EdgeSPtr edge_from)
  {
    CGAL_precondition(edge_into && edge_from);
    CGAL_precondition(edge_into != edge_from);
    CGAL_precondition(HdsUtils::getSheet(edge_into) && HdsUtils::getSheet(edge_from));

    if (HdsUtils::getSheet(edge_into) == HdsUtils::getSheet(edge_from)) {
      return;
    }

    skel_result_->mergeSheets(HdsUtils::getSheet(edge_into), HdsUtils::getSheet(edge_from));
    CGAL_assertion(bool(HdsUtils::getSheet(edge_into)));
    HdsUtils::setSheet(edge_from, HdsUtils::getSheet(edge_into));
  }

  /**
    * Store within each facet the coefficients of the plane at t=0
    */
  void cacheBasePlanes(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    for (const FacetSPtr& facet : polyhedron->facets()) {
      Plane3SPtr plane = facet->getPlane();
      CGAL_SS3_DEBUG_SPTR(plane);
      CGAL_assertion(KernelWrapper::hasNormalizedPlane(plane));
      HdsUtils::setBasePlane(facet, plane);
    }
  }

  /**
    * Split all vertices with degree > 3 and
    * initializes the data variables of all edges and vertices.
    */
  bool init(const PolyhedronSPtr& polyhedron,
            AbstractVertexSplitterSPtr vertex_splitter)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    bool result = true;

    CGAL_SS3_CORE_TRACE("Input: " << polyhedron->vertices().size() << " NV " << polyhedron->facets().size() << " NF");

    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      CGAL_precondition(vertex->degree() >= 3);
      if (!vertex->hasData()) {
        SkelVertexData::create(vertex);
      }
    }

    // vertex splitting
    std::list<VertexSPtr> vertices_tosplit;

    for (const VertexSPtr& vertex : polyhedron->vertices()) {
#ifndef CGAL_SS3_NO_SKELETON_DS
      NodeSPtr node = createNode(vertex);
      CGAL_SS3_DEBUG_SPTR(node);
      skel_result_->addNode(node);
#endif

      if (vertex->degree() > 3) {
        vertices_tosplit.push_back(vertex);
      }
    }

    CGAL_SS3_CORE_TRACE_V(2, vertices_tosplit.size() << " vertices to split");
    CGAL_SS3_CORE_TRACE_V(2, "Using " << vertex_splitter->toString() << " to split vertices.");

    for (const VertexSPtr& vertex : vertices_tosplit) {
      CGAL_SS3_SPLITTER_TRACE_V(8, "Generic split vertex:\n" << vertex->toString());
      if (vertex->degree() > 15) {
        CGAL_SS3_SPLITTER_TRACE_V(1, "Warning: degree of vertex is so high that even a combinatorial split with eary exit will take forever.");
      }
      vertex_splitter->splitVertex(vertex);
    }

#ifdef CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("results/post_split.obj", polyhedron, false /*do not triangulate*/);
#endif

    CGAL_postcondition_code(for (auto v : polyhedron->vertices()))
    CGAL_postcondition(v->getID() != -1);
    CGAL_postcondition_code(for (auto e : polyhedron->edges()))
    CGAL_postcondition(e->getID() != -1);
    CGAL_postcondition_code(for (auto f : polyhedron->facets()))
    CGAL_postcondition(f->getID() != -1);

#ifndef CGAL_SS3_NO_SKELETON_DS
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      CGAL_assertion(vertex->degree() == 3);
      ArcSPtr arc = createArc(vertex);
      CGAL_SS3_DEBUG_SPTR(arc);
      skel_result_->addArc(arc);
    }
#endif

    for (const EdgeSPtr& edge : polyhedron->edges()) {
      if (!edge->hasData()) {
        SkelEdgeData::create(edge);
      }

#ifndef CGAL_SS3_NO_SKELETON_DS
      SheetSPtr sheet = createSheet(edge);
      CGAL_SS3_DEBUG_SPTR(sheet);

      // Squatting the arc type, this isn't really an arc, but a contour edge.
      // Contours are useful to get a closed polygon (with holes) to draw sheets.
      ArcSPtr contour = std::make_shared<Arc>(HdsUtils::getNode(edge->getVertexSrc()),
                                              HdsUtils::getNode(edge->getVertexDst()));
      sheet->addContour(contour);

      skel_result_->addSheet(sheet);
#endif
    }

    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (!facet->hasData()) {
        SkelFacetData::create(facet);
      }
    }

    return result;
  }

  static bool checkBisectorsV2(const EdgeSPtr& edge,
                               const Point3SPtr& point,
                               const FT& event_time)
  {
    std::optional<FT> vanish_time = HdsUtils::getVanishTime(edge);
    Point3SPtr o_src = Transformation::offsetPointFromBase(edge->getVertexSrc(), event_time);
    if (vanish_time.has_value() && event_time == *vanish_time) {
      return (*point == *o_src);
    }

    Point3SPtr o_dst = Transformation::offsetPointFromBase(edge->getVertexDst(), event_time);
    CGAL_assertion(*o_src != *o_dst);
    CGAL_assertion(CGAL::collinear(*o_src, *point, *o_dst));

    return CGAL::collinear_are_ordered_along_line(*o_src, *point, *o_dst);
  }

  /**
  * Returns the intersection time of the 4 shifting planes.
  */
  static Point3SPtr intersectionPointOffsetPlanes(const FacetSPtr& facet_0,
                                                  const FacetSPtr& facet_1,
                                                  const FacetSPtr& facet_2,
                                                  const FacetSPtr& facet_3)
  {
    CGAL_SS3_DEBUG_SPTR(facet_0);
    CGAL_SS3_DEBUG_SPTR(facet_1);
    CGAL_SS3_DEBUG_SPTR(facet_2);
    CGAL_SS3_DEBUG_SPTR(facet_3);

    Plane3SPtr plane_0 = HdsUtils::getBasePlane(facet_0);
    Plane3SPtr plane_1 = HdsUtils::getBasePlane(facet_1);
    Plane3SPtr plane_2 = HdsUtils::getBasePlane(facet_2);
    Plane3SPtr plane_3 = HdsUtils::getBasePlane(facet_3);

    const FT& speed_0 = HdsUtils::getSpeed(facet_0);
    const FT& speed_1 = HdsUtils::getSpeed(facet_1);
    const FT& speed_2 = HdsUtils::getSpeed(facet_2);
    const FT& speed_3 = HdsUtils::getSpeed(facet_3);

    return GeomUtils::intersectionPointOffsetPlanes(plane_0, speed_0, plane_1, speed_1,
                                                    plane_2, speed_2, plane_3, speed_3);
  }

  /**
    * Returns the intersection time of the 4 shifting planes.
    */
  FT intersectionTimeOffsetPlanes(const FacetSPtr& facet_0,
                                  const FacetSPtr& facet_1,
                                  const FacetSPtr& facet_2,
                                  const FacetSPtr& facet_3)
  {
    CGAL_SS3_DEBUG_SPTR(facet_0);
    CGAL_SS3_DEBUG_SPTR(facet_1);
    CGAL_SS3_DEBUG_SPTR(facet_2);
    CGAL_SS3_DEBUG_SPTR(facet_3);

    Plane3SPtr plane_0 = HdsUtils::getBasePlane(facet_0);
    Plane3SPtr plane_1 = HdsUtils::getBasePlane(facet_1);
    Plane3SPtr plane_2 = HdsUtils::getBasePlane(facet_2);
    Plane3SPtr plane_3 = HdsUtils::getBasePlane(facet_3);

    const FT& speed_0 = HdsUtils::getSpeed(facet_0);
    const FT& speed_1 = HdsUtils::getSpeed(facet_1);
    const FT& speed_2 = HdsUtils::getSpeed(facet_2);
    const FT& speed_3 = HdsUtils::getSpeed(facet_3);

    return GeomUtils::intersectionTimeOffsetPlanes(plane_0, speed_0, plane_1, speed_1,
                                                   plane_2, speed_2, plane_3, speed_3);
  }

  Point3SPtr vanishesAtPoint(const EdgeSPtr& edge)
  {
    CGAL_SS3_CORE_TRACE_V(16, "vanishesAtPoint " << edge->toString());

    CGAL_SS3_DEBUG_SPTR(edge);

    FacetSPtr facetL = edge->getFacetL();
    FacetSPtr facetR = edge->getFacetR();
    CGAL_SS3_CORE_TRACE_V(16, "facetL: " << facetL->getID());
    CGAL_SS3_CORE_TRACE_V(16, "facetR: " << facetR->getID());
    CGAL_assertion(facetL && facetR && facetL != facetR);

    FacetSPtr facetP = edge->prev(facetL)->other(facetL);
    FacetSPtr facetN = edge->next(facetL)->other(facetL);
    CGAL_SS3_CORE_TRACE_V(16, "facetP: " << facetP->getID());
    CGAL_SS3_CORE_TRACE_V(16, "facetN: " << facetN->getID());
    CGAL_assertion(facetP && facetP != facetL && facetP != facetR);
    CGAL_assertion(facetN && facetN != facetL && facetN != facetR && facetN != facetP);

    return intersectionPointOffsetPlanes(facetL, facetP, facetR, facetN);
  }

  /**
    * Returns the time at which the edge will vanish.
    */
  std::optional<FT> vanishesAtTime(const EdgeSPtr& edge,
                                   const std::optional<FT>& time_past_bound = std::nullopt,
                                   const std::optional<FT>& time_future_bound = std::nullopt)
  {
    CGAL_SS3_CORE_TRACE_V(16, "Computing vanish time of " << edge->toString());

    CGAL_SS3_DEBUG_SPTR(edge);

    FacetSPtr facetL = edge->getFacetL();
    FacetSPtr facetR = edge->getFacetR();
    CGAL_SS3_CORE_TRACE_V(16, "facetL: " << facetL->getID());
    CGAL_SS3_CORE_TRACE_V(16, "facetR: " << facetR->getID());
    CGAL_assertion(facetL && facetR && facetL != facetR);

    FacetSPtr facetP = edge->prev(facetL)->other(facetL);
    FacetSPtr facetN = edge->next(facetL)->other(facetL);
    CGAL_SS3_CORE_TRACE_V(16, "facetP: " << facetP->getID());
    CGAL_SS3_CORE_TRACE_V(16, "facetN: " << facetN->getID());
    CGAL_assertion(facetP && facetP != facetL && facetP != facetR);
    CGAL_assertion(facetN && facetN != facetL && facetN != facetR && facetN != facetP);

    const FT vanish_time = intersectionTimeOffsetPlanes(facetL, facetP, facetR, facetN);

    if (time_past_bound && vanish_time >= *time_past_bound) {
      CGAL_SS3_TRAITS_TRACE("Vanish event is strictly in the past");
      HdsUtils::setVanishTime(edge, std::nullopt);
      return { };
    }

    if (time_future_bound && vanish_time < *time_future_bound) {
      CGAL_SS3_TRAITS_TRACE("Vanish event is too far in the future");
      HdsUtils::setVanishTime(edge, std::nullopt);
      return { };
    }

    HdsUtils::setVanishTime(edge, vanish_time);
    return vanish_time;
  }

  // for pierce events
  std::optional<FT> crashAtTime(const VertexSPtr& vertex, const FacetSPtr& facet,
                                const std::optional<FT>& time_past_bound = std::nullopt,
                                const std::optional<FT>& time_future_bound = std::nullopt)
  {
    CGAL_SS3_CORE_TRACE_V(16, "-- Crash At Time (Pierce event)\n  " << vertex->toString() << "\n  " << facet->toString());

    CGAL_assertion(vertex->facets().size() >= 3);
    std::array<FacetSPtr, 3> fs;
    for (int i = 0; i < 3; ++i) {
      FacetWPtr wf = *(std::next(vertex->facets().begin(), i));
      CGAL_assertion(!wf.expired());
      fs[i] = wf.lock();
    }

    FT event_time = intersectionTimeOffsetPlanes(facet, fs[0], fs[1], fs[2]);

    if (time_past_bound && event_time >= *time_past_bound) {
      CGAL_SS3_TRAITS_TRACE("Contact event is strictly in the past");
      return { };
    }

    if (time_future_bound && event_time < *time_future_bound) {
      CGAL_SS3_TRAITS_TRACE("Contact event is too far in the future");
      return { };
    }

    CGAL_SS3_CORE_TRACE_V(16, "Tentative event @ " << event_time);

    return event_time;
  }

  std::optional<FT> crashAtTime(const EdgeSPtr& edge_1, const EdgeSPtr& edge_2,
                                const std::optional<FT>& time_past_bound = std::nullopt,
                                const std::optional<FT>& time_future_bound = std::nullopt)
  {
    CGAL_SS3_CORE_TRACE_V(16, "-- Crash At Time\n  " << edge_1->toString() << "\n  " << edge_2->toString());

    CGAL_SS3_DEBUG_SPTR(edge_1);
    CGAL_SS3_DEBUG_SPTR(edge_2);

    FacetSPtr facet_l1 = edge_1->getFacetL();
    FacetSPtr facet_r1 = edge_1->getFacetR();
    FacetSPtr facet_l2 = edge_2->getFacetL();
    FacetSPtr facet_r2 = edge_2->getFacetR();

    CGAL_SS3_CORE_TRACE_V(16, "Facet L1 = " << facet_l1->getID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet R1 = " << facet_r1->getID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet L2 = " << facet_l2->getID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet R2 = " << facet_r2->getID());

    // It is pointless to check for contact events that are farther in time than any of the two
    // involved edges as either they will be gone, or new events will be computed
    std::optional<FT> tight_future_bound = time_future_bound;
    for (const EdgeSPtr& edge : {edge_1, edge_2}) {
      const std::optional<FT>& vanish_time = HdsUtils::getVanishTime(edge);
      // In the current context, the future time bound is the last save event,
      // and vanish event times are only stored if they are earlier.
      CGAL_assertion(!vanish_time.has_value() || !time_future_bound.has_value() || *vanish_time > *time_future_bound);
      if (vanish_time.has_value()) {
        if (tight_future_bound.has_value()) {
          if (*vanish_time > *tight_future_bound) {
            tight_future_bound = vanish_time;
          }
        } else {
          tight_future_bound = vanish_time;
        }
      }
    }

    FT event_time = intersectionTimeOffsetPlanes(facet_l1, facet_r1, facet_l2, facet_r2);

    if (time_past_bound && event_time >= *time_past_bound) {
      CGAL_SS3_TRAITS_TRACE("Contact event is strictly in the past");
      return { };
    }

    if (tight_future_bound && event_time < *tight_future_bound) {
      CGAL_SS3_TRAITS_TRACE("Contact event is too far in the future");
      return { };
    }

    CGAL_SS3_CORE_TRACE_V(16, "Tentative event @ " << event_time);

    Point3SPtr point = intersectionPointOffsetPlanes(facet_l1, facet_r1, facet_l2, facet_r2);

    return event_time;
  }

  /**
    * Returns `true` if the event is in the past
    */
  static bool isEventInThePast(const AbstractEventSPtr& event,
                               const FT& current_time)
  {
    CGAL_SS3_DEBUG_SPTR(event);
    CGAL_precondition(event->isValid());
    return event->getTime() >= current_time;
  }

  /**
    * Returns `true` if the neighborhood of an event has changed.
    */
  static bool isEventObsolete(const AbstractEventSPtr& event)
  {
    CGAL_SS3_DEBUG_SPTR(event);
    CGAL_precondition(event->isValid());
    return event->isObsolete();
  }

  static bool isActualVertexEvent(const VertexEventSPtr& event)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "#######  Tentative Vertex Event  #######");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->getTime();

    VertexSPtr vertex_1 = event->getVertex1();
    VertexSPtr vertex_2 = event->getVertex2();
    FacetSPtr facet_1 = event->getFacet1();
    FacetSPtr facet_2 = event->getFacet2();

    // @todo avoid all this duplication...
    EdgeSPtr edge_11 = EdgeSPtr();
    for (EdgeWPtr edge_1_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge_1 = edge_1_wptr.lock()) {
        FacetSPtr facet_1l = edge_1->getFacetL();
        FacetSPtr facet_1r = edge_1->getFacetR();
        if ((facet_1l == facet_1 && facet_1r != facet_2) ||
            (facet_1r == facet_1 && facet_1l != facet_2)) {
          edge_11 = edge_1;
        }
      }
    }

    EdgeSPtr edge_22 = EdgeSPtr();
    for (EdgeWPtr edge_2_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge_2 = edge_2_wptr.lock()) {
        FacetSPtr facet_2l = edge_2->getFacetL();
        FacetSPtr facet_2r = edge_2->getFacetR();
        if ((facet_2l == facet_2 && facet_2r != facet_1) ||
            (facet_2r == facet_2 && facet_2l != facet_1)) {
          edge_22 = edge_2;
        }
      }
    }

    bool conv_split_event = false;
    FacetSPtr facet_1b = facet_2->next(vertex_1);
    FacetSPtr facet_2b = facet_1->next(vertex_2);
    EdgeSPtr edge_cur = edge_11->next(facet_1b);
    while (edge_cur != edge_11) {
      if ((edge_cur->getFacetL() == facet_1b && edge_cur->getFacetR() == facet_2b) ||
          (edge_cur->getFacetR() == facet_1b && edge_cur->getFacetL() == facet_2b)) {
        conv_split_event = true;
        break;
      }
      edge_cur = edge_cur->next(facet_1b);
    }
    if (conv_split_event) {
      CGAL_SS3_CORE_TRACE_V(8, "Vertex event: convex split detected");
      return false;
    }

    // Bisector check
    Point3SPtr point = intersectionPointOffsetPlanes(edge_11->getFacetL(), edge_11->getFacetR(),
                                                     edge_22->getFacetL(), edge_22->getFacetR());

    if (!checkBisectorsV2(edge_11, point, event_time) ||
        !checkBisectorsV2(edge_22, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Vertex event: bisector check failure");
      return false;
    }

    event->setPoint(point);

    CGAL_SS3_CORE_TRACE_V(8, "Vertex event: accepted");
    return true;
  }

  static bool isActualFlipVertexEvent(const FlipVertexEventSPtr& event)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "#####  Tentative Flip Vertex Event  ####");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->getTime();
    VertexSPtr vertex_1 = event->getVertex1();
    VertexSPtr vertex_2 = event->getVertex2();
    FacetSPtr facet_1 = event->getFacet1();
    FacetSPtr facet_2 = event->getFacet2();

    // convex split event checks
    EdgeSPtr edge_11 = EdgeSPtr();
    for (EdgeWPtr edge_1_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge_1 = edge_1_wptr.lock()) {
        FacetSPtr facet_1l = edge_1->getFacetL();
        FacetSPtr facet_1r = edge_1->getFacetR();
        if ((facet_1l == facet_1 && facet_1r != facet_2) ||
            (facet_1r == facet_1 && facet_1l != facet_2)) {
          edge_11 = edge_1;
        }
      }
    }

    EdgeSPtr edge_22 = EdgeSPtr();
    for (EdgeWPtr edge_2_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge_2 = edge_2_wptr.lock()) {
        FacetSPtr facet_2l = edge_2->getFacetL();
        FacetSPtr facet_2r = edge_2->getFacetR();
        if ((facet_2l == facet_2 && facet_2r != facet_1) ||
            (facet_2r == facet_2 && facet_2l != facet_1)) {
          edge_22 = edge_2;
        }
      }
    }

    bool conv_split_event = false;
    FacetSPtr facet_1b = facet_2->next(vertex_1);
    FacetSPtr facet_2b = facet_2->next(vertex_2);
    EdgeSPtr edge_cur = edge_11->next(facet_1b);
    while (edge_cur != edge_11) {
      if ((edge_cur->getFacetL() == facet_1b && edge_cur->getFacetR() == facet_2b) ||
          (edge_cur->getFacetR() == facet_1b && edge_cur->getFacetL() == facet_2b)) {
        conv_split_event = true;
        break;
      }
      edge_cur = edge_cur->next(facet_1b);
    }
    if (conv_split_event) {
      CGAL_SS3_CORE_TRACE_V(8, "Flip vertex event: convex split detected");
      return false;
    }

    // Bisector check
    Point3SPtr point = intersectionPointOffsetPlanes(edge_11->getFacetL(), edge_11->getFacetR(),
                                                     edge_22->getFacetL(), edge_22->getFacetR());

    if (!checkBisectorsV2(edge_11, point, event_time) ||
        !checkBisectorsV2(edge_22, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Flip Vertex event: bisector check failure");
      return false;
    }

    event->setPoint(point);

    CGAL_SS3_CORE_TRACE_V(8, "Flip vertex event: accepted");
    return true;
  }

  static bool isActualSurfaceEvent(const SurfaceEventSPtr& event)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "######  Tentative Surface Event  #######");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->getTime();
    EdgeSPtr edge_1 = event->getEdge1();
    EdgeSPtr edge_2 = event->getEdge2();
    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();

    // convex split checks
    bool conv_split_event = false;
    std::list<EdgeSPtr> common_edges = edge_1->getFacetL()->findEdges(edge_1->getFacetR());
    for (const EdgeSPtr& edge : common_edges) {
      if (edge == edge_1) {
        continue;
      }
      FacetSPtr facet_src = edge->getFacetSrc();
      FacetSPtr facet_dst = edge->getFacetDst();
      if (facet_1_src == edge_2->getFacetL() ||
          facet_1_dst == edge_2->getFacetL()) {
        if (facet_src == edge_2->getFacetR() ||
            facet_dst == edge_2->getFacetR()) {
          conv_split_event = true;
          break;
        }
      } else if (facet_1_src == edge_2->getFacetR() ||
                 facet_1_dst == edge_2->getFacetR()) {
        if (facet_src == edge_2->getFacetL() ||
            facet_dst == edge_2->getFacetL()) {
          conv_split_event = true;
          break;
        }
      }
    }

    if (conv_split_event) {
      CGAL_SS3_CORE_TRACE_V(8, "Surface event: convex split detected");
      return false;
    }

    // Bisector check
    Point3SPtr point = intersectionPointOffsetPlanes(edge_1->getFacetL(), edge_1->getFacetR(),
                                                     edge_2->getFacetL(), edge_2->getFacetR());

    if (!checkBisectorsV2(edge_1, point, event_time) ||
        !checkBisectorsV2(edge_2, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Surface event: bisector check failure");
      return false;
    }

    event->setPoint(point);

    CGAL_SS3_CORE_TRACE_V(8, "Surface event: accepted");
    return true;
  }

  static bool isActualPolyhedronSplitEvent(const PolyhedronSplitEventSPtr& event,
                                           const FT& current_time)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "####  Tentative Split Merge Event  #####");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->getTime();

    // Bisector check
    EdgeSPtr edge_1 = event->getEdge1();
    EdgeSPtr edge_2 = event->getEdge2();

    Point3SPtr point = intersectionPointOffsetPlanes(edge_1->getFacetL(), edge_1->getFacetR(),
                                                     edge_2->getFacetL(), edge_2->getFacetR());

    if (!checkBisectorsV2(edge_1, point, event_time) ||
        !checkBisectorsV2(edge_2, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Polyhedron split event: bisector check failure");
      return false;
    }

    event->setPoint(point);

    // @speed is_degenerate(4 planes)? But it would be a sure filter failure, so, costly...
    FT shift = event_time - current_time;
    Segment3SPtr e1o = Transformation::shiftEdge(event->getEdge1(), shift);
    if (!e1o->is_degenerate()) {
      CGAL_SS3_CORE_TRACE_V(8, "Polyhedron split event: bisector check failure");
      return false;
    }

    CGAL_SS3_CORE_TRACE_V(8, "Polyhedron split event accepted");
    return true;
  }

  static bool isActualSplitMergeEvent(const SplitMergeEventSPtr& event)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "####  Tentative Split Merge Event  #####");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->getTime();
    VertexSPtr vertex_1 = event->getVertex1();
    VertexSPtr vertex_2 = event->getVertex2();
    FacetSPtr facet_1 = event->getFacet1();
    FacetSPtr facet_2 = event->getFacet2();

    // convex split checks
    EdgeSPtr edge_11 = EdgeSPtr();
    for (EdgeWPtr edge_1_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge_1 = edge_1_wptr.lock()) {
        FacetSPtr facet_1l = edge_1->getFacetL();
        FacetSPtr facet_1r = edge_1->getFacetR();
        if ((facet_1l == facet_1 && facet_1r != facet_2) ||
            (facet_1r == facet_1 && facet_1l != facet_2)) {
          edge_11 = edge_1;
        }
      }
    }

    EdgeSPtr edge_22 = EdgeSPtr();
    for (EdgeWPtr edge_2_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge_2 = edge_2_wptr.lock()) {
        FacetSPtr facet_2l = edge_2->getFacetL();
        FacetSPtr facet_2r = edge_2->getFacetR();
        if ((facet_2l == facet_2 && facet_2r != facet_1) ||
            (facet_2r == facet_2 && facet_2l != facet_1)) {
          edge_22 = edge_2;
        }
      }
    }

    bool conv_split_event = false;
    FacetSPtr facet_1b = facet_2->next(vertex_1);
    FacetSPtr facet_2b = facet_1->next(vertex_2);
    if (facet_2b == facet_2) {
      facet_2b = facet_2b->next(vertex_2);
    }

    EdgeSPtr edge_cur = edge_11->next(facet_1b);
    while (edge_cur != edge_11) {
      if ((edge_cur->getFacetL() == facet_1b && edge_cur->getFacetR() == facet_2b) ||
          (edge_cur->getFacetR() == facet_1b && edge_cur->getFacetL() == facet_2b)) {
        conv_split_event = true;
        break;
      }
      edge_cur = edge_cur->next(facet_1b);
    }
    if (!conv_split_event) {
      CGAL_SS3_CORE_TRACE_V(8, "Split merge event: Convex split event detected");
      return false;
    }

    // Bisector check
    Point3SPtr point = intersectionPointOffsetPlanes(edge_11->getFacetL(), edge_11->getFacetR(),
                                                     edge_22->getFacetL(), edge_22->getFacetR());

    if (!checkBisectorsV2(edge_11, point, event_time) ||
        !checkBisectorsV2(edge_22, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Split merge event: bisector check failure");
      return false;
    }

    event->setPoint(point);

    CGAL_SS3_CORE_TRACE_V(8, "Split merge event accepted");
    return true;
  }

  static bool isActualEdgeSplitEvent(const EdgeSplitEventSPtr& event)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "#####  Tentative Edge Split Event  #####");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->getTime();
    EdgeSPtr edge_1 = event->getEdge1(); // @todo const& all of this
    EdgeSPtr edge_2 = event->getEdge2();

    // Bisector check
    Point3SPtr point = intersectionPointOffsetPlanes(edge_1->getFacetL(), edge_1->getFacetR(),
                                                     edge_2->getFacetL(), edge_2->getFacetR());
    CGAL_SS3_DEBUG_SPTR(point);

    if (!checkBisectorsV2(edge_1, point, event_time) ||
        !checkBisectorsV2(edge_2, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Edge split rejected at pop time");
      return false;
    }

    event->setPoint(point);

    CGAL_SS3_CORE_TRACE_V(8, "Edge split event accepted");
    return true;
  }

  static bool isActualPierceEvent(const PierceEventSPtr& event,
                                  const FT& current_time,
                                  const std::optional<FT>& time_future_bound)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "######  Tentative Pierce Event  ########");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    VertexSPtr pv = event->getVertex();
    FacetSPtr pf = event->getFacet();

    // Filter #1
    //
    // This combinatorics filter is here because some events such as surface events
    // can increase the (combinatorial) separation between a vertex and a facet, so a pierce
    // event can be revealed without anything changing around the vertex.
    //
    // This is problematic if we are using local updates because we (currently) only check
    // modified vertices. We could increase the range of vertices being considered in local
    // updates, but it's simpler to delay the check until the event has become the next event
    // to treat.
    bool has_edge_to_facet = false;
    for (EdgeWPtr edge_wptr : pv->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        FacetSPtr facet_src = edge->getFacetSrc();
        FacetSPtr facet_dst = edge->getFacetDst();
        if (pf == facet_src || pf == facet_dst) {
          has_edge_to_facet = true;
          break;
        }
      }
    }
    if (has_edge_to_facet) {
      CGAL_SS3_CORE_TRACE_V(8, "Pierce event rejected at pop time (a)");
      return false;
    }

    // Filter #2
    //
    // Only done here because pierce events are asymetrical (vertex <-> facet)
    // and if it were done at collect time, we would have to seek pierce events
    // when facets get modified - and not only when reflex vertices are modified -,
    // and that seems more costly (and complicated).
    if (time_future_bound.has_value()) {
      CGAL::Bbox_3 b1;
      b1 += pv->getPoint()->bbox();
      b1 += HdsUtils::getFinalPoint(pv, *time_future_bound)->bbox();

      CGAL::Bbox_3 b2;
      for (const VertexSPtr& v : pf->vertices()) {
        b2 += v->getPoint()->bbox();
        b2 += HdsUtils::getFinalPoint(v, *time_future_bound)->bbox();
      }

      if (!CGAL::do_overlap(b1, b2)) {
        CGAL_SS3_CORE_TRACE_V(8, "Pierce event rejected at pop time (b)");
        return false;
      }
    }

    // Compute the event position
    CGAL_assertion(pv->facets().size() >= 3);
    std::array<FacetSPtr, 3> fs;
    for (int i = 0; i < 3; ++i) {
      FacetWPtr wf = *(std::next(pv->facets().begin(), i));
      CGAL_assertion(!wf.expired());
      fs[i] = wf.lock();
    }

    event->setPoint(intersectionPointOffsetPlanes(pf, fs[0], fs[1], fs[2]));

    // Filter #3
    //
    // Filter if the event point is on an edge (and a fortiori on a vertex)
    // as it will be a different kind of event
    Point3SPtr point = event->getPoint();
    FacetSPtr facet_clone = pf->clone();

    FT shift = event->getTime() - current_time;
    const FT& speed = HdsUtils::getSpeed(pf);
    Plane3SPtr offset_plane = GeomUtils::offsetPlane(pf->getPlane(), shift*speed);
    facet_clone->setPlane(offset_plane);

    // abusing the fact that vertices will have the same order in both facets
    typename std::list<VertexSPtr>::iterator it_v = pf->vertices().begin();
    typename std::list<VertexSPtr>::iterator it_v_offset = facet_clone->vertices().begin();
    while (it_v != pf->vertices().end()) {
      VertexSPtr vertex = *it_v++;
      VertexSPtr offset_vertex = *it_v_offset++;
      Point3SPtr point_offset = Transformation::shiftPoint(vertex, shift);
      offset_vertex->setPoint(point_offset);
    }

#ifdef CGAL_SLS3_NEW_IS_INSIDE
    if (!SelfIntersection::isInsideWithRayShootingV2(point, facet_clone)) {
      CGAL_SS3_CORE_TRACE_V(8, "Pierce event rejected at pop time (c)");
      return false;
    }
#else
    if (!SelfIntersection::isInsideWithRayShooting(point, facet_clone)) {
      CGAL_SS3_CORE_TRACE_V(8, "Pierce event rejected at pop time (c)");
      return false;
    }

    bool boundary_rejection = false;
    for (const EdgeSPtr& edge : facet_clone->edges()) {
      Segment3SPtr seg = KernelFactory::createSegment3(edge->getVertexSrc()->getPoint(),
                                                       edge->getVertexDst()->getPoint());
      if (!seg || seg->is_degenerate()) {
        continue;
      }

      if (seg->has_on(*point)) {
        boundary_rejection = true;
        break;
      }
    }

    if (boundary_rejection) {
      CGAL_SS3_CORE_TRACE_V(8, "Pierce event rejected at pop time (d)");
      return false;
    }
#endif

    CGAL_SS3_CORE_TRACE_V(8, "Pierce event accepted");
    return true;
  }

  /**
    * Some combinatorial or geometric checks are very expensive to perform,
    * so delay them until the event is the best in the queue.
    * The gain is that the event can be invalidated combinatorially by the time
    * it gets popped.
    */
  static bool isActualEvent(const AbstractEventSPtr& event,
                            const FT& current_time,
                            const std::optional<FT>& time_future_bound)
  {
    CGAL_SS3_DEBUG_SPTR(event);
    CGAL_precondition(event->isValid());

    bool result = true;

    if (event->getType() == AbstractEvent::VERTEX_EVENT) {
      result = isActualVertexEvent(std::dynamic_pointer_cast<VertexEvent>(event));
    } else if (event->getType() == AbstractEvent::FLIP_VERTEX_EVENT) {
      result = isActualFlipVertexEvent(std::dynamic_pointer_cast<FlipVertexEvent>(event));
    } else if (event->getType() == AbstractEvent::SURFACE_EVENT) {
      result = isActualSurfaceEvent(std::dynamic_pointer_cast<SurfaceEvent>(event));
    } else if (event->getType() == AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
      result = isActualPolyhedronSplitEvent(std::dynamic_pointer_cast<PolyhedronSplitEvent>(event), current_time);
    } else if (event->getType() == AbstractEvent::SPLIT_MERGE_EVENT) {
      result = isActualSplitMergeEvent(std::dynamic_pointer_cast<SplitMergeEvent>(event));
    } else if (event->getType() == AbstractEvent::EDGE_SPLIT_EVENT) {
      result = isActualEdgeSplitEvent(std::dynamic_pointer_cast<EdgeSplitEvent>(event));
    } else if (event->getType() == AbstractEvent::PIERCE_EVENT) {
      result = isActualPierceEvent(std::dynamic_pointer_cast<PierceEvent>(event), current_time, time_future_bound);
    }

    return result;
  }

  /**
    * Vanish events.
    */
  void collectVanishEvents(const std::list<EdgeSPtr>& edges,
                           const PolyhedronSPtr& /*polyhedron*/,
                           const FT& current_time,
                           const std::optional<FT>& time_future_bound,
                           PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Vanish Events [" << current_time << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    for (const EdgeSPtr& edge : edges) {
      CGAL_SS3_DEBUG_SPTR(edge);

      VertexSPtr vertex_src = edge->getVertexSrc();
      VertexSPtr vertex_dst = edge->getVertexDst();
      if (vertex_src->getPoint() == vertex_dst->getPoint()) {
        HdsUtils::setVanishTime(edge, std::nullopt);
        continue;
      }

      std::optional<FT> event_time = vanishesAtTime(edge, current_time, time_future_bound);
      if (!event_time) {
        continue;
      }

      CGAL_assertion(*event_time < current_time);
      CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

      VanishEventSPtr event = VanishEvent::create();
      event->setTime(*event_time);
      event->setPoint(vanishesAtPoint(edge));
      event->setEdge(edge);
      queue.push(event);
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Vanish Events in: " << timer.time());
#endif
  }

  void collectVanishEvents(const PolyhedronSPtr& polyhedron,
                           const FT& current_time,
                           const std::optional<FT>& time_future_bound,
                           PQ& queue)
  {
    return collectVanishEvents(polyhedron->edges(), polyhedron, current_time, time_future_bound, queue);
  }

  /**
    * Two vertices crash into each other.
    */
  void collectVertexEvents(const std::list<VertexSPtr>& vertices,
                           const PolyhedronSPtr& /*polyhedron*/,
                           const bool use_canonical_event_reps,
                           const FT& current_time,
                           const std::optional<FT>& time_future_bound,
                           PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Vertex Events [" << current_time << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    for (const VertexSPtr& vertex_1 : vertices) {
      CGAL_SS3_DEBUG_SPTR(vertex_1);

      if (HdsUtils::isConvex(vertex_1)) {
        continue;
      }

      std::set<VertexSPtr> vertices_2;
      for (FacetWPtr facet_wptr : vertex_1->facets()) {
        if (FacetSPtr facet = facet_wptr.lock()) {
            vertices_2.insert(facet->vertices().begin(), facet->vertices().end());
            // @speed sort vertices_2 as to get the first vertex after start_e->dst(facet)
        }
      }

      for (const VertexSPtr& vertex_2 : vertices_2) {
        CGAL_SS3_DEBUG_SPTR(vertex_2);
        if (vertex_1 == vertex_2) {
          continue;
        }
#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        if (use_canonical_event_reps) {
          CGAL_assertion(vertex_1->getID() != -1 && vertex_2->getID() != -1);
          if (vertex_1->getID() > vertex_2->getID()) {
            continue;
          }
        }
#endif
        if (vertex_1->getPoint() == vertex_2->getPoint()) {
          continue;
        }
        if (vertex_1->findEdge(vertex_2)) {
          // edge event
          continue;
        }
        if (HdsUtils::isConvex(vertex_2)) {
          continue;
        }

        // Subtlety here: this event is not symmetrical because the two chosen edges
        // incident to vertex_1 and vertex_2 depend on the respective order of the vertices.
        // Afterwards, we will check the validity of a potential intersection point
        // with respect to these edges, but not the other potential pair.
        // However, one pair could be valid while the other one is not.
        // Thus, whether we look at vertex_1-vertex_2 or vertex_2-vertex_1, we could
        // get told that the event exists, or that it does not.
        // But, this is not a real inconsistency: if the event exists for the current order
        // but does not for the other one, it's because there is another event (e.g. a simple
        // edge event) that prevents this event from actually existing.
        // As such, it does not really matter that we do not see that the event does not in fact
        // exist, because even if it gets put in the queue, it will be invalidated
        // by the nearer events.
        //
        // The point of the swap() below is to ensure consistency whether we are filling
        // a global queue or a local queue, because otherwise we can get an error in due to the
        // asymmetry: for example, the event does not exist in the local queue (as vertex_2 -
        // vertex_1), but exists in the global queue (as vertex_1 - vertex_2).
        // It is a false positive in the consistency check, but still, might as well ensure
        // consistency.
        VertexSPtr v1 = vertex_1;
        VertexSPtr v2 = vertex_2;
        if (v1->getID() > v2->getID()) {
          std::swap(v1, v2);
        }

        FacetSPtr facet_1;
        FacetSPtr facet_2;
        int num_equal_facets = 0;
        for (FacetWPtr facet_1_wptr : v1->facets()) {
          if (FacetSPtr f1 = facet_1_wptr.lock()) {
            for (FacetWPtr facet_2_wptr : v2->facets()) {
              if (FacetSPtr f2 = facet_2_wptr.lock()) {
                if (f1 == f2) {
                  if (num_equal_facets == 0) {
                    facet_1 = f1;
                  } else {
                    facet_2 = f2;
                  }
                  ++num_equal_facets;
                }
              }
            }
          }
        }
        if (num_equal_facets != 2) {
          continue;
        }
        if (facet_1->next(v1) != facet_2) {
          FacetSPtr facet_tmp = facet_1;
          facet_1 = facet_2;
          facet_2 = facet_tmp;
        }
        if (v1->next(facet_1)->next(facet_1) == v2 ||
            v1->next(facet_2)->next(facet_2) == v2 ||
            v1->prev(facet_1)->prev(facet_1) == v2 ||
            v1->prev(facet_2)->prev(facet_2) == v2) {
          // edge merge event
          continue;
        }

        EdgeSPtr edge_11 = EdgeSPtr();
        EdgeSPtr edge_12 = EdgeSPtr();
        for (EdgeWPtr edge_1_wptr : v1->edges()) {
          if (EdgeSPtr edge_1 = edge_1_wptr.lock()) {
            FacetSPtr facet_1l = edge_1->getFacetL();
            FacetSPtr facet_1r = edge_1->getFacetR();
            if ((facet_1l == facet_1 && facet_1r != facet_2) ||
                (facet_1r == facet_1 && facet_1l != facet_2)) {
              edge_11 = edge_1;
            } else if ((facet_1l == facet_2 && facet_1r != facet_1) ||
                       (facet_1r == facet_2 && facet_1l != facet_1)) {
              edge_12 = edge_1;
            }
          }
        }
        EdgeSPtr edge_21 = EdgeSPtr();
        EdgeSPtr edge_22 = EdgeSPtr();
        for (EdgeWPtr edge_2_wptr : v2->edges()) {
          if (EdgeSPtr edge_2 = edge_2_wptr.lock()) {
            FacetSPtr facet_2l = edge_2->getFacetL();
            FacetSPtr facet_2r = edge_2->getFacetR();
            if ((facet_2l == facet_1 && facet_2r != facet_2) ||
                (facet_2r == facet_1 && facet_2l != facet_2)) {
              edge_21 = edge_2;
            } else if ((facet_2l == facet_2 && facet_2r != facet_1) ||
                       (facet_2r == facet_2 && facet_2l != facet_1)) {
              edge_22 = edge_2;
            }
          }
        }
        if (!((edge_11->next(v1) == edge_12 && edge_22->next(v2) == edge_21) ||
              (edge_12->next(v1) == edge_11 && edge_21->next(v2) == edge_22))) {
          // flip vertex event
          continue;
        }

        // convex split event checks are performed at pop time - see isActualVertexEvent()

        std::optional<FT> event_time = crashAtTime(edge_11, edge_22, current_time, time_future_bound);
        if (!event_time) {
          continue;
        }

        CGAL_assertion(*event_time < current_time);
        CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

        VertexEventSPtr event = VertexEvent::create();
        event->setTime(*event_time);
        event->setVertex1(v1);
        event->setVertex2(v2);
        event->setFacet1(facet_1);
        event->setFacet2(facet_2);
        queue.push(event);
      }
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Vertex Events in: " << timer.time());
#endif
  }

  void collectVertexEvents(const PolyhedronSPtr& polyhedron,
                           const FT& current_time,
                           const std::optional<FT>& time_future_bound,
                           PQ& queue)
  {
    return collectVertexEvents(polyhedron->vertices(), polyhedron, true /*use canonical reps*/,
                               current_time, time_future_bound, queue);
  }

  /**
    * Flip vertex event
    */
  void collectFlipVertexEvents(const std::list<VertexSPtr>& vertices,
                               const PolyhedronSPtr& /*polyhedron*/,
                               const bool use_canonical_event_reps,
                               const FT& current_time,
                               const std::optional<FT>& time_future_bound,
                               PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Flip Vertex Events [" << current_time << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    for (const VertexSPtr& vertex_1 : vertices) {
      CGAL_assertion(vertex_1->getID() != -1);

      if (HdsUtils::isConvex(vertex_1)) {
        continue;
      }

      std::set<VertexSPtr> vertices_2;
      for (FacetWPtr facet_wptr : vertex_1->facets()) {
        if (FacetSPtr facet = facet_wptr.lock()) {
          vertices_2.insert(facet->vertices().begin(), facet->vertices().end());
        }
      }

      for (const VertexSPtr& vertex_2 : vertices_2) {
        CGAL_assertion(vertex_2->getID() != -1);

        if (vertex_1 == vertex_2) {
          continue;
        }
#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        if (use_canonical_event_reps) {
          if (vertex_1->getID() > vertex_2->getID()) {
            continue;
          }
        }
#endif
        if (vertex_1->getPoint() == vertex_2->getPoint()) {
          continue;
        }
        if (vertex_1->findEdge(vertex_2)) {
          // edge event
          continue;
        }
        if (HdsUtils::isConvex(vertex_2)) {
          continue;
        }

        // See comment in the first instance of this swap
        VertexSPtr v1 = vertex_1;
        VertexSPtr v2 = vertex_2;
        if (v1->getID() > v2->getID()) {
          std::swap(v1, v2);
        }

        FacetSPtr facet_1;
        FacetSPtr facet_2;
        int num_equal_facets = 0;
        for (FacetWPtr facet_1_wptr : v1->facets()) {
          if (FacetSPtr f1 = facet_1_wptr.lock()) {
            for (FacetWPtr facet_2_wptr : v2->facets()) {
              if (FacetSPtr f2 = facet_2_wptr.lock()) {
                if (f1 == f2) {
                  if (num_equal_facets == 0) {
                    facet_1 = f1;
                  } else {
                    facet_2 = f2;
                  }
                  ++num_equal_facets;
                }
              }
            }
          }
        }
        if (num_equal_facets != 2) {
          continue;
        }
#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        if (use_canonical_event_reps) {
          if (facet_1->getID() > facet_2->getID()) {
            continue;
          }
        }
#endif
        if (facet_1->next(v1) != facet_2) {
          FacetSPtr facet_tmp = facet_1;
          facet_1 = facet_2;
          facet_2 = facet_tmp;
        }
        if (v1->next(facet_1)->next(facet_1) == v2 ||
            v1->next(facet_2)->next(facet_2) == v2 ||
            v1->prev(facet_1)->prev(facet_1) == v2 ||
            v1->prev(facet_2)->prev(facet_2) == v2) {
          // edge merge event
          continue;
        }

        EdgeSPtr edge_11 = EdgeSPtr();
        EdgeSPtr edge_12 = EdgeSPtr();
        for (EdgeWPtr edge_1_wptr : v1->edges()) {
          if (EdgeSPtr edge_1 = edge_1_wptr.lock()) {
            FacetSPtr facet_1l = edge_1->getFacetL();
            FacetSPtr facet_1r = edge_1->getFacetR();
            if ((facet_1l == facet_1 && facet_1r != facet_2) ||
                (facet_1r == facet_1 && facet_1l != facet_2)) {
              edge_11 = edge_1;
            } else if ((facet_1l == facet_2 && facet_1r != facet_1) ||
                       (facet_1r == facet_2 && facet_1l != facet_1)) {
              edge_12 = edge_1;
            }
          }
        }
        EdgeSPtr edge_21 = EdgeSPtr();
        EdgeSPtr edge_22 = EdgeSPtr();
        for (EdgeWPtr edge_2_wptr : v2->edges()) {
          if (EdgeSPtr edge_2 = edge_2_wptr.lock()) {
            FacetSPtr facet_2l = edge_2->getFacetL();
            FacetSPtr facet_2r = edge_2->getFacetR();
            if ((facet_2l == facet_1 && facet_2r != facet_2) ||
                (facet_2r == facet_1 && facet_2l != facet_2)) {
              edge_21 = edge_2;
            } else if ((facet_2l == facet_2 && facet_2r != facet_1) ||
                       (facet_2r == facet_2 && facet_2l != facet_1)) {
              edge_22 = edge_2;
            }
          }
        }
        if (!(edge_12->next(v1) == edge_11 && edge_22->next(v2) == edge_21)) {
          // vertex event
          continue;
        }

        // convex split event checks are performed at pop time - see isActualFlipVertexEvent()

        std::optional<FT> event_time = crashAtTime(edge_11, edge_22, current_time, time_future_bound);
        if (!event_time) {
          continue;
        }

        CGAL_assertion(*event_time < current_time);
        CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

        FlipVertexEventSPtr event = FlipVertexEvent::create();
        event->setTime(*event_time);
        event->setVertex1(v1);
        event->setVertex2(v2);
        event->setFacet1(facet_1);
        event->setFacet2(facet_2);
        queue.push(event);
      }
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Flip Vertex Events in: " << timer.time());
#endif
  }

  void collectFlipVertexEvents(const PolyhedronSPtr& polyhedron,
                               const FT& current_time,
                               const std::optional<FT>& time_future_bound,
                               PQ& queue)
  {
    return collectFlipVertexEvents(polyhedron->vertices(), polyhedron, true /*use canonical reps*/,
                                   current_time, time_future_bound, queue);
  }

  /**
    * Split event on the surface.
    * Edges do not need to be reflex.
    */
  void collectSurfaceEvent(const EdgeSPtr& edge_1,
                           const EdgeSPtr& edge_2,
                           const PolyhedronSPtr& /*polyhedron*/,
                           const FT& current_time,
                           const std::optional<FT>& time_future_bound,
                           PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(8, ">>> Collect Surface Event [\n  " << edge_1->toString() << "\n  " << edge_2->toString() << "]");

    CGAL_SS3_DEBUG_SPTR(edge_1);
    CGAL_SS3_DEBUG_SPTR(edge_2);

    if (edge_1 == edge_2) {
      return;
    }

    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();

    if (edge_1->getFacetL() == edge_2->getFacetL() ||
        edge_1->getFacetL() == edge_2->getFacetR() ||
        edge_1->getFacetR() == edge_2->getFacetL() ||
        edge_1->getFacetR() == edge_2->getFacetR()) {
      // on same facet
      return;
    }

    if (edge_1->getVertexSrc()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
        edge_1->getVertexSrc()->getPoint() == edge_2->getVertexDst()->getPoint() ||
        edge_1->getVertexDst()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
        edge_1->getVertexDst()->getPoint() == edge_2->getVertexDst()->getPoint()) {
      // share a vertex
      return;
    }

    // vertex of edge_1 splits edge_2
    if (!((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() != facet_1_dst) ||
          (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() != facet_1_src) ||
          (edge_2->getFacetR() == facet_1_src && edge_2->getFacetL() != facet_1_dst) ||
          (edge_2->getFacetR() == facet_1_dst && edge_2->getFacetL() != facet_1_src))) {
      // no surface event
      return;
    }

    FacetSPtr facet_2_src = edge_2->getFacetSrc();
    FacetSPtr facet_2_dst = edge_2->getFacetDst();
    if ((edge_1->getFacetL() == facet_2_src && edge_1->getFacetR() != facet_2_dst) ||
        (edge_1->getFacetL() == facet_2_dst && edge_1->getFacetR() != facet_2_src) ||
        (edge_1->getFacetR() == facet_2_src && edge_1->getFacetL() != facet_2_dst) ||
        (edge_1->getFacetR() == facet_2_dst && edge_1->getFacetL() != facet_2_src)) {
      // flip vertex event
      return;
    }

    if (edge_1->getVertexSrc()->findEdge(edge_2->getVertexSrc()) ||
        edge_1->getVertexSrc()->findEdge(edge_2->getVertexDst()) ||
        edge_1->getVertexDst()->findEdge(edge_2->getVertexSrc()) ||
        edge_1->getVertexDst()->findEdge(edge_2->getVertexDst()) ) {
      // edge event (when a pyramid grows outwards)
      // a surface split is not possible with only one edge in between
      return;
    }

    if ((edge_1->getFacetL() == facet_2_src && facet_1_src == edge_2->getFacetL()) ||
        (edge_1->getFacetL() == facet_2_dst && facet_1_src == edge_2->getFacetR()) ||
        (edge_1->getFacetR() == facet_2_src && facet_1_dst == edge_2->getFacetL()) ||
        (edge_1->getFacetR() == facet_2_dst && facet_1_dst == edge_2->getFacetR()) ||
        (edge_1->getFacetR() == facet_2_src && facet_1_src == edge_2->getFacetR()) ||
        (edge_1->getFacetR() == facet_2_dst && facet_1_src == edge_2->getFacetL()) ||
        (edge_1->getFacetL() == facet_2_src && facet_1_dst == edge_2->getFacetR()) ||
        (edge_1->getFacetL() == facet_2_dst && facet_1_dst == edge_2->getFacetL())) {
      // vertex event
      return;
    }

    // convex split event are performed at pop time - see isActualSurfaceEvent()

    // let's just check if bboxes overlap first
    if (time_future_bound.has_value()) {
      CGAL::Bbox_3 b1;
      b1 += edge_1->getVertexSrc()->getPoint()->bbox();
      b1 += edge_1->getVertexDst()->getPoint()->bbox();
      b1 += HdsUtils::getFinalPoint(edge_1->getVertexSrc(), *time_future_bound)->bbox();
      b1 += HdsUtils::getFinalPoint(edge_1->getVertexDst(), *time_future_bound)->bbox();

      CGAL::Bbox_3 b2;
      b2 += edge_2->getVertexSrc()->getPoint()->bbox();
      b2 += edge_2->getVertexDst()->getPoint()->bbox();
      b2 += HdsUtils::getFinalPoint(edge_2->getVertexSrc(), *time_future_bound)->bbox();
      b2 += HdsUtils::getFinalPoint(edge_2->getVertexDst(), *time_future_bound)->bbox();

      if (!CGAL::do_overlap(b1, b2)) {
        CGAL_SS3_CORE_TRACE_V(32, "Filtered possible surface event candidates\n\t" << edge_1->toString() << "\n\t"
                                                                                   << edge_2->toString());
        return;
      } else {
        CGAL_SS3_CORE_TRACE_V(32, "Checking possible surface event event\n\t" << edge_1->toString() << "\n\t"
                                                                              << edge_2->toString());
      }
    }

    std::optional<FT> event_time = crashAtTime(edge_1, edge_2, current_time, time_future_bound);
    if (!event_time) {
      return;
    }

    CGAL_assertion(*event_time < current_time);
    CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

    SurfaceEventSPtr event = SurfaceEvent::create();
    event->setTime(*event_time);
    event->setEdge1(edge_1);
    event->setEdge2(edge_2);
    queue.push(event);
  }

  void collectSurfaceEvents(const std::list<EdgeSPtr>& edges,
                            const PolyhedronSPtr& /*polyhedron*/,
                            const FT& current_time,
                            const std::optional<FT>& time_future_bound,
                            PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Surface Events [" << current_time << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    unsigned int tested_candidates = 0;
    unsigned int total_candidates = 0;
#endif

    for (const EdgeSPtr& edge_1 : edges) {
      CGAL_SS3_DEBUG_SPTR(edge_1);

      FacetSPtr facet_1_src = edge_1->getFacetSrc();
      FacetSPtr facet_1_dst = edge_1->getFacetDst();
      std::list<EdgeSPtr> edges_2;
      edges_2.insert(edges_2.end(), facet_1_src->edges().begin(), facet_1_src->edges().end());
      edges_2.insert(edges_2.end(), facet_1_dst->edges().begin(), facet_1_dst->edges().end());

      // not outside of the loop just because maybe one day this will be called
      // as the first collect function with an initial bound that gets updated...
      CGAL::Bbox_3 b1;
      if (time_future_bound.has_value()) {
        b1 += edge_1->getVertexSrc()->getPoint()->bbox();
        b1 += edge_1->getVertexDst()->getPoint()->bbox();
        b1 += HdsUtils::getFinalPoint(edge_1->getVertexSrc(), *time_future_bound)->bbox();
        b1 += HdsUtils::getFinalPoint(edge_1->getVertexDst(), *time_future_bound)->bbox();
      }

      for (const EdgeSPtr& edge_2 : edges_2) {
        CGAL_SS3_DEBUG_SPTR(edge_2);

        if (edge_1 == edge_2) {
          continue;
        }

        CGAL_SS3_CORE_TRACE_V(64, "Possible surface event:\n\t" << edge_1->toString() << "\n\t"
                                                                << edge_2->toString());

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
        ++total_candidates;
#endif

        if (edge_1->getFacetL() == edge_2->getFacetL() ||
            edge_1->getFacetL() == edge_2->getFacetR() ||
            edge_1->getFacetR() == edge_2->getFacetL() ||
            edge_1->getFacetR() == edge_2->getFacetR()) {
          // on same facet
          continue;
        }

        if (edge_1->getVertexSrc()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
            edge_1->getVertexSrc()->getPoint() == edge_2->getVertexDst()->getPoint() ||
            edge_1->getVertexDst()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
            edge_1->getVertexDst()->getPoint() == edge_2->getVertexDst()->getPoint()) {
          // share a vertex
          continue;
        }

        // vertex of edge_1 splits edge_2
        if (!((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() != facet_1_dst) ||
              (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() != facet_1_src) ||
              (edge_2->getFacetR() == facet_1_src && edge_2->getFacetL() != facet_1_dst) ||
              (edge_2->getFacetR() == facet_1_dst && edge_2->getFacetL() != facet_1_src))) {
          // no surface event
          continue;
        }

        FacetSPtr facet_2_src = edge_2->getFacetSrc();
        FacetSPtr facet_2_dst = edge_2->getFacetDst();
        if ((edge_1->getFacetL() == facet_2_src && edge_1->getFacetR() != facet_2_dst) ||
            (edge_1->getFacetL() == facet_2_dst && edge_1->getFacetR() != facet_2_src) ||
            (edge_1->getFacetR() == facet_2_src && edge_1->getFacetL() != facet_2_dst) ||
            (edge_1->getFacetR() == facet_2_dst && edge_1->getFacetL() != facet_2_src)) {
          // flip vertex event
          continue;
        }

        if (edge_1->getVertexSrc()->findEdge(edge_2->getVertexSrc()) ||
            edge_1->getVertexSrc()->findEdge(edge_2->getVertexDst()) ||
            edge_1->getVertexDst()->findEdge(edge_2->getVertexSrc()) ||
            edge_1->getVertexDst()->findEdge(edge_2->getVertexDst()) ) {
          // edge event (when a pyramid grows outwards)
          // a surface split is not possible with only one edge in between
          continue;
        }

        if ((edge_1->getFacetL() == facet_2_src && facet_1_src == edge_2->getFacetL()) ||
            (edge_1->getFacetL() == facet_2_dst && facet_1_src == edge_2->getFacetR()) ||
            (edge_1->getFacetR() == facet_2_src && facet_1_dst == edge_2->getFacetL()) ||
            (edge_1->getFacetR() == facet_2_dst && facet_1_dst == edge_2->getFacetR()) ||
            (edge_1->getFacetR() == facet_2_src && facet_1_src == edge_2->getFacetR()) ||
            (edge_1->getFacetR() == facet_2_dst && facet_1_src == edge_2->getFacetL()) ||
            (edge_1->getFacetL() == facet_2_src && facet_1_dst == edge_2->getFacetR()) ||
            (edge_1->getFacetL() == facet_2_dst && facet_1_dst == edge_2->getFacetL())) {
          // vertex event
          continue;
        }

        // convex split event checks are performed at pop time - see isActualSurfaceEvent()

        // let's just check if bboxes overlap first
        if (time_future_bound.has_value()) {
          CGAL::Bbox_3 b2;
          b2 += edge_2->getVertexSrc()->getPoint()->bbox();
          b2 += edge_2->getVertexDst()->getPoint()->bbox();
          b2 += HdsUtils::getFinalPoint(edge_2->getVertexSrc(), *time_future_bound)->bbox();
          b2 += HdsUtils::getFinalPoint(edge_2->getVertexDst(), *time_future_bound)->bbox();

          if (!CGAL::do_overlap(b1, b2)) {
            CGAL_SS3_CORE_TRACE_V(64, "Filtered possible surface event candidates\n\t" << edge_1->toString() << "\n\t"
                                                                                        << edge_2->toString());
            continue;
          }
        }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
        ++tested_candidates;
#endif

        std::optional<FT> event_time = crashAtTime(edge_1, edge_2, current_time, time_future_bound);
        if (!event_time) {
          continue;
        }

        CGAL_assertion(*event_time < current_time);
        CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

        SurfaceEventSPtr event = SurfaceEvent::create();
        event->setTime(*event_time);
        event->setEdge1(edge_1);
        event->setEdge2(edge_2);
        queue.push(event);
      }
    }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    CGAL_SS3_CORE_TRACE_V(16, "  Tested: " << tested_candidates << " / " << total_candidates << " candidates");
#endif

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Surface Events in: " << timer.time());
#endif
  }

  void collectSurfaceEvents(const PolyhedronSPtr& polyhedron,
                            const FT& current_time,
                            const std::optional<FT>& time_future_bound,
                            PQ& queue)
  {
    return collectSurfaceEvents(polyhedron->edges(), polyhedron,
                                current_time, time_future_bound, queue);
  }

  /**
    * This event occurs when two edges collide.
    * The first edge is always reflex.
    */
  void collectPolyhedronSplitEvent(const EdgeSPtr& edge_1,
                                   const EdgeSPtr& edge_2,
                                   const PolyhedronSPtr& /*polyhedron*/,
                                   const FT& current_time,
                                   const std::optional<FT>& time_future_bound,
                                   PQ& queue)
  {
    CGAL_SS3_DEBUG_SPTR(edge_1);
    CGAL_SS3_DEBUG_SPTR(edge_2);

    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();

    if (edge_1->getVertexSrc()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
        edge_1->getVertexSrc()->getPoint() == edge_2->getVertexDst()->getPoint() ||
        edge_1->getVertexDst()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
        edge_1->getVertexDst()->getPoint() == edge_2->getVertexDst()->getPoint()) {
      // share a vertex
      return;
    }
    if (!((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() == facet_1_dst) ||
          (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() == facet_1_src))) {
      // no polyhedron split event
      return;
    }
    if (edge_1->getVertexSrc()->findEdge(edge_2->getVertexSrc()) ||
        edge_1->getVertexSrc()->findEdge(edge_2->getVertexDst()) ||
        edge_1->getVertexDst()->findEdge(edge_2->getVertexSrc()) ||
        edge_1->getVertexDst()->findEdge(edge_2->getVertexDst())) {
      // does not work when there is only one edge in between
      return;
    }

    std::optional<FT> event_time = crashAtTime(edge_1, edge_2, current_time, time_future_bound);
    if (!event_time) {
      return;
    }

    CGAL_assertion(*event_time < current_time);
    CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

    // degeneracy check is performed at pop time - see isActualPolyhedronSplitEvent()

    PolyhedronSplitEventSPtr event = PolyhedronSplitEvent::create();
    event->setTime(*event_time);
    event->setEdge1(edge_1);
    event->setEdge2(edge_2);
    queue.push(event);
  }

  void collectPolyhedronSplitEvents(const std::list<EdgeSPtr>& edges,
                                    const PolyhedronSPtr& polyhedron,
                                    const FT& current_time,
                                    const std::optional<FT>& time_future_bound,
                                    PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Polyhedron Split Events [" << current_time << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    for (const EdgeSPtr& edge_1 : edges) {
      CGAL_SS3_DEBUG_SPTR(edge_1);

      if (!HdsUtils::isReflex(edge_1)) {
        continue;
      }

      FacetSPtr facet_1_src = edge_1->getFacetSrc();
      CGAL_assertion(bool(facet_1_src));
      for (const EdgeSPtr& edge_2 : facet_1_src->edges()) {
        collectPolyhedronSplitEvent(edge_1, edge_2, polyhedron,
                                    current_time, time_future_bound,
                                    queue);
      }
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Polyhedron Split Events in: " << timer.time());
#endif
  }

  void collectPolyhedronSplitEvents(const PolyhedronSPtr& polyhedron,
                                    const FT& current_time,
                                    const std::optional<FT>& time_future_bound,
                                    PQ& queue)
  {
    return collectPolyhedronSplitEvents(polyhedron->edges(), polyhedron,
                                        current_time, time_future_bound, queue);
  }

  /**
    * Split Merge event
    */
  void collectSplitMergeEvents(const std::list<VertexSPtr>& vertices,
                               const PolyhedronSPtr& /*polyhedron*/,
                               const bool use_canonical_event_reps,
                               const FT& current_time,
                               const std::optional<FT>& time_future_bound,
                               PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Split Merge Events [" << current_time << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    for (const VertexSPtr& vertex_1 : vertices) {
      CGAL_SS3_DEBUG_SPTR(vertex_1);

      if (HdsUtils::isConvex(vertex_1)) {
        continue;
      }

      std::set<VertexSPtr> vertices_2;
      for (FacetWPtr facet_wptr : vertex_1->facets()) {
        if (FacetSPtr facet = facet_wptr.lock()) {
          vertices_2.insert(facet->vertices().begin(), facet->vertices().end());
        }
    }

# ifdef CGAL_SS3_RUN_TIMERS
    // CGAL_SS3_CORE_TRACE("filled vertices_2 after " << timer.time() << "s");
# endif

    for (const VertexSPtr& vertex_2 : vertices_2) {
        CGAL_SS3_DEBUG_SPTR(vertex_2);

        if (vertex_1 == vertex_2) {
          continue;
        }
#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        if (use_canonical_event_reps) {
          if (vertex_1->getID() > vertex_2->getID()) {
            continue;
          }
        }
#endif
        if (vertex_1->getPoint() == vertex_2->getPoint()) {
          continue;
        }
        if (vertex_1->findEdge(vertex_2)) {
          // edge event
          continue;
        }
        if (HdsUtils::isConvex(vertex_2)) {
          continue;
        }

        // See comment in the first instance of this swap
        VertexSPtr v1 = vertex_1;
        VertexSPtr v2 = vertex_2;
        if (v1->getID() > v2->getID()) {
          std::swap(v1, v2);
        }

        FacetSPtr facet_1;
        FacetSPtr facet_2;
        int num_equal_facets = 0;
        for (FacetWPtr facet_1_wptr : v1->facets()) {
          if (FacetSPtr f1 = facet_1_wptr.lock()) {
            for (FacetWPtr facet_2_wptr : v2->facets()) {
              if (FacetSPtr f2 = facet_2_wptr.lock()) {
                if (f1 == f2) {
                  if (num_equal_facets == 0) {
                    facet_1 = f1;
                  } else {
                    facet_2 = f2;
                  }
                  ++num_equal_facets;
                }
              }
            }
          }
        }
        if (num_equal_facets != 2) {
          continue;
        }
        if (facet_1->next(v1) != facet_2) {
          FacetSPtr facet_tmp = facet_1;
          facet_1 = facet_2;
          facet_2 = facet_tmp;
        }
        CGAL_assertion(facet_1->next(v1) == facet_2);
        if (v1->next(facet_1)->next(facet_1) == v2 ||
            v1->next(facet_2)->next(facet_2) == v2 ||
            v1->prev(facet_1)->prev(facet_1) == v2 ||
            v1->prev(facet_2)->prev(facet_2) == v2) {
          // edge merge event
          continue;
        }

        EdgeSPtr edge_11 = EdgeSPtr();
        EdgeSPtr edge_12 = EdgeSPtr();
        for (EdgeWPtr edge_1_wptr : v1->edges()) {
          if (EdgeSPtr edge_1 = edge_1_wptr.lock()) {
            FacetSPtr facet_1l = edge_1->getFacetL();
            FacetSPtr facet_1r = edge_1->getFacetR();
            if ((facet_1l == facet_1 && facet_1r != facet_2) ||
                (facet_1r == facet_1 && facet_1l != facet_2)) {
              edge_11 = edge_1;
            } else if ((facet_1l == facet_2 && facet_1r != facet_1) ||
                       (facet_1r == facet_2 && facet_1l != facet_1)) {
              edge_12 = edge_1;
            }
          }
        }
        EdgeSPtr edge_21 = EdgeSPtr();
        EdgeSPtr edge_22 = EdgeSPtr();
        for (EdgeWPtr edge_2_wptr : v2->edges()) {
          if (EdgeSPtr edge_2 = edge_2_wptr.lock()) {
            FacetSPtr facet_2l = edge_2->getFacetL();
            FacetSPtr facet_2r = edge_2->getFacetR();
            if ((facet_2l == facet_1 && facet_2r != facet_2) ||
                (facet_2r == facet_1 && facet_2l != facet_2)) {
              edge_21 = edge_2;
            } else if ((facet_2l == facet_2 && facet_2r != facet_1) ||
                       (facet_2r == facet_2 && facet_2l != facet_1)) {
              edge_22 = edge_2;
            }
          }
        }

        // convex split event checks are performed at pop time - see isActualSplitMergeEvent()

        std::optional<FT> event_time = crashAtTime(edge_11, edge_22, current_time, time_future_bound);
        if (!event_time) {
          continue;
        }

        CGAL_assertion(*event_time < current_time);
        CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

        SplitMergeEventSPtr event = SplitMergeEvent::create();
        event->setTime(*event_time);
        event->setVertex1(v1);
        event->setVertex2(v2);
        event->setFacet1(facet_1);
        event->setFacet2(facet_2);
        queue.push(event);
      }
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Split Merge Events in: " << timer.time());
#endif
  }

  void collectSplitMergeEvents(const PolyhedronSPtr& polyhedron,
                               const FT& current_time,
                               const std::optional<FT>& time_future_bound,
                               PQ& queue)
  {
    return collectSplitMergeEvents(polyhedron->vertices(), polyhedron, true /*use canonical reps*/,
                                   current_time, time_future_bound, queue);
  }

  /**
    * This event occurs when two edges collide.
    * The first edge is always reflex.
    */
  void collectEdgeSplitEvents(const std::list<EdgeSPtr>& edges_1,
                              const std::list<EdgeSPtr>& edges_2,
                              const PolyhedronSPtr& /*polyhedron*/,
                              const bool use_canonical_event_reps,
                              const FT& current_time,
                              const std::optional<FT>& time_future_bound,
                              PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Edge Split Events [" << current_time << "]");

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    unsigned int tested_candidates = 0;
    unsigned int total_candidates = 0;
#endif

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    std::list<EdgeSPtr> edges_reflex_1, edges_reflex_2;
    auto fill_reflex_edges = [&](const std::list<EdgeSPtr>& edges,
                                 std::list<EdgeSPtr>& edges_reflex) {
      for (const EdgeSPtr& edge : edges) {
        if (HdsUtils::isReflex(edge)) {
          edges_reflex.push_back(edge);
        }
      }
    };

    fill_reflex_edges(edges_1, edges_reflex_1);

    if (edges_reflex_1.empty()) {
#ifdef CGAL_SS3_RUN_TIMERS
      timer.stop();
      CGAL_SS3_CORE_TRACE_V(4, "  Sought Edge Split Events in: " << timer.time());
#endif
      return;
    }

    fill_reflex_edges(edges_2, edges_reflex_2);

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL_SS3_CORE_TRACE_V(16, "  Collect reflex edges: " << timer.time());
#endif

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    CGAL_SS3_CORE_TRACE_V(16, "  " << edges_reflex_1.size() << " reflex edges [1]");
    CGAL_SS3_CORE_TRACE_V(16, "  " << edges_reflex_2.size() << " reflex edges [2]");
#endif

    // maybe we will use canonical reps with non symmetrical ranges one day, but not for now
    CGAL_warning(!use_canonical_event_reps || edges_reflex_1 == edges_reflex_2);

    for (const EdgeSPtr& edge_1 : edges_reflex_1) {
      CGAL_SS3_DEBUG_SPTR(edge_1);

      FacetSPtr facet_1_src = edge_1->getFacetSrc();
      FacetSPtr facet_1_dst = edge_1->getFacetDst();

      CGAL::Bbox_3 b1;
      if (time_future_bound.has_value()) {
        b1 += edge_1->getVertexSrc()->getPoint()->bbox();
        b1 += edge_1->getVertexDst()->getPoint()->bbox();
        b1 += HdsUtils::getFinalPoint(edge_1->getVertexSrc(), *time_future_bound)->bbox();
        b1 += HdsUtils::getFinalPoint(edge_1->getVertexDst(), *time_future_bound)->bbox();
      }

      for (const EdgeSPtr& edge_2 : edges_reflex_2) {
        CGAL_SS3_DEBUG_SPTR(edge_2);

        CGAL_SS3_CORE_TRACE_V(64, "Possible edge split event\n\t" << edge_1->toString() << "\n\t"
                                                                  << edge_2->toString());

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
        ++total_candidates;
#endif

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        if (use_canonical_event_reps) {
          if (edge_1->getID() > edge_2->getID()) {
            continue;
          }
        }
#else
        if (edge_1 == edge_2) {
          continue;
        }
#endif

        if (edge_1->getFacetL() == edge_2->getFacetL() ||
            edge_1->getFacetL() == edge_2->getFacetR() ||
            edge_1->getFacetR() == edge_2->getFacetL() ||
            edge_1->getFacetR() == edge_2->getFacetR()) {
          // on same facet
          continue;
        }
        if (edge_1->getVertexSrc()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
            edge_1->getVertexSrc()->getPoint() == edge_2->getVertexDst()->getPoint() ||
            edge_1->getVertexDst()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
            edge_1->getVertexDst()->getPoint() == edge_2->getVertexDst()->getPoint()) {
          // share a vertex
          continue;
        }

        if (((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() == facet_1_dst) ||
             (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() == facet_1_src))) {
          // polyhedron split event
          continue;
        }
        FacetSPtr facet_2_src = edge_2->getFacetSrc();
        FacetSPtr facet_2_dst = edge_2->getFacetDst();
        if ((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() != facet_1_dst) ||
            (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() != facet_1_src) ||
            (edge_2->getFacetR() == facet_1_src && edge_2->getFacetL() != facet_1_dst) ||
            (edge_2->getFacetR() == facet_1_dst && edge_2->getFacetL() != facet_1_src)) {
          // surface event
          continue;
        }
        if ((edge_1->getFacetL() == facet_2_src && edge_1->getFacetR() != facet_2_dst) ||
            (edge_1->getFacetL() == facet_2_dst && edge_1->getFacetR() != facet_2_src) ||
            (edge_1->getFacetR() == facet_2_src && edge_1->getFacetL() != facet_2_dst) ||
            (edge_1->getFacetR() == facet_2_dst && edge_1->getFacetL() != facet_2_src)) {
          // surface event
          continue;
        }

        // Build 2 pairs of quads from each edge + shifted edge and check intersections.
        // - The quad is planar because a shifting edge moves on a plane
        // - The shifted edge could be different due to other events... but then this other event
        //   will be treated first and the event will be invalid

        // let's just check if bboxes overlap first
        if (time_future_bound.has_value()) {
          CGAL::Bbox_3 b2;
          b2 += edge_2->getVertexSrc()->getPoint()->bbox();
          b2 += edge_2->getVertexDst()->getPoint()->bbox();
          b2 += HdsUtils::getFinalPoint(edge_2->getVertexSrc(), *time_future_bound)->bbox();
          b2 += HdsUtils::getFinalPoint(edge_2->getVertexDst(), *time_future_bound)->bbox();

          if (!CGAL::do_overlap(b1, b2)) {
            CGAL_SS3_CORE_TRACE_V(64, "Filtered edge split candidates\n\t" << edge_1->toString() << "\n\t"
                                                                           << edge_2->toString() << "\n\t"
                                                                           << b1 << "\n\t" << b2);
            continue;
          }
        }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
        ++tested_candidates;
#endif

        std::optional<FT> event_time = crashAtTime(edge_1, edge_2, current_time, time_future_bound);
        if (!event_time) {
          continue;
        }

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        // @fixme below needs to be double checked
        if (use_canonical_event_reps) {
          FT shift = *event_time - current_time;
          Segment3SPtr e1o = Transformation::shiftEdge(edge_1, shift);
          if (e1o->is_degenerate()) {
            // polyhedron split
            continue;
          }
          Segment3SPtr e2o = Transformation::shiftEdge(edge_2, shift);
          if (e2o->is_degenerate()) {
            // polyhedron split
            continue;
          }
        }
#endif

        EdgeSplitEventSPtr event = EdgeSplitEvent::create();
        event->setTime(*event_time);
        event->setEdge1(edge_1);
        event->setEdge2(edge_2);
        queue.push(event);
      }
    }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    CGAL_SS3_CORE_TRACE_V(16, "  Tested: " << tested_candidates << " / " << total_candidates << " candidates");
#endif

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Edge Split Events in: " << timer.time());
#endif
  }

  void collectEdgeSplitEvents(const PolyhedronSPtr& polyhedron,
                              const FT& current_time,
                              const std::optional<FT>& time_future_bound,
                              PQ& queue)
  {
    return collectEdgeSplitEvents(polyhedron->edges(), polyhedron->edges(), polyhedron, true /*use canonical reps*/,
                                  current_time, time_future_bound, queue);
  }

  // The function below is meant to be used as the callback of the box_d spatial searching:
  // it does not perform any box-box filtering
  void collectEdgeSplitEvent(const EdgeSPtr& edge_1,
                             const EdgeSPtr& edge_2,
                             const PolyhedronSPtr& /*polyhedron*/,
                             const bool use_canonical_event_reps,
                             const FT& current_time,
                             const std::optional<FT>& time_future_bound,
                             PQ& queue)
  {
    CGAL::Real_timer timer;
    timer.start();

    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();

    CGAL_SS3_DEBUG_SPTR(edge_1);
    CGAL_SS3_DEBUG_SPTR(edge_2);

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
    if (use_canonical_event_reps) {
      if (edge_1->getID() > edge_2->getID()) {
        return;
      }
    }
#else
    if (edge_1 == edge_2) {
      return;
    }
#endif

    if (edge_1->getFacetL() == edge_2->getFacetL() ||
        edge_1->getFacetL() == edge_2->getFacetR() ||
        edge_1->getFacetR() == edge_2->getFacetL() ||
        edge_1->getFacetR() == edge_2->getFacetR()) {
      // on same facet
      return;
    }
    if (edge_1->getVertexSrc()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
        edge_1->getVertexSrc()->getPoint() == edge_2->getVertexDst()->getPoint() ||
        edge_1->getVertexDst()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
        edge_1->getVertexDst()->getPoint() == edge_2->getVertexDst()->getPoint()) {
      // share a vertex
      return;
    }

    if (((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() == facet_1_dst) ||
         (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() == facet_1_src))) {
        // polyhedron split event
        return;
    }
    FacetSPtr facet_2_src = edge_2->getFacetSrc();
    FacetSPtr facet_2_dst = edge_2->getFacetDst();
    if ((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() != facet_1_dst) ||
        (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() != facet_1_src) ||
        (edge_2->getFacetR() == facet_1_src && edge_2->getFacetL() != facet_1_dst) ||
        (edge_2->getFacetR() == facet_1_dst && edge_2->getFacetL() != facet_1_src)) {
      // surface event
      return;
    }
    if ((edge_1->getFacetL() == facet_2_src && edge_1->getFacetR() != facet_2_dst) ||
        (edge_1->getFacetL() == facet_2_dst && edge_1->getFacetR() != facet_2_src) ||
        (edge_1->getFacetR() == facet_2_src && edge_1->getFacetL() != facet_2_dst) ||
        (edge_1->getFacetR() == facet_2_dst && edge_1->getFacetL() != facet_2_src)) {
      // surface event
        return;
    }

    std::optional<FT> event_time = crashAtTime(edge_1, edge_2, current_time, time_future_bound);
    if (!event_time) {
      return;
    }

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
    // @fixme below needs to be double checked
    if (use_canonical_event_reps) {
      FT shift = *event_time - current_time;
      Segment3SPtr e1o = Transformation::shiftEdge(edge_1, shift);
      if (e1o->is_degenerate()) {
        // polyhedron split
        return;
      }
      Segment3SPtr e2o = Transformation::shiftEdge(edge_2, shift);
      if (e2o->is_degenerate()) {
        // polyhedron split
        return;
      }
    }
#endif

    EdgeSplitEventSPtr event = EdgeSplitEvent::create();
    event->setTime(*event_time);
    event->setEdge1(edge_1);
    event->setEdge2(edge_2);
    queue.push(event);
  }

  void collectEdgeSplitEventsWithBoxD(const std::list<EdgeSPtr>& edges_1,
                                      const std::list<EdgeSPtr>& edges_2,
                                      const PolyhedronSPtr& polyhedron,
                                      const bool use_canonical_event_reps,
                                      const FT& current_time,
                                      const std::optional<FT>& time_future_bound,
                                      PQ& queue)
  {
    CGAL_precondition(time_future_bound.has_value());

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    using Box = CGAL::Box_intersection_d::Box_with_handle_d<double, 3, EdgeSPtr>;

    auto fill_boxes = [&](const std::list<EdgeSPtr>& edges, std::vector<Box>& boxes)
    {
      for (const EdgeSPtr& edge : edges) {
        if (!HdsUtils::isReflex(edge)) {
          continue;
        }

        CGAL::Bbox_3 b = edge->getVertexSrc()->getPoint()->bbox();
        b += edge->getVertexDst()->getPoint()->bbox();
        b += HdsUtils::getFinalPoint(edge->getVertexSrc(), *time_future_bound)->bbox();
        b += HdsUtils::getFinalPoint(edge->getVertexDst(), *time_future_bound)->bbox();

        boxes.emplace_back(b, edge);
      }
    };

    std::vector<Box> boxes_1, boxes_2;
    fill_boxes(edges_1, boxes_1);
    fill_boxes(edges_2, boxes_2);

    // std::cout << boxes_1.size() << " boxes (1)" << std::endl;
    // std::cout << boxes_2.size() << " boxes (2)" << std::endl;

    // int callback_count = 0;

    auto callback = [&](const Box& box_a, const Box& box_b)
    {
      // ++callback_count;
      EdgeSPtr edge_1 = box_a.handle();
      EdgeSPtr edge_2 = box_b.handle();
      return collectEdgeSplitEvent(edge_1, edge_2, polyhedron, use_canonical_event_reps,
                                   current_time, time_future_bound, queue);
    };

    // std::cout << "Queue before: " << queue.size() << std::endl;
    CGAL::box_intersection_d(boxes_1.begin(), boxes_1.end(), boxes_2.begin(), boxes_2.end(), callback);
    // std::cout << "Callbacks: " << callback_count << std::endl;
    // std::cout << "Queue after: " << queue.size() << std::endl;

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Edge Split Events in: " << timer.time());
#endif
  }

  void collectEdgeSplitEventsWithBoxD(const PolyhedronSPtr& polyhedron,
                                      const FT& current_time,
                                      const std::optional<FT>& time_future_bound,
                                      PQ& queue)
  {
    CGAL_precondition(time_future_bound.has_value());

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    using Box = CGAL::Box_intersection_d::Box_with_handle_d<double, 3, EdgeSPtr>;

    std::vector<Box> boxes;

    for (const EdgeSPtr& edge : polyhedron->edges()) {
      if (!HdsUtils::isReflex(edge)) {
        continue;
      }

      CGAL::Bbox_3 b = edge->getVertexSrc()->getPoint()->bbox();
      b += edge->getVertexDst()->getPoint()->bbox();
      b += HdsUtils::getFinalPoint(edge->getVertexSrc(), *time_future_bound)->bbox();
      b += HdsUtils::getFinalPoint(edge->getVertexDst(), *time_future_bound)->bbox();

      boxes.emplace_back(b, edge);
    }

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL_SS3_CORE_TRACE_V(16, "  Built boxes: " << timer.time());
#endif

    auto callback = [&](const Box& box_a, const Box& box_b)
    {
      EdgeSPtr edge_1 = box_a.handle();
      EdgeSPtr edge_2 = box_b.handle();

      // no canonical reps because box_d already returns a single instance for each pair
      return collectEdgeSplitEvent(edge_1, edge_2, polyhedron, false /*use_canonical_event_reps*/,
                                   current_time, time_future_bound, queue);
    };

    CGAL::box_self_intersection_d(boxes.begin(), boxes.end(), callback);

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Edge Split Events in: " << timer.time());
#endif
  }

  /**
    * A reflex vertex reaches a facet.
    */
  // @speed should keep a (sorted) list of vertex ---> facet of known contact events
  // to filter how far we need to look and also avoid checking again multiple times
  // something like SLS2
  void collectPierceEvents(const std::list<VertexSPtr>& vertices,
                           const std::list<FacetSPtr>& facets,
                           const PolyhedronSPtr& /*polyhedron*/,
                           const FT& current_time,
                           const std::optional<FT>& time_future_bound,
                           PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Pierce Events [" << current_time << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    unsigned int pierce_vertex_counter = 0;
    unsigned int total_candidates = 0;
    unsigned int tested_candidates = 0;
#endif

    for (const VertexSPtr& vertex : vertices) {
      CGAL_assertion(vertex->getID() != -1);

      if (HdsUtils::isReflex(vertex)) {
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
        ++pierce_vertex_counter;
#endif

        for (const FacetSPtr& facet : facets) {
          // @todo contains_vertex is redundant with has_edge_to_facet?
          bool contains_vertex = false;
          for (const VertexSPtr& vertex_2 : facet->vertices()) {
            if (vertex_2->getPoint() == vertex->getPoint()) {
              contains_vertex = true;
              break;
            }
          }
          if (contains_vertex) {
            continue;
          }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
          ++total_candidates;
#endif
          // # `has_edge_to_facet` check is done at pop time - see isActualPierceEvent()

          // # Shrinking, so the vertex must be on the backside of the plane
          if (KernelWrapper::side(facet->getPlane(), vertex->getPoint()) > 0) {
            continue;
          }

          // # If the facet is so far that even when shifting point and plane by the maximal
          // time, the vertex would not cross it, then we are done
          if (time_future_bound.has_value()) {
            Point3SPtr shifted_pt = HdsUtils::getFinalPoint(vertex, *time_future_bound);
            Plane3SPtr shifted_plane = HdsUtils::getFinalPlane(facet, *time_future_bound);

            if (KernelWrapper::side(shifted_plane, shifted_pt) < 0) {
              continue;
            }

            // # It would be tempting here to use a bbox around the vertex and facet's trajectories,
            // but if we did, then we would have to detect potential piercing events involving
            // a facet after each modification of the facet, and that seems costly
            // (even if potential piercing vertices were cached...)

            // # If both move in the same direction, we cannot pierce
            Vector_3 d (*(vertex->getPoint()), *shifted_pt);
            Vector_3 n = facet->getPlane()->orthogonal_vector();
            if (n * d < 0) { // negative because the facet moves in a direction opposite of the normal
              continue;
            }
          }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
          ++tested_candidates;
#endif

          std::optional<FT> event_time = crashAtTime(vertex, facet, current_time, time_future_bound);
          if (!event_time) {
            continue;
          }

          CGAL_assertion(*event_time < current_time);
          CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

          // actual intersection checks are performed at pop time - see isActualPierceEvent()

          PierceEventSPtr event = PierceEvent::create();
          event->setTime(*event_time);
          event->setFacet(facet);
          event->setVertex(vertex);
          queue.push(event);
        }
      }
    }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    CGAL_SS3_CORE_TRACE_V(16, "  " << pierce_vertex_counter << " reflex vertices");
    CGAL_SS3_CORE_TRACE_V(16, "  Tested: " << tested_candidates << " / " << total_candidates << " candidates");
#endif

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Pierce Events in: " << timer.time());
#endif
  }

  void collectPierceEvents(const PolyhedronSPtr& polyhedron,
                           const FT& current_time,
                           const std::optional<FT>& time_future_bound,
                           PQ& queue)
  {
    return collectPierceEvents(polyhedron->vertices(), polyhedron->facets(), polyhedron,
                               current_time, time_future_bound, queue);
  }

  /**
    * Collect all events for a subset of the polyhedron
    */
  void collectLocalEvents(const PolyhedronSPtr& polyhedron,
                          const FT& current_time,
                          const std::optional<FT>& time_future_bound,
                          PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(2, "collectLocalEvents(" << current_time << ")");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    CGAL_SS3_CORE_TRACE_V(16, "Past bound = " << current_time);
    CGAL_SS3_CORE_TRACE_IF(time_future_bound, 4, "Initial future bound = " << *time_future_bound);

    // == VANISH EVENTS ==
    {
      std::list<EdgeSPtr> local_edges(post_op_edges_.begin(), post_op_edges_.end());

      CGAL_SS3_CORE_TRACE_V(8, "Local Edges for Vanish Events (" << local_edges.size() << ")");
      CGAL_SS3_CORE_TRACE_CODE(for (const EdgeSPtr& e : local_edges))
      CGAL_SS3_CORE_TRACE_V(8, "\t" << e->toString());

      collectVanishEvents(local_edges, polyhedron, current_time, time_future_bound, queue);
    }

    // == VERTEX-VERTEX EVENTS ==
    {
#if 1
      // for these three events, we need to look farther than the vertices that were involved
      // in the event because an event might add new edges to an edge cycle and so what was
      // before maybe a simple edge merge event now becomes a vertex event...
      std::list<VertexSPtr> local_vertices_VV(post_op_vertices_VV_.begin(),
                                              post_op_vertices_VV_.end());

      CGAL_SS3_CORE_TRACE_V(8, "Local Vertices for Vertex-Vertex Events (" << local_vertices_VV.size() << "):")
      CGAL_SS3_CORE_TRACE_CODE(std::stringstream ss;)
      CGAL_SS3_CORE_TRACE_CODE(for (const VertexSPtr& v : local_vertices_VV))
      CGAL_SS3_CORE_TRACE_CODE(ss << " " << v->getID();)
      CGAL_SS3_CORE_TRACE_V(8, ss.str());

      const bool use_canonical_reps = false;
#else
      const std::list<VertexSPtr>& local_vertices_VV = polyhedron->vertices();
      const bool use_canonical_reps = true;
#endif

      collectVertexEvents(local_vertices_VV, polyhedron, use_canonical_reps,
                          current_time, time_future_bound, queue);
      collectFlipVertexEvents(local_vertices_VV, polyhedron, use_canonical_reps,
                              current_time, time_future_bound, queue);
      collectSplitMergeEvents(local_vertices_VV, polyhedron, use_canonical_reps,
                              current_time, time_future_bound, queue);
    }

    // == POLYHEDRON SPLIT EVENTS ==
    {
#if 1
      // Polyhedron Split Events are not symmetrical, so we need to check both modified edges
      // as the first edge parameter, but also grab all the cases where a modified
      // edge is the second edge.

      std::list<EdgeSPtr> local_edges_EE(post_op_edges_.begin(), post_op_edges_.end());

      CGAL_SS3_CORE_TRACE_V(8, "Local Edges for Polyhedron Events (" << local_edges_EE.size() << ")");
      CGAL_SS3_CORE_TRACE_CODE(for (const EdgeSPtr& e : local_edges_EE))
      CGAL_SS3_CORE_TRACE_V(8, "\t" << e->toString());

      // this is the modified edges as 'edge_1'
      collectPolyhedronSplitEvents(local_edges_EE, polyhedron,
                                   current_time, time_future_bound, queue);

      // this is the modified edges as 'edge_2'
      // since we know that for a polyhedron split event, edge_2's LR facets
      // are the src and dst of edge_1, we can build a subset of all edges for edge_1s
      for (const EdgeSPtr& edge_2 : local_edges_EE) {
        std::set<EdgeSPtr> edges_1;
        for (const FacetSPtr& facet_2 : {edge_2->getFacetL(), edge_2->getFacetR()}) {
          for (const VertexSPtr& vertex_1 : facet_2->vertices()) {
            // keep the edge incident to 'vertex_1' which is not incident to 'facet_2'
            EdgeSPtr edge_1;
            for (EdgeWPtr edge_wptr : vertex_1->edges()) {
              if (EdgeSPtr edge = edge_wptr.lock()) {
                if (edge->getFacetL() != facet_2 && edge->getFacetR() != facet_2) {
                  edge_1 = edge;
                  break;
                }
              }
            }

            CGAL_SS3_DEBUG_SPTR(edge_1);
            edges_1.insert(edge_1);
          }
        }

        for (const EdgeSPtr& edge_1 : edges_1) {
          collectPolyhedronSplitEvent(edge_1, edge_2, polyhedron,
                                      current_time, time_future_bound, queue);
        }
      }


#else
      collectPolyhedronSplitEvents(polyhedron->edges(), polyhedron,
                                   current_time, time_future_bound, queue);
#endif
    }

    // == PIERCE EVENTS ==
    {
#if 1
      // For Pierce events, we add new concave corners but we must also add the other endpoints
      // of modified edges because a Pierce event might be created at an unmodified vertex
      // after an incident edge got disconnected from a facet by an event
      //
      // However, this does not need to be done after all events (e.g. a triangle event
      // cannot disconnect an edge from another facet so this is a pointless addition,
      // but surface events can disconnect)
      //
      // @speed also, for these vertices, we should not check ALL the extremities and ALL the faces,
      // but just the vertex and the facet which got disconnected?
      // something like extra_combinations_for_pierce_events...
      std::list<VertexSPtr> local_vertices_VF(post_op_vertices_pierce_.begin(), post_op_vertices_pierce_.end());

      CGAL_SS3_CORE_TRACE_V(8, "Local Vertices for Pierce Events (" << local_vertices_VF.size() << "):");
      CGAL_SS3_CORE_TRACE_CODE(std::stringstream ss;)
      CGAL_SS3_CORE_TRACE_CODE(for (const VertexSPtr& v : local_vertices_VF))
      CGAL_SS3_CORE_TRACE_CODE(ss << " " << v->getID());
      CGAL_SS3_CORE_TRACE_V(8, ss.str());
#else
      std::list<VertexSPtr> local_vertices_VF = polyhedron->vertices();
#endif

      // Pierce events must check with all faces, no choice about this
      collectPierceEvents(local_vertices_VF, polyhedron->facets(), polyhedron,
                          current_time, time_future_bound, queue);
    }

    // == SURFACE EVENTS ==
    {
      // Surface events are not symmetrical, so we need to check both modified edges
      // as the first edge parameter, but also grab all the cases where a modified
      // edge is the second edge.
#if 1
      std::list<EdgeSPtr> local_edges_EE(post_op_edges_.begin(), post_op_edges_.end());

      CGAL_SS3_CORE_TRACE_V(8, "Local Edges for Surface Events (" << local_edges_EE.size() << ")");
      CGAL_SS3_CORE_TRACE_CODE(for (const EdgeSPtr& e : local_edges_EE))
      CGAL_SS3_CORE_TRACE_V(8, "\t" << e->toString());

      // this is the modified edges as 'edge_1'
      collectSurfaceEvents(local_edges_EE, polyhedron,
                            current_time, time_future_bound, queue);

      // this is the modified edges as 'edge_2'
      for (const EdgeSPtr& edge_2 : local_edges_EE) {
        std::set<EdgeSPtr> edges_1;
        for (const FacetSPtr& facet_2 : {edge_2->getFacetL(), edge_2->getFacetR()}) {
          for (const VertexSPtr& vertex_1 : facet_2->vertices()) {
            // keep the edge incident to 'vertex_1' which is not incident to 'facet_2'
            EdgeSPtr edge_1;
            for (EdgeWPtr edge_wptr : vertex_1->edges()) {
              if (EdgeSPtr edge = edge_wptr.lock()) {
                if (edge->getFacetL() != facet_2 && edge->getFacetR() != facet_2) {
                  edge_1 = edge;
                  break;
                }
              }
            }

            CGAL_SS3_DEBUG_SPTR(edge_1);
            edges_1.insert(edge_1);
          }
        }

        for (const EdgeSPtr& edge_1 : edges_1) {
          collectSurfaceEvent(edge_1, edge_2, polyhedron,
                              current_time, time_future_bound, queue);
        }
      }
#else
      collectSurfaceEvents(polyhedron->edges(), polyhedron,
                            current_time, time_future_bound, queue);
#endif
    }

    // == EDGE SPLIT EVENTS ==
    {
#if 1
      // During collection of edge split events, filtering is only based on incident faces.
      // Thus, regardless of the N-1 event, two unmodified edges cannot have a change
      // of incident faces (since otherwise they would be modified by definition),
      // and thus it is enough to look at modified edges to get the new events.
      //
      // Note that even events that split and create multiple CCs (like an edge merge)
      // do not create new facets, they only create new edge cycles within the same facet.
      std::list<EdgeSPtr> local_edges_EE(post_op_edges_.begin(), post_op_edges_.end());

      CGAL_SS3_CORE_TRACE_V(8, "Local Edges for Edge Split Events (" << local_edges_EE.size() << ")");
      CGAL_SS3_CORE_TRACE_CODE(for (const EdgeSPtr& e : local_edges_EE))
      CGAL_SS3_CORE_TRACE_V(8, "\t" << e->toString());

      const bool use_canonical_reps = false;
#else
      std::list<EdgeSPtr> local_edges_EE = polyhedron->edges();
      const bool use_canonical_reps = true;
#endif

      // not worth the effort? 'local_edges_EE' is always very small
#ifdef CGAL_SS3_DETECT_EDGE_SPLIT_EVENTS_WITH_BOX_D
      if (time_future_bound.has_value()) {
        collectEdgeSplitEventsWithBoxD(local_edges_EE, polyhedron->edges(), polyhedron,
                                       false /*use_canonical_reps*/, // box_d returns a canonical pair
                                       current_time, time_future_bound, queue);
      } else
#endif
      {
        collectEdgeSplitEvents(local_edges_EE, polyhedron->edges(), polyhedron, use_canonical_reps,
                                current_time, time_future_bound, queue);
      }
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(2, "Sought All Local Events in: " << timer.time());
#endif

#ifdef CGAL_SS3_DEBUG_PRINT_QUEUE
    printQueue(queue);
#endif

    CGAL_postcondition(checkQueueCorrectness(queue, polyhedron, current_time, time_future_bound));
  }

  /**
    * Collect all events for the polyhedron
    */
  void collectEvents(const PolyhedronSPtr& polyhedron,
                     const FT& current_time,
                     const std::optional<FT>& time_future_bound,
                     PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(2, "Collecting events in [" << current_time
                          << " | " << (time_future_bound.has_value() ? IO::StringFactory::fromDouble(CGAL::to_double(*time_future_bound)) : "") << "]");

    AbstractEventSPtr result = AbstractEventSPtr();
    if (!polyhedron || polyhedron->facets().size() == 0) {
      return;
    }

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    CGAL_SS3_CORE_TRACE_V(8, "Past bound = " << current_time);
    CGAL_SS3_CORE_TRACE_IF(time_future_bound, 8, "Initial future bound = " << *time_future_bound);

    // --- Vanish Events
    collectVanishEvents(polyhedron, current_time, time_future_bound, queue);

    CGAL_assertion_code(for (const EdgeSPtr& edge : polyhedron->edges()))
    CGAL_assertion_code(HdsUtils::getVanishTime(edge);) // check is_known within skeledgedata

    // --- Contact Event
    collectVertexEvents(polyhedron, current_time, time_future_bound, queue);
    collectFlipVertexEvents(polyhedron, current_time, time_future_bound, queue);
    collectPolyhedronSplitEvents(polyhedron, current_time, time_future_bound, queue);
    collectSplitMergeEvents(polyhedron, current_time, time_future_bound, queue);

    // the next event types are particularly slow, so reduce the bound by doing them last
    // so other events lower the bound
    collectPierceEvents(polyhedron, current_time, time_future_bound, queue);
    collectSurfaceEvents(polyhedron, current_time, time_future_bound, queue);

#ifdef CGAL_SS3_DETECT_EDGE_SPLIT_EVENTS_WITH_BOX_D
    if (time_future_bound.has_value()) {
      collectEdgeSplitEventsWithBoxD(polyhedron, current_time, time_future_bound, queue);
    } else
#endif
    {
      collectEdgeSplitEvents(polyhedron, current_time, time_future_bound, queue);
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(2, "Sought All Events in: " << timer.time());
#endif

#ifdef CGAL_SS3_DEBUG_PRINT_QUEUE
    printQueue(queue);
#endif
  }

  void printQueue(const PQ& queue)
  {
    CGAL_SS3_CORE_TRACE("-------------------------------------------------");
    CGAL_SS3_CORE_TRACE("--- Event queue (size = " << queue.size() << "; iter = " << step_id_ << ") ---");
    CGAL_SS3_CORE_TRACE("-------------------------------------------------");

    PQ duplicate_queue = queue;
    while (!duplicate_queue.empty()) {
      AbstractEventSPtr event = duplicate_queue.top();
      CGAL_SS3_CORE_TRACE("Event E" << event->getID()
                            << " T" << event->getType()
                           << " @ " << event->getTime());
      if (event->isValid()) {
        if (isEventObsolete(event)) {
          CGAL_SS3_CORE_TRACE("  Event is obsolete");
        } else {
          CGAL_SS3_CORE_TRACE(event->toString());
        }
      } else {
        CGAL_SS3_CORE_TRACE(" Event is invalid");
      }
      duplicate_queue.pop();
    }

    CGAL_SS3_CORE_TRACE_CODE(std::stringstream ss;)
    CGAL_SS3_CORE_TRACE_CODE(ss << "Saves:";)
    CGAL_SS3_CORE_TRACE_CODE(for (const FT& save_time : save_times_))
    CGAL_SS3_CORE_TRACE_CODE(ss << " " << save_time;)
    CGAL_SS3_CORE_TRACE(ss.str());
    CGAL_SS3_CORE_TRACE("-------------------------------------------------");
    CGAL_SS3_CORE_TRACE("-------------------------------------------------");
  }

  /**
    * Checks if the queue updated with local events is displaying the same events
    * as if it had been built from scratch.
    */
  // this function checks the correctness of the local queue by computing the queue from scratch
  // and checking that the first valid event is the same for both queues
  bool checkQueueCorrectness(const PQ& queue,
                             const PolyhedronSPtr& polyhedron,
                             const FT& current_time,
                             const std::optional<FT>& time_future_bound)
  {
    CGAL_SS3_CORE_TRACE("Checking queue correctness...");
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    // Compute a queue from scratch using collectEvents()
    PQ queue_from_scratch;
    collectEvents(polyhedron, current_time, time_future_bound, queue_from_scratch);

    // Duplicate the queue since we need to pop events
    PQ duplicate_queue = queue;

    // Tops should be the same
    auto purge_top = [&](PQ& q)
    {
      while (!q.empty()) {
        AbstractEventSPtr event = q.top();
        if (!event->isValid() ||
            isEventInThePast(event, current_time) ||
            isEventObsolete(event) ||
            !isActualEvent(event, current_time, time_future_bound)) {
          q.pop();
        } else {
          break;
        }
      }
    };

    auto is_same_event = [](const AbstractEventSPtr& event_1,
                            const AbstractEventSPtr& event_2)
    {
      if (event_1->getType() != event_2->getType()) {
        return false;
      }

      switch (event_1->getType()) {
        case AbstractEvent::SAVE_EVENT: {
          auto save_event_1 = std::dynamic_pointer_cast<SaveEvent>(event_1);
          auto save_event_2 = std::dynamic_pointer_cast<SaveEvent>(event_2);
          return *save_event_1 == *save_event_2;
        }
        case AbstractEvent::CONST_TIME_EVENT: {
          auto const_event_1 = std::dynamic_pointer_cast<ConstTimeEvent>(event_1);
          auto const_event_2 = std::dynamic_pointer_cast<ConstTimeEvent>(event_2);
          return *const_event_1 == *const_event_2;
        }
        case AbstractEvent::VANISH_EVENT: {
          auto vanish_event_1 = std::dynamic_pointer_cast<VanishEvent>(event_1);
          auto vanish_event_2 = std::dynamic_pointer_cast<VanishEvent>(event_2);
          return *vanish_event_1 == *vanish_event_2;
        }
        case AbstractEvent::EDGE_EVENT: {
          auto edge_event_1 = std::dynamic_pointer_cast<EdgeEvent>(event_1);
          auto edge_event_2 = std::dynamic_pointer_cast<EdgeEvent>(event_2);
          return *edge_event_1 == *edge_event_2;
        }
        case AbstractEvent::EDGE_MERGE_EVENT: {
          auto edge_merge_event_1 = std::dynamic_pointer_cast<EdgeMergeEvent>(event_1);
          auto edge_merge_event_2 = std::dynamic_pointer_cast<EdgeMergeEvent>(event_2);
          return *edge_merge_event_1 == *edge_merge_event_2;
        }
        case AbstractEvent::TRIANGLE_EVENT: {
          auto triangle_event_1 = std::dynamic_pointer_cast<TriangleEvent>(event_1);
          auto triangle_event_2 = std::dynamic_pointer_cast<TriangleEvent>(event_2);
          return *triangle_event_1 == *triangle_event_2;
        }
        case AbstractEvent::DBL_EDGE_MERGE_EVENT: {
          auto dbl_edge_merge_event_1 = std::dynamic_pointer_cast<DblEdgeMergeEvent>(event_1);
          auto dbl_edge_merge_event_2 = std::dynamic_pointer_cast<DblEdgeMergeEvent>(event_2);
          return *dbl_edge_merge_event_1 == *dbl_edge_merge_event_2;
        }
        case AbstractEvent::DBL_TRIANGLE_EVENT: {
          auto dbl_triangle_event_1 = std::dynamic_pointer_cast<DblTriangleEvent>(event_1);
          auto dbl_triangle_event_2 = std::dynamic_pointer_cast<DblTriangleEvent>(event_2);
          return *dbl_triangle_event_1 == *dbl_triangle_event_2;
        }
        case AbstractEvent::TETRAHEDRON_EVENT: {
          auto tetrahedron_event_1 = std::dynamic_pointer_cast<TetrahedronEvent>(event_1);
          auto tetrahedron_event_2 = std::dynamic_pointer_cast<TetrahedronEvent>(event_2);
          return *tetrahedron_event_1 == *tetrahedron_event_2;
        }
        case AbstractEvent::VERTEX_EVENT: {
          auto vertex_event_1 = std::dynamic_pointer_cast<VertexEvent>(event_1);
          auto vertex_event_2 = std::dynamic_pointer_cast<VertexEvent>(event_2);
          return *vertex_event_1 == *vertex_event_2;
        }
        case AbstractEvent::FLIP_VERTEX_EVENT: {
          auto flip_vertex_event_1 = std::dynamic_pointer_cast<FlipVertexEvent>(event_1);
          auto flip_vertex_event_2 = std::dynamic_pointer_cast<FlipVertexEvent>(event_2);
          return *flip_vertex_event_1 == *flip_vertex_event_2;
        }
        case AbstractEvent::SURFACE_EVENT: {
          auto surface_event_1 = std::dynamic_pointer_cast<SurfaceEvent>(event_1);
          auto surface_event_2 = std::dynamic_pointer_cast<SurfaceEvent>(event_2);
          return *surface_event_1 == *surface_event_2;
        }
        case AbstractEvent::POLYHEDRON_SPLIT_EVENT: {
          auto polyhedron_split_event_1 = std::dynamic_pointer_cast<PolyhedronSplitEvent>(event_1);
          auto polyhedron_split_event_2 = std::dynamic_pointer_cast<PolyhedronSplitEvent>(event_2);
          return *polyhedron_split_event_1 == *polyhedron_split_event_2;
        }
        case AbstractEvent::EDGE_SPLIT_EVENT: {
          auto edge_event_1 = std::dynamic_pointer_cast<EdgeSplitEvent>(event_1);
          auto edge_event_2 = std::dynamic_pointer_cast<EdgeSplitEvent>(event_2);
          return *edge_event_1 == *edge_event_2;
        }
        case AbstractEvent::SPLIT_MERGE_EVENT: {
          auto split_merge_event_1 = std::dynamic_pointer_cast<SplitMergeEvent>(event_1);
          auto split_merge_event_2 = std::dynamic_pointer_cast<SplitMergeEvent>(event_2);
          return *split_merge_event_1 == *split_merge_event_2;
        }
        case AbstractEvent::PIERCE_EVENT: {
          auto pierce_event_1 = std::dynamic_pointer_cast<PierceEvent>(event_1);
          auto pierce_event_2 = std::dynamic_pointer_cast<PierceEvent>(event_2);
          return *pierce_event_1 == *pierce_event_2;
        }
        default:
          return false;
      }
    };

    // Get rid of the fake events at the top
    purge_top(queue_from_scratch);
    purge_top(duplicate_queue);

    if (!queue_from_scratch.empty() && !duplicate_queue.empty()) {
      AbstractEventSPtr event_scratch = queue_from_scratch.top();
      AbstractEventSPtr event_duplicate = duplicate_queue.top();
      CGAL_SS3_CORE_TRACE("First event from scratch: " << event_scratch->toString());
      CGAL_SS3_CORE_TRACE("First event from duplicate: " << event_duplicate->toString());

      if (!is_same_event(event_scratch, event_duplicate)) {
        CGAL_SS3_CORE_TRACE("Error: top events differ");
        CGAL_assertion(false);
        return false;
      }
    }

    // We must find all valid 'from-scratch' events in the duplicate queue
    for (;;) {
      purge_top(queue_from_scratch);
      if (queue_from_scratch.empty()) {
        break;
      }

      AbstractEventSPtr event_scratch = queue_from_scratch.top();
      queue_from_scratch.pop();

      CGAL_SS3_CORE_TRACE("Seek event @ " << event_scratch->getTime() << " Type " << event_scratch->getType());

      CGAL_SS3_CORE_TRACE("Event E" << event_scratch->getID()
                          << " T" << event_scratch->getType()
                          << " @ " << event_scratch->getTime());
      CGAL_assertion(event_scratch->isValid() && !isEventObsolete(event_scratch));
      CGAL_SS3_CORE_TRACE(event_scratch->toString());

      // Find the event in the duplicate queue
      bool found = false;
      duplicate_queue = queue;
      while (!duplicate_queue.empty()) {
        AbstractEventSPtr event_duplicate = duplicate_queue.top();
        if (event_duplicate->isValid()) {
          if (is_same_event(event_scratch, event_duplicate)) {
            found = true;
          }

          if (found) {
            if (isEventObsolete(event_duplicate)) {
                CGAL_SS3_CORE_TRACE("Warning: Found event in duplicate queue but it's marked as obsolete");
                CGAL_SS3_CORE_TRACE("Event: " << event_duplicate->toString());
            }
          }
        }

        if (found) {
          break;
        }

        duplicate_queue.pop();
      }

      if (!found) {
          CGAL_SS3_CORE_TRACE("Error: could not find event in duplicate queue");
          CGAL_SS3_CORE_TRACE("Event: " << event_scratch->toString());
          return false;
      }
    }

    CGAL_SS3_CORE_TRACE("OK: found all 'from-scratch' events!");

    return true;
  }

  /**
    * Determines the next event.
    */
  AbstractEventSPtr nextEvent(PQ& queue,
                              const PolyhedronSPtr& polyhedron,
                              const FT& current_time)
  {
    // We could "break" directly in run(), but like this, the save events
    // that are farther than the last event are still processed (arbitrarily decision).
    if (polyhedron->empty()) {
      queue = { };
      // do not return here as to allow save events to exist after the last 'real' event
    }

    AbstractEventSPtr event;
    FT time = 0;

    // If a save event is closest, delay building it in case a const event is even closer
    bool is_save_event = false;

    while (!queue.empty() || !save_times_.empty()) {
      // purge the queue lazily as to avoid wasting time if we stop on the last save event
      if (save_times_.empty()) {
        event = queue.top();
        time = event->getTime();
        is_save_event = false;
      } else {
        // If we have upcoming save events, compare
        const FT& next_save_time = save_times_.front();

        if (queue.empty()) {
          // queue is empty but save_times_ cannot be empty as well
          is_save_event = true;
          time = save_times_.front();
        } else {
          // neither queue nor save_times_ are empty, take the earliest
          if (next_save_time > queue.top()->getTime()) { // save is strictly earlier
            is_save_event = true;
            time = next_save_time;
          } else {
            // save times exist, but are farther in the future than the next event
            event = queue.top();
            time = event->getTime();
          }
        }
      }

      // Tentative next event is not a save event, test its sanity.
      // Do this here because we don't want a const event to get created
      // because it's before a zombie event.
      if (!is_save_event) {
        if (!event->isValid()) {
          CGAL_SS3_CORE_TRACE_V(16, "Skipping invalid event E" << event->getID());
          event = {};
          queue.pop();
          continue;
        }

        if (isEventInThePast(event, current_time)) {
          CGAL_SS3_CORE_TRACE_V(16, "Skipping event-in-the-past E" << event->getID());
          event = {};
          queue.pop();
          continue;
        }

        // This "isObsolete()" function is only a sufficient condition: if the neighborhoods
        // have changed, then the event should be discarded.
        // "IsActual...Event()" takes care of the other checks.
        if (isEventObsolete(event)) {
          CGAL_SS3_CORE_TRACE_V(16, "Skipping obsolete event E" << event->getID());
          event = {};
          queue.pop();
          continue;
        }
      }

      break;
    }

    if (queue.empty() && save_times_.empty()) {
      return { };
    }

    CGAL_assertion(is_save_event || bool(event));
    CGAL_assertion(time != 0);

    // Check if the next const event would be (strictly) closer
    FT pulse = Configuration::getInstance()->getDouble("Algorithm", "const_offset");
    if (pulse != 0) {
      FT next_const_time = floor(CGAL::to_double(current_time / pulse) + 1.0) * pulse;
      if (current_time > next_const_time && next_const_time > time) {
        is_save_event = false;
        ConstTimeEventSPtr const_event = ConstTimeEvent::create();
        const_event->setTime(next_const_time);
        return const_event;
      }
    }

    if (is_save_event) {
      save_times_.erase(save_times_.begin());
      SaveEventSPtr save_event = SaveEvent::create();
      save_event->setTime(time);
      return save_event;
    }

    CGAL_assertion(bool(event));

    // if here, the topmost is neither a const nor save event
    queue.pop();

    return event;
  }

  void addEvent(const AbstractEventSPtr& event)
  {
    typename std::list<AbstractEventSPtr>::iterator it = events_.insert(events_.end(), event);
    event->setListIt(it);
  }

  bool removeEvent(const AbstractEventSPtr& event)
  {
    bool result = false;
    events_.erase(event->getListIt());
    event->setListIt(typename std::list<AbstractEventSPtr>::iterator());
    return result;
  }

  void shiftToEventTime(const PolyhedronSPtr& polyhedron,
                        const FT& current_time,
                        const FT& target_time)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    const FT shift = target_time - current_time;
    CGAL_precondition(!is_zero(shift));

    Transformation::shiftFacets(polyhedron, shift);

#ifdef CGAL_SS3_DUMP_FILES
    // below will have degeneracies since we have not yet treated the event
    static int shift_id = -1;
    IO::OBJFile::save("results/shift_" + std::to_string(++shift_id) + ".obj",
                      polyhedron,
                      false /*do not triangulate*/);
#endif
  }

  bool savePolyhedron(PolyhedronSPtr polyhedron,
                      const FT& current_time)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    bool result = true;

    std::stringstream ss_filename, ss_filename_triangulated, ss_filename_exact;
    ss_filename.precision(17);
    ss_filename_triangulated.precision(17);
    ss_filename_exact.precision(17);
    ss_filename << save_path_.string() << "/time_" << current_time << ".obj";
    ss_filename_triangulated << save_path_.string() << "/time_" << current_time << "_triangulated.obj";
    ss_filename_exact << save_path_.string() << "/time_" << current_time << "_exact.obj";

    result = (IO::OBJFile::save(ss_filename.str(), polyhedron,
                                false /*do_triangulate*/,
                                true /*convert_to_double*/) && result);
    result = (IO::OBJFile::save(ss_filename_triangulated.str(), polyhedron,
                                true /*do_triangulate*/,
                                true /*convert_to_double*/) && result);
    result = (IO::OBJFile::save(ss_filename_exact.str(), polyhedron,
                                true /*do_triangulate*/,
                                false /*convert_to_double*/) && result);

    return result;
  }

  bool saveSkeleton(const StraightSkeletonSPtr& skeleton,
                    const FT& current_time)
  {
    CGAL_SS3_DEBUG_SPTR(skeleton);

    bool result = true;

    std::stringstream ss_filename, ss_filename_triangulated, ss_filename_exact;
    ss_filename << save_path_.string() << "/skel_time_" << current_time << ".obj";
    ss_filename_exact << save_path_.string() << "/skel_time_" << current_time << "_exact.obj";

    result = (IO::OBJFile::save(ss_filename.str(), skeleton,
                                true /*convert_to_double*/) && result);
    result = (IO::OBJFile::save(ss_filename_exact.str(), skeleton,
                                false /*convert_to_double*/) && result);

    return result;
  }

  EventStatus handleSaveEvent(const SaveEventSPtr& event,
                              const FT& current_time,
                              const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#########  Handle Save Event  ##########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->getTime();
    shiftToEventTime(polyhedron, current_time, event_time);

    bool res = true;

#ifdef CGAL_SS3_DUMP_FILES
    res = savePolyhedron(polyhedron, event_time);
    // res = saveSkeleton(skel_result_, event_time) && res;
#endif

    if (res) {
      addEvent(event);
    }

    return (res ? EventStatus::EVENT_HANDLED : EventStatus::EVENT_NOT_HANDLED);
  }

  // This 'handle' is in fact more akin to a collect, but the interesting point
  // is that it happens after pop time
  EventStatus handleConstTimeEvent(ConstTimeEventSPtr event,
                                   const FT& current_time,
                                   const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#########  Handle Const Event  #########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->getTime();
    shiftToEventTime(polyhedron, current_time, event_time);

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  // This 'handle' is in fact more akin to a collect, but the interesting point
  // is that it happens after pop time.
  // This function does NOT shift the polyhedron, it is for the "real" handler
  // to deal with it.
  // This function might not do anything, for example if the vanish event is in fact
  // escalated as a contact event.
  EventStatus handleVanishEvent(const VanishEventSPtr& event,
                                const FT& current_time,
                                const std::optional<FT>& time_future_bound,
                                const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "####### Tentative Vanish Event  ########");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    EdgeSPtr edge = event->getEdge();
    const FT& event_time = event->getTime();

    event->setPoint(vanishesAtPoint(edge)); // @tmp somewhere proper
    Point3SPtr point = event->getPoint();

    // @todo would nice:
    // - to avoid redundant checks (e.g. isTriangle multiple times)
    // - to avoid code duplication with collectXYZEvents()
    //
    // To avoid redundant code, it could be re-ordered by descending amount of collapsing edges (6 --> 1)?
    // Might not speed things up if majority of events are a small number of edge collapses.
    // Anyway, vanish events are cheap to collect, and cheap to process.

    // The fake "for (;;)"" infinite loops are only there as to be able to keep the 'break's
    // such that the code is exactly identical to that of the respective collectXYZEvent() functions

    // Edge Event
    for (;;)
    {
      FacetSPtr facet_l = edge->getFacetL();
      FacetSPtr facet_r = edge->getFacetR();
      if (HdsUtils::isTriangle(facet_l, edge) ||
          HdsUtils::isTriangle(facet_r, edge)) {
        // triangle event
        CGAL_SS3_CORE_TRACE_V(8, "Not an Edge Event (one incident facet is a triangle)");
        break;
      }

      FacetSPtr facet_src = edge->getFacetSrc();
      FacetSPtr facet_dst = edge->getFacetDst();

      // This does not work when there is more than one edge between both facets.
      // EdgeSPtr edge_2 = facet_src->findEdge(facet_dst);
      std::list<EdgeSPtr> edges_2 = facet_src->findEdges(facet_dst); // @todo shouldn't this check also happen in other events?...

      bool split_event = false;
      for (const EdgeSPtr& edge_2 : edges_2) {
        FacetSPtr facet_l2 = edge_2->getFacetL();
        FacetSPtr facet_r2 = edge_2->getFacetR();
        FacetSPtr facet_2_src = edge_2->getFacetSrc();
        FacetSPtr facet_2_dst = edge_2->getFacetDst();

        Plane3SPtr plane_l2 = facet_l2->getPlane();
        CGAL_assertion_code(const FT& l2a = plane_l2->a();)
        CGAL_assertion_code(const FT& l2b = plane_l2->b();)
        CGAL_assertion_code(const FT& l2c = plane_l2->c();)
        CGAL_assertion_code(const FT& l2d = plane_l2->d();)
        CGAL_assertion_code(const FT& speed_l2 = HdsUtils::getSpeed(facet_l2);)
        CGAL_assertion_code(FT lt2 = (l2a * point->x() + l2b * point->y() + l2c * point->z() + l2d) / speed_l2);

        CGAL_assertion_code(Plane3SPtr plane_r2 = facet_r2->getPlane();)
        CGAL_assertion_code(const FT& r2a = plane_r2->a();)
        CGAL_assertion_code(const FT& r2b = plane_r2->b();)
        CGAL_assertion_code(const FT& r2c = plane_r2->c();)
        CGAL_assertion_code(const FT& r2d = plane_r2->d();)
        CGAL_assertion_code(const FT& speed_r2 = HdsUtils::getSpeed(facet_r2);)
        CGAL_assertion_code(FT rt2 = (r2a * point->x() + r2b * point->y() + r2c * point->z() + r2d) / speed_r2);

        CGAL_assertion(lt2 == rt2);
        CGAL_assertion(!is_positive(lt2));

        if (!checkBisectorsV2(edge_2, point, event_time)) {
          continue;
        }

        split_event = true;
        break;
      }

      if (split_event) {
        CGAL_SS3_CORE_TRACE_V(8, "Not an Edge Event (Split)");
        break;
      }

      // edge merge event
      EdgeSPtr edge_prev = edge->prev(facet_l);
      EdgeSPtr edge_next = edge->next(facet_l)->next(facet_l);
      if (edge_prev->hasSameFacets(edge_next)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not an Edge Event (Edge Merge #1)");
        break;
      }
      edge_prev = edge->prev(facet_l)->prev(facet_l);
      edge_next = edge->next(facet_l);
      if (edge_prev->hasSameFacets(edge_next)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not an Edge Event (Edge Merge #2)");
        break;
      }
      edge_prev = edge->prev(facet_r);
      edge_next = edge->next(facet_r)->next(facet_r);
      if (edge_prev->hasSameFacets(edge_next)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not an Edge Event (Edge Merge #3)");
        break;
      }
      edge_prev = edge->prev(facet_r)->prev(facet_r);
      edge_next = edge->next(facet_r);
      if (edge_prev->hasSameFacets(edge_next)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not an Edge Event (Edge Merge #4)");
        break;
      }

      // if here, it's an edge event
      EdgeEventSPtr edge_event = EdgeEvent::create();
      edge_event->setTime(event_time);
      edge_event->setPoint(point);
      edge_event->setEdge(edge);

      return handleEdgeEvent(edge_event, current_time, time_future_bound, polyhedron);
    }

    // EdgeMergeEvent
    for (;;)
    {
      FacetSPtr facet_l = edge->getFacetL();
      FacetSPtr facet_r = edge->getFacetR();
      if (HdsUtils::isTriangle(facet_l, edge) ||
          HdsUtils::isTriangle(facet_r, edge)) {
        // triangle event
        CGAL_SS3_CORE_TRACE_V(8, "Not an EdgeMerge Event (one incident facet is a triangle)");
        break;
      }

      FacetSPtr facet_other = edge->getFacetL();
      EdgeSPtr edge_next = edge->next(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->prev(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->next(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->prev(facet_other);
      if (edge_next == edge) {
        // dbl edge merge event
        CGAL_SS3_CORE_TRACE_V(8, "Not an EdgeMerge Event (Dbl #1)");
        break;
      }

      facet_other = edge->getFacetR();
      edge_next = edge->prev(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->next(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->prev(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->next(facet_other);
      if (edge_next == edge) {
        // dbl edge merge event
        CGAL_SS3_CORE_TRACE_V(8, "Not an EdgeMerge Event (Dbl #2)");
        break;
      }

      FacetSPtr facet = FacetSPtr();
      EdgeSPtr edge_1 = EdgeSPtr();
      EdgeSPtr edge_2 = EdgeSPtr();

      // @todo do we still need to test other combinations if a previous one matched?
      EdgeSPtr edge_prev = edge->prev(facet_l);
      edge_next = edge->next(facet_l)->next(facet_l);
      if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
        facet = facet_l;
        edge_1 = edge_prev;
        edge_2 = edge_next;
      }

      edge_prev = edge->prev(facet_l)->prev(facet_l);
      edge_next = edge->next(facet_l);
      if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
        facet = facet_l;
        edge_1 = edge_prev;
        edge_2 = edge_next;
      }

      edge_prev = edge->prev(facet_r);
      edge_next = edge->next(facet_r)->next(facet_r);
      if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
        facet = facet_r;
        edge_1 = edge_prev;
        edge_2 = edge_next;
      }

      edge_prev = edge->prev(facet_r)->prev(facet_r);
      edge_next = edge->next(facet_r);
      if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
        facet = facet_r;
        edge_1 = edge_prev;
        edge_2 = edge_next;
      }

      if (!(facet && edge_1 && edge_2)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not an EdgeMerge Event (Neighborhood)");
        break;
      }

      // if here, it's an edge merge event
      EdgeMergeEventSPtr edge_merge_event = EdgeMergeEvent::create();
      edge_merge_event->setTime(event_time);
      edge_merge_event->setPoint(point);
      edge_merge_event->setFacet(facet);
      edge_merge_event->setEdge1(edge_1);
      edge_merge_event->setEdge2(edge_2);

      return handleEdgeMergeEvent(edge_merge_event, current_time, time_future_bound, polyhedron);
    }

    // TriangleEvent
    for (;;)
    {
      if (HdsUtils::isTetrahedron(edge)) {
        // tetrahedron event
        CGAL_SS3_CORE_TRACE_V(8, "Not a Triangle Event (Tetrahedron)");
        break;
      }

      FacetSPtr facet;
      if (HdsUtils::isTriangle(edge->getFacetL(), edge)) {
        facet = edge->getFacetL();
      } else if (HdsUtils::isTriangle(edge->getFacetR(), edge)) {
        facet = edge->getFacetR();
      } else {
        CGAL_SS3_CORE_TRACE_V(8, "Not a Triangle Event (not triangle)");
        break;
      }

      bool dbl_triangle_event = false;
      EdgeSPtr edge_tmp = edge;
      for (unsigned int i = 0; i < 3; ++i) {
        FacetSPtr facet_tmp_l = edge_tmp->getFacetL();
        FacetSPtr facet_tmp_r = edge_tmp->getFacetR();
        if (facet_tmp_l && facet_tmp_r) {
          if (HdsUtils::isTriangle(facet_tmp_l, edge_tmp) &&
              HdsUtils::isTriangle(facet_tmp_r, edge_tmp)) {
            dbl_triangle_event = true;
            break;
          }
        }
        edge_tmp = edge_tmp->next(facet);
      }
      if (dbl_triangle_event) {
        CGAL_SS3_CORE_TRACE_V(8, "Not a Triangle Event (2 triangles)");
        break;
      }

      // if here, it's a triangle event
      TriangleEventSPtr triangle_event = TriangleEvent::create();
      triangle_event->setTime(event_time);
      triangle_event->setPoint(point);
      triangle_event->setFacet(facet);
      triangle_event->setEdgeBegin(edge);

      return handleTriangleEvent(triangle_event, current_time, time_future_bound, polyhedron);
    }

    // DblEdgeMergeEvent
    for (;;)
    {
      // At pop time, edge is degenerate, but by default isReflex() does a symbolic
      // shift into the future. However, here we want to know if the edge was reflex
      // BEFORE it became degenerate.
      if (!HdsUtils::isReflex(edge, false /*shift into the past*/)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not a DblEdgeMerge Event (not reflex)");
        break;
      }

      bool is_dbl_edge_merge_event = false;
      FacetSPtr facet_1 = FacetSPtr();
      EdgeSPtr edge_11 = EdgeSPtr();
      EdgeSPtr edge_12 = EdgeSPtr();
      FacetSPtr facet_2 = FacetSPtr();
      EdgeSPtr edge_21 = EdgeSPtr();
      EdgeSPtr edge_22 = EdgeSPtr();
      FacetSPtr facet_other = edge->getFacetL();
      EdgeSPtr edge_next = edge->next(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->prev(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->next(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->prev(facet_other);
      if (edge_next == edge) {
        is_dbl_edge_merge_event = true;
        facet_1 = edge->getFacetL();
        edge_11 = edge->prev(facet_1);
        edge_12 = edge->next(facet_1)->next(facet_1);
        facet_2 = edge->getFacetR();
        edge_21 = edge->prev(facet_2);
        edge_22 = edge->next(facet_2)->next(facet_2);
      }

      facet_other = edge->getFacetR();
      edge_next = edge->prev(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->next(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->prev(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->next(facet_other);
      if (edge_next == edge) {
        is_dbl_edge_merge_event = true;
        facet_1 = edge->getFacetR();
        edge_11 = edge->prev(facet_1)->prev(facet_1);
        edge_12 = edge->next(facet_1);
        facet_2 = edge->getFacetL();
        edge_21 = edge->prev(facet_2)->prev(facet_2);
        edge_22 = edge->next(facet_2);
      }

      if (edge_11 == edge_12 || edge_21 == edge_22) {
        // double triangle event
        CGAL_SS3_CORE_TRACE_V(8, "Not a DblEdgeMerge Event (DblTriangle)");
        break;
      }

      if (!is_dbl_edge_merge_event) {
        CGAL_SS3_CORE_TRACE_V(8, "Not a DblEdgeMerge Event");
        break;
      }

      // if here, it's a double edge merge event
      DblEdgeMergeEventSPtr dbl_edge_merge_event = DblEdgeMergeEvent::create();
      dbl_edge_merge_event->setTime(event_time);
      dbl_edge_merge_event->setPoint(point);
      dbl_edge_merge_event->setFacet1(facet_1);
      dbl_edge_merge_event->setEdge11(edge_11);
      dbl_edge_merge_event->setEdge12(edge_12);
      dbl_edge_merge_event->setFacet2(facet_2);
      dbl_edge_merge_event->setEdge21(edge_21);
      dbl_edge_merge_event->setEdge22(edge_22);

      return handleDblEdgeMergeEvent(dbl_edge_merge_event, current_time, time_future_bound, polyhedron);
    }

    // DblTriangleEvent
    for (;;)
    {
      if (HdsUtils::isTetrahedron(edge)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not a DblTriangleMerge Event (Tetrahedron)");
        break;
      }
      FacetSPtr facet_l = edge->getFacetL();
      FacetSPtr facet_r = edge->getFacetR();
      if (!facet_l || !facet_r) {
        CGAL_SS3_CORE_TRACE_V(8, "Not a DblTriangleMerge Event (neighborhood)");
        break;
      }
      if (!(HdsUtils::isTriangle(facet_l, edge) && HdsUtils::isTriangle(facet_r, edge))) {
        CGAL_SS3_CORE_TRACE_V(8, "Not a DblTriangleMerge Event (not triangles)");
        break;
      }

      // if here, it's a double triangle event
      DblTriangleEventSPtr dbl_triangle_event = DblTriangleEvent::create();
      dbl_triangle_event->setTime(event_time);
      dbl_triangle_event->setPoint(point);
      dbl_triangle_event->setEdge(edge);

      return handleDblTriangleEvent(dbl_triangle_event, current_time, time_future_bound, polyhedron);
    }

    // TetrahedronEvent
    for (;;)
    {
      if (!HdsUtils::isTetrahedron(edge)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not a Tetrahedron Event");
        break;
      }

      // if here, it's a tetrahedron event
      TetrahedronEventSPtr tetrahedron_event = TetrahedronEvent::create();
      tetrahedron_event->setTime(event_time);
      tetrahedron_event->setPoint(point);
      tetrahedron_event->setEdgeBegin(edge);

      return handleTetrahedronEvent(tetrahedron_event, current_time, time_future_bound, polyhedron);
    }

    return EventStatus::NON_EVENT;
  }

  EventStatus handleEdgeEvent(const EdgeEventSPtr& event,
                              const FT& current_time,
                              const std::optional<FT>& time_future_bound,
                              const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#########  Handle Edge Event  ##########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    EdgeSPtr edge = event->getEdge();
    VertexSPtr vertex_src = edge->getVertexSrc();
    VertexSPtr vertex_dst = edge->getVertexDst();

    CGAL_SS3_CORE_TRACE_V(4,"Edge:\n" << edge->toString());

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    HdsUtils::getArc(vertex_src)->closeArc(node);
    HdsUtils::getArc(vertex_dst)->closeArc(node);

    HdsUtils::getSheet(edge)->addNode(node);
    HdsUtils::clearSheet(edge);
#endif

    // grab vertices and facets counter clockwise around node
    std::array<VertexSPtr, 4> vertices;
    vertices[0] = vertex_src->prev(edge->getFacetL());
    vertices[1] = vertex_src->next(edge->getFacetR());
    vertices[2] = vertex_dst->prev(edge->getFacetR());
    vertices[3] = vertex_dst->next(edge->getFacetL());

    std::array<FacetSPtr, 4> facets;
    facets[1] = edge->getFacetR();
    facets[3] = edge->getFacetL();
    facets[0] = facets[3]->next(vertex_src);
    facets[2] = facets[1]->next(vertex_dst);

    // check if edge should be flipped
    bool flip_edge = true;

    // Try without a flip first. In practice, it is quite rare so try it first:
    // - checking for self-intersections exits as soon as one is found, so it is faster
    //   when it fails than when it succeeeds
    // - as soon as we fail for the unflipped one, we know that we have to flip and don't need
    //   to check for self-intersections.
    bool not_flipped_valid = false;
    {
      CGAL_SS3_CORE_TRACE_V(16, "== Trying WITHOUT a flip ==");

      /*
      NO FLIP:

            V0                       V3
            \                      /
            E0          F3       E3
              \                  /
      F0     src ------------ dst      F2
              /                 \
            E1         F1       E2
          /                      \
        V1                       V2

        dst of Ei = Vi
      */

      VertexSPtr vertex_src_clone = vertex_src->clone();
      VertexSPtr vertex_dst_clone = vertex_dst->clone();
      EdgeSPtr edge_no_flip = Edge::create(vertex_src_clone, vertex_dst_clone);
      std::array<EdgeSPtr, 4> edges;
      edges[0] = Edge::create(vertex_src_clone, vertices[0]->clone());
      edges[1] = Edge::create(vertex_src_clone, vertices[1]->clone());
      edges[2] = Edge::create(vertex_dst_clone, vertices[2]->clone());
      edges[3] = Edge::create(vertex_dst_clone, vertices[3]->clone());
      std::vector<FacetSPtr> facets_clone(4);
      for (unsigned int i = 0; i < 4; ++i) {
        facets_clone[i] = Facet::create(); // combinatorics built manually below
        facets_clone[i]->setPlane(facets[i]->getPlane());
        CGAL_assertion(facets[i]->hasData());
        facets_clone[i]->setData(facets[i]->getData()->clone(facets[i]));
      }

      edge_no_flip->setFacetR(facets_clone[1]);
      edge_no_flip->setFacetL(facets_clone[3]);

      // In edge events, we have unbounded faces (degree 1 vertices)
      facets_clone[3]->addEdge(edge_no_flip);
      facets_clone[1]->addEdge(edge_no_flip);
      for (unsigned int i = 0; i < 4; ++i) {
        edges[i]->setFacetR(facets_clone[(i+3)%4]);
        edges[i]->setFacetL(facets_clone[i]);
        facets_clone[(i+3)%4]->addEdge(edges[i]);
        facets_clone[i]->addEdge(edges[i]);
      }

      PolyhedronSPtr polyhedron_no_flip = Polyhedron::create(facets_clone);
      not_flipped_valid =
          (Transformation::shiftFacetsDegree1(polyhedron_no_flip, -1.0) &&
           !SelfIntersection::hasSelfIntersectingSurface(polyhedron_no_flip));
    }

    CGAL_SS3_CORE_TRACE_V(16, "not_flipped_valid = " << not_flipped_valid);

    bool flipped_valid = false;
#if 1
    // if one fails, the other one must be valid, and we can avoid testing for self-intersections
    if (!not_flipped_valid) {
      flipped_valid = true;
    } else
#endif
    {
      CGAL_SS3_CORE_TRACE_V(16, "== Trying WITH a flip ==");

      /*
      FLIP:

                            V0         V3
                F0          \   F3   /
                            E0     E3
                              \   /
              src ------------ dst
            /   \
          E1     E2
          /   F1   \        F2
        V1         V2

        dst of Ei = Vi
      */

      VertexSPtr vertex_src_clone = vertex_src->clone();
      VertexSPtr vertex_dst_clone = vertex_dst->clone();
      EdgeSPtr edge_flipped = Edge::create(vertex_src_clone, vertex_dst_clone);
      std::array<EdgeSPtr, 4> edges;
      edges[0] = Edge::create(vertex_dst_clone, vertices[0]->clone());
      edges[1] = Edge::create(vertex_src_clone, vertices[1]->clone());
      edges[2] = Edge::create(vertex_src_clone, vertices[2]->clone());
      edges[3] = Edge::create(vertex_dst_clone, vertices[3]->clone());
      std::vector<FacetSPtr> facets_clone(4);
      for (unsigned int i = 0; i < 4; ++i) {
        facets_clone[i] = Facet::create(); // combinatorics built manually below
        facets_clone[i]->setPlane(facets[i]->getPlane());
        CGAL_assertion(facets[i]->hasData());
        facets_clone[i]->setData(facets[i]->getData()->clone(facets[i]));
      }

      edge_flipped->setFacetR(facets_clone[2]);
      edge_flipped->setFacetL(facets_clone[0]);

      facets_clone[0]->addEdge(edge_flipped);
      facets_clone[2]->addEdge(edge_flipped);
      for (unsigned int i = 0; i < 4; ++i) {
        edges[i]->setFacetR(facets_clone[(i+3)%4]);
        edges[i]->setFacetL(facets_clone[i]);
        facets_clone[(i+3)%4]->addEdge(edges[i]);
        facets_clone[i]->addEdge(edges[i]);
      }

      PolyhedronSPtr polyhedron_flipped = Polyhedron::create(facets_clone);
      flipped_valid =
          (Transformation::shiftFacetsDegree1(polyhedron_flipped, -1.0) &&
           !SelfIntersection::hasSelfIntersectingSurface(polyhedron_flipped));
    }

    CGAL_SS3_CORE_TRACE_V(16, "flipped_valid = " << flipped_valid);

    if (flipped_valid && !not_flipped_valid) {
      flip_edge = true;
    } else if (not_flipped_valid && !flipped_valid) {
      flip_edge = false;
    } else if (flipped_valid && not_flipped_valid) {
      Vector3SPtr n0 = KernelFactory::createVector3(facets[0]->getPlane());
      Vector3SPtr n2 = KernelFactory::createVector3(facets[2]->getPlane());
      Vector3SPtr n1 = KernelFactory::createVector3(facets[1]->getPlane());
      Vector3SPtr n3 = KernelFactory::createVector3(facets[3]->getPlane());
      CGAL::Comparison_result dac = CGAL::compare_angle(*n0, *n2, *n1, *n3);

      if (edge_event_ == 0) {
        // convex
        flip_edge = (dac != CGAL::LARGER); // (angle_flipped <= angle_no_flip);
      } else if (edge_event_ == 1) {
        // reflex
        // choose edge that moves faster
        flip_edge = (dac != CGAL::SMALLER); // (angle_flipped >= angle_no_flip);
      } else {
        CGAL_assertion(edge_event_ == 2);
        // flip_when_possible
        flip_edge = true;
      }
    } else {
      throw std::runtime_error("Error: not able to handle EdgeEvent (2).");
    }

    std::array<EdgeSPtr, 4> edges;
    edges[0] = edge->next(vertex_src);
    edges[1] = edges[0]->next(vertex_src);
    edges[2] = edge->next(vertex_dst);
    edges[3] = edges[2]->next(vertex_dst);

    CGAL_SS3_CORE_TRACE_V(16, "flip_edge = " << flip_edge);

    if (flip_edge) {
      facets[3]->removeVertex(vertex_src);
      facets[2]->addVertex(vertex_src);
      facets[1]->removeVertex(vertex_dst);
      facets[0]->addVertex(vertex_dst);

      facets[1]->removeEdge(edge);
      facets[3]->removeEdge(edge);
      edge->setFacetL(facets[0]);
      edge->setFacetR(facets[2]);
      facets[0]->addEdge(edge);
      facets[2]->addEdge(edge);

      if (edges[0]->getVertexSrc() == vertex_src) {
        edges[0]->replaceVertexSrc(vertex_dst);
      } else if (edges[0]->getVertexDst() == vertex_src) {
        edges[0]->replaceVertexDst(vertex_dst);
      }
      if (edges[2]->getVertexSrc() == vertex_dst) {
        edges[2]->replaceVertexSrc(vertex_src);
      } else if (edges[2]->getVertexDst() == vertex_dst) {
        edges[2]->replaceVertexDst(vertex_src);
      }

      if (time_future_bound.has_value()) {
        HdsUtils::setFinalPoint(vertex_src, nullptr);
        HdsUtils::setFinalPoint(vertex_dst, nullptr);
      }

      post_op_vertices_VV_ = {{ vertex_src, vertex_dst }};

      // now, in the facets that have grown in size, we also need to collect
      // vertices that might create new events because they are now separated-enough
      // from other vertices in the cycle
      post_op_vertices_VV_.insert(vertex_src->prev(facets[0]));
      post_op_vertices_VV_.insert(vertex_dst->next(facets[0]));
      post_op_vertices_VV_.insert(vertex_src->next(facets[2]));
      post_op_vertices_VV_.insert(vertex_dst->prev(facets[2]));
      CGAL_assertion(post_op_vertices_VV_.size() == 6);

      // if there was no flip, then (unmodified) vertices should not have new events
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    for (int i = 0; i < 4; ++i) {
      HdsUtils::getSheet(edges[i])->addNode(node);
    }

    // update arcs and sheets
    for (VertexSPtr vertex_ext : { vertex_src, vertex_dst }) {
      HdsUtils::setNode(vertex_ext, node);

      ArcSPtr arc = createArc(vertex_ext);
      skel_result_->addArc(arc);

      // 'edge' is common to both extremities, and is lacking a new sheet at this point
      for (EdgeWPtr inc_edge_w : vertex_ext->edges()) {
        if (EdgeSPtr inc_edge = inc_edge_w.lock()) {
          if (inc_edge != edge) {
            HdsUtils::getSheet(inc_edge)->addArc(arc);
          }
        }
      }
    }

    // this sets up the two incident arcs, the nodes, etc.
    SheetSPtr sheet = createSheet(edge);
    skel_result_->addSheet(sheet);
#endif

    // Gather relevant elements for local queue updates
    post_op_vertices_ = {{ vertex_src, vertex_dst }};
    post_op_edges_ = {{ edge,
                        edge->next(vertex_src),
                        edge->prev(vertex_src),
                        edge->next(vertex_dst),
                        edge->prev(vertex_dst) }};
    post_op_facets_ = {{ facets[0], facets[1], facets[2], facets[3] }};
    CGAL_postcondition(post_op_vertices_.size() == 2 &&
                       post_op_edges_.size() == 5 &&
                       post_op_facets_.size() == 4);

    for (const EdgeSPtr& poe : post_op_edges_) {
        post_op_vertices_pierce_.insert(poe->getVertexSrc());
        post_op_vertices_pierce_.insert(poe->getVertexDst());
    }
    CGAL_postcondition(post_op_vertices_pierce_.size() == 6);

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handleEdgeMergeEvent(const EdgeMergeEventSPtr& event,
                                   const FT& current_time,
                                   const std::optional<FT>& time_future_bound,
                                   const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Edge Merge Event  #######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    FacetSPtr facet = event->getFacet();
    EdgeSPtr edge_1 = event->getEdge1();
    EdgeSPtr edge_2 = event->getEdge2();
    FacetSPtr facet_l = edge_1->getFacetL();
    FacetSPtr facet_r = edge_1->getFacetR();
    EdgeSPtr edge_toremove_1 = edge_1->next(facet);
    EdgeSPtr edge_toremove_2 = edge_toremove_1->next(facet);

    CGAL_assertion(edge_2 == edge_toremove_2->next(facet));

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    HdsUtils::getArc(edge_toremove_1->src(facet))->closeArc(node);
    HdsUtils::getArc(edge_toremove_1->dst(facet))->closeArc(node);
    HdsUtils::getArc(edge_toremove_2->dst(facet))->closeArc(node);

    HdsUtils::getSheet(edge_toremove_1)->addNode(node);
    HdsUtils::getSheet(edge_toremove_2)->addNode(node);
    HdsUtils::getSheet(edge_1)->addNode(node);
    // no need for edge_2, since edge_2's and edge_1's sheets are merged

    mergeSheets(edge_1, edge_2);
#endif

    VertexSPtr vertex = edge_toremove_1->dst(facet);
    vertex->setPoint(point); // @fixme unnecessary
    VertexSPtr vertex_1 = edge_1->dst(facet);
    VertexSPtr vertex_2 = edge_2->src(facet);
    EdgeSPtr edge_b = edge_toremove_1->prev(edge_toremove_1->other(facet));
    EdgeSPtr edge_b1 = edge_1->prev(edge_1->other(facet));
    EdgeSPtr edge_b2 = edge_2->next(edge_2->other(facet));
    facet->removeVertex(vertex);
    edge_1->other(facet)->addVertex(vertex);
    if (edge_b1->getVertexSrc() == vertex_1) {
      edge_b1->replaceVertexSrc(vertex);
    } else {
      edge_b1->replaceVertexDst(vertex);
    }
    if (edge_b2->getVertexSrc() == vertex_2) {
      edge_b2->replaceVertexSrc(vertex);
    } else {
      edge_b2->replaceVertexDst(vertex);
    }
    if (edge_1->getVertexDst() == vertex_1) {
      if (edge_2->getVertexSrc() == vertex_2) {
        edge_1->replaceVertexDst(edge_2->getVertexDst());
      } else {
        edge_1->replaceVertexDst(edge_2->getVertexSrc());
      }
    } else {
      if (edge_2->getVertexSrc() == vertex_2) {
        edge_1->replaceVertexSrc(edge_2->getVertexDst());
      } else {
        edge_1->replaceVertexSrc(edge_2->getVertexSrc());
      }
    }
    edge_toremove_1->getFacetL()->removeEdge(edge_toremove_1);
    edge_toremove_1->getFacetR()->removeEdge(edge_toremove_1);
    polyhedron->removeEdge(edge_toremove_1);
    edge_toremove_2->getFacetL()->removeEdge(edge_toremove_2);
    edge_toremove_2->getFacetR()->removeEdge(edge_toremove_2);
    polyhedron->removeEdge(edge_toremove_2);
    edge_2->getFacetL()->removeEdge(edge_2);
    edge_2->getFacetR()->removeEdge(edge_2);
    polyhedron->removeEdge(edge_2);
    for (auto it_f = vertex_1->facets().begin(); it_f != vertex_1->facets().end(); ) { // no C++11
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        facet->removeVertex(vertex_1);
      }
    }
    polyhedron->removeVertex(vertex_1);
    for (auto it_f = vertex_2->facets().begin(); it_f != vertex_2->facets().end(); ) { // no C++11
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        facet->removeVertex(vertex_2);
      }
    }
    polyhedron->removeVertex(vertex_2);

    if (time_future_bound.has_value()) {
      HdsUtils::setFinalPoint(vertex, nullptr);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    HdsUtils::setNode(vertex, node);

    ArcSPtr arc = createArc(vertex);
    skel_result_->addArc(arc);

    for (EdgeWPtr edge_w : vertex->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        HdsUtils::getSheet(edge)->addNode(node);
        HdsUtils::getSheet(edge)->addArc(arc);
      }
    }
#endif

    // Gather relevant elements for local queue updates
    post_op_vertices_ = {{ vertex }};
    post_op_edges_ = {{ edge_1, edge_b1, edge_b, edge_b2 }};
    post_op_facets_ = {{ facet_l, facet_r }};
    for (FacetWPtr wf : vertex->facets()) { post_op_facets_.insert(wf.lock()); }
    CGAL_postcondition(post_op_vertices_.size() == 1 && post_op_edges_.size() == 4 && post_op_facets_.size() == 4);

    // faces are smaller so nothing from unmodified vertices
    post_op_vertices_VV_ = {{ vertex }};

    CGAL_assertion(!HdsUtils::isReflex(vertex)); // just to see the configurations where this could not be the case
    post_op_vertices_pierce_.clear();

    // since all faces are getting smaller, we don't need to check unmodified edges
    post_op_edges_edgesplit_ = {{ edge_1, edge_b1, edge_b, edge_b2 }};
    CGAL_postcondition(post_op_edges_edgesplit_.size() == 4);

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handleTriangleEvent(const TriangleEventSPtr& event,
                                  const FT& current_time,
                                  const std::optional<FT>& time_future_bound,
                                  const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#######  Handle Triangle Event  ########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    VertexSPtr vertices[3];
    event->getVertices(vertices);

    FacetSPtr facet = event->getFacet();

    CGAL_SS3_CORE_TRACE_V(4, "Facet: " << facet->getID());
    CGAL_SS3_CORE_TRACE_V(4, "VS:\n" << vertices[0]->toString() << "\n"
                                     << vertices[1]->toString() << "\n"
                                     << vertices[2]->toString());

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    for (unsigned int i = 0; i < 3; ++i) {
      HdsUtils::getArc(vertices[i])->closeArc(node);
    }
    EdgeSPtr edges[3];
    event->getEdges(edges);
    for (unsigned int i = 0; i < 3; ++i) {
      HdsUtils::getSheet(edges[i])->addNode(node);
    }
#endif

    if (facet->vertices().size() == 3) {
      polyhedron->removeFacet(facet);
    }
    for (unsigned int i = 0; i < 3; ++i) {
      EdgeSPtr edge = vertices[i]->findEdge(vertices[(i+1)%3]);
      if (edge->getFacetL()) {
        edge->getFacetL()->removeEdge(edge);
      }
      if (edge->getFacetR()) {
        edge->getFacetR()->removeEdge(edge);
      }
      polyhedron->removeEdge(edge);
    }
    VertexSPtr new_vertex = Vertex::create(point);
    SkelVertexData::create(new_vertex);
    for (unsigned int i = 0; i < 3; ++i) {
      EdgeSPtr edge = vertices[i]->edges().front().lock();
      if (edge->getVertexSrc() == vertices[i]) {
        edge->replaceVertexSrc(new_vertex);
      } else if (edge->getVertexDst() == vertices[i]) {
        edge->replaceVertexDst(new_vertex);
      }
      edge->getFacetL()->removeVertex(vertices[i]);
      edge->getFacetR()->removeVertex(vertices[i]);
      facet->removeVertex(vertices[i]);
      polyhedron->removeVertex(vertices[i]);
      if (!edge->getFacetL()->containsVertex(new_vertex)) {
        edge->getFacetL()->addVertex(new_vertex);
      }
      if (!edge->getFacetR()->containsVertex(new_vertex)) {
        edge->getFacetR()->addVertex(new_vertex);
      }
    }
    polyhedron->addVertex(new_vertex);

    if (time_future_bound.has_value()) {
      HdsUtils::setFinalPoint(new_vertex, nullptr);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    HdsUtils::setNode(new_vertex, node);

    ArcSPtr arc = createArc(new_vertex);
    skel_result_->addArc(arc);

    for (EdgeWPtr edge_w : new_vertex->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        HdsUtils::getSheet(edge)->addNode(node);
        HdsUtils::getSheet(edge)->addArc(arc);
      }
    }

    CGAL_postcondition(arc->sheets().size() == 3);
#endif

    // Gather relevant elements for local queue updates
    post_op_vertices_ = {{ new_vertex }};
    for (EdgeWPtr we : new_vertex->edges()) {
      post_op_edges_.insert(EdgeSPtr(we.lock()));
    }
    for (FacetWPtr wf : new_vertex->facets()) {
      post_op_facets_.insert(wf.lock());
    }
    CGAL_postcondition(post_op_vertices_.size() == 1 && post_op_edges_.size() == 3 && post_op_facets_.size() == 3);

    // faces are smaller so nothing from unmodified vertices
    post_op_vertices_VV_ = {{ new_vertex }};

    // Assume a triangle event with a facet whose incident edges are all reflex.
    // At a reflex edge, an epsilon shift is a growth of the facet.
    // If all edges are reflex, the facet will grow and there could not have been
    // a triangle event. So at least one edge is convex and the resulting vertex
    // is not reflex. Since it is not reflex, we don't need to look at its pierce events.
    //
    // @todo checking for reflexness is (relatively) so cheap that we should insert it
    // anyway, in case the reasoning above is wrong!
    CGAL_assertion(!HdsUtils::isReflex(new_vertex));
    post_op_vertices_pierce_.clear();

    // faces are getting smaller so no need to check unmodified edges
    post_op_edges_edgesplit_ = post_op_edges_;
    CGAL_postcondition(post_op_edges_edgesplit_.size() == 3);

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handleDblEdgeMergeEvent(const DblEdgeMergeEventSPtr& event,
                                      const FT& current_time,
                                      const std::optional<FT>& /*time_future_bound*/,
                                      const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#######  Handle Dbl Edge Event  ########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    EdgeSPtr edge_11 = event->getEdge11();
    EdgeSPtr edge_12 = event->getEdge12();
    EdgeSPtr edge_21 = event->getEdge21();
    EdgeSPtr edge_22 = event->getEdge22();

    VertexSPtr vertices[4];
    event->getVertices(vertices);
    EdgeSPtr edges[4];
    event->getEdges(edges);

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    for (unsigned int i = 0; i < 4; ++i) {
      HdsUtils::getArc(vertices[i])->closeArc(node);
    }
    for (unsigned int i = 0; i < 4; ++i) {
      HdsUtils::getSheet(edges[i])->addNode(node);
    }

    HdsUtils::getSheet(edge_11)->addNode(node);
    HdsUtils::getSheet(edge_21)->addNode(node);
    // no need for _2 because we merge sheets

    mergeSheets(edge_11, edge_12);
    mergeSheets(edge_21, edge_22);
#endif

    for (unsigned int i = 0; i < 4; ++i) {
      EdgeSPtr edge = edges[i];
      edge->getFacetL()->removeEdge(edge);
      edge->getFacetR()->removeEdge(edge);
      polyhedron->removeEdge(edge);
    }

    if (edge_11->getVertexDst() == vertices[0]) {
      if (edge_12->getVertexSrc() == vertices[2]) {
        edge_11->replaceVertexDst(edge_12->getVertexDst());
      } else {
        edge_11->replaceVertexDst(edge_12->getVertexSrc());
      }
    } else {
      if (edge_12->getVertexSrc() == vertices[2]) {
        edge_11->replaceVertexSrc(edge_12->getVertexDst());
      } else {
        edge_11->replaceVertexSrc(edge_12->getVertexSrc());
      }
    }
    edge_12->getFacetL()->removeEdge(edge_12);
    edge_12->getFacetR()->removeEdge(edge_12);
    polyhedron->removeEdge(edge_12);
    if (edge_21->getVertexDst() == vertices[1]) {
      if (edge_22->getVertexSrc() == vertices[3]) {
        edge_21->replaceVertexDst(edge_22->getVertexDst());
      } else {
        edge_21->replaceVertexDst(edge_22->getVertexSrc());
      }
    } else {
      if (edge_22->getVertexSrc() == vertices[3]) {
        edge_21->replaceVertexSrc(edge_22->getVertexDst());
      } else {
        edge_21->replaceVertexSrc(edge_22->getVertexSrc());
      }
    }
    edge_22->getFacetL()->removeEdge(edge_22);
    edge_22->getFacetR()->removeEdge(edge_22);
    polyhedron->removeEdge(edge_22);
    for (unsigned int i = 0; i < 4; ++i) {
      VertexSPtr vertex = vertices[i];
      for (auto it_f = vertex->facets().begin(); it_f != vertex->facets().end(); ) { // no C++11
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          facet->removeVertex(vertex);
        }
      }
      polyhedron->removeVertex(vertex);
    }

    // Gather relevant elements for local queue updates
    post_op_vertices_.clear();
    post_op_edges_ = {{ edge_11, edge_21 }};
    post_op_facets_ = {{ edge_11->getFacetL(),
                         edge_11->getFacetR(),
                         edge_21->getFacetL(),
                         edge_21->getFacetR() }};
    CGAL_postcondition(post_op_vertices_.empty() && post_op_edges_.size() == 2 && post_op_facets_.size() == 4);

    // faces are smaller so nothing from unmodified vertices
    post_op_vertices_VV_.clear();

    // no new vertices & only reducing the size of facets so no edge disconnection
    post_op_vertices_pierce_.clear();

    CGAL_assertion(!HdsUtils::isReflex(edge_11));
    CGAL_assertion(!HdsUtils::isReflex(edge_21));
    post_op_edges_edgesplit_.clear();

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handleDblTriangleEvent(const DblTriangleEventSPtr& event,
                                     const FT& current_time,
                                     const std::optional<FT>& /*time_future_bound*/,
                                     const PolyhedronSPtr& polyhedron) {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#####  Handle Dbl Triangle Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    EdgeSPtr edge = event->getEdge();

    VertexSPtr vertices[4];
    event->getVertices(vertices);

    EdgeSPtr edges[5];
    event->getEdges(edges);

    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    VertexSPtr vertex_l = edge->next(facet_l)->dst(facet_l);
    VertexSPtr vertex_r = edge->next(facet_r)->dst(facet_r);
    EdgeSPtr edge_l = edge->next(facet_l);
    FacetSPtr facet_ll = edge_l->other(facet_l);
    edge_l = edge_l->prev(facet_ll);
    EdgeSPtr edge_r = edge->next(facet_r);
    FacetSPtr facet_rr = edge_r->other(facet_r);
    edge_r = edge_r->prev(facet_rr);

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    for (unsigned int i = 0; i < 4; ++i) {
      HdsUtils::getArc(vertices[i])->closeArc(node);
    }
    for (unsigned int i = 0; i < 5; ++i) {
      HdsUtils::getSheet(edges[i])->addNode(node);
    }

    HdsUtils::getSheet(edge_l)->addNode(node);

    mergeSheets(edge_l, edge_r);
#endif

    if (facet_l->edges().size() == 3) {
      polyhedron->removeFacet(facet_l);
    }
    if (facet_r->edges().size() == 3) {
      polyhedron->removeFacet(facet_r);
    }
    for (unsigned int i = 0; i < 5; ++i) {
      EdgeSPtr edge = edges[i];
      if (edge->getFacetL()) {
        edge->getFacetL()->removeEdge(edge);
      }
      if (edge->getFacetR()) {
        edge->getFacetR()->removeEdge(edge);
      }
      polyhedron->removeEdge(edge);
    }

    if (edge_l->getVertexSrc() == vertex_l) {
      if (edge_r->getVertexDst() == vertex_r) {
        edge_l->replaceVertexSrc(edge_r->getVertexSrc());
      } else {
        edge_l->replaceVertexSrc(edge_r->getVertexDst());
      }
    } else {
      if (edge_r->getVertexDst() == vertex_r) {
        edge_l->replaceVertexDst(edge_r->getVertexSrc());
      } else {
        edge_l->replaceVertexDst(edge_r->getVertexDst());
      }
    }
    edge_r->getFacetL()->removeEdge(edge_r);
    edge_r->getFacetR()->removeEdge(edge_r);
    polyhedron->removeEdge(edge_r);

    for (unsigned int i = 0; i < 4; ++i) {
      VertexSPtr vertex = vertices[i];
      for (auto it_f = vertex->facets().begin(); it_f != vertex->facets().end(); ) { // no C++11
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          facet->removeVertex(vertex);
        }
      }
      polyhedron->removeVertex(vertex);
    }

    // Gather relevant elements for local queue updates
    post_op_vertices_.clear();
    post_op_edges_ = {{ edge_l }};
    post_op_facets_ = {{ facet_ll, facet_rr }};
    CGAL_postcondition(post_op_vertices_.empty() && post_op_edges_.size() == 1 && post_op_facets_.size() == 2);

    // faces are smaller so nothing from unmodified vertices
    post_op_vertices_VV_.clear();

    // no new vertices & only reducing the size of facets so no edge disconnection
    post_op_vertices_pierce_.clear();

    // faces are getting smaller so no need to check unmodified edges
    post_op_edges_edgesplit_ = {{ edge_l }};

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handleTetrahedronEvent(const TetrahedronEventSPtr& event,
                                     const FT& current_time,
                                     const std::optional<FT>& /*time_future_bound*/,
                                     const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Tetrahedron Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    VertexSPtr vertices[4];
    event->getVertices(vertices);

    EdgeSPtr edges[6];
    event->getEdges(edges);

    FacetSPtr facets[4];
    event->getFacets(facets);

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    for (unsigned int i = 0; i < 4; ++i) {
      HdsUtils::getArc(vertices[i])->closeArc(node);
    }
    for (unsigned int i = 0; i < 6; ++i) {
      HdsUtils::getSheet(edges[i])->addNode(node);
    }
#endif

    for (unsigned int i = 0; i < 4; ++i) {
      if (facets[i]->vertices().size() == 3) {
        polyhedron->removeFacet(facets[i]);
      }
    }
    for (unsigned int i = 0; i < 6; ++i) {
      EdgeSPtr edge = edges[i];
      if (edge->getFacetL()) {
        edge->getFacetL()->removeEdge(edge);
      }
      if (edge->getFacetR()) {
        edge->getFacetR()->removeEdge(edge);
      }
      polyhedron->removeEdge(edge);
    }
    for (unsigned int i = 0; i < 4; ++i) {
      VertexSPtr vertex = vertices[i];
      for (auto it_f = vertex->facets().begin(); it_f != vertex->facets().end(); ) { // no C++11
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          facet->removeVertex(vertex);
        }
      }
      polyhedron->removeVertex(vertex);
    }

    // Gather relevant elements for local queue updates
    post_op_vertices_.clear();
    post_op_edges_.clear();
    post_op_facets_.clear();
    CGAL_postcondition(post_op_vertices_.empty() && post_op_edges_.empty() && post_op_facets_.empty());

    // no new vertices & only reducing the size of facets so no edge disconnection
    post_op_vertices_VV_.clear();
    post_op_vertices_pierce_.clear();
    post_op_edges_edgesplit_.clear();

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handleVertexEvent(const VertexEventSPtr& event,
                                const FT& current_time,
                                const std::optional<FT>& time_future_bound,
                                const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "########  Handle Vertex Event  #########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    VertexSPtr vertex_1 = event->getVertex1();
    VertexSPtr vertex_2 = event->getVertex2();
    FacetSPtr facet_1 = event->getFacet1();
    FacetSPtr facet_2 = event->getFacet2();

    EdgeSPtr edge_tomerge_1 = EdgeSPtr();
    EdgeSPtr edge_11 = EdgeSPtr();
    EdgeSPtr edge_12 = EdgeSPtr();
    EdgeSPtr edge_tomerge_2 = EdgeSPtr();
    EdgeSPtr edge_21 = EdgeSPtr();
    EdgeSPtr edge_22 = EdgeSPtr();
    for (EdgeWPtr edge_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->getFacetL() == facet_1 && edge->getFacetR() == facet_2) ||
            (edge->getFacetL() == facet_2 && edge->getFacetR() == facet_1)) {
          edge_tomerge_1 = edge;
          continue;
        }
        if (edge->getFacetL() == facet_1 || edge->getFacetR() == facet_1) {
          edge_11 = edge;
        }
        if (edge->getFacetL() == facet_2 || edge->getFacetR() == facet_2) {
          edge_12 = edge;
        }
      }
    }
    for (EdgeWPtr edge_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->getFacetL() == facet_1 && edge->getFacetR() == facet_2) ||
            (edge->getFacetL() == facet_2 && edge->getFacetR() == facet_1)) {
          edge_tomerge_2 = edge;
          continue;
        }
        if (edge->getFacetL() == facet_1 || edge->getFacetR() == facet_1) {
          edge_21 = edge;
        }
        if (edge->getFacetL() == facet_2 || edge->getFacetR() == facet_2) {
          edge_22 = edge;
        }
      }
    }
    FacetSPtr facet_1b = edge_11->getFacetL();
    if (facet_1b == facet_1 || facet_1b == facet_2) {
      facet_1b = edge_11->getFacetR();
    }
    FacetSPtr facet_2b = edge_21->getFacetL();
    if (facet_2b == facet_1 || facet_2b == facet_2) {
      facet_2b = edge_21->getFacetR();
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    HdsUtils::getArc(vertex_1)->closeArc(node);
    HdsUtils::getArc(vertex_2)->closeArc(node);

    HdsUtils::getSheet(edge_11)->addNode(node);
    HdsUtils::getSheet(edge_12)->addNode(node);
    HdsUtils::getSheet(edge_21)->addNode(node);
    HdsUtils::getSheet(edge_22)->addNode(node);

    HdsUtils::getSheet(edge_tomerge_1)->addNode(node);
    // no need it to edge_tomerge_2's sheet because that sheet gets merged with edge_tomerge_1's

    mergeSheets(edge_tomerge_1, edge_tomerge_2);
#endif

    if (edge_tomerge_1->getVertexSrc() == vertex_1) {
      if (edge_tomerge_2->getVertexSrc() == vertex_2) {
        edge_tomerge_1->replaceVertexSrc(edge_tomerge_2->getVertexDst());
      } else {
        edge_tomerge_1->replaceVertexSrc(edge_tomerge_2->getVertexSrc());
      }
    } else {
      if (edge_tomerge_2->getVertexSrc() == vertex_2) {
        edge_tomerge_1->replaceVertexDst(edge_tomerge_2->getVertexDst());
      } else {
        edge_tomerge_1->replaceVertexDst(edge_tomerge_2->getVertexSrc());
      }
    }
    facet_1->removeVertex(vertex_2);
    facet_2->removeVertex(vertex_1);
    facet_1b->addVertex(vertex_2);
    facet_2b->addVertex(vertex_1);
    edge_tomerge_2->replaceVertexSrc(vertex_1);
    edge_tomerge_2->replaceVertexDst(vertex_2);
    edge_tomerge_2->replaceFacetL(facet_1b);
    edge_tomerge_2->replaceFacetR(facet_2b);
    if (edge_12->getVertexSrc() == vertex_1) {
      edge_12->replaceVertexSrc(vertex_2);
    } else {
      edge_12->replaceVertexDst(vertex_2);
    }
    if (edge_21->getVertexSrc() == vertex_2) {
      edge_21->replaceVertexSrc(vertex_1);
    } else {
      edge_21->replaceVertexDst(vertex_1);
    }

    if (time_future_bound.has_value()) {
      HdsUtils::setFinalPoint(vertex_1, nullptr);
      HdsUtils::setFinalPoint(vertex_2, nullptr);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    HdsUtils::clearSheet(edge_tomerge_2);

    HdsUtils::setNode(vertex_1, node);
    HdsUtils::setNode(vertex_2, node);

    ArcSPtr arc_1 = createArc(vertex_1);
    skel_result_->addArc(arc_1);
    ArcSPtr arc_2 = createArc(vertex_2);
    skel_result_->addArc(arc_2);

    HdsUtils::getSheet(edge_11)->addArc(arc_1);
    HdsUtils::getSheet(edge_21)->addArc(arc_1);
    HdsUtils::getSheet(edge_12)->addArc(arc_2);
    HdsUtils::getSheet(edge_22)->addArc(arc_2);

    SheetSPtr sheet = createSheet(edge_tomerge_2);
    skel_result_->addSheet(sheet);

    CGAL_postcondition(arc_1->sheets().size() == 3);
    CGAL_postcondition(arc_2->sheets().size() == 3);
#endif

    // Gather relevant elements for local queue updates
    post_op_vertices_ = {{ vertex_1, vertex_2 }};
    post_op_edges_ = {{ edge_tomerge_1, edge_12, edge_21, edge_tomerge_2, edge_11, edge_22 }};
    for (FacetWPtr wf : vertex_1->facets()) { post_op_facets_.insert(wf.lock()); }
    for (FacetWPtr wf : vertex_2->facets()) { post_op_facets_.insert(wf.lock()); }
    CGAL_postcondition(post_op_vertices_.size() == 2 && post_op_edges_.size() == 6 && post_op_facets_.size() == 4);

    // facets incident to 'edge_merge_2' have grown
    post_op_vertices_VV_ = {{ vertex_1, vertex_2 }};
    post_op_vertices_VV_.insert(vertex_1->prev(facet_1b));
    post_op_vertices_VV_.insert(vertex_2->next(facet_1b));
    // post_op_vertices_VV_.insert(vertex_1->next(facet_2b)); redundant?
    // post_op_vertices_VV_.insert(vertex_2->prev(facet_2b));

    // probably could improve this but _vertex events are rare_, so we are rarely
    // looking for new pierce events after a vertex event so it doesn't matter much
    for (const EdgeSPtr& poe : post_op_edges_) {
      post_op_vertices_pierce_.insert(poe->getVertexSrc());
      post_op_vertices_pierce_.insert(poe->getVertexDst());
    }

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handleFlipVertexEvent(const FlipVertexEventSPtr& event,
                                    const FT& current_time,
                                    const std::optional<FT>& time_future_bound,
                                    const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Flip Vertex Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    VertexSPtr vertex_1 = event->getVertex1();
    VertexSPtr vertex_2 = event->getVertex2();
    FacetSPtr facet_1 = event->getFacet1();
    FacetSPtr facet_2 = event->getFacet2();

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    HdsUtils::getArc(vertex_1)->closeArc(node);
    HdsUtils::getArc(vertex_2)->closeArc(node);

    for (VertexSPtr vertex : { vertex_1, vertex_2 }) {
      for (EdgeWPtr edge_w : vertex->edges()) {
        if (EdgeSPtr edge = edge_w.lock()) {
          HdsUtils::getSheet(edge)->addNode(node);
        }
      }
    }
#endif

    EdgeSPtr edge_1;
    for (EdgeWPtr edge_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->getFacetL() == facet_1 && edge->getFacetR() == facet_2) ||
            (edge->getFacetL() == facet_2 && edge->getFacetR() == facet_1)) {
          edge_1 = edge;
          break;
        }
      }
    }
    EdgeSPtr edge_2;
    for (EdgeWPtr edge_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->getFacetL() == facet_1 && edge->getFacetR() == facet_2) ||
            (edge->getFacetL() == facet_2 && edge->getFacetR() == facet_1)) {
          edge_2 = edge;
          break;
        }
      }
    }

    if (edge_1->getVertexSrc() == vertex_1) {
      edge_1->replaceVertexSrc(vertex_2);
    } else if (edge_1->getVertexDst() == vertex_1) {
      edge_1->replaceVertexDst(vertex_2);
    }
    if (edge_2->getVertexSrc() == vertex_2) {
      edge_2->replaceVertexSrc(vertex_1);
    } else if (edge_2->getVertexDst() == vertex_2) {
      edge_2->replaceVertexDst(vertex_1);
    }

    if (time_future_bound.has_value()) {
      HdsUtils::setFinalPoint(vertex_1, nullptr);
      HdsUtils::setFinalPoint(vertex_2, nullptr);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    HdsUtils::setNode(vertex_1, node);
    HdsUtils::setNode(vertex_2, node);

    ArcSPtr arc_1 = createArc(vertex_1);
    skel_result_->addArc(arc_1);
    ArcSPtr arc_2 = createArc(vertex_2);
    skel_result_->addArc(arc_2);

    for (VertexSPtr vertex : { vertex_1, vertex_2 }) {
      for (EdgeWPtr edge_w : vertex->edges()) {
        if (EdgeSPtr edge = edge_w.lock()) {
          HdsUtils::getSheet(edge)->addArc(HdsUtils::getArc(vertex));
        }
      }
    }
#endif

    // Gather relevant elements for local queue updates
    post_op_vertices_ = {{ vertex_1, vertex_2 }};
    for (EdgeWPtr we : vertex_1->edges()) { post_op_edges_.insert(EdgeSPtr(we.lock())); }
    for (EdgeWPtr we : vertex_2->edges()) { post_op_edges_.insert(EdgeSPtr(we.lock())); }
    for (FacetWPtr wf : vertex_1->facets()) { post_op_facets_.insert(wf.lock()); }
    for (FacetWPtr wf : vertex_2->facets()) { post_op_facets_.insert(wf.lock()); }
    CGAL_postcondition(post_op_vertices_.size() == 2 && post_op_edges_.size() == 6 && post_op_facets_.size() == 4);

    // facets did not grow
    post_op_vertices_VV_ = {{ vertex_1, vertex_2 }};

    // probably could improve this but flip vertex events are uncommon, so we are rarely
    // looking for new pierce events after a vertex event, and it does not matter much
    for (const EdgeSPtr& poe : post_op_edges_) {
      post_op_vertices_pierce_.insert(poe->getVertexSrc());
      post_op_vertices_pierce_.insert(poe->getVertexDst());
    }

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handleSurfaceEvent(const SurfaceEventSPtr& event,
                                 const FT& current_time,
                                 const std::optional<FT>& time_future_bound,
                                 const PolyhedronSPtr& polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#######  Handle Surface Event  #########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    CGAL_SS3_CORE_TRACE_V(4, "Edge A = " << event->getEdge1()->toString());
    CGAL_SS3_CORE_TRACE_V(4, "Edge B = " << event->getEdge2()->toString());

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    EdgeSPtr edge_1 = event->getEdge1();
    EdgeSPtr edge_2 = event->getEdge2();
    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();

    VertexSPtr vertex = VertexSPtr();
    EdgeSPtr edge_b1 = EdgeSPtr();
    EdgeSPtr edge_b2 = EdgeSPtr();
    if (edge_2->getFacetL() == facet_1_src) {
      vertex = edge_1->getVertexSrc();
      edge_b1 = edge_1->prev(edge_1->getFacetL());
      edge_b2 = edge_1->next(edge_1->getFacetR());
    } else if (edge_2->getFacetR() == facet_1_src) {
      vertex = edge_1->getVertexSrc();
      edge_b1 = edge_1->next(edge_1->getFacetR());
      edge_b2 = edge_1->prev(edge_1->getFacetL());
    } else if (edge_2->getFacetL() == facet_1_dst) {
      vertex = edge_1->getVertexDst();
      edge_b1 = edge_1->prev(edge_1->getFacetR());
      edge_b2 = edge_1->next(edge_1->getFacetL());
    } else if (edge_2->getFacetR() == facet_1_dst) {
      vertex = edge_1->getVertexDst();
      edge_b1 = edge_1->next(edge_1->getFacetL());
      edge_b2 = edge_1->prev(edge_1->getFacetR());
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    if (facet_1_src == edge_2->getFacetL() || facet_1_src == edge_2->getFacetR()) {
      HdsUtils::getArc(edge_1->getVertexSrc())->closeArc(node);
    }
    if (facet_1_dst == edge_2->getFacetL() || facet_1_dst == edge_2->getFacetR()) {
      HdsUtils::getArc(edge_1->getVertexDst())->closeArc(node);
    }

    for (EdgeWPtr edge_w : vertex->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        HdsUtils::getSheet(edge)->addNode(node);
      }
    }

    HdsUtils::getSheet(edge_2)->addNode(node);
#endif

    VertexSPtr vertex_21 = Vertex::create(point);
    VertexSPtr vertex_22 = Vertex::create(point);
    SkelVertexData::create(vertex_21);
    SkelVertexData::create(vertex_22);
    polyhedron->addVertex(vertex_21);
    polyhedron->addVertex(vertex_22);
    if (edge_b1->getVertexSrc() == vertex) {
      edge_b1->replaceVertexSrc(vertex_21);
    } else if (edge_b1->getVertexDst() == vertex) {
      edge_b1->replaceVertexDst(vertex_21);
    }
    if (edge_b2->getVertexSrc() == vertex) {
      edge_b2->replaceVertexSrc(vertex_22);
    } else if (edge_b2->getVertexDst() == vertex) {
      edge_b2->replaceVertexDst(vertex_22);
    }
    edge_b1->getFacetL()->addVertex(vertex_21);
    edge_b1->getFacetR()->addVertex(vertex_21);
    edge_b2->getFacetL()->addVertex(vertex_22);
    edge_b2->getFacetR()->addVertex(vertex_22);

    EdgeSPtr edge_tmp = edge_2->split(vertex);
    EdgeSPtr edge_21 = edge_2->split(vertex_21);
    EdgeSPtr edge_22 = edge_tmp;
    edge_tmp = edge_22->split(vertex_22);

    SkelEdgeData::create(edge_tmp);
    SkelEdgeData::create(edge_21);
    SkelEdgeData::create(edge_22);

    if (edge_2->getFacetL() == facet_1_src ||
        edge_2->getFacetL() == facet_1_dst) {
      edge_2->getFacetL()->removeVertex(vertex);
      if (vertex == edge_1->getVertexSrc()) {
        edge_21->replaceFacetL(edge_1->getFacetL());
        edge_22->replaceFacetL(edge_1->getFacetR());
      } else if (vertex == edge_1->getVertexDst()) {
        edge_21->replaceFacetL(edge_1->getFacetR());
        edge_22->replaceFacetL(edge_1->getFacetL());
      }
    } else {
      edge_2->getFacetR()->removeVertex(vertex);
      if (vertex == edge_1->getVertexSrc()) {
        edge_21->replaceFacetR(edge_1->getFacetR());
        edge_22->replaceFacetR(edge_1->getFacetL());
      } else if (vertex == edge_1->getVertexDst()) {
        edge_21->replaceFacetR(edge_1->getFacetL());
        edge_22->replaceFacetR(edge_1->getFacetR());
      }
    }

    if (time_future_bound.has_value()) {
      HdsUtils::setFinalPoint(vertex, nullptr);
      HdsUtils::setFinalPoint(vertex_21, nullptr);
      HdsUtils::setFinalPoint(vertex_22, nullptr);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    HdsUtils::setSheet(edge_tmp, HdsUtils::getSheet(edge_2));

    HdsUtils::setNode(vertex, node);
    ArcSPtr arc = createArc(vertex);
    skel_result_->addArc(arc);
    HdsUtils::getSheet(edge_1)->addArc(arc); // other 2 sheets are edge_21/edge_22's, created below

    HdsUtils::setNode(vertex_21, node);
    ArcSPtr arc_21 = createArc(vertex_21);
    skel_result_->addArc(arc_21);
    HdsUtils::getSheet(edge_b1)->addArc(arc_21);
    HdsUtils::getSheet(edge_2)->addArc(arc_21); // third sheet is edge_21's, created below

    HdsUtils::setNode(vertex_22, node);
    ArcSPtr arc_22 = createArc(vertex_22);
    skel_result_->addArc(arc_22);
    HdsUtils::getSheet(edge_b2)->addArc(arc_22);
    HdsUtils::getSheet(edge_tmp)->addArc(arc_22); // third sheet is edge_22's, created below

    SheetSPtr sheet_21 = createSheet(edge_21);
    skel_result_->addSheet(sheet_21);

    SheetSPtr sheet_22 = createSheet(edge_22);
    skel_result_->addSheet(sheet_22);
#endif

    // Gather relevant elements for local queue updates
    post_op_vertices_ = {{ vertex, vertex_21, vertex_22 }};
    for (const VertexSPtr& v : post_op_vertices_) {
      for (EdgeWPtr edge_w : v->edges()) {
        if (EdgeSPtr edge = edge_w.lock()) {
          post_op_edges_.insert(edge);
        }
      }
    }

    post_op_facets_ = {{ edge_1->getFacetL(), edge_1->getFacetR(),
                         edge_2->getFacetL(), edge_2->getFacetR() }};

    // not sure if we need that much stuff
    post_op_vertices_VV_ = {{ vertex, vertex_21, vertex_22 }};
    for (const VertexSPtr& v : { vertex, vertex_21, vertex_22 }) {
      for (EdgeWPtr we : v->edges()) {
        EdgeSPtr e = we.lock();
        post_op_vertices_VV_.insert(e->other(v));
      }
    }

    CGAL_postcondition(post_op_vertices_.size() == 3 &&
                       post_op_edges_.size() == 7 &&
                       post_op_facets_.size() == 4);

    // @speed the vertex at the top of the reflex edge is not a modified vertex, so the only
    // event that could appear is with the facet it just go disconnected from (i.e.,
    // edge_2->other_face)?
    post_op_vertices_pierce_ = {{ edge_1->getVertexSrc(), edge_1->getVertexDst(), vertex_21, vertex_22 }};
    CGAL_postcondition(post_op_vertices_pierce_.size() == 4);

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handlePolyhedronSplitEvent(const PolyhedronSplitEventSPtr& event,
                                         const FT& current_time,
                                         const std::optional<FT>& time_future_bound,
                                         const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "####  Handle Polyhedron Split Event  ###");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    EdgeSPtr edge_1 = event->getEdge1();
    EdgeSPtr edge_2 = event->getEdge2();
    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();

    CGAL_SS3_CORE_TRACE_V(4, "Edge 1 = " << edge_1->toString());
    CGAL_SS3_CORE_TRACE_V(4, "Edge 2 = " << edge_2->toString());

    VertexSPtr vertex_l = VertexSPtr();
    VertexSPtr vertex_r = VertexSPtr();
    if (edge_2->getFacetL() == facet_1_src &&
        edge_2->getFacetR() == facet_1_dst) {
      vertex_l = edge_1->getVertexSrc();
      vertex_r = edge_1->getVertexDst();
    }
    if (edge_2->getFacetL() == facet_1_dst &&
        edge_2->getFacetR() == facet_1_src) {
      vertex_l = edge_1->getVertexDst();
      vertex_r = edge_1->getVertexSrc();
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    if (facet_1_src == edge_2->getFacetL() || facet_1_src == edge_2->getFacetR()) {
      HdsUtils::getArc(edge_1->getVertexSrc())->closeArc(node);
    }
    if (facet_1_dst == edge_2->getFacetL() || facet_1_dst == edge_2->getFacetR()) {
      HdsUtils::getArc(edge_1->getVertexDst())->closeArc(node);
    }

    std::array<EdgeSPtr, 4> edges;
    edges[0] = edge_1->next(vertex_l);
    edges[1] = edges[0]->next(vertex_l);
    edges[2] = edge_1->next(vertex_r);
    edges[3] = edges[2]->next(vertex_r);

    HdsUtils::getSheet(edge_1)->addNode(node);
    for (int i = 0; i < 4; ++i) {
      HdsUtils::getSheet(edges[i])->addNode(node);
    }
    HdsUtils::getSheet(edge_2)->addNode(node);
#endif

    EdgeSPtr edge_22;
    if (edge_2->getFacetL() == facet_1_src &&
        edge_2->getFacetR() == facet_1_dst) {
      EdgeSPtr edge_l = edge_1->prev(edge_1->getVertexDst());
      if (edge_l->getVertexSrc() == vertex_r) {
        edge_l->replaceVertexSrc(vertex_l);
      } else if (edge_l->getVertexDst() == vertex_r) {
        edge_l->replaceVertexDst(vertex_l);
      }
      EdgeSPtr edge_r = edge_1->prev(edge_1->getVertexSrc());
      if (edge_r->getVertexSrc() == vertex_l) {
        edge_r->replaceVertexSrc(vertex_r);
      } else if (edge_r->getVertexDst() == vertex_l) {
        edge_r->replaceVertexDst(vertex_r);
      }

      edge_1->getFacetR()->removeVertex(vertex_l);
      edge_2->getFacetR()->addVertex(vertex_l);
      edge_1->getFacetL()->removeVertex(vertex_r);
      edge_2->getFacetL()->addVertex(vertex_r);

      edge_22 = edge_1;
      edge_22->replaceVertexDst(edge_2->getVertexDst());
      edge_2->replaceVertexDst(vertex_l);
      edge_22->replaceVertexSrc(vertex_r);

      edge_22->replaceFacetL(edge_2->getFacetL());
      edge_22->replaceFacetR(edge_2->getFacetR());
    }
    // @todo just else it...
    if (edge_2->getFacetL() == facet_1_dst &&
        edge_2->getFacetR() == facet_1_src) {
      vertex_l = edge_1->getVertexDst();
      vertex_r = edge_1->getVertexSrc();

      EdgeSPtr edge_l = edge_1->next(edge_1->getVertexSrc());
      if (edge_l->getVertexSrc() == vertex_r) {
        edge_l->replaceVertexSrc(vertex_l);
      } else if (edge_l->getVertexDst() == vertex_r) {
        edge_l->replaceVertexDst(vertex_l);
      }
      EdgeSPtr edge_r = edge_1->next(edge_1->getVertexDst());
      if (edge_r->getVertexSrc() == vertex_l) {
        edge_r->replaceVertexSrc(vertex_r);
      } else if (edge_r->getVertexDst() == vertex_l) {
        edge_r->replaceVertexDst(vertex_r);
      }

      edge_1->getFacetR()->removeVertex(vertex_l);
      edge_2->getFacetR()->addVertex(vertex_l);
      edge_1->getFacetL()->removeVertex(vertex_r);
      edge_2->getFacetL()->addVertex(vertex_r);

      edge_22 = edge_1;
      edge_22->replaceVertexDst(edge_2->getVertexDst());
      edge_2->replaceVertexDst(vertex_r);
      edge_22->replaceVertexSrc(vertex_l);

      edge_22->replaceFacetL(edge_2->getFacetL());
      edge_22->replaceFacetR(edge_2->getFacetR());
    }

    if (time_future_bound.has_value()) {
      HdsUtils::setFinalPoint(vertex_l, nullptr);
      HdsUtils::setFinalPoint(vertex_r, nullptr);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    HdsUtils::setNode(vertex_l, node);
    HdsUtils::setNode(vertex_r, node);

    ArcSPtr arc_l = createArc(vertex_l);
    skel_result_->addArc(arc_l);
    ArcSPtr arc_r = createArc(vertex_r);
    skel_result_->addArc(arc_r);

    HdsUtils::setSheet(edge_22, HdsUtils::getSheet(edge_2));

    HdsUtils::getSheet(edge_2)->addArc(arc_l);
    HdsUtils::getSheet(edge_2)->addArc(arc_r);

    for (EdgeWPtr edge_w : vertex_l->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        if (HdsUtils::getSheet(edge) != HdsUtils::getSheet(edge_2)) {
          HdsUtils::getSheet(edge)->addArc(arc_l);
        }
      }
    }

    for (EdgeWPtr edge_w : vertex_r->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        if (HdsUtils::getSheet(edge) != HdsUtils::getSheet(edge_2)) {
          HdsUtils::getSheet(edge)->addArc(arc_r);
        }
      }
    }
#endif

    // Gather relevant elements for local queue updates
    post_op_vertices_ = {{ vertex_l, vertex_r }};
    for (const VertexSPtr& v : post_op_vertices_) {
      for (EdgeWPtr we : v->edges()) {
        post_op_edges_.insert(EdgeSPtr(we.lock()));
      }
    }
    for (FacetWPtr wf : vertex_l->facets()) { post_op_facets_.insert(wf.lock()); }
    for (FacetWPtr wf : vertex_r->facets()) { post_op_facets_.insert(wf.lock()); }
    CGAL_postcondition(post_op_vertices_.size() == 2 && post_op_edges_.size() == 6 && post_op_facets_.size() == 4);

    post_op_vertices_VV_ = {{ vertex_l, vertex_r }};

    // probably could improve this but flip vertex events are rare_, so we are rarely
    // looking for new pierce events after a vertex event so it doesn't matter much
    post_op_vertices_pierce_ = {{ vertex_l, vertex_r }};

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handleSplitMergeEvent(const SplitMergeEventSPtr& event,
                                    const FT& current_time,
                                    const std::optional<FT>& time_future_bound,
                                    const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Split Merge Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    VertexSPtr vertex_1 = event->getVertex1();
    VertexSPtr vertex_2 = event->getVertex2();
    FacetSPtr facet_1 = event->getFacet1();
    FacetSPtr facet_2 = event->getFacet2();

    EdgeSPtr edge_tomerge_1 = EdgeSPtr();
    EdgeSPtr edge_11 = EdgeSPtr();
    EdgeSPtr edge_12 = EdgeSPtr();
    EdgeSPtr edge_tomerge_2 = EdgeSPtr();
    EdgeSPtr edge_21 = EdgeSPtr();
    EdgeSPtr edge_22 = EdgeSPtr();
    for (EdgeWPtr edge_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->getFacetL() == facet_1 && edge->getFacetR() == facet_2) ||
            (edge->getFacetL() == facet_2 && edge->getFacetR() == facet_1)) {
          edge_tomerge_1 = edge;
          continue;
        }
        if (edge->getFacetL() == facet_1 || edge->getFacetR() == facet_1) {
          edge_11 = edge;
        }
        if (edge->getFacetL() == facet_2 || edge->getFacetR() == facet_2) {
          edge_12 = edge;
        }
      }
    }
    for (EdgeWPtr edge_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->getFacetL() == facet_1 && edge->getFacetR() == facet_2) ||
            (edge->getFacetL() == facet_2 && edge->getFacetR() == facet_1)) {
          edge_tomerge_2 = edge;
          continue;
        }
        if (edge->getFacetL() == facet_1 || edge->getFacetR() == facet_1) {
          edge_21 = edge;
        }
        if (edge->getFacetL() == facet_2 || edge->getFacetR() == facet_2) {
          edge_22 = edge;
        }
      }
    }
    FacetSPtr facet_1b = edge_11->getFacetL();
    if (facet_1b == facet_1 || facet_1b == facet_2) {
      facet_1b = edge_11->getFacetR();
    }
    FacetSPtr facet_2b = edge_21->getFacetL();
    if (facet_2b == facet_1 || facet_2b == facet_2) {
      facet_2b = edge_21->getFacetR();
    }
    EdgeSPtr edge_tosplit = EdgeSPtr();
    // edge_tosplit = facet_1b->findEdge(facet_2b);
    EdgeSPtr edge_cur = edge_11->next(facet_1b);
    while (edge_cur != edge_11) {
      if ((edge_cur->getFacetL() == facet_1b && edge_cur->getFacetR() == facet_2b) ||
          (edge_cur->getFacetR() == facet_1b && edge_cur->getFacetL() == facet_2b)) {
        edge_tosplit = edge_cur;
        break;
      }
      edge_cur = edge_cur->next(facet_1b);
    }
    CGAL_SS3_DEBUG_SPTR(edge_tosplit);

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    HdsUtils::getArc(vertex_1)->closeArc(node);
    HdsUtils::getArc(vertex_2)->closeArc(node);

    for (EdgeSPtr e : {edge_11, edge_12, edge_21, edge_22}) {
      HdsUtils::getSheet(e)->addNode(node);
    }
    HdsUtils::getSheet(edge_tosplit)->addNode(node);
    HdsUtils::getSheet(edge_tomerge_1)->addNode(node);

    mergeSheets(edge_tomerge_1, edge_tomerge_2);
#endif

    if (edge_tomerge_1->getVertexSrc() == vertex_1) {
      if (edge_tomerge_2->getVertexSrc() == vertex_2) {
        edge_tomerge_1->replaceVertexSrc(edge_tomerge_2->getVertexDst());
      } else {
        edge_tomerge_1->replaceVertexSrc(edge_tomerge_2->getVertexSrc());
      }
    } else {
      if (edge_tomerge_2->getVertexSrc() == vertex_2) {
        edge_tomerge_1->replaceVertexDst(edge_tomerge_2->getVertexDst());
      } else {
        edge_tomerge_1->replaceVertexDst(edge_tomerge_2->getVertexSrc());
      }
    }
    if (edge_12->getVertexDst() == vertex_1) {
      edge_12->replaceVertexDst(vertex_2);
    } else {
      edge_12->replaceVertexSrc(vertex_2);
    }
    if (edge_21->getVertexDst() == vertex_2) {
      edge_21->replaceVertexDst(vertex_1);
    } else {
      edge_21->replaceVertexSrc(vertex_1);
    }
    facet_1->removeVertex(vertex_2);
    facet_2->removeVertex(vertex_1);
    facet_1b->addVertex(vertex_2);
    facet_2b->addVertex(vertex_1);
    if (edge_tosplit->getFacetL() == facet_1b &&
        edge_tosplit->getFacetR() == facet_2b) {
      edge_tomerge_2->replaceVertexSrc(edge_tosplit->getVertexSrc());
      edge_tomerge_2->replaceVertexDst(vertex_2);
      edge_tosplit->replaceVertexSrc(vertex_1);
    } else if (edge_tosplit->getFacetL() == facet_2b &&
               edge_tosplit->getFacetR() == facet_1b) {
      edge_tomerge_2->replaceVertexDst(edge_tosplit->getVertexDst());
      edge_tomerge_2->replaceVertexSrc(vertex_2);
      edge_tosplit->replaceVertexDst(vertex_1);
    }
    edge_tomerge_2->replaceFacetL(edge_tosplit->getFacetL());
    edge_tomerge_2->replaceFacetR(edge_tosplit->getFacetR());

    if (time_future_bound.has_value()) {
      HdsUtils::setFinalPoint(vertex_1, nullptr);
      HdsUtils::setFinalPoint(vertex_2, nullptr);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    HdsUtils::setNode(vertex_1, node);
    HdsUtils::setNode(vertex_2, node);

    ArcSPtr arc_1 = createArc(vertex_1);
    skel_result_->addArc(arc_1);
    ArcSPtr arc_2 = createArc(vertex_2);
    skel_result_->addArc(arc_2);

    HdsUtils::setSheet(edge_tomerge_2, HdsUtils::getSheet(edge_tosplit));

    for (VertexSPtr vertex : { vertex_1, vertex_2 }) {
      for (EdgeWPtr edge_w : vertex->edges()) {
        if (EdgeSPtr edge = edge_w.lock()) {
          HdsUtils::getSheet(edge)->addArc(HdsUtils::getArc(vertex));
        }
      }
    }
#endif

    // Gather relevant elements for local queue updates
    post_op_vertices_ = {{ vertex_1, vertex_2 }};
    for (const VertexSPtr& v : post_op_vertices_) {
      for (EdgeWPtr we : v->edges()) {
        post_op_edges_.insert(EdgeSPtr(we.lock()));
      }
    }
    post_op_edges_.insert(edge_tomerge_1);
    post_op_edges_.insert(edge_tomerge_2);
    for (FacetWPtr wf : vertex_1->facets()) { post_op_facets_.insert(wf.lock()); }
    for (FacetWPtr wf : vertex_2->facets()) { post_op_facets_.insert(wf.lock()); }
    CGAL_postcondition(post_op_vertices_.size() == 2 && post_op_edges_.size() == 7 && post_op_facets_.size() == 4);

    post_op_vertices_VV_ = {{ vertex_1, vertex_2 }};

    // and all faces are getting smaller so shouldn't there be a need to check disconnections
    // @todo actually assert that all faces involved get smaller
    CGAL_assertion(!HdsUtils::isReflex(vertex_1));
    CGAL_assertion(!HdsUtils::isReflex(vertex_2));
    post_op_vertices_pierce_.clear();

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handleEdgeSplitEvent(const EdgeSplitEventSPtr& event,
                                   const FT& current_time,
                                   const std::optional<FT>& time_future_bound,
                                   const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Edge Split Event  #######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    // this must be done **BEFORE** the shift because after the shift there is
    // an intersection between the two edges, and so the orientation becomes coplanar
    int orientation = KernelWrapper::orientation(HdsUtils::line(event->getEdge1()),
                                                 HdsUtils::line(event->getEdge2()));
    CGAL_assertion(orientation != 0);

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    EdgeSPtr edge_1 = event->getEdge1();
    EdgeSPtr edge_2 = event->getEdge2();

    CGAL_SS3_CORE_TRACE("edge_1 = " << event->getEdge1()->toString());
    CGAL_SS3_CORE_TRACE("edge_2 = " << event->getEdge2()->toString());

    FacetSPtr facet_l1 = edge_1->getFacetL();
    FacetSPtr facet_r1 = edge_1->getFacetR();
    FacetSPtr facet_l2 = edge_2->getFacetL();
    FacetSPtr facet_r2 = edge_2->getFacetR();

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();

    if (facet_1_src == facet_l2 || facet_1_src == facet_r2) {
      HdsUtils::getArc(edge_1->getVertexSrc())->closeArc(node);
    }
    if (facet_1_dst == facet_l2 || facet_1_dst == facet_r2) {
      HdsUtils::getArc(edge_1->getVertexDst())->closeArc(node);
    }

    HdsUtils::getSheet(edge_1)->addNode(node);
    HdsUtils::getSheet(edge_2)->addNode(node);
#endif

    std::array<VertexSPtr, 4> vertices;
    for (unsigned int i = 0; i < 4; ++i) {
      vertices[i] = Vertex::create(point);
      SkelVertexData::create(vertices[i]);
      polyhedron->addVertex(vertices[i]);
    }
    std::array<EdgeSPtr, 4> edges;
    for (unsigned int i = 0; i < 4; ++i) {
      edges[i] = Edge::create(vertices[i], vertices[(i+1)%4]);
      SkelEdgeData::create(edges[i]);
      polyhedron->addEdge(edges[i]);
    }

    if (orientation > 0) {
      edges[0]->setFacetL(edge_2->getFacetR());
      edges[0]->setFacetR(edge_1->getFacetR());
      edges[1]->setFacetL(edge_2->getFacetL());
      edges[1]->setFacetR(edge_1->getFacetR());
      edges[2]->setFacetL(edge_2->getFacetL());
      edges[2]->setFacetR(edge_1->getFacetL());
      edges[3]->setFacetL(edge_2->getFacetR());
      edges[3]->setFacetR(edge_1->getFacetL());
    } else {
      edges[0]->setFacetL(edge_1->getFacetL());
      edges[0]->setFacetR(edge_2->getFacetL());
      edges[1]->setFacetL(edge_1->getFacetL());
      edges[1]->setFacetR(edge_2->getFacetR());
      edges[2]->setFacetL(edge_1->getFacetR());
      edges[2]->setFacetR(edge_2->getFacetR());
      edges[3]->setFacetL(edge_1->getFacetR());
      edges[3]->setFacetR(edge_2->getFacetL());
    }
    EdgeSPtr edge_12 = edge_1->split(vertices[2]);
    EdgeSPtr edge_22 = edge_2->split(vertices[3]);
    SkelEdgeData::create(edge_12);
    SkelEdgeData::create(edge_22);
    edge_1->replaceVertexDst(vertices[0]);
    edge_2->replaceVertexDst(vertices[1]);
    for (unsigned int i = 0; i < 4; ++i) {
      // adds the vertices also
      edges[i]->getFacetL()->addEdge(edges[i]);
      edges[i]->getFacetR()->addEdge(edges[i]);
    }

    if (time_future_bound.has_value()) {
      for (std::size_t i=0; i<4; ++i) {
        HdsUtils::setFinalPoint(vertices[i], nullptr);
      }
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    HdsUtils::setSheet(edge_12, HdsUtils::getSheet(edge_1));
    HdsUtils::setSheet(edge_22, HdsUtils::getSheet(edge_2));

    for (unsigned int i = 0; i < 4; ++i) {
      SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(vertices[i]->getData());
      vertex_data->setNode(node);
      ArcSPtr arc = createArc(vertices[i]);
      skel_result_->addArc(arc);

      for (EdgeWPtr ew : vertices[i]->edges()) {
        if (EdgeSPtr e = ew.lock()) {
          // could abuse the fact that new edges have no data or sheet yet, but this is clearer
          if (std::find(edges.begin(), edges.end(), e) == edges.end()) {
            HdsUtils::getSheet(e)->addArc(arc); // other 2 are added when sheets are created below
          }
        }
      }
    }
    for (unsigned int i = 0; i < 4; ++i) {
      SheetSPtr sheet = createSheet(edges[i]);
      skel_result_->addSheet(sheet);
    }
#endif

    // Gather relevant elements for local queue updates
    post_op_vertices_ = {{ vertices[0], vertices[1], vertices[2], vertices[3] }};
    for (const VertexSPtr& v : post_op_vertices_) {
      for (EdgeWPtr we : v->edges()) {
        post_op_edges_.insert(EdgeSPtr(we.lock()));
      }
    }
    post_op_facets_ = {{ facet_l1, facet_r1, facet_l2, facet_r2 }};
    CGAL_postcondition(post_op_vertices_.size() == 4 && post_op_edges_.size() == 8 && post_op_facets_.size() == 4);

    // @todo do it a bit more elegantly
    post_op_vertices_VV_ = {{ vertices[0], vertices[1], vertices[2], vertices[3] }};
    for (const VertexSPtr& v : { vertices[0], vertices[1], vertices[2], vertices[3] }) {
      for (EdgeWPtr we : v->edges()) {
        EdgeSPtr e = we.lock();
        post_op_vertices_VV_.insert(e->other(v));
      }
    }

    // facets grow so we also need to check vertices that are extremities of edges
    // being subdivided
    for (const EdgeSPtr& poe : post_op_edges_) {
      post_op_vertices_pierce_.insert(poe->getVertexSrc());
      post_op_vertices_pierce_.insert(poe->getVertexDst());
    }
    CGAL_postcondition(post_op_vertices_pierce_.size() == 8);

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handlePierceEvent(const PierceEventSPtr& event,
                                const FT& current_time,
                                const std::optional<FT>& time_future_bound,
                                const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "########  Handle Pierce Event  #########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->getTime();
    Point3SPtr point = event->getPoint();

    shiftToEventTime(polyhedron, current_time, event_time);

    VertexSPtr vertex = event->getVertex();
    FacetSPtr facet = event->getFacet();

    CGAL_SS3_CORE_TRACE("V: " << event->getVertex()->toString());
    CGAL_SS3_CORE_TRACE("F: " << event->getFacet()->toString());

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->setTime(event_time);
    node->setPoint(point);
    skel_result_->addNode(node);

    HdsUtils::getArc(vertex)->closeArc(node);

    for (EdgeWPtr edge_w : vertex->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        HdsUtils::getSheet(edge)->addNode(node);
      }
    }
#endif

    // the 3 new vertices cannot be reflex, but since we grow faces,
    // we need to check the other extremities of the edges
    for (EdgeWPtr ew : vertex->edges()) {
      EdgeSPtr edge = ew.lock();
      if (edge->getVertexSrc() != vertex) {
        post_op_vertices_pierce_.insert(edge->getVertexSrc());
      } else {
        post_op_vertices_pierce_.insert(edge->getVertexDst());
      }
    }

    CGAL_postcondition(post_op_vertices_pierce_.size() == 3);


    std::array<FacetSPtr, 3> facets;
    std::array<EdgeSPtr, 3> edges;
    EdgeSPtr edge = vertex->firstEdge();
    for (unsigned int i = 0; i < 3; ++i) {
      edges[i] = edge;
      if (edge->getVertexSrc() == vertex) {
        facets[i] = edge->getFacetL();
      } else if (edge->getVertexDst() == vertex) {
        facets[i] = edge->getFacetR();
      }
      edge = edge->next(vertex);
    }

    std::array<VertexSPtr, 3> vertices;
    for (unsigned int i = 0; i < 3; ++i) {
      vertices[i] = Vertex::create(point);
      SkelVertexData::create(vertices[i]);
      facet->addVertex(vertices[i]);
      polyhedron->addVertex(vertices[i]);
    }
    for (unsigned int i = 0; i < 3; ++i) {
      EdgeSPtr edge = edges[i];
      if (edge->getVertexSrc() == vertex) {
        edge->replaceVertexSrc(vertices[i]);
      } else if (edge->getVertexDst() == vertex) {
        edge->replaceVertexDst(vertices[i]);
      }
      facets[i]->removeVertex(vertex);
      facets[i]->addVertex(vertices[i]);
      facets[(i+2)%3]->addVertex(vertices[i]);
    }
    vertex->facets().clear();
    vertex->edges().clear();
    polyhedron->removeVertex(vertex);
    std::array<EdgeSPtr, 3> new_edges;
    for (unsigned int i = 0; i < 3; ++i) {
      new_edges[i] = Edge::create(vertices[i], vertices[(i+1)%3]);
      SkelEdgeData::create(new_edges[i]);
      new_edges[i]->setFacetL(facet);
      new_edges[i]->setFacetR(facets[i]);
      facet->addEdge(new_edges[i]);
      facets[i]->addEdge(new_edges[i]);
      polyhedron->addEdge(new_edges[i]);
    }

    if (time_future_bound.has_value()) {
      for (std::size_t i=0; i<3; ++i) {
        HdsUtils::setFinalPoint(vertices[i], nullptr);
      }
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    for (unsigned int i = 0; i < 3; ++i) {
      HdsUtils::setNode(vertices[i], node);
      ArcSPtr arc = createArc(vertices[i]);
      skel_result_->addArc(arc);
      HdsUtils::getSheet(edges[i])->addArc(arc); // other 2 sheets are the new edges', created below
    }
    for (unsigned int i = 0; i < 3; ++i) {
      SheetSPtr sheet = createSheet(new_edges[i]);
      skel_result_->addSheet(sheet);
    }
#endif

    // Gather relevant elements for local queue updates
    post_op_vertices_ = {{ vertices[0], vertices[1], vertices[2] }};
    for (const VertexSPtr& v : post_op_vertices_) {
      for (EdgeWPtr we : v->edges()) {
        if (EdgeSPtr e = we.lock()) {
          post_op_edges_.insert(e);
        }
      }
    }
    for (unsigned int i = 0; i < 3; ++i) {
      post_op_facets_.insert(new_edges[i]->getFacetL());
      post_op_facets_.insert(new_edges[i]->getFacetR());
    }
    CGAL_postcondition(post_op_vertices_.size() == 3 && post_op_edges_.size() == 6 && post_op_facets_.size() == 4);

    post_op_vertices_VV_ = {{ vertices[0], vertices[1], vertices[2] }};
    for (const VertexSPtr& v : { vertices[0], vertices[1], vertices[2] }) {
      for (EdgeWPtr we : v->edges()) {
        if (EdgeSPtr e = we.lock()) {
          post_op_vertices_VV_.insert(e->other(v));
        }
      }
    }

    // the edge with 'facet offset' is convex so these vertices are not interesting
    CGAL_postcondition(!HdsUtils::isReflex(vertices[0]));
    CGAL_postcondition(!HdsUtils::isReflex(vertices[1]));
    CGAL_postcondition(!HdsUtils::isReflex(vertices[2]));

    addEvent(event);

    return EventStatus::EVENT_HANDLED;
  }

  EventStatus handleEvent(const AbstractEventSPtr& event,
                          const FT& current_time,
                          const std::optional<FT>& time_future_bound,
                          const PolyhedronSPtr& polyhedron)
  {
    EventStatus result = EventStatus::NON_EVENT;

    // The point of these checks being here and not in the collect is two-fold:
    // - Some of them perform geometric or combinatorial checks that would be difficult
    //   to perform at collect time *when using local updates*. For example, to know
    //   if we do have a pierce event, we need to check if the vertex and the facet do not share
    //   an edge. But combinatorial information can evolve with events between the moment
    //   we first try to collect (locally) the event, and the moment where we pop it. So we
    // - We often waste a lot of time performing checks at collect time. Taking the same pierce
    //   event type, we need to check if the movement of the vertex would actually touch the
    //   polygon facet. This is an expensive check and the facet could change quite a bit
    //   between collect time and event time. Thus, it's cheaper to delay the check,
    //   even if that means computing the event time.
    if (!isActualEvent(event, current_time, time_future_bound)) {
      return result;
    }

    if (event->getType() == AbstractEvent::SAVE_EVENT) {
      result = handleSaveEvent(std::dynamic_pointer_cast<SaveEvent>(event),
                               current_time, polyhedron);
    } else if (event->getType() == AbstractEvent::CONST_TIME_EVENT) {
      result = handleConstTimeEvent(std::dynamic_pointer_cast<ConstTimeEvent>(event),
                                    current_time, polyhedron);
    } else if (event->getType() == AbstractEvent::VANISH_EVENT) {
      result = handleVanishEvent(std::dynamic_pointer_cast<VanishEvent>(event),
                                 current_time, time_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::VERTEX_EVENT) {
      result = handleVertexEvent(std::dynamic_pointer_cast<VertexEvent>(event),
                                 current_time, time_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::FLIP_VERTEX_EVENT) {
      result = handleFlipVertexEvent(std::dynamic_pointer_cast<FlipVertexEvent>(event),
                                     current_time, time_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::SURFACE_EVENT) {
      result = handleSurfaceEvent(std::dynamic_pointer_cast<SurfaceEvent>(event),
                                  current_time, time_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
      result = handlePolyhedronSplitEvent(std::dynamic_pointer_cast<PolyhedronSplitEvent>(event),
                                          current_time, time_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::SPLIT_MERGE_EVENT) {
      result = handleSplitMergeEvent(std::dynamic_pointer_cast<SplitMergeEvent>(event),
                                     current_time, time_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::EDGE_SPLIT_EVENT) {
      result = handleEdgeSplitEvent(std::dynamic_pointer_cast<EdgeSplitEvent>(event),
                                    current_time, time_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::PIERCE_EVENT) {
      result = handlePierceEvent(std::dynamic_pointer_cast<PierceEvent>(event),
                                 current_time, time_future_bound, polyhedron);
    } else {
      CGAL_SS3_CORE_TRACE("Error: Cannot handle event of type " << event->getType());
      CGAL_assertion(false);
      result = EventStatus::EVENT_NOT_HANDLED;
    }

    CGAL_postcondition_code(for (const VertexSPtr& v : polyhedron->vertices()))
    CGAL_postcondition(v->getID() != -1);
    CGAL_postcondition_code(for (const EdgeSPtr& e : polyhedron->edges()))
    CGAL_postcondition(e->getID() != -1);
    CGAL_postcondition_code(for (const FacetSPtr& f : polyhedron->facets()))
    CGAL_postcondition(f->getID() != -1);

#ifndef CGAL_SS3_NO_SKELETON_DS
    CGAL_postcondition_code(for (const VertexSPtr& v : polyhedron->vertices()) {)
    CGAL_postcondition(HdsUtils::getNode(v) != NodeSPtr());
    CGAL_postcondition(HdsUtils::getArc(v) != ArcSPtr());
    CGAL_postcondition_code(})
    CGAL_postcondition_code(for (const EdgeSPtr& e : polyhedron->edges()))
    CGAL_postcondition(HdsUtils::getSheet(e) != SheetSPtr());
#endif

    CGAL_SS3_CORE_TRACE_V(4, "-- Finished handling Event --");
    return result;
  }

  StraightSkeletonSPtr getResult() const
  {
#ifndef CGAL_SS3_NO_SKELETON_DS
    CGAL_SS3_CORE_TRACE("Warning: no skeleton to return as it was not built");
#endif
    return this->skel_result_;
  }

  std::list<AbstractEventSPtr>& events()
  {
    return this->events_;
  }

  int countEvents(int type) const
  {
    int result = 0;
    typename std::list<AbstractEventSPtr>::const_iterator it_e = events_.begin();
    while (it_e != events_.end()) {
      AbstractEventSPtr event = *it_e++;
      if (event->getType() == type) {
        result += 1;
      }
    }
    return result;
  }

  std::string eventSummary() const
  {
    std::stringstream sstr;
    sstr << "Events: " << events_.size() << std::endl;
    sstr << "    ConstTimeEvents:       " << countEvents(AbstractEvent::CONST_TIME_EVENT) << std::endl;
    sstr << "    SaveEvents:            " << countEvents(AbstractEvent::SAVE_EVENT) << std::endl;
    sstr << "  VanishEvents:" << std::endl;
    sstr << "    Generic VanishEvents:  " << countEvents(AbstractEvent::VANISH_EVENT) << std::endl;
    sstr << "    EdgeEvents:            " << countEvents(AbstractEvent::EDGE_EVENT) << std::endl;
    sstr << "    EdgeMergeEvents:       " << countEvents(AbstractEvent::EDGE_MERGE_EVENT) << std::endl;
    sstr << "    TriangleEvents:        " << countEvents(AbstractEvent::TRIANGLE_EVENT) << std::endl;
    sstr << "    DblEdgeMergeEvents:    " << countEvents(AbstractEvent::DBL_EDGE_MERGE_EVENT) << std::endl;
    sstr << "    DblTriangleEvents:     " << countEvents(AbstractEvent::DBL_TRIANGLE_EVENT) << std::endl;
    sstr << "    TetrahedronEvents:     " << countEvents(AbstractEvent::TETRAHEDRON_EVENT) << std::endl;
    sstr << "  ContactEvents:" << std::endl;
    sstr << "    VertexEvents:          " << countEvents(AbstractEvent::VERTEX_EVENT) << std::endl;
    sstr << "    FlipVertexEvents:      " << countEvents(AbstractEvent::FLIP_VERTEX_EVENT) << std::endl;
    sstr << "    SurfaceEvents:         " << countEvents(AbstractEvent::SURFACE_EVENT) << std::endl;
    sstr << "    PolyhedronSplitEvents: " << countEvents(AbstractEvent::POLYHEDRON_SPLIT_EVENT) << std::endl;
    sstr << "    SplitMergeEvents:      " << countEvents(AbstractEvent::SPLIT_MERGE_EVENT) << std::endl;
    sstr << "    EdgeSplitEvents:       " << countEvents(AbstractEvent::EDGE_SPLIT_EVENT) << std::endl;
    sstr << "    PierceEvents:          " << countEvents(AbstractEvent::PIERCE_EVENT) << std::endl;
    sstr << ")" << std::endl;
    return sstr.str();

  }

private:
  PolyhedronSPtr polyhedron_;
  AbstractVertexSplitterSPtr vertex_splitter_;
  int edge_event_;

  Base_mesh_offset_visitor<Traits>* visitor_ = nullptr;

  std::vector<FT> save_times_;
  std::filesystem::path save_path_;

  std::list<AbstractEventSPtr> events_;
  StraightSkeletonSPtr skel_result_;

  int step_id_;

  std::set<VertexSPtr> post_op_vertices_;
  std::set<EdgeSPtr> post_op_edges_;
  std::set<FacetSPtr> post_op_facets_;
  std::set<VertexSPtr> post_op_vertices_VV_;
  std::set<VertexSPtr> post_op_vertices_pierce_;
  std::set<EdgeSPtr> post_op_edges_edgesplit_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_STRAIGHT_SKELETON_BUILDER_H */
