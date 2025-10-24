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

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_STRAIGHT_SKELETON_BUILDER_3_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_STRAIGHT_SKELETON_BUILDER_3_H

// @fixme yesterday:
// - perturbation depth-limit without hardcoded digit size values
// - queue correctness assertion failures (362913? doublebox_n? verworrtakelt_n?)
// - 57419 (perturbation failure)
// - 168077 Expr: !Self_intersection::has_self_intersecting_surface(polyhedron_)
// - 80084 Expr: !Self_intersection::has_self_intersecting_surface(polyhedron_)
// - 67817 Expr: !Self_intersection::has_self_intersecting_surface(polyhedron_)
// - 100638 polyhedron->is_consistent()
// - 55583 Expr: !Self_intersection::has_self_intersecting_surface(polyhedron_)
// - 153956 Expr: !Self_intersection::has_self_intersecting_surface(polyhedron_)
// - 500087

// @fixme:
// - skeleton arc initial directions are wrong
// - handle save times at event time: merge vertices with equal position in output

// @fixme later:
// - more combinatorial checks should happen at pop time (?)
// - double check self-intersection detection on boundaries
// - re-using items in event handlers goes against detecting obsolete events with expired pointers

// @fixme latest:
// - Fix simultaneous events still happening sometimes (likely the same event multiple times
//   since we don't check the queue before pushing)
// - EPECK -> EPICK embedding could create self-intersections

// @speed
// - If an edge is growing, there is no point computing its vanish event
// - Avoid recomputing the vanish time if the edge has not been modified
// - Don't actually shift at all intermediate steps: we can do everything with base planes,
//   including at places where we use shifted point positions
// - Use CGAL_SS3_NEWER_VV_VERTEX_2_DETECTION and CGAL_SS3_VV_VERTEX_2_WALK_FACES_FOR_DETECTION
//   in vertex-vertex events (check the graveyard)
// - For contact events: exit early if the 4 planes are clearly not intersecting (diametral spheres
//   around the edges of size [something]?)
// - Avoid recomputing the denominator when the point is required; store event time as a quotient
// - In crash_time() for pierce events, we could tighten the future bound using the *farthest* vanish
//   event of any edge of the facet (probably not worth it, though)

// @todo: cleaning
// - clean up all the code related to local queue updates (horrible variable names, duplicates, etc.)
// - improve the visitor & give access to the skeleton data structure

// @todo
// - add tests; doc figures
// - Do not accept non manifold inputs, non-triangulated inputs (clarify doc)
// - Do not triangulate outputs, use the code from remesh_planar_faces() for not simply connected faces

// @todo later:
// - if checking perturbation fails because of self-intersections, use a smaller epsilon
// - perform facet merging using CGAL's region growing and remesh_planar_faces() (?)
// - re-enable the option to translate and scale (?)
// - tolerate non-triangulated inputs (?)
// - check for overly shared objects, redundant function calls (plane normalization, for example)

// @todo latest:
// - get rid of all the shared ptr stuff, we only need to zombie the elements and use IDs
// - use traits' functors
// - get rid of the exact construction requirement? At least if we do not have to split high-degree
//   vertices, it should be possible, but that would require writing filtered predicates like SLS2's.
// - lighter & faster polyhedron data structures
// - write a sanitization algorithm without perturbation, something akin to: apply_rand_plane_tilts_V3(p, eps=0),
//   which only ensures that points are on supporting planes, but does NOT perturb the planes.
// - splitting high degree vertices in reasonable time

// ----

/*
  As to not waste energy building the skeleton if we do not care about it.
*/
// #define CGAL_SS3_NO_SKELETON_DS

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
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/Straight_skeleton_3.h>
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
  using Abstract_event_sptr = std::shared_ptr<Abstract_event<Traits> >;

  virtual bool go_further(int, PolyhedronSPtr, FT) = 0;
  virtual void before_event(PolyhedronSPtr, FT, Abstract_event_sptr) = 0;
  virtual void on_save_event(PolyhedronSPtr, FT) = 0;
  virtual void after_event(PolyhedronSPtr, FT) = 0;
};

template <typename Traits>
struct Default_mesh_offset_visitor
  : public Base_mesh_offset_visitor<Traits>
{
  using FT = typename Traits::FT;
  using PolyhedronSPtr = std::shared_ptr<HDS::Polyhedron<Traits> >;
  using Abstract_event_sptr = std::shared_ptr<Abstract_event<Traits> >;

  bool go_further(int, PolyhedronSPtr, FT) override { return true; }
  void before_event(PolyhedronSPtr, FT, Abstract_event_sptr) override { }
  void on_save_event(PolyhedronSPtr, FT) override { }
  void after_event(PolyhedronSPtr, FT) override { }
};

template <typename Traits>
class Straight_skeleton_builder_3
{
  using Straight_skeleton_builder_sptr = std::shared_ptr<Straight_skeleton_builder_3<Traits> >;

private:
  // Geometry
  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;
  using Segment_3 = typename Traits::Segment_3;
  using Vector_3 = typename Traits::Vector_3;
  using Line_3 = typename Traits::Line_3;
  using Plane_3 = typename Traits::Plane_3;

private:
  // Polyhedron Data Structure
  using Polyhedron = HDS::Polyhedron<Traits>;
  using PolyhedronWPtr = typename Polyhedron::PolyhedronWPtr;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Vertex = typename Polyhedron::Vertex;
  using VertexWPtr = typename Polyhedron::VertexWPtr;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using Vertex_data = typename Polyhedron::Vertex_data;
  using VertexDataSPtr = typename Polyhedron::VertexDataSPtr;
  using Skeleton_vertex_data = typename Polyhedron::Skeleton_vertex_data;
  using SkelVertexDataSPtr = typename Polyhedron::SkelVertexDataSPtr;
  using Edge = typename Polyhedron::Edge;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using Edge_data = typename Polyhedron::Edge_data;
  using EdgeDataSPtr = typename Polyhedron::EdgeDataSPtr;
  using Skeleton_edge_data = typename Polyhedron::Skeleton_edge_data;
  using SkelEdgeDataSPtr = typename Polyhedron::SkelEdgeDataSPtr;
  using Facet = typename Polyhedron::Facet;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;
  using Facet_data = typename Polyhedron::Facet_data;
  using FacetDataSPtr = typename Polyhedron::FacetDataSPtr;
  using Skeleton_facet_data = typename Polyhedron::Skeleton_facet_data;
  using SkelFacetDataSPtr = typename Polyhedron::SkelFacetDataSPtr;

private:
  // Straight Skeleton Data Structure
  using Straight_skeleton_3 = CGAL::Straight_skeleton_3<Traits>;
  using StraightSkeletonWPtr = typename Straight_skeleton_3::StraightSkeletonWPtr;
  using StraightSkeletonSPtr = typename Straight_skeleton_3::StraightSkeletonSPtr;

  using Node = typename Straight_skeleton_3::Node;
  using NodeSPtr = typename Straight_skeleton_3::NodeSPtr;
  using Arc = typename Straight_skeleton_3::Arc;
  using ArcWPtr = typename Straight_skeleton_3::ArcWPtr;
  using ArcSPtr = typename Straight_skeleton_3::ArcSPtr;
  using Sheet = typename Straight_skeleton_3::Sheet;
  using SheetWPtr = typename Straight_skeleton_3::SheetWPtr;
  using SheetSPtr = typename Straight_skeleton_3::SheetSPtr;

private:
  // Vertex Splitters
  using Abstract_vertex_splitter = algorithm::Abstract_vertex_splitter<Traits>;
  using Abstract_vertex_splitter_sptr = std::shared_ptr<Abstract_vertex_splitter>;
  using Combi_vertex_splitter = algorithm::Combi_vertex_splitter<Traits>;
  using Convex_vertex_splitter = algorithm::Convex_vertex_splitter<Traits>;

private:
  // Events
  using Abstract_event = algorithm::Abstract_event<Traits>;
  using Abstract_event_sptr = std::shared_ptr<Abstract_event>;

  using Const_time_event = algorithm::Const_time_event<Traits>;
  using Const_time_event_sptr = std::shared_ptr<Const_time_event>;
  using Save_event = algorithm::Save_event<Traits>;
  using Save_event_sptr = std::shared_ptr<Save_event>;

  using Vanish_event = algorithm::Vanish_event<Traits>;
  using Vanish_event_sptr = std::shared_ptr<Vanish_event>;
  using Edge_event = algorithm::Edge_event<Traits>;
  using Edge_event_sptr = std::shared_ptr<Edge_event>;
  using Edge_merge_event = algorithm::Edge_merge_event<Traits>;
  using Edge_merge_event_sptr = std::shared_ptr<Edge_merge_event>;
  using Triangle_event = algorithm::Triangle_event<Traits>;
  using Triangle_event_sptr = std::shared_ptr<Triangle_event>;
  using Dbl_edge_merge_event = algorithm::Dbl_edge_merge_event<Traits>;
  using Dbl_edge_merge_event_sptr = std::shared_ptr<Dbl_edge_merge_event>;
  using Dbl_triangle_event = algorithm::Dbl_triangle_event<Traits>;
  using Dbl_triangle_event_sptr = std::shared_ptr<Dbl_triangle_event>;
  using Tetrahedron_event = algorithm::Tetrahedron_event<Traits>;
  using Tetrahedron_event_sptr = std::shared_ptr<Tetrahedron_event>;

  using Vertex_event = algorithm::Vertex_event<Traits>;
  using Vertex_event_sptr = std::shared_ptr<Vertex_event>;
  using Flip_vertex_event = algorithm::Flip_vertex_event<Traits>;
  using Flip_vertex_event_sptr = std::shared_ptr<Flip_vertex_event>;
  using Surface_event = algorithm::Surface_event<Traits>;
  using Surface_event_sptr = std::shared_ptr<Surface_event>;
  using Polyhedron_split_event = algorithm::Polyhedron_split_event<Traits>;
  using Polyhedron_split_event_sptr = std::shared_ptr<Polyhedron_split_event>;
  using Split_merge_event = algorithm::Split_merge_event<Traits>;
  using Split_merge_event_sptr = std::shared_ptr<Split_merge_event>;
  using Edge_split_event = algorithm::Edge_split_event<Traits>;
  using Edge_split_event_sptr = std::shared_ptr<Edge_split_event>;
  using Pierce_event = algorithm::Pierce_event<Traits>;
  using Pierce_event_sptr = std::shared_ptr<Pierce_event>;

  enum class Event_status {
    NON_EVENT = 0,
    EVENT_HANDLED,
    EVENT_NOT_HANDLED
  };

private:
  using Kernel_wrapper = kernel::Kernel_wrapper<Traits>;
  using Geom_utils = algorithm::Geom_utils<Traits>;
  using Hds_utils = algorithm::Hds_utils<Traits>;
  using Transformation = algorithm::Polyhedron_transformation<Traits>;
  using Perturbation = algorithm::Polyhedron_perturbation<Traits>;
  using Self_intersection = algorithm::Self_intersection<Traits>;

private:
  using PQ = std::priority_queue<Abstract_event_sptr,
                                 std::vector<Abstract_event_sptr>,
                                 Abstract_event_compare<Traits> >;

public:
  Straight_skeleton_builder_3(PolyhedronSPtr polyhedron)
    : polyhedron_(polyhedron),
      save_path_(std::filesystem::current_path()),
      skeleton_(Straight_skeleton_3::create())
  {
    init_vertex_splitter();
    init_edge_event();
  }

  Straight_skeleton_builder_3(PolyhedronSPtr polyhedron,
                              const std::vector<FT>& save_times,
                              const std::filesystem::path& save_path)
    : polyhedron_(polyhedron),
      save_times_(save_times), // intentional copy
      save_path_(save_path),
      skeleton_(Straight_skeleton_3::create())
  {
    std::sort(save_times_.begin(), save_times_.end(),
              [](const FT& a, const FT& b) { return CGAL::abs(a) < CGAL::abs(b); });

    init_vertex_splitter();
    init_edge_event();
  }

  ~Straight_skeleton_builder_3()
  {
    polyhedron_.reset();
    vertex_splitter_.reset();
    events_.clear();
    skeleton_.reset();
  }

  static Straight_skeleton_builder_sptr create(PolyhedronSPtr polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    return std::make_shared<Straight_skeleton_builder_3>(polyhedron);
  }

  static Straight_skeleton_builder_sptr create(PolyhedronSPtr polyhedron,
                                               const std::vector<FT>& save_times)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    return std::make_shared<Straight_skeleton_builder_3>(polyhedron, save_times);
  }

  static Straight_skeleton_builder_sptr create(PolyhedronSPtr polyhedron,
                                               const std::vector<FT>& save_times,
                                               const std::filesystem::path& save_path)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    return std::make_shared<Straight_skeleton_builder_3>(polyhedron, save_times, save_path);
  }

  void setVisitor(Base_mesh_offset_visitor<Traits>* visitor)
  {
    visitor_ = visitor;
  }

  void init_vertex_splitter()
  {
    ConfigurationSPtr config = Configuration::get_instance();
    std::string s_vertex_splitter;
    if (config->is_loaded()) {
      s_vertex_splitter = config->get_string("Algorithm", "vertex_splitter");
      if (s_vertex_splitter.compare("Combi_vertex_splitter") == 0) {
        vertex_splitter_ = Combi_vertex_splitter::create();
      } else if (s_vertex_splitter.compare("Convex_vertex_splitter") == 0) {
        vertex_splitter_ = Convex_vertex_splitter::create();
      } else {
        CGAL_SS3_SPLITTER_TRACE("Warning: option '" << s_vertex_splitter << "' not found.");
        CGAL_SS3_SPLITTER_TRACE("Using 'Combi_vertex_splitter'.");
        vertex_splitter_ = Combi_vertex_splitter::create();
      }
    } else {
      vertex_splitter_ = Combi_vertex_splitter::create();
    }
  }

  void init_edge_event()
  {
    ConfigurationSPtr config = Configuration::get_instance();
    std::string s_edge_event;
    if (config->is_loaded()) {
      s_edge_event = Configuration::get_instance()->get_string("Algorithm", "edge_event");
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

    CGAL_SS3_DEBUG_SPTR(polyhedron_);
    CGAL_assertion(polyhedron_->is_consistent());

#if 1//def CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("results/input.obj", polyhedron_, false /*do not triangulate*/);
#endif

    CGAL_SS3_CORE_TRACE_V(1, polyhedron_->vertices().size() << " NV " << polyhedron_->facets().size() << " NF");

    CGAL_assertion(Perturbation::do_all_plane_pairs_intersect(polyhedron_));
    CGAL_assertion(Perturbation::do_all_plane_triplets_intersect(polyhedron_));
    CGAL_assertion(!Self_intersection::has_self_intersecting_surface(polyhedron_));

    skeleton_->set_polyhedron(polyhedron_); // skeleton's polyhedron is fixed
    PolyhedronSPtr polyhedron = polyhedron_->clone();

    // store base plane coefficients
    cache_base_planes(polyhedron);

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
      const auto pl = facet->get_plane();
      const auto normal = pl.orthogonal_vector();
      // DEBUG_PRINT("SP X " << CGAL::scalar_product(*normal, Vector_3(1,0,0)));
      // DEBUG_PRINT("SP Y " << CGAL::scalar_product(*normal, Vector_3(0,1,0)));
      // DEBUG_PRINT("SP Z " << CGAL::scalar_product(*normal, Vector_3(0,0,1)));
      if (CGAL::abs(CGAL::abs(CGAL::scalar_product(normal, Vector_3(1,0,0))) - 1) < 1e-3)
        speed = x_speed;
      if (CGAL::abs(CGAL::abs(CGAL::scalar_product(normal, Vector_3(0,1,0))) - 1) < 1e-3)
        speed = y_speed;
      if (CGAL::abs(CGAL::abs(CGAL::scalar_product(normal, Vector_3(0,0,1))) - 1) < 1e-3)
        speed = z_speed;

      Hds_utils::set_speed(facet, speed);
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
      ConfigurationSPtr config = Configuration::get_instance();
      if (config->is_loaded()) {
        if ((config->contains("Algorithm", "stop_after_last_save_event") &&
              config->get_Boolean("Algorithm", "stop_after_last_save_event"))) {
          time_future_bound = save_times_.back();
        }
      }
    }

    step_id_ = -1;
    FT current_time = 0;
    FT upcoming_event_time;

    CGAL_assertion_code(const bool is_emptiness_expected = save_times_.empty();)

    PQ queue;
    collect_events(polyhedron, current_time, time_future_bound, queue);

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
      CGAL_assertion(facet->get_plane().a() == Hds_utils::get_base_plane(facet).a());
      CGAL_assertion(facet->get_plane().b() == Hds_utils::get_base_plane(facet).b());
      CGAL_assertion(facet->get_plane().c() == Hds_utils::get_base_plane(facet).c());
      CGAL_assertion_code(FT speed = Hds_utils::get_speed(facet);)
      CGAL_assertion(facet->get_plane().d() == Hds_utils::get_base_plane(facet).d() - speed * current_time);
      CGAL_assertion_code(})

      Abstract_event_sptr event = nextEvent(queue, polyhedron, current_time);
      if (!event) {
        CGAL_SS3_CORE_TRACE_V(2, "No more events to treat");
        break;
      }

      CGAL_SS3_CORE_TRACE_V(2, "popped E" << event->get_ID() << " Type [" << event->getType() << "]");

      static int event_id = -1;
      CGAL_SS3_CORE_TRACE_V(2, "--> Accepted event #" << ++event_id << " " << event->to_string() << " --");

      upcoming_event_time = event->time();

      // the next event should be at a time that is further away than the current one
      CGAL_assertion(upcoming_event_time < current_time);

#ifdef CGAL_SS3_RUN_TIMERS
      CGAL_SS3_CORE_TRACE_V(2, "current elapsed time: " << timer.time());
#endif

      if (visitor_) {
        visitor_->before_event(polyhedron, current_time, event);
      }

      // Event treatment
      Event_status es = handle_event(event, current_time, time_future_bound, polyhedron);
      CGAL_assertion(es != Event_status::EVENT_NOT_HANDLED);
      if (es == Event_status::NON_EVENT) {
        continue;
      }

      current_time = upcoming_event_time;

#ifdef CGAL_SS3_DUMP_FILES
      IO::OBJFile::save("results/event_" + std::to_string(event_id) + ".obj", polyhedron, false /*do_triangulate*/);
      IO::OBJFile::save("results/event_" + std::to_string(event_id) + "_triangulated.obj", polyhedron);
#endif

      if (visitor_ && event->getType() == Abstract_event::SAVE_EVENT) {
        visitor_->on_save_event(polyhedron, current_time);
      }

      CGAL_SS3_CORE_TRACE_V(2, skeleton_->to_string());

#ifdef CGAL_SS3_DUMP_FILES
      // Dump skeleton nodes in an .xyz file
      std::ofstream nodes_out("final_nodes.xyz");
      nodes_out.precision(17);
      for (NodeSPtr node : skeleton_->nodes()) {
        nodes_out << node->point() << "\n";
      }
      nodes_out.close();

      // Dump skeleton arcs as CGAL polylines
      std::ofstream arcs_out("final_arcs.polylines.txt");
      arcs_out.precision(17);
      for (ArcSPtr arc : skeleton_->arcs()) {
        arcs_out << "2 ";
        arcs_out << arc->get_node_src()->point() << " ";
        if (arc->has_node_dst()) {
          arcs_out << arc->get_node_dst()->point() << "\n";
        } else {
          const Point_3& src_pt = arc->get_node_src()->point();
          const Vector_3& dir = arc->get_direction();
          constexpr double ray_length = 0.1; // @todo relative value
          Point_3 ray_pt = src_pt + ray_length * dir;
          arcs_out << ray_pt << "\n";
        }
      }
      arcs_out.close();
#endif

      CGAL_postcondition(polyhedron->is_consistent());
      CGAL_postcondition(skeleton_->is_consistent());

      if (visitor_) {
        visitor_->after_event(polyhedron, current_time);
      }

      // If we are only interested in specific times, there is no point going further
      if (event->getType() == Abstract_event::SAVE_EVENT && save_times_.empty()) {
        ConfigurationSPtr config = Configuration::get_instance();
        if (config->is_loaded() &&
            config->contains("Algorithm", "stop_after_last_save_event") &&
            config->get_Boolean("Algorithm", "stop_after_last_save_event")) {
          break;
        }
      }

      // Update the event priority queue
      collect_local_events(polyhedron, current_time, time_future_bound, queue);

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

    CGAL_SS3_CORE_TRACE_V(2, events_summary());


    CGAL_assertion(skeleton_->is_consistent(false /*is_partial*/));
#ifdef CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("final_skeleton.obj", skeleton_, true /*convert_to_double*/);
#endif

    return true;
  }

  /**
    * Creates a new node for the vertex data.
    * Used by init(...) only.
    */
  static NodeSPtr create_node(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    NodeSPtr result = Node::create();
    result->set_time(0);
    result->set_point(vertex->point());
    Hds_utils::set_node(vertex, result);
    return result;
  }

  /**
    * Creates a new arc for the vertex data.
    * The node of the vertex data has to be set before.
    */
  static ArcSPtr create_arc(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->degree() == 3);
    ArcSPtr result = ArcSPtr();
    CGAL_precondition(vertex->has_data());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());

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

    const Plane_3& plane_1 = facets[0]->get_plane();
    const Plane_3& plane_2 = facets[1]->get_plane();
    const Plane_3& plane_3 = facets[2]->get_plane();
    const FT& speed_1 = Hds_utils::get_speed(facets[0]);
    const FT& speed_2 = Hds_utils::get_speed(facets[1]);
    const FT& speed_3 = Hds_utils::get_speed(facets[2]);
    const Vector_3 n_1 = plane_1.orthogonal_vector();
    const Vector_3 n_2 = plane_2.orthogonal_vector();
    const Vector_3 n_3 = plane_3.orthogonal_vector();
    const Vector_3 direction = Vector_3 { speed_1*n_1 + speed_2*n_2 + speed_3*n_3 };

    result = Arc::create(data->get_node(), direction);
    data->set_arc(result);
    CGAL_SS3_DEBUG_SPTR(result);

    return result;
  }

  /**
    * Creates a new sheet for the edge data.
    * The arcs of the vertices have to be set before.
    */
  static SheetSPtr create_sheet(const EdgeSPtr& edge)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    CGAL_precondition(edge->get_vertex_src() != edge->get_vertex_dst());

    SheetSPtr result = SheetSPtr();

    FacetSPtr facet_l = edge->get_facet_L();
    FacetSPtr facet_r = edge->get_facet_R();
    CGAL_SS3_DEBUG_SPTR(facet_l);
    CGAL_SS3_DEBUG_SPTR(facet_r);

    const Plane_3& plane_l = facet_l->get_plane();
    const Plane_3& plane_r = facet_r->get_plane();
    FacetSPtr facet_b = Hds_utils::get_facet_origin(facet_l);
    FT speed_l = Hds_utils::get_speed(facet_l);
    FacetSPtr facet_f = Hds_utils::get_facet_origin(facet_r);
    FT speed_r = Hds_utils::get_speed(facet_r);

    Plane_3 plane_sheet;
    if (speed_l == speed_r) {
      if (Hds_utils::is_reflex(edge)) {
        plane_sheet = Kernel_wrapper::bisector(plane_l.opposite(), plane_r);
      } else {
        plane_sheet = Kernel_wrapper::bisector(plane_l, plane_r.opposite());
      }
    } else {
      std::optional<Line_3> line = Kernel_wrapper::intersection(plane_l, plane_r);
      Plane_3 offset_l = Geom_utils::offset_plane(plane_l, -speed_l);
      Plane_3 offset_r = Geom_utils::offset_plane(plane_r, -speed_r);
      std::optional<Line_3> line_offset = Kernel_wrapper::intersection(offset_l, offset_r);
      Point_3 point_1 = line->point();
      Vector_3 direction = line->to_vector();
      Point_3 point_2 = point_1 + direction;
      Point_3 point_3 = line_offset->point();
      if (Hds_utils::is_reflex(edge)) {
        plane_sheet = Plane_3 { point_3, point_2, point_1 };
      } else {
        plane_sheet = Plane_3 { point_1, point_2, point_3 };
      }
    }

    result = Sheet::create();
    result->set_plane(plane_sheet);
    result->set_facet_B(facet_b);
    result->set_facet_F(facet_f);

    set_sheet(edge, result);

    NodeSPtr node_src = Hds_utils::get_node(edge->get_vertex_src());
    NodeSPtr node_dst = Hds_utils::get_node(edge->get_vertex_dst());
    result->add_node(node_src);
    if (node_src != node_dst) {
      result->add_node(node_dst);
    }
    result->add_arc(Hds_utils::get_arc(edge->get_vertex_src()));
    result->add_arc(Hds_utils::get_arc(edge->get_vertex_dst()));

    return result;
  }

  static void set_sheet(EdgeSPtr edge,
                        SheetSPtr sheet)
  {
    Hds_utils::set_sheet(edge, sheet);
    sheet->add_edge(edge);
  }

  void merge_sheets(EdgeSPtr edge_into,
                    EdgeSPtr edge_from)
  {
    CGAL_precondition(edge_into && edge_from);
    CGAL_precondition(edge_into != edge_from);

    SheetSPtr sheet_into = Hds_utils::get_sheet(edge_into);
    SheetSPtr sheet_from = Hds_utils::get_sheet(edge_from);
    CGAL_precondition(sheet_into && sheet_from);

    std::cout << "merge E" << edge_from->get_ID() << " S" << sheet_from->get_ID() << std::endl;
    std::cout << "into E" << edge_into->get_ID() << " S" << sheet_into->get_ID() << std::endl;

    if (sheet_into == sheet_from) {
      return;
    }

    std::list<EdgeWPtr> edges_from = sheet_from->edges(); // intentional copy
    CGAL_assertion(!edges_from.empty());

    skeleton_->merge_sheets(sheet_into, sheet_from);
    CGAL_SS3_DEBUG_SPTR(sheet_into);

    for (EdgeWPtr edge_wptr : edges_from) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        std::cout << "Set Sheet of E" << edge->get_ID() << std::endl;
        set_sheet(edge, sheet_into);
      }
    }

    CGAL_postcondition(Hds_utils::get_sheet(edge_from) == sheet_into);
  }

  /**
    * Store within each facet the coefficients of the plane at t=0
    */
  void cache_base_planes(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    for (const FacetSPtr& facet : polyhedron->facets()) {
      const Plane_3& plane = facet->get_plane();
      CGAL_SS3_DEBUG_SPTR(plane);
      CGAL_assertion(Kernel_wrapper::has_normalized_plane(plane));
      Hds_utils::set_base_plane(facet, plane);
    }
  }

  /**
    * Split all vertices with degree > 3 and
    * initializes the data variables of all edges and vertices.
    */
  bool init(const PolyhedronSPtr& polyhedron,
            Abstract_vertex_splitter_sptr vertex_splitter)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    bool result = true;

    CGAL_SS3_CORE_TRACE("Input: " << polyhedron->vertices().size() << " NV " << polyhedron->facets().size() << " NF");

    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      CGAL_precondition(vertex->degree() >= 3);
      if (!vertex->has_data()) {
        Skeleton_vertex_data::create(vertex);
      }
    }

    // vertex splitting
    std::list<VertexSPtr> vertices_tosplit;

    for (const VertexSPtr& vertex : polyhedron->vertices()) {
#ifndef CGAL_SS3_NO_SKELETON_DS
      NodeSPtr node = create_node(vertex);
      CGAL_SS3_DEBUG_SPTR(node);
      skeleton_->add_node(node);
#endif

      if (vertex->degree() > 3) {
        vertices_tosplit.push_back(vertex);
      }
    }

    CGAL_SS3_CORE_TRACE_V(2, vertices_tosplit.size() << " vertices to split");
    CGAL_SS3_CORE_TRACE_V(2, "Using " << vertex_splitter->to_string() << " to split vertices.");

    for (const VertexSPtr& vertex : vertices_tosplit) {
      CGAL_SS3_SPLITTER_TRACE_V(8, "Generic split vertex:\n" << vertex->to_string());
      if (vertex->degree() > 15) {
        CGAL_SS3_SPLITTER_TRACE_V(1, "Warning: degree of vertex (" << vertex->degree() << ") at "
                                        << vertex->point() << " is extremely high!");
      }
      vertex_splitter->split_vertex(vertex);
    }

#ifdef CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("results/post_split.obj", polyhedron, false /*do not triangulate*/);
#endif

    CGAL_postcondition_code(for (auto v : polyhedron->vertices()))
    CGAL_postcondition(v->get_ID() != -1);
    CGAL_postcondition_code(for (auto e : polyhedron->edges()))
    CGAL_postcondition(e->get_ID() != -1);
    CGAL_postcondition_code(for (auto f : polyhedron->facets()))
    CGAL_postcondition(f->get_ID() != -1);

#ifndef CGAL_SS3_NO_SKELETON_DS
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      CGAL_assertion(vertex->degree() == 3);
      ArcSPtr arc = create_arc(vertex);
      CGAL_SS3_DEBUG_SPTR(arc);
      skeleton_->add_arc(arc);
    }
#endif

    for (const EdgeSPtr& edge : polyhedron->edges()) {
      if (!edge->has_data()) {
        Skeleton_edge_data::create(edge);
      }

#ifndef CGAL_SS3_NO_SKELETON_DS
      SheetSPtr sheet = create_sheet(edge);
      CGAL_SS3_DEBUG_SPTR(sheet);

      // Squatting the arc type, this isn't really an arc, but a contour edge.
      // Contours are useful to get a closed polygon (with holes) to draw sheets.
      ArcSPtr contour = std::make_shared<Arc>(Hds_utils::get_node(edge->get_vertex_src()),
                                              Hds_utils::get_node(edge->get_vertex_dst()));
      sheet->add_contour(contour);

      skeleton_->add_sheet(sheet);
#endif
    }

    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (!facet->has_data()) {
        Skeleton_facet_data::create(facet);
      }
    }

    return result;
  }

  static bool check_bisectors_V2(const EdgeSPtr& edge,
                                 const Point_3& point,
                                 const FT& event_time)
  {
    std::optional<FT> vanish_time = Hds_utils::get_vanish_time(edge);
    Point_3 o_src = Transformation::offset_point_from_base(edge->get_vertex_src(), event_time);
    if (vanish_time.has_value() && event_time == *vanish_time) {
      return (point == o_src);
    }

    Point_3 o_dst = Transformation::offset_point_from_base(edge->get_vertex_dst(), event_time);
    CGAL_assertion(o_src != o_dst);
    CGAL_assertion(CGAL::collinear(o_src, point, o_dst));

    return CGAL::collinear_are_ordered_along_line(o_src, point, o_dst);
  }

  /**
  * Returns the intersection time of the 4 shifting planes.
  */
  static Point_3 intersection_point_offset_planes(const FacetSPtr& facet_0,
                                                  const FacetSPtr& facet_1,
                                                  const FacetSPtr& facet_2,
                                                  const FacetSPtr& facet_3)
  {
    CGAL_SS3_DEBUG_SPTR(facet_0);
    CGAL_SS3_DEBUG_SPTR(facet_1);
    CGAL_SS3_DEBUG_SPTR(facet_2);
    CGAL_SS3_DEBUG_SPTR(facet_3);

    const Plane_3& plane_0 = Hds_utils::get_base_plane(facet_0);
    const Plane_3& plane_1 = Hds_utils::get_base_plane(facet_1);
    const Plane_3& plane_2 = Hds_utils::get_base_plane(facet_2);
    const Plane_3& plane_3 = Hds_utils::get_base_plane(facet_3);
    const FT& speed_0 = Hds_utils::get_speed(facet_0);
    const FT& speed_1 = Hds_utils::get_speed(facet_1);
    const FT& speed_2 = Hds_utils::get_speed(facet_2);
    const FT& speed_3 = Hds_utils::get_speed(facet_3);

    return Geom_utils::intersection_point_offset_planes(plane_0, speed_0, plane_1, speed_1,
                                                        plane_2, speed_2, plane_3, speed_3);
  }

  /**
    * Returns the intersection time of the 4 shifting planes.
    */
  FT intersection_time_offset_planes(const FacetSPtr& facet_0,
                                     const FacetSPtr& facet_1,
                                     const FacetSPtr& facet_2,
                                     const FacetSPtr& facet_3)
  {
    CGAL_SS3_DEBUG_SPTR(facet_0);
    CGAL_SS3_DEBUG_SPTR(facet_1);
    CGAL_SS3_DEBUG_SPTR(facet_2);
    CGAL_SS3_DEBUG_SPTR(facet_3);

    const Plane_3& plane_0 = Hds_utils::get_base_plane(facet_0);
    const Plane_3& plane_1 = Hds_utils::get_base_plane(facet_1);
    const Plane_3& plane_2 = Hds_utils::get_base_plane(facet_2);
    const Plane_3& plane_3 = Hds_utils::get_base_plane(facet_3);
    const FT& speed_0 = Hds_utils::get_speed(facet_0);
    const FT& speed_1 = Hds_utils::get_speed(facet_1);
    const FT& speed_2 = Hds_utils::get_speed(facet_2);
    const FT& speed_3 = Hds_utils::get_speed(facet_3);

    return Geom_utils::intersection_time_offset_planes(plane_0, speed_0, plane_1, speed_1,
                                                       plane_2, speed_2, plane_3, speed_3);
  }

  Point_3 vanish_point(const EdgeSPtr& edge)
  {
    CGAL_SS3_CORE_TRACE_V(16, "VanishesAt() for " << edge->to_string());

    CGAL_SS3_DEBUG_SPTR(edge);

    FacetSPtr facetL = edge->get_facet_L();
    FacetSPtr facetR = edge->get_facet_R();
    CGAL_SS3_CORE_TRACE_V(16, "facetL: " << facetL->get_ID());
    CGAL_SS3_CORE_TRACE_V(16, "facetR: " << facetR->get_ID());
    CGAL_assertion(facetL && facetR && facetL != facetR);

    FacetSPtr facetP = edge->prev(facetL)->other(facetL);
    FacetSPtr facetN = edge->next(facetL)->other(facetL);
    CGAL_SS3_CORE_TRACE_V(16, "facetP: " << facetP->get_ID());
    CGAL_SS3_CORE_TRACE_V(16, "facetN: " << facetN->get_ID());
    CGAL_assertion(facetP && facetP != facetL && facetP != facetR);
    CGAL_assertion(facetN && facetN != facetL && facetN != facetR && facetN != facetP);

    return intersection_point_offset_planes(facetL, facetP, facetR, facetN);
  }

  /**
    * Returns the time at which the edge will vanish.
    */
  std::optional<FT> vanish_time(const EdgeSPtr& edge,
                                const std::optional<FT>& time_past_bound = std::nullopt,
                                const std::optional<FT>& time_future_bound = std::nullopt)
  {
    CGAL_SS3_CORE_TRACE_V(16, "Computing vanish time of " << edge->to_string());

    CGAL_SS3_DEBUG_SPTR(edge);

    FacetSPtr facetL = edge->get_facet_L();
    FacetSPtr facetR = edge->get_facet_R();
    CGAL_SS3_CORE_TRACE_V(16, "facetL: " << facetL->get_ID());
    CGAL_SS3_CORE_TRACE_V(16, "facetR: " << facetR->get_ID());
    CGAL_assertion(facetL && facetR && facetL != facetR);

    FacetSPtr facetP = edge->prev(facetL)->other(facetL);
    FacetSPtr facetN = edge->next(facetL)->other(facetL);
    CGAL_SS3_CORE_TRACE_V(16, "facetP: " << facetP->get_ID());
    CGAL_SS3_CORE_TRACE_V(16, "facetN: " << facetN->get_ID());
    CGAL_assertion(facetP && facetP != facetL && facetP != facetR);
    CGAL_assertion(facetN && facetN != facetL && facetN != facetR && facetN != facetP);

    const FT vanish_time = intersection_time_offset_planes(facetL, facetP, facetR, facetN);

    if (time_past_bound && vanish_time >= *time_past_bound) {
      CGAL_SS3_TRAITS_TRACE("Vanish event is strictly in the past");
      Hds_utils::set_vanish_time(edge, std::nullopt);
      return { };
    }

    if (time_future_bound && vanish_time < *time_future_bound) {
      CGAL_SS3_TRAITS_TRACE("Vanish event is too far in the future");
      Hds_utils::set_vanish_time(edge, std::nullopt);
      return { };
    }

    Hds_utils::set_vanish_time(edge, vanish_time);
    return vanish_time;
  }

  // for pierce events only
  std::optional<FT> crash_time(const VertexSPtr& vertex, const FacetSPtr& facet,
                               const std::optional<FT>& time_past_bound = std::nullopt,
                               const std::optional<FT>& time_future_bound = std::nullopt)
  {
    CGAL_SS3_CORE_TRACE_V(16, "-- Crash At Time (Pierce event)\n  " << vertex->to_string() << "\n  " << facet->to_string());

    CGAL_assertion(vertex->facets().size() >= 3);
    std::array<FacetSPtr, 3> fs;
    for (int i = 0; i < 3; ++i) {
      FacetWPtr wf = *(std::next(vertex->facets().begin(), i));
      CGAL_assertion(!wf.expired());
      fs[i] = wf.lock();
    }

    FT event_time = intersection_time_offset_planes(facet, fs[0], fs[1], fs[2]);

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

  std::optional<FT> crash_time(const EdgeSPtr& edge_1, const EdgeSPtr& edge_2,
                               const std::optional<FT>& time_past_bound = std::nullopt,
                               const std::optional<FT>& time_future_bound = std::nullopt)
  {
    CGAL_SS3_CORE_TRACE_V(16, "-- Crash At Time\n  " << edge_1->to_string() << "\n  " << edge_2->to_string());

    CGAL_SS3_DEBUG_SPTR(edge_1);
    CGAL_SS3_DEBUG_SPTR(edge_2);

    FacetSPtr facet_l1 = edge_1->get_facet_L();
    FacetSPtr facet_r1 = edge_1->get_facet_R();
    FacetSPtr facet_l2 = edge_2->get_facet_L();
    FacetSPtr facet_r2 = edge_2->get_facet_R();

    CGAL_SS3_CORE_TRACE_V(16, "Facet L1 = " << facet_l1->get_ID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet R1 = " << facet_r1->get_ID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet L2 = " << facet_l2->get_ID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet R2 = " << facet_r2->get_ID());

    // It is pointless to check for contact events that are farther in time than any of the two
    // involved edges as either they will be gone, or new events will be computed
    std::optional<FT> tight_future_bound = time_future_bound;
    for (const EdgeSPtr& edge : {edge_1, edge_2}) {
      const std::optional<FT>& vanish_time = Hds_utils::get_vanish_time(edge);
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

    FT event_time = intersection_time_offset_planes(facet_l1, facet_r1, facet_l2, facet_r2);

    if (time_past_bound && event_time >= *time_past_bound) {
      CGAL_SS3_TRAITS_TRACE("Contact event is strictly in the past");
      return { };
    }

    if (tight_future_bound && event_time < *tight_future_bound) {
      CGAL_SS3_TRAITS_TRACE("Contact event is too far in the future");
      return { };
    }

    CGAL_SS3_CORE_TRACE_V(16, "Tentative event @ " << event_time);
    return event_time;
  }

  /**
    * Returns `true` if the event is in the past
    */
  static bool is_event_in_the_past(const Abstract_event_sptr& event,
                                   const FT& current_time)
  {
    CGAL_SS3_DEBUG_SPTR(event);
    CGAL_precondition(event->is_valid());
    return event->time() >= current_time;
  }

  /**
    * Returns `true` if the neighborhood of an event has changed.
    */
  static bool is_event_obsolete(const Abstract_event_sptr& event)
  {
    CGAL_SS3_DEBUG_SPTR(event);
    CGAL_precondition(event->is_valid());
    return event->is_obsolete();
  }

  static bool is_actual_vertex_event(const Vertex_event_sptr& event)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "#######  Tentative Vertex Event  #######");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->time();

    VertexSPtr vertex_1 = event->get_vertex_1();
    VertexSPtr vertex_2 = event->get_vertex_2();
    FacetSPtr facet_1 = event->get_facet_1();
    FacetSPtr facet_2 = event->get_facet_2();

    // @todo avoid all this duplication...
    EdgeSPtr edge_11 = EdgeSPtr();
    for (EdgeWPtr edge_1_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge_1 = edge_1_wptr.lock()) {
        FacetSPtr facet_1l = edge_1->get_facet_L();
        FacetSPtr facet_1r = edge_1->get_facet_R();
        if ((facet_1l == facet_1 && facet_1r != facet_2) ||
            (facet_1r == facet_1 && facet_1l != facet_2)) {
          edge_11 = edge_1;
        }
      }
    }

    EdgeSPtr edge_22 = EdgeSPtr();
    for (EdgeWPtr edge_2_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge_2 = edge_2_wptr.lock()) {
        FacetSPtr facet_2l = edge_2->get_facet_L();
        FacetSPtr facet_2r = edge_2->get_facet_R();
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
      if ((edge_cur->get_facet_L() == facet_1b && edge_cur->get_facet_R() == facet_2b) ||
          (edge_cur->get_facet_R() == facet_1b && edge_cur->get_facet_L() == facet_2b)) {
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
    Point_3 point = intersection_point_offset_planes(edge_11->get_facet_L(),
                                                     edge_11->get_facet_R(),
                                                     edge_22->get_facet_L(),
                                                     edge_22->get_facet_R());

    if (!check_bisectors_V2(edge_11, point, event_time) ||
        !check_bisectors_V2(edge_22, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Vertex event: bisector check failure");
      return false;
    }

    event->set_point(point);

    CGAL_SS3_CORE_TRACE_V(8, "Vertex event: accepted");
    return true;
  }

  static bool is_actual_flip_vertex_event(const Flip_vertex_event_sptr& event)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "#####  Tentative Flip Vertex Event  ####");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->time();
    VertexSPtr vertex_1 = event->get_vertex_1();
    VertexSPtr vertex_2 = event->get_vertex_2();
    FacetSPtr facet_1 = event->get_facet_1();
    FacetSPtr facet_2 = event->get_facet_2();

    // convex split event checks
    EdgeSPtr edge_11 = EdgeSPtr();
    for (EdgeWPtr edge_1_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge_1 = edge_1_wptr.lock()) {
        FacetSPtr facet_1l = edge_1->get_facet_L();
        FacetSPtr facet_1r = edge_1->get_facet_R();
        if ((facet_1l == facet_1 && facet_1r != facet_2) ||
            (facet_1r == facet_1 && facet_1l != facet_2)) {
          edge_11 = edge_1;
        }
      }
    }

    EdgeSPtr edge_22 = EdgeSPtr();
    for (EdgeWPtr edge_2_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge_2 = edge_2_wptr.lock()) {
        FacetSPtr facet_2l = edge_2->get_facet_L();
        FacetSPtr facet_2r = edge_2->get_facet_R();
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
      if ((edge_cur->get_facet_L() == facet_1b && edge_cur->get_facet_R() == facet_2b) ||
          (edge_cur->get_facet_R() == facet_1b && edge_cur->get_facet_L() == facet_2b)) {
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
    Point_3 point = intersection_point_offset_planes(edge_11->get_facet_L(),
                                                     edge_11->get_facet_R(),
                                                     edge_22->get_facet_L(),
                                                     edge_22->get_facet_R());

    if (!check_bisectors_V2(edge_11, point, event_time) ||
        !check_bisectors_V2(edge_22, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Flip Vertex event: bisector check failure");
      return false;
    }

    event->set_point(point);

    CGAL_SS3_CORE_TRACE_V(8, "Flip vertex event: accepted");
    return true;
  }

  static bool is_actual_surface_event(const Surface_event_sptr& event)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "######  Tentative Surface Event  #######");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->time();
    EdgeSPtr edge_1 = event->get_edge_1();
    EdgeSPtr edge_2 = event->get_edge_2();
    FacetSPtr facet_1_src = edge_1->get_facet_src();
    FacetSPtr facet_1_dst = edge_1->get_facet_dst();

    // convex split checks
    bool conv_split_event = false;
    std::list<EdgeSPtr> common_edges = edge_1->get_facet_L()->find_edges(edge_1->get_facet_R());
    for (const EdgeSPtr& edge : common_edges) {
      if (edge == edge_1) {
        continue;
      }
      FacetSPtr facet_src = edge->get_facet_src();
      FacetSPtr facet_dst = edge->get_facet_dst();
      if (facet_1_src == edge_2->get_facet_L() ||
          facet_1_dst == edge_2->get_facet_L()) {
        if (facet_src == edge_2->get_facet_R() ||
            facet_dst == edge_2->get_facet_R()) {
          conv_split_event = true;
          break;
        }
      } else if (facet_1_src == edge_2->get_facet_R() ||
                 facet_1_dst == edge_2->get_facet_R()) {
        if (facet_src == edge_2->get_facet_L() ||
            facet_dst == edge_2->get_facet_L()) {
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
    Point_3 point = intersection_point_offset_planes(edge_1->get_facet_L(),
                                                     edge_1->get_facet_R(),
                                                     edge_2->get_facet_L(),
                                                     edge_2->get_facet_R());

    if (!check_bisectors_V2(edge_1, point, event_time) ||
        !check_bisectors_V2(edge_2, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Surface event: bisector check failure");
      return false;
    }

    event->set_point(point);

    CGAL_SS3_CORE_TRACE_V(8, "Surface event: accepted");
    return true;
  }

  static bool is_actual_polyhedron_split_event(const Polyhedron_split_event_sptr& event,
                                               const FT& current_time)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "####  Tentative Split Merge Event  #####");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->time();

    // Bisector check
    EdgeSPtr edge_1 = event->get_edge_1();
    EdgeSPtr edge_2 = event->get_edge_2();

    Point_3 point = intersection_point_offset_planes(edge_1->get_facet_L(),
                                                     edge_1->get_facet_R(),
                                                     edge_2->get_facet_L(),
                                                     edge_2->get_facet_R());

    if (!check_bisectors_V2(edge_1, point, event_time) ||
        !check_bisectors_V2(edge_2, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Polyhedron split event: bisector check failure");
      return false;
    }

    event->set_point(point);

    // @speed is_degenerate(4 planes)? But it would be a sure filter failure, so, costly...
    FT shift = event_time - current_time;
    Segment_3 e1o = Transformation::shift_edge(event->get_edge_1(), shift);
    if (!e1o.is_degenerate()) {
      CGAL_SS3_CORE_TRACE_V(8, "Polyhedron split event: bisector check failure");
      return false;
    }

    CGAL_SS3_CORE_TRACE_V(8, "Polyhedron split event accepted");
    return true;
  }

  static bool is_actual_split_merge_event(const Split_merge_event_sptr& event)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "####  Tentative Split Merge Event  #####");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->time();
    VertexSPtr vertex_1 = event->get_vertex_1();
    VertexSPtr vertex_2 = event->get_vertex_2();
    FacetSPtr facet_1 = event->get_facet_1();
    FacetSPtr facet_2 = event->get_facet_2();

    // convex split checks
    EdgeSPtr edge_11 = EdgeSPtr();
    for (EdgeWPtr edge_1_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge_1 = edge_1_wptr.lock()) {
        FacetSPtr facet_1l = edge_1->get_facet_L();
        FacetSPtr facet_1r = edge_1->get_facet_R();
        if ((facet_1l == facet_1 && facet_1r != facet_2) ||
            (facet_1r == facet_1 && facet_1l != facet_2)) {
          edge_11 = edge_1;
        }
      }
    }

    EdgeSPtr edge_22 = EdgeSPtr();
    for (EdgeWPtr edge_2_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge_2 = edge_2_wptr.lock()) {
        FacetSPtr facet_2l = edge_2->get_facet_L();
        FacetSPtr facet_2r = edge_2->get_facet_R();
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
      if ((edge_cur->get_facet_L() == facet_1b && edge_cur->get_facet_R() == facet_2b) ||
          (edge_cur->get_facet_R() == facet_1b && edge_cur->get_facet_L() == facet_2b)) {
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
    Point_3 point = intersection_point_offset_planes(edge_11->get_facet_L(),
                                                     edge_11->get_facet_R(),
                                                     edge_22->get_facet_L(),
                                                     edge_22->get_facet_R());

    if (!check_bisectors_V2(edge_11, point, event_time) ||
        !check_bisectors_V2(edge_22, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Split merge event: bisector check failure");
      return false;
    }

    event->set_point(point);

    CGAL_SS3_CORE_TRACE_V(8, "Split merge event accepted");
    return true;
  }

  static bool is_actual_edge_split_event(const Edge_split_event_sptr& event)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "#####  Tentative Edge Split Event  #####");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->time();
    EdgeSPtr edge_1 = event->get_edge_1(); // @todo const& all of this
    EdgeSPtr edge_2 = event->get_edge_2();

    // Bisector check
    Point_3 point = intersection_point_offset_planes(edge_1->get_facet_L(),
                                                     edge_1->get_facet_R(),
                                                     edge_2->get_facet_L(),
                                                     edge_2->get_facet_R());
    CGAL_SS3_DEBUG_SPTR(point);

    if (!check_bisectors_V2(edge_1, point, event_time) ||
        !check_bisectors_V2(edge_2, point, event_time)) {
      CGAL_SS3_CORE_TRACE_V(8, "Edge split rejected at pop time");
      return false;
    }

    event->set_point(point);

    CGAL_SS3_CORE_TRACE_V(8, "Edge split event accepted");
    return true;
  }

  static bool is_actual_pierce_event(const Pierce_event_sptr& event,
                                     const FT& current_time,
                                     const std::optional<FT>& time_future_bound)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "######  Tentative Pierce Event  ########");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    VertexSPtr pv = event->get_vertex();
    FacetSPtr pf = event->get_facet();

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
        FacetSPtr facet_src = edge->get_facet_src();
        FacetSPtr facet_dst = edge->get_facet_dst();
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
      b1 += pv->point().bbox();
      b1 += Hds_utils::get_final_point(pv, *time_future_bound).bbox();

      CGAL::Bbox_3 b2;
      for (const VertexSPtr& v : pf->vertices()) {
        b2 += v->point().bbox();
        b2 += Hds_utils::get_final_point(v, *time_future_bound).bbox();
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

    event->set_point(intersection_point_offset_planes(pf, fs[0], fs[1], fs[2]));

    // Filter #3
    //
    // Filter if the event point is on an edge (and a fortiori on a vertex)
    // as it will be a different kind of event
    const Point_3& point = event->point();
    FacetSPtr facet_clone = pf->clone();

    const FT shift = event->time() - current_time;
    const FT& speed = Hds_utils::get_speed(pf);
    const Plane_3 offset_plane = Geom_utils::offset_plane(pf->get_plane(), shift*speed);
    facet_clone->set_plane(offset_plane);

    // abusing the fact that vertices will have the same order in both facets
    typename std::list<VertexSPtr>::iterator it_v = pf->vertices().begin();
    typename std::list<VertexSPtr>::iterator it_v_offset = facet_clone->vertices().begin();
    while (it_v != pf->vertices().end()) {
      VertexSPtr vertex = *it_v++;
      VertexSPtr offset_vertex = *it_v_offset++;
      Point_3 point_offset = Transformation::shift_point(vertex, shift);
      offset_vertex->set_point(point_offset);
    }

#ifdef CGAL_SLS3_NEW_IS_INSIDE
    if (!Self_intersection::is_inside_using_ray_shooting_V2(point, facet_clone)) {
      CGAL_SS3_CORE_TRACE_V(8, "Pierce event rejected at pop time (c)");
      return false;
    }
#else
    if (!Self_intersection::is_inside_using_ray_shooting(point, facet_clone)) {
      CGAL_SS3_CORE_TRACE_V(8, "Pierce event rejected at pop time (c)");
      return false;
    }

    bool boundary_rejection = false;
    for (const EdgeSPtr& edge : facet_clone->edges()) {
      Segment_3 seg { edge->get_vertex_src()->point(),
                      edge->get_vertex_dst()->point() };
      if (seg.is_degenerate()) {
        continue;
      }

      if (seg.has_on(point)) {
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
  static bool is_actual_event(const Abstract_event_sptr& event,
                              const FT& current_time,
                              const std::optional<FT>& time_future_bound)
  {
    CGAL_SS3_DEBUG_SPTR(event);
    CGAL_precondition(event->is_valid());

    bool result = true;

    if (event->getType() == Abstract_event::VERTEX_EVENT) {
      result = is_actual_vertex_event(std::dynamic_pointer_cast<Vertex_event>(event));
    } else if (event->getType() == Abstract_event::FLIP_VERTEX_EVENT) {
      result = is_actual_flip_vertex_event(std::dynamic_pointer_cast<Flip_vertex_event>(event));
    } else if (event->getType() == Abstract_event::SURFACE_EVENT) {
      result = is_actual_surface_event(std::dynamic_pointer_cast<Surface_event>(event));
    } else if (event->getType() == Abstract_event::POLYHEDRON_SPLIT_EVENT) {
      result = is_actual_polyhedron_split_event(std::dynamic_pointer_cast<Polyhedron_split_event>(event), current_time);
    } else if (event->getType() == Abstract_event::SPLIT_MERGE_EVENT) {
      result = is_actual_split_merge_event(std::dynamic_pointer_cast<Split_merge_event>(event));
    } else if (event->getType() == Abstract_event::EDGE_SPLIT_EVENT) {
      result = is_actual_edge_split_event(std::dynamic_pointer_cast<Edge_split_event>(event));
    } else if (event->getType() == Abstract_event::PIERCE_EVENT) {
      result = is_actual_pierce_event(std::dynamic_pointer_cast<Pierce_event>(event), current_time, time_future_bound);
    }

    return result;
  }

  /**
    * Vanish events.
    */
  void collect_vanish_events(const std::list<EdgeSPtr>& edges,
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

      VertexSPtr vertex_src = edge->get_vertex_src();
      VertexSPtr vertex_dst = edge->get_vertex_dst();
      if (vertex_src->point() == vertex_dst->point()) {
        Hds_utils::set_vanish_time(edge, std::nullopt);
        continue;
      }

      std::optional<FT> event_time = vanish_time(edge, current_time, time_future_bound);
      if (!event_time) {
        continue;
      }

      CGAL_assertion(*event_time < current_time);
      CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

      Vanish_event_sptr event = Vanish_event::create();
      event->set_time(*event_time);
      event->set_point(vanish_point(edge));
      event->set_edge(edge);
      queue.push(event);
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Vanish Events in: " << timer.time());
#endif
  }

  void collect_vanish_events(const PolyhedronSPtr& polyhedron,
                             const FT& current_time,
                             const std::optional<FT>& time_future_bound,
                             PQ& queue)
  {
    return collect_vanish_events(polyhedron->edges(), polyhedron, current_time, time_future_bound, queue);
  }

  /**
    * Two vertices crash into each other.
    */
  void collect_vertex_events(const std::list<VertexSPtr>& vertices,
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

      if (Hds_utils::is_convex(vertex_1)) {
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
          CGAL_assertion(vertex_1->get_ID() != -1 && vertex_2->get_ID() != -1);
          if (vertex_1->get_ID() > vertex_2->get_ID()) {
            continue;
          }
        }
#endif
        if (vertex_1->point() == vertex_2->point()) {
          continue;
        }
        if (vertex_1->find_edge(vertex_2)) {
          // edge event
          continue;
        }
        if (Hds_utils::is_convex(vertex_2)) {
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
        if (v1->get_ID() > v2->get_ID()) {
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
            FacetSPtr facet_1l = edge_1->get_facet_L();
            FacetSPtr facet_1r = edge_1->get_facet_R();
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
            FacetSPtr facet_2l = edge_2->get_facet_L();
            FacetSPtr facet_2r = edge_2->get_facet_R();
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

        // convex split event checks are performed at pop time - see is_actual_vertex_event()

        std::optional<FT> event_time = crash_time(edge_11, edge_22, current_time, time_future_bound);
        if (!event_time) {
          continue;
        }

        CGAL_assertion(*event_time < current_time);
        CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

        Vertex_event_sptr event = Vertex_event::create();
        event->set_time(*event_time);
        event->set_vertex_1(v1);
        event->set_vertex_2(v2);
        event->set_facet_1(facet_1);
        event->set_facet_2(facet_2);
        queue.push(event);
      }
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Vertex Events in: " << timer.time());
#endif
  }

  void collect_vertex_events(const PolyhedronSPtr& polyhedron,
                             const FT& current_time,
                             const std::optional<FT>& time_future_bound,
                             PQ& queue)
  {
    return collect_vertex_events(polyhedron->vertices(), polyhedron, true /*use canonical reps*/,
                               current_time, time_future_bound, queue);
  }

  /**
    * Flip vertex event
    */
  void collect_flip_vertex_events(const std::list<VertexSPtr>& vertices,
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
      CGAL_assertion(vertex_1->get_ID() != -1);

      if (Hds_utils::is_convex(vertex_1)) {
        continue;
      }

      std::set<VertexSPtr> vertices_2;
      for (FacetWPtr facet_wptr : vertex_1->facets()) {
        if (FacetSPtr facet = facet_wptr.lock()) {
          vertices_2.insert(facet->vertices().begin(), facet->vertices().end());
        }
      }

      for (const VertexSPtr& vertex_2 : vertices_2) {
        CGAL_assertion(vertex_2->get_ID() != -1);

        if (vertex_1 == vertex_2) {
          continue;
        }
#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        if (use_canonical_event_reps) {
          if (vertex_1->get_ID() > vertex_2->get_ID()) {
            continue;
          }
        }
#endif
        if (vertex_1->point() == vertex_2->point()) {
          continue;
        }
        if (vertex_1->find_edge(vertex_2)) {
          // edge event
          continue;
        }
        if (Hds_utils::is_convex(vertex_2)) {
          continue;
        }

        // See comment in the first instance of this swap
        VertexSPtr v1 = vertex_1;
        VertexSPtr v2 = vertex_2;
        if (v1->get_ID() > v2->get_ID()) {
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
          if (facet_1->get_ID() > facet_2->get_ID()) {
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
            FacetSPtr facet_1l = edge_1->get_facet_L();
            FacetSPtr facet_1r = edge_1->get_facet_R();
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
            FacetSPtr facet_2l = edge_2->get_facet_L();
            FacetSPtr facet_2r = edge_2->get_facet_R();
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

        // convex split event checks are performed at pop time - see is_actual_flip_vertex_event()

        std::optional<FT> event_time = crash_time(edge_11, edge_22, current_time, time_future_bound);
        if (!event_time) {
          continue;
        }

        CGAL_assertion(*event_time < current_time);
        CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

        Flip_vertex_event_sptr event = Flip_vertex_event::create();
        event->set_time(*event_time);
        event->set_vertex_1(v1);
        event->set_vertex_2(v2);
        event->set_facet_1(facet_1);
        event->set_facet_2(facet_2);
        queue.push(event);
      }
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Flip Vertex Events in: " << timer.time());
#endif
  }

  void collect_flip_vertex_events(const PolyhedronSPtr& polyhedron,
                                  const FT& current_time,
                                  const std::optional<FT>& time_future_bound,
                                  PQ& queue)
  {
    return collect_flip_vertex_events(polyhedron->vertices(), polyhedron, true /*use canonical reps*/,
                                   current_time, time_future_bound, queue);
  }

  /**
    * Split event on the surface.
    * Edges do not need to be reflex.
    */
  void collect_surface_event(const EdgeSPtr& edge_1,
                             const EdgeSPtr& edge_2,
                             const PolyhedronSPtr& /*polyhedron*/,
                             const FT& current_time,
                             const std::optional<FT>& time_future_bound,
                             PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(8, ">>> Collect Surface Event [\n  " << edge_1->to_string() << "\n  " << edge_2->to_string() << "]");

    CGAL_SS3_DEBUG_SPTR(edge_1);
    CGAL_SS3_DEBUG_SPTR(edge_2);

    if (edge_1 == edge_2) {
      return;
    }

    FacetSPtr facet_1_src = edge_1->get_facet_src();
    FacetSPtr facet_1_dst = edge_1->get_facet_dst();

    if (edge_1->get_facet_L() == edge_2->get_facet_L() ||
        edge_1->get_facet_L() == edge_2->get_facet_R() ||
        edge_1->get_facet_R() == edge_2->get_facet_L() ||
        edge_1->get_facet_R() == edge_2->get_facet_R()) {
      // on same facet
      return;
    }

    if (edge_1->get_vertex_src()->point() == edge_2->get_vertex_src()->point() ||
        edge_1->get_vertex_src()->point() == edge_2->get_vertex_dst()->point() ||
        edge_1->get_vertex_dst()->point() == edge_2->get_vertex_src()->point() ||
        edge_1->get_vertex_dst()->point() == edge_2->get_vertex_dst()->point()) {
      // share a vertex
      return;
    }

    // vertex of edge_1 splits edge_2
    if (!((edge_2->get_facet_L() == facet_1_src && edge_2->get_facet_R() != facet_1_dst) ||
          (edge_2->get_facet_L() == facet_1_dst && edge_2->get_facet_R() != facet_1_src) ||
          (edge_2->get_facet_R() == facet_1_src && edge_2->get_facet_L() != facet_1_dst) ||
          (edge_2->get_facet_R() == facet_1_dst && edge_2->get_facet_L() != facet_1_src))) {
      // no surface event
      return;
    }

    FacetSPtr facet_2_src = edge_2->get_facet_src();
    FacetSPtr facet_2_dst = edge_2->get_facet_dst();
    if ((edge_1->get_facet_L() == facet_2_src && edge_1->get_facet_R() != facet_2_dst) ||
        (edge_1->get_facet_L() == facet_2_dst && edge_1->get_facet_R() != facet_2_src) ||
        (edge_1->get_facet_R() == facet_2_src && edge_1->get_facet_L() != facet_2_dst) ||
        (edge_1->get_facet_R() == facet_2_dst && edge_1->get_facet_L() != facet_2_src)) {
      // flip vertex event
      return;
    }

    if (edge_1->get_vertex_src()->find_edge(edge_2->get_vertex_src()) ||
        edge_1->get_vertex_src()->find_edge(edge_2->get_vertex_dst()) ||
        edge_1->get_vertex_dst()->find_edge(edge_2->get_vertex_src()) ||
        edge_1->get_vertex_dst()->find_edge(edge_2->get_vertex_dst()) ) {
      // edge event (when a pyramid grows outwards)
      // a surface split is not possible with only one edge in between
      return;
    }

    if ((edge_1->get_facet_L() == facet_2_src && facet_1_src == edge_2->get_facet_L()) ||
        (edge_1->get_facet_L() == facet_2_dst && facet_1_src == edge_2->get_facet_R()) ||
        (edge_1->get_facet_R() == facet_2_src && facet_1_dst == edge_2->get_facet_L()) ||
        (edge_1->get_facet_R() == facet_2_dst && facet_1_dst == edge_2->get_facet_R()) ||
        (edge_1->get_facet_R() == facet_2_src && facet_1_src == edge_2->get_facet_R()) ||
        (edge_1->get_facet_R() == facet_2_dst && facet_1_src == edge_2->get_facet_L()) ||
        (edge_1->get_facet_L() == facet_2_src && facet_1_dst == edge_2->get_facet_R()) ||
        (edge_1->get_facet_L() == facet_2_dst && facet_1_dst == edge_2->get_facet_L())) {
      // vertex event
      return;
    }

    // convex split event are performed at pop time - see is_actual_surface_event()

    // let's just check if bboxes overlap first
    if (time_future_bound.has_value()) {
      CGAL::Bbox_3 b1;
      b1 += edge_1->get_vertex_src()->point().bbox();
      b1 += edge_1->get_vertex_dst()->point().bbox();
      b1 += Hds_utils::get_final_point(edge_1->get_vertex_src(), *time_future_bound).bbox();
      b1 += Hds_utils::get_final_point(edge_1->get_vertex_dst(), *time_future_bound).bbox();

      CGAL::Bbox_3 b2;
      b2 += edge_2->get_vertex_src()->point().bbox();
      b2 += edge_2->get_vertex_dst()->point().bbox();
      b2 += Hds_utils::get_final_point(edge_2->get_vertex_src(), *time_future_bound).bbox();
      b2 += Hds_utils::get_final_point(edge_2->get_vertex_dst(), *time_future_bound).bbox();

      if (!CGAL::do_overlap(b1, b2)) {
        CGAL_SS3_CORE_TRACE_V(32, "Filtered possible surface event candidates\n\t" << edge_1->to_string() << "\n\t"
                                                                                   << edge_2->to_string());
        return;
      } else {
        CGAL_SS3_CORE_TRACE_V(32, "Checking possible surface event event\n\t" << edge_1->to_string() << "\n\t"
                                                                              << edge_2->to_string());
      }
    }

    std::optional<FT> event_time = crash_time(edge_1, edge_2, current_time, time_future_bound);
    if (!event_time) {
      return;
    }

    CGAL_assertion(*event_time < current_time);
    CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

    Surface_event_sptr event = Surface_event::create();
    event->set_time(*event_time);
    event->set_edge_1(edge_1);
    event->set_edge_2(edge_2);
    queue.push(event);
  }

  void collect_surface_events(const std::list<EdgeSPtr>& edges,
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

      FacetSPtr facet_1_src = edge_1->get_facet_src();
      FacetSPtr facet_1_dst = edge_1->get_facet_dst();
      std::list<EdgeSPtr> edges_2;
      edges_2.insert(edges_2.end(), facet_1_src->edges().begin(), facet_1_src->edges().end());
      edges_2.insert(edges_2.end(), facet_1_dst->edges().begin(), facet_1_dst->edges().end());

      // not outside of the loop just because maybe one day this will be called
      // as the first collect function with an initial bound that gets updated...
      CGAL::Bbox_3 b1;
      if (time_future_bound.has_value()) {
        b1 += edge_1->get_vertex_src()->point().bbox();
        b1 += edge_1->get_vertex_dst()->point().bbox();
        b1 += Hds_utils::get_final_point(edge_1->get_vertex_src(), *time_future_bound).bbox();
        b1 += Hds_utils::get_final_point(edge_1->get_vertex_dst(), *time_future_bound).bbox();
      }

      for (const EdgeSPtr& edge_2 : edges_2) {
        CGAL_SS3_DEBUG_SPTR(edge_2);

        if (edge_1 == edge_2) {
          continue;
        }

        CGAL_SS3_CORE_TRACE_V(64, "Possible surface event:\n\t" << edge_1->to_string() << "\n\t"
                                                                << edge_2->to_string());

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
        ++total_candidates;
#endif

        if (edge_1->get_facet_L() == edge_2->get_facet_L() ||
            edge_1->get_facet_L() == edge_2->get_facet_R() ||
            edge_1->get_facet_R() == edge_2->get_facet_L() ||
            edge_1->get_facet_R() == edge_2->get_facet_R()) {
          // on same facet
          continue;
        }

        if (edge_1->get_vertex_src()->point() == edge_2->get_vertex_src()->point() ||
            edge_1->get_vertex_src()->point() == edge_2->get_vertex_dst()->point() ||
            edge_1->get_vertex_dst()->point() == edge_2->get_vertex_src()->point() ||
            edge_1->get_vertex_dst()->point() == edge_2->get_vertex_dst()->point()) {
          // share a vertex
          continue;
        }

        // vertex of edge_1 splits edge_2
        if (!((edge_2->get_facet_L() == facet_1_src && edge_2->get_facet_R() != facet_1_dst) ||
              (edge_2->get_facet_L() == facet_1_dst && edge_2->get_facet_R() != facet_1_src) ||
              (edge_2->get_facet_R() == facet_1_src && edge_2->get_facet_L() != facet_1_dst) ||
              (edge_2->get_facet_R() == facet_1_dst && edge_2->get_facet_L() != facet_1_src))) {
          // no surface event
          continue;
        }

        FacetSPtr facet_2_src = edge_2->get_facet_src();
        FacetSPtr facet_2_dst = edge_2->get_facet_dst();
        if ((edge_1->get_facet_L() == facet_2_src && edge_1->get_facet_R() != facet_2_dst) ||
            (edge_1->get_facet_L() == facet_2_dst && edge_1->get_facet_R() != facet_2_src) ||
            (edge_1->get_facet_R() == facet_2_src && edge_1->get_facet_L() != facet_2_dst) ||
            (edge_1->get_facet_R() == facet_2_dst && edge_1->get_facet_L() != facet_2_src)) {
          // flip vertex event
          continue;
        }

        if (edge_1->get_vertex_src()->find_edge(edge_2->get_vertex_src()) ||
            edge_1->get_vertex_src()->find_edge(edge_2->get_vertex_dst()) ||
            edge_1->get_vertex_dst()->find_edge(edge_2->get_vertex_src()) ||
            edge_1->get_vertex_dst()->find_edge(edge_2->get_vertex_dst()) ) {
          // edge event (when a pyramid grows outwards)
          // a surface split is not possible with only one edge in between
          continue;
        }

        if ((edge_1->get_facet_L() == facet_2_src && facet_1_src == edge_2->get_facet_L()) ||
            (edge_1->get_facet_L() == facet_2_dst && facet_1_src == edge_2->get_facet_R()) ||
            (edge_1->get_facet_R() == facet_2_src && facet_1_dst == edge_2->get_facet_L()) ||
            (edge_1->get_facet_R() == facet_2_dst && facet_1_dst == edge_2->get_facet_R()) ||
            (edge_1->get_facet_R() == facet_2_src && facet_1_src == edge_2->get_facet_R()) ||
            (edge_1->get_facet_R() == facet_2_dst && facet_1_src == edge_2->get_facet_L()) ||
            (edge_1->get_facet_L() == facet_2_src && facet_1_dst == edge_2->get_facet_R()) ||
            (edge_1->get_facet_L() == facet_2_dst && facet_1_dst == edge_2->get_facet_L())) {
          // vertex event
          continue;
        }

        // convex split event checks are performed at pop time - see is_actual_surface_event()

        // let's just check if bboxes overlap first
        if (time_future_bound.has_value()) {
          CGAL::Bbox_3 b2;
          b2 += edge_2->get_vertex_src()->point().bbox();
          b2 += edge_2->get_vertex_dst()->point().bbox();
          b2 += Hds_utils::get_final_point(edge_2->get_vertex_src(), *time_future_bound).bbox();
          b2 += Hds_utils::get_final_point(edge_2->get_vertex_dst(), *time_future_bound).bbox();

          if (!CGAL::do_overlap(b1, b2)) {
            CGAL_SS3_CORE_TRACE_V(64, "Filtered possible surface event candidates\n\t" << edge_1->to_string() << "\n\t"
                                                                                        << edge_2->to_string());
            continue;
          }
        }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
        ++tested_candidates;
#endif

        std::optional<FT> event_time = crash_time(edge_1, edge_2, current_time, time_future_bound);
        if (!event_time) {
          continue;
        }

        CGAL_assertion(*event_time < current_time);
        CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

        Surface_event_sptr event = Surface_event::create();
        event->set_time(*event_time);
        event->set_edge_1(edge_1);
        event->set_edge_2(edge_2);
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

  void collect_surface_events(const PolyhedronSPtr& polyhedron,
                              const FT& current_time,
                              const std::optional<FT>& time_future_bound,
                              PQ& queue)
  {
    return collect_surface_events(polyhedron->edges(), polyhedron,
                                current_time, time_future_bound, queue);
  }

  /**
    * This event occurs when two edges collide.
    * The first edge is always reflex.
    */
  void collect_polyhedron_split_event(const EdgeSPtr& edge_1,
                                      const EdgeSPtr& edge_2,
                                      const PolyhedronSPtr& /*polyhedron*/,
                                      const FT& current_time,
                                      const std::optional<FT>& time_future_bound,
                                      PQ& queue)
  {
    CGAL_SS3_DEBUG_SPTR(edge_1);
    CGAL_SS3_DEBUG_SPTR(edge_2);

    FacetSPtr facet_1_src = edge_1->get_facet_src();
    FacetSPtr facet_1_dst = edge_1->get_facet_dst();

    if (edge_1->get_vertex_src()->point() == edge_2->get_vertex_src()->point() ||
        edge_1->get_vertex_src()->point() == edge_2->get_vertex_dst()->point() ||
        edge_1->get_vertex_dst()->point() == edge_2->get_vertex_src()->point() ||
        edge_1->get_vertex_dst()->point() == edge_2->get_vertex_dst()->point()) {
      // share a vertex
      return;
    }
    if (!((edge_2->get_facet_L() == facet_1_src && edge_2->get_facet_R() == facet_1_dst) ||
          (edge_2->get_facet_L() == facet_1_dst && edge_2->get_facet_R() == facet_1_src))) {
      // no polyhedron split event
      return;
    }
    if (edge_1->get_vertex_src()->find_edge(edge_2->get_vertex_src()) ||
        edge_1->get_vertex_src()->find_edge(edge_2->get_vertex_dst()) ||
        edge_1->get_vertex_dst()->find_edge(edge_2->get_vertex_src()) ||
        edge_1->get_vertex_dst()->find_edge(edge_2->get_vertex_dst())) {
      // does not work when there is only one edge in between
      return;
    }

    std::optional<FT> event_time = crash_time(edge_1, edge_2, current_time, time_future_bound);
    if (!event_time) {
      return;
    }

    CGAL_assertion(*event_time < current_time);
    CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

    // degeneracy check is performed at pop time - see is_actual_polyhedron_split_event()

    Polyhedron_split_event_sptr event = Polyhedron_split_event::create();
    event->set_time(*event_time);
    event->set_edge_1(edge_1);
    event->set_edge_2(edge_2);
    queue.push(event);
  }

  void collect_polyhedron_split_events(const std::list<EdgeSPtr>& edges,
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

      if (!Hds_utils::is_reflex(edge_1)) {
        continue;
      }

      FacetSPtr facet_1_src = edge_1->get_facet_src();
      CGAL_SS3_DEBUG_SPTR(facet_1_src);
      for (const EdgeSPtr& edge_2 : facet_1_src->edges()) {
        collect_polyhedron_split_event(edge_1, edge_2, polyhedron,
                                    current_time, time_future_bound,
                                    queue);
      }
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Polyhedron Split Events in: " << timer.time());
#endif
  }

  void collect_polyhedron_split_events(const PolyhedronSPtr& polyhedron,
                                       const FT& current_time,
                                       const std::optional<FT>& time_future_bound,
                                       PQ& queue)
  {
    return collect_polyhedron_split_events(polyhedron->edges(), polyhedron,
                                        current_time, time_future_bound, queue);
  }

  /**
    * Split Merge event
    */
  void collect_split_merge_events(const std::list<VertexSPtr>& vertices,
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

      if (Hds_utils::is_convex(vertex_1)) {
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
          if (vertex_1->get_ID() > vertex_2->get_ID()) {
            continue;
          }
        }
#endif
        if (vertex_1->point() == vertex_2->point()) {
          continue;
        }
        if (vertex_1->find_edge(vertex_2)) {
          // edge event
          continue;
        }
        if (Hds_utils::is_convex(vertex_2)) {
          continue;
        }

        // See comment in the first instance of this swap
        VertexSPtr v1 = vertex_1;
        VertexSPtr v2 = vertex_2;
        if (v1->get_ID() > v2->get_ID()) {
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
            FacetSPtr facet_1l = edge_1->get_facet_L();
            FacetSPtr facet_1r = edge_1->get_facet_R();
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
            FacetSPtr facet_2l = edge_2->get_facet_L();
            FacetSPtr facet_2r = edge_2->get_facet_R();
            if ((facet_2l == facet_1 && facet_2r != facet_2) ||
                (facet_2r == facet_1 && facet_2l != facet_2)) {
              edge_21 = edge_2;
            } else if ((facet_2l == facet_2 && facet_2r != facet_1) ||
                       (facet_2r == facet_2 && facet_2l != facet_1)) {
              edge_22 = edge_2;
            }
          }
        }

        // convex split event checks are performed at pop time - see is_actual_split_merge_event()

        std::optional<FT> event_time = crash_time(edge_11, edge_22, current_time, time_future_bound);
        if (!event_time) {
          continue;
        }

        CGAL_assertion(*event_time < current_time);
        CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

        Split_merge_event_sptr event = Split_merge_event::create();
        event->set_time(*event_time);
        event->set_vertex_1(v1);
        event->set_vertex_2(v2);
        event->set_facet_1(facet_1);
        event->set_facet_2(facet_2);
        queue.push(event);
      }
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Split Merge Events in: " << timer.time());
#endif
  }

  void collect_split_merge_events(const PolyhedronSPtr& polyhedron,
                                  const FT& current_time,
                                  const std::optional<FT>& time_future_bound,
                                  PQ& queue)
  {
    return collect_split_merge_events(polyhedron->vertices(), polyhedron, true /*use canonical reps*/,
                                   current_time, time_future_bound, queue);
  }

  /**
    * This event occurs when two edges collide.
    * The first edge is always reflex.
    */
  void collect_edge_split_events(const std::list<EdgeSPtr>& edges_1,
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
        if (Hds_utils::is_reflex(edge)) {
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

      FacetSPtr facet_1_src = edge_1->get_facet_src();
      FacetSPtr facet_1_dst = edge_1->get_facet_dst();

      CGAL::Bbox_3 b1;
      if (time_future_bound.has_value()) {
        b1 += edge_1->get_vertex_src()->point().bbox();
        b1 += edge_1->get_vertex_dst()->point().bbox();
        b1 += Hds_utils::get_final_point(edge_1->get_vertex_src(), *time_future_bound).bbox();
        b1 += Hds_utils::get_final_point(edge_1->get_vertex_dst(), *time_future_bound).bbox();
      }

      for (const EdgeSPtr& edge_2 : edges_reflex_2) {
        CGAL_SS3_DEBUG_SPTR(edge_2);

        CGAL_SS3_CORE_TRACE_V(64, "Possible edge split event\n\t" << edge_1->to_string() << "\n\t"
                                                                  << edge_2->to_string());

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
        ++total_candidates;
#endif

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        if (use_canonical_event_reps) {
          if (edge_1->get_ID() > edge_2->get_ID()) {
            continue;
          }
        }
#else
        if (edge_1 == edge_2) {
          continue;
        }
#endif

        if (edge_1->get_facet_L() == edge_2->get_facet_L() ||
            edge_1->get_facet_L() == edge_2->get_facet_R() ||
            edge_1->get_facet_R() == edge_2->get_facet_L() ||
            edge_1->get_facet_R() == edge_2->get_facet_R()) {
          // on same facet
          continue;
        }
        if (edge_1->get_vertex_src()->point() == edge_2->get_vertex_src()->point() ||
            edge_1->get_vertex_src()->point() == edge_2->get_vertex_dst()->point() ||
            edge_1->get_vertex_dst()->point() == edge_2->get_vertex_src()->point() ||
            edge_1->get_vertex_dst()->point() == edge_2->get_vertex_dst()->point()) {
          // share a vertex
          continue;
        }

        if (((edge_2->get_facet_L() == facet_1_src && edge_2->get_facet_R() == facet_1_dst) ||
             (edge_2->get_facet_L() == facet_1_dst && edge_2->get_facet_R() == facet_1_src))) {
          // polyhedron split event
          continue;
        }
        FacetSPtr facet_2_src = edge_2->get_facet_src();
        FacetSPtr facet_2_dst = edge_2->get_facet_dst();
        if ((edge_2->get_facet_L() == facet_1_src && edge_2->get_facet_R() != facet_1_dst) ||
            (edge_2->get_facet_L() == facet_1_dst && edge_2->get_facet_R() != facet_1_src) ||
            (edge_2->get_facet_R() == facet_1_src && edge_2->get_facet_L() != facet_1_dst) ||
            (edge_2->get_facet_R() == facet_1_dst && edge_2->get_facet_L() != facet_1_src)) {
          // surface event
          continue;
        }
        if ((edge_1->get_facet_L() == facet_2_src && edge_1->get_facet_R() != facet_2_dst) ||
            (edge_1->get_facet_L() == facet_2_dst && edge_1->get_facet_R() != facet_2_src) ||
            (edge_1->get_facet_R() == facet_2_src && edge_1->get_facet_L() != facet_2_dst) ||
            (edge_1->get_facet_R() == facet_2_dst && edge_1->get_facet_L() != facet_2_src)) {
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
          b2 += edge_2->get_vertex_src()->point().bbox();
          b2 += edge_2->get_vertex_dst()->point().bbox();
          b2 += Hds_utils::get_final_point(edge_2->get_vertex_src(), *time_future_bound).bbox();
          b2 += Hds_utils::get_final_point(edge_2->get_vertex_dst(), *time_future_bound).bbox();

          if (!CGAL::do_overlap(b1, b2)) {
            CGAL_SS3_CORE_TRACE_V(64, "Filtered edge split candidates\n\t" << edge_1->to_string() << "\n\t"
                                                                           << edge_2->to_string() << "\n\t"
                                                                           << b1 << "\n\t" << b2);
            continue;
          }
        }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
        ++tested_candidates;
#endif

        std::optional<FT> event_time = crash_time(edge_1, edge_2, current_time, time_future_bound);
        if (!event_time) {
          continue;
        }

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        // @fixme below needs to be double checked
        if (use_canonical_event_reps) {
          const FT shift = *event_time - current_time;
          const Segment_3 e1o = Transformation::shift_edge(edge_1, shift);
          if (e1o.is_degenerate()) {
            // polyhedron split
            continue;
          }
          const Segment_3 e2o = Transformation::shift_edge(edge_2, shift);
          if (e2o.is_degenerate()) {
            // polyhedron split
            continue;
          }
        }
#endif

        Edge_split_event_sptr event = Edge_split_event::create();
        event->set_time(*event_time);
        event->set_edge_1(edge_1);
        event->set_edge_2(edge_2);
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

  void collect_edge_split_events(const PolyhedronSPtr& polyhedron,
                                 const FT& current_time,
                                 const std::optional<FT>& time_future_bound,
                                 PQ& queue)
  {
    return collect_edge_split_events(polyhedron->edges(), polyhedron->edges(), polyhedron, true /*use canonical reps*/,
                                  current_time, time_future_bound, queue);
  }

  // The function below is meant to be used as the callback of the box_d spatial searching:
  // it does not perform any box-box filtering
  void collect_edge_split_event(const EdgeSPtr& edge_1,
                                const EdgeSPtr& edge_2,
                                const PolyhedronSPtr& /*polyhedron*/,
                                const bool use_canonical_event_reps,
                                const FT& current_time,
                                const std::optional<FT>& time_future_bound,
                                PQ& queue)
  {
    CGAL::Real_timer timer;
    timer.start();

    FacetSPtr facet_1_src = edge_1->get_facet_src();
    FacetSPtr facet_1_dst = edge_1->get_facet_dst();

    CGAL_SS3_DEBUG_SPTR(edge_1);
    CGAL_SS3_DEBUG_SPTR(edge_2);

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
    if (use_canonical_event_reps) {
      if (edge_1->get_ID() > edge_2->get_ID()) {
        return;
      }
    }
#else
    if (edge_1 == edge_2) {
      return;
    }
#endif

    if (edge_1->get_facet_L() == edge_2->get_facet_L() ||
        edge_1->get_facet_L() == edge_2->get_facet_R() ||
        edge_1->get_facet_R() == edge_2->get_facet_L() ||
        edge_1->get_facet_R() == edge_2->get_facet_R()) {
      // on same facet
      return;
    }
    if (edge_1->get_vertex_src()->point() == edge_2->get_vertex_src()->point() ||
        edge_1->get_vertex_src()->point() == edge_2->get_vertex_dst()->point() ||
        edge_1->get_vertex_dst()->point() == edge_2->get_vertex_src()->point() ||
        edge_1->get_vertex_dst()->point() == edge_2->get_vertex_dst()->point()) {
      // share a vertex
      return;
    }

    if (((edge_2->get_facet_L() == facet_1_src && edge_2->get_facet_R() == facet_1_dst) ||
         (edge_2->get_facet_L() == facet_1_dst && edge_2->get_facet_R() == facet_1_src))) {
        // polyhedron split event
        return;
    }
    FacetSPtr facet_2_src = edge_2->get_facet_src();
    FacetSPtr facet_2_dst = edge_2->get_facet_dst();
    if ((edge_2->get_facet_L() == facet_1_src && edge_2->get_facet_R() != facet_1_dst) ||
        (edge_2->get_facet_L() == facet_1_dst && edge_2->get_facet_R() != facet_1_src) ||
        (edge_2->get_facet_R() == facet_1_src && edge_2->get_facet_L() != facet_1_dst) ||
        (edge_2->get_facet_R() == facet_1_dst && edge_2->get_facet_L() != facet_1_src)) {
      // surface event
      return;
    }
    if ((edge_1->get_facet_L() == facet_2_src && edge_1->get_facet_R() != facet_2_dst) ||
        (edge_1->get_facet_L() == facet_2_dst && edge_1->get_facet_R() != facet_2_src) ||
        (edge_1->get_facet_R() == facet_2_src && edge_1->get_facet_L() != facet_2_dst) ||
        (edge_1->get_facet_R() == facet_2_dst && edge_1->get_facet_L() != facet_2_src)) {
      // surface event
        return;
    }

    std::optional<FT> event_time = crash_time(edge_1, edge_2, current_time, time_future_bound);
    if (!event_time) {
      return;
    }

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
    // @fixme below needs to be double checked
    if (use_canonical_event_reps) {
      FT shift = *event_time - current_time;
      const Segment_3 e1o = Transformation::shift_edge(edge_1, shift);
      if (e1o.is_degenerate()) {
        // polyhedron split
        return;
      }
      const Segment_3 e2o = Transformation::shift_edge(edge_2, shift);
      if (e2o.is_degenerate()) {
        // polyhedron split
        return;
      }
    }
#endif

    Edge_split_event_sptr event = Edge_split_event::create();
    event->set_time(*event_time);
    event->set_edge_1(edge_1);
    event->set_edge_2(edge_2);
    queue.push(event);
  }

  void collect_edge_split_events_with_boxd(const std::list<EdgeSPtr>& edges_1,
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
        if (!Hds_utils::is_reflex(edge)) {
          continue;
        }

        CGAL::Bbox_3 b = edge->get_vertex_src()->point().bbox();
        b += edge->get_vertex_dst()->point().bbox();
        b += Hds_utils::get_final_point(edge->get_vertex_src(), *time_future_bound).bbox();
        b += Hds_utils::get_final_point(edge->get_vertex_dst(), *time_future_bound).bbox();

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
      return collect_edge_split_event(edge_1, edge_2, polyhedron, use_canonical_event_reps,
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

  void collect_edge_split_events_with_boxd(const PolyhedronSPtr& polyhedron,
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
      if (!Hds_utils::is_reflex(edge)) {
        continue;
      }

      CGAL::Bbox_3 b = edge->get_vertex_src()->point().bbox();
      b += edge->get_vertex_dst()->point().bbox();
      b += Hds_utils::get_final_point(edge->get_vertex_src(), *time_future_bound).bbox();
      b += Hds_utils::get_final_point(edge->get_vertex_dst(), *time_future_bound).bbox();

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
      return collect_edge_split_event(edge_1, edge_2, polyhedron, false /*use_canonical_event_reps*/,
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
  void collect_pierce_events(const std::list<VertexSPtr>& vertices,
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
      CGAL_assertion(vertex->get_ID() != -1);

      if (Hds_utils::is_reflex(vertex)) {
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
        ++pierce_vertex_counter;
#endif

        for (const FacetSPtr& facet : facets) {
          // @todo contains_vertex is redundant with has_edge_to_facet?
          bool contains_vertex = false;
          for (const VertexSPtr& vertex_2 : facet->vertices()) {
            if (vertex_2->point() == vertex->point()) {
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
          // # `has_edge_to_facet` check is done at pop time - see is_actual_pierce_event()

          // # Shrinking, so the vertex must be on the backside of the plane
          if (Kernel_wrapper::side(facet->get_plane(), vertex->point()) > 0) {
            continue;
          }

          // # If the facet is so far that even when shifting point and plane by the maximal
          // time, the vertex would not cross it, then we are done
          if (time_future_bound.has_value()) {
            Point_3 shifted_pt = Hds_utils::get_final_point(vertex, *time_future_bound);
            Plane_3 shifted_plane = Hds_utils::get_final_plane(facet, *time_future_bound);

            if (Kernel_wrapper::side(shifted_plane, shifted_pt) < 0) {
              continue;
            }

            // # It would be tempting here to use a bbox around the vertex and facet's trajectories,
            // but if we did, then we would have to detect potential piercing events involving
            // a facet after each modification of the facet, and that seems costly
            // (even if potential piercing vertices were cached...)

            // # If both move in the same direction, we cannot pierce
            const Vector_3 d (vertex->point(), shifted_pt);
            const Vector_3 n = facet->get_plane().orthogonal_vector();
            if (n * d < 0) { // negative because the facet moves in a direction opposite of the normal
              continue;
            }
          }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
          ++tested_candidates;
#endif

          std::optional<FT> event_time = crash_time(vertex, facet, current_time, time_future_bound);
          if (!event_time) {
            continue;
          }

          CGAL_assertion(*event_time < current_time);
          CGAL_assertion(!time_future_bound.has_value() || *event_time >= *time_future_bound);

          // actual intersection checks are performed at pop time - see is_actual_pierce_event()

          Pierce_event_sptr event = Pierce_event::create();
          event->set_time(*event_time);
          event->set_facet(facet);
          event->set_vertex(vertex);
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

  void collect_pierce_events(const PolyhedronSPtr& polyhedron,
                             const FT& current_time,
                             const std::optional<FT>& time_future_bound,
                             PQ& queue)
  {
    return collect_pierce_events(polyhedron->vertices(), polyhedron->facets(), polyhedron,
                               current_time, time_future_bound, queue);
  }

  /**
    * Collect all events for a subset of the polyhedron
    */
  void collect_local_events(const PolyhedronSPtr& polyhedron,
                            const FT& current_time,
                            const std::optional<FT>& time_future_bound,
                            PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(2, "collect_local_events(" << current_time << ")");

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
      CGAL_SS3_CORE_TRACE_V(8, "\t" << e->to_string());

      collect_vanish_events(local_edges, polyhedron, current_time, time_future_bound, queue);
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
      CGAL_SS3_CORE_TRACE_CODE(ss << " " << v->get_ID();)
      CGAL_SS3_CORE_TRACE_V(8, ss.str());

      const bool use_canonical_reps = false;
#else
      const std::list<VertexSPtr>& local_vertices_VV = polyhedron->vertices();
      const bool use_canonical_reps = true;
#endif

      collect_vertex_events(local_vertices_VV, polyhedron, use_canonical_reps,
                          current_time, time_future_bound, queue);
      collect_flip_vertex_events(local_vertices_VV, polyhedron, use_canonical_reps,
                              current_time, time_future_bound, queue);
      collect_split_merge_events(local_vertices_VV, polyhedron, use_canonical_reps,
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
      CGAL_SS3_CORE_TRACE_V(8, "\t" << e->to_string());

      // this is the modified edges as 'edge_1'
      collect_polyhedron_split_events(local_edges_EE, polyhedron,
                                   current_time, time_future_bound, queue);

      // this is the modified edges as 'edge_2'
      // since we know that for a polyhedron split event, edge_2's LR facets
      // are the src and dst of edge_1, we can build a subset of all edges for edge_1s
      for (const EdgeSPtr& edge_2 : local_edges_EE) {
        std::set<EdgeSPtr> edges_1;
        for (const FacetSPtr& facet_2 : {edge_2->get_facet_L(), edge_2->get_facet_R()}) {
          for (const VertexSPtr& vertex_1 : facet_2->vertices()) {
            // keep the edge incident to 'vertex_1' which is not incident to 'facet_2'
            EdgeSPtr edge_1;
            for (EdgeWPtr edge_wptr : vertex_1->edges()) {
              if (EdgeSPtr edge = edge_wptr.lock()) {
                if (edge->get_facet_L() != facet_2 && edge->get_facet_R() != facet_2) {
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
          collect_polyhedron_split_event(edge_1, edge_2, polyhedron,
                                      current_time, time_future_bound, queue);
        }
      }


#else
      collect_polyhedron_split_events(polyhedron->edges(), polyhedron,
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
      CGAL_SS3_CORE_TRACE_CODE(ss << " " << v->get_ID());
      CGAL_SS3_CORE_TRACE_V(8, ss.str());
#else
      std::list<VertexSPtr> local_vertices_VF = polyhedron->vertices();
#endif

      // Pierce events must check with all faces, no choice about this
      collect_pierce_events(local_vertices_VF, polyhedron->facets(), polyhedron,
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
      CGAL_SS3_CORE_TRACE_V(8, "\t" << e->to_string());

      // this is the modified edges as 'edge_1'
      collect_surface_events(local_edges_EE, polyhedron,
                            current_time, time_future_bound, queue);

      // this is the modified edges as 'edge_2'
      for (const EdgeSPtr& edge_2 : local_edges_EE) {
        std::set<EdgeSPtr> edges_1;
        for (const FacetSPtr& facet_2 : {edge_2->get_facet_L(), edge_2->get_facet_R()}) {
          for (const VertexSPtr& vertex_1 : facet_2->vertices()) {
            // keep the edge incident to 'vertex_1' which is not incident to 'facet_2'
            EdgeSPtr edge_1;
            for (EdgeWPtr edge_wptr : vertex_1->edges()) {
              if (EdgeSPtr edge = edge_wptr.lock()) {
                if (edge->get_facet_L() != facet_2 && edge->get_facet_R() != facet_2) {
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
          collect_surface_event(edge_1, edge_2, polyhedron,
                              current_time, time_future_bound, queue);
        }
      }
#else
      collect_surface_events(polyhedron->edges(), polyhedron,
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
      CGAL_SS3_CORE_TRACE_V(8, "\t" << e->to_string());

      const bool use_canonical_reps = false;
#else
      std::list<EdgeSPtr> local_edges_EE = polyhedron->edges();
      const bool use_canonical_reps = true;
#endif

      // not worth the effort? 'local_edges_EE' is always very small
#ifdef CGAL_SS3_DETECT_EDGE_SPLIT_EVENTS_WITH_BOX_D
      if (time_future_bound.has_value()) {
        collect_edge_split_events_with_boxd(local_edges_EE, polyhedron->edges(), polyhedron,
                                       false /*use_canonical_reps*/, // box_d returns a canonical pair
                                       current_time, time_future_bound, queue);
      } else
#endif
      {
        collect_edge_split_events(local_edges_EE, polyhedron->edges(), polyhedron, use_canonical_reps,
                                current_time, time_future_bound, queue);
      }
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(2, "Sought All Local Events in: " << timer.time());
#endif

#ifdef CGAL_SS3_DEBUG_PRINT_QUEUE
    print_queue(queue);
#endif

    // CGAL_postcondition(check_queue_correctness(queue, polyhedron, current_time, time_future_bound));
  }

  /**
    * Collect all events for the polyhedron
    */
  void collect_events(const PolyhedronSPtr& polyhedron,
                      const FT& current_time,
                      const std::optional<FT>& time_future_bound,
                      PQ& queue)
  {
    CGAL_SS3_CORE_TRACE_V(2, "Collecting events in [" << current_time
                          << " | " << (time_future_bound.has_value() ? IO::String_factory::fromDouble(CGAL::to_double(*time_future_bound)) : "") << "]");

    Abstract_event_sptr result = Abstract_event_sptr();
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
    collect_vanish_events(polyhedron, current_time, time_future_bound, queue);

    CGAL_assertion_code(for (const EdgeSPtr& edge : polyhedron->edges()))
    CGAL_assertion_code(Hds_utils::get_vanish_time(edge);) // check is_known within skeledgedata

    // --- Contact Event
    collect_vertex_events(polyhedron, current_time, time_future_bound, queue);
    collect_flip_vertex_events(polyhedron, current_time, time_future_bound, queue);
    collect_polyhedron_split_events(polyhedron, current_time, time_future_bound, queue);
    collect_split_merge_events(polyhedron, current_time, time_future_bound, queue);

    // the next event types are particularly slow, so reduce the bound by doing them last
    // so other events lower the bound
    collect_pierce_events(polyhedron, current_time, time_future_bound, queue);
    collect_surface_events(polyhedron, current_time, time_future_bound, queue);

#ifdef CGAL_SS3_DETECT_EDGE_SPLIT_EVENTS_WITH_BOX_D
    if (time_future_bound.has_value()) {
      collect_edge_split_events_with_boxd(polyhedron, current_time, time_future_bound, queue);
    } else
#endif
    {
      collect_edge_split_events(polyhedron, current_time, time_future_bound, queue);
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(2, "Sought All Events in: " << timer.time());
#endif

#ifdef CGAL_SS3_DEBUG_PRINT_QUEUE
    print_queue(queue);
#endif
  }

  void print_queue(const PQ& queue)
  {
    CGAL_SS3_CORE_TRACE("-------------------------------------------------");
    CGAL_SS3_CORE_TRACE("--- Event queue (size = " << queue.size() << "; iter = " << step_id_ << ") ---");
    CGAL_SS3_CORE_TRACE("-------------------------------------------------");

    PQ duplicate_queue = queue;
    while (!duplicate_queue.empty()) {
      Abstract_event_sptr event = duplicate_queue.top();
      CGAL_SS3_CORE_TRACE("Event E" << event->get_ID()
                            << " T" << event->getType()
                           << " @ " << event->time());
      if (event->is_valid()) {
        if (is_event_obsolete(event)) {
          CGAL_SS3_CORE_TRACE("  Event is obsolete");
        } else {
          CGAL_SS3_CORE_TRACE(event->to_string());
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
  bool check_queue_correctness(const PQ& queue,
                               const PolyhedronSPtr& polyhedron,
                               const FT& current_time,
                               const std::optional<FT>& time_future_bound)
  {
    CGAL_SS3_CORE_TRACE("Checking queue correctness...");
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    // Compute a queue from scratch using collect_events()
    PQ queue_from_scratch;
    collect_events(polyhedron, current_time, time_future_bound, queue_from_scratch);

    // Duplicate the queue since we need to pop events
    PQ duplicate_queue = queue;

    // Tops should be the same
    auto purge_top = [&](PQ& q)
    {
      while (!q.empty()) {
        Abstract_event_sptr event = q.top();
        if (!event->is_valid() ||
            is_event_in_the_past(event, current_time) ||
            is_event_obsolete(event) ||
            !is_actual_event(event, current_time, time_future_bound)) {
          q.pop();
        } else {
          break;
        }
      }
    };

    auto is_same_event = [](const Abstract_event_sptr& event_1,
                            const Abstract_event_sptr& event_2)
    {
      if (event_1->getType() != event_2->getType()) {
        return false;
      }

      switch (event_1->getType()) {
        case Abstract_event::SAVE_EVENT: {
          auto save_event_1 = std::dynamic_pointer_cast<Save_event>(event_1);
          auto save_event_2 = std::dynamic_pointer_cast<Save_event>(event_2);
          return *save_event_1 == *save_event_2;
        }
        case Abstract_event::CONST_TIME_EVENT: {
          auto const_event_1 = std::dynamic_pointer_cast<Const_time_event>(event_1);
          auto const_event_2 = std::dynamic_pointer_cast<Const_time_event>(event_2);
          return *const_event_1 == *const_event_2;
        }
        case Abstract_event::VANISH_EVENT: {
          auto vanish_event_1 = std::dynamic_pointer_cast<Vanish_event>(event_1);
          auto vanish_event_2 = std::dynamic_pointer_cast<Vanish_event>(event_2);
          return *vanish_event_1 == *vanish_event_2;
        }
        case Abstract_event::EDGE_EVENT: {
          auto edge_event_1 = std::dynamic_pointer_cast<Edge_event>(event_1);
          auto edge_event_2 = std::dynamic_pointer_cast<Edge_event>(event_2);
          return *edge_event_1 == *edge_event_2;
        }
        case Abstract_event::EDGE_MERGE_EVENT: {
          auto edge_merge_event_1 = std::dynamic_pointer_cast<Edge_merge_event>(event_1);
          auto edge_merge_event_2 = std::dynamic_pointer_cast<Edge_merge_event>(event_2);
          return *edge_merge_event_1 == *edge_merge_event_2;
        }
        case Abstract_event::TRIANGLE_EVENT: {
          auto triangle_event_1 = std::dynamic_pointer_cast<Triangle_event>(event_1);
          auto triangle_event_2 = std::dynamic_pointer_cast<Triangle_event>(event_2);
          return *triangle_event_1 == *triangle_event_2;
        }
        case Abstract_event::DBL_EDGE_MERGE_EVENT: {
          auto dbl_edge_merge_event_1 = std::dynamic_pointer_cast<Dbl_edge_merge_event>(event_1);
          auto dbl_edge_merge_event_2 = std::dynamic_pointer_cast<Dbl_edge_merge_event>(event_2);
          return *dbl_edge_merge_event_1 == *dbl_edge_merge_event_2;
        }
        case Abstract_event::DBL_TRIANGLE_EVENT: {
          auto dbl_triangle_event_1 = std::dynamic_pointer_cast<Dbl_triangle_event>(event_1);
          auto dbl_triangle_event_2 = std::dynamic_pointer_cast<Dbl_triangle_event>(event_2);
          return *dbl_triangle_event_1 == *dbl_triangle_event_2;
        }
        case Abstract_event::TETRAHEDRON_EVENT: {
          auto tetrahedron_event_1 = std::dynamic_pointer_cast<Tetrahedron_event>(event_1);
          auto tetrahedron_event_2 = std::dynamic_pointer_cast<Tetrahedron_event>(event_2);
          return *tetrahedron_event_1 == *tetrahedron_event_2;
        }
        case Abstract_event::VERTEX_EVENT: {
          auto vertex_event_1 = std::dynamic_pointer_cast<Vertex_event>(event_1);
          auto vertex_event_2 = std::dynamic_pointer_cast<Vertex_event>(event_2);
          return *vertex_event_1 == *vertex_event_2;
        }
        case Abstract_event::FLIP_VERTEX_EVENT: {
          auto flip_vertex_event_1 = std::dynamic_pointer_cast<Flip_vertex_event>(event_1);
          auto flip_vertex_event_2 = std::dynamic_pointer_cast<Flip_vertex_event>(event_2);
          return *flip_vertex_event_1 == *flip_vertex_event_2;
        }
        case Abstract_event::SURFACE_EVENT: {
          auto surface_event_1 = std::dynamic_pointer_cast<Surface_event>(event_1);
          auto surface_event_2 = std::dynamic_pointer_cast<Surface_event>(event_2);
          return *surface_event_1 == *surface_event_2;
        }
        case Abstract_event::POLYHEDRON_SPLIT_EVENT: {
          auto polyhedron_split_event_1 = std::dynamic_pointer_cast<Polyhedron_split_event>(event_1);
          auto polyhedron_split_event_2 = std::dynamic_pointer_cast<Polyhedron_split_event>(event_2);
          return *polyhedron_split_event_1 == *polyhedron_split_event_2;
        }
        case Abstract_event::EDGE_SPLIT_EVENT: {
          auto edge_event_1 = std::dynamic_pointer_cast<Edge_split_event>(event_1);
          auto edge_event_2 = std::dynamic_pointer_cast<Edge_split_event>(event_2);
          return *edge_event_1 == *edge_event_2;
        }
        case Abstract_event::SPLIT_MERGE_EVENT: {
          auto split_merge_event_1 = std::dynamic_pointer_cast<Split_merge_event>(event_1);
          auto split_merge_event_2 = std::dynamic_pointer_cast<Split_merge_event>(event_2);
          return *split_merge_event_1 == *split_merge_event_2;
        }
        case Abstract_event::PIERCE_EVENT: {
          auto pierce_event_1 = std::dynamic_pointer_cast<Pierce_event>(event_1);
          auto pierce_event_2 = std::dynamic_pointer_cast<Pierce_event>(event_2);
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
      Abstract_event_sptr event_scratch = queue_from_scratch.top();
      Abstract_event_sptr event_duplicate = duplicate_queue.top();
      CGAL_SS3_CORE_TRACE("First event from scratch: " << event_scratch->to_string());
      CGAL_SS3_CORE_TRACE("First event from duplicate: " << event_duplicate->to_string());

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

      Abstract_event_sptr event_scratch = queue_from_scratch.top();
      queue_from_scratch.pop();

      CGAL_SS3_CORE_TRACE("Seek event @ " << event_scratch->time() << " Type " << event_scratch->getType());

      CGAL_SS3_CORE_TRACE("Event E" << event_scratch->get_ID()
                          << " T" << event_scratch->getType()
                          << " @ " << event_scratch->time());
      CGAL_assertion(event_scratch->is_valid() && !is_event_obsolete(event_scratch));
      CGAL_SS3_CORE_TRACE(event_scratch->to_string());

      // Find the event in the duplicate queue
      bool found = false;
      duplicate_queue = queue;
      while (!duplicate_queue.empty()) {
        Abstract_event_sptr event_duplicate = duplicate_queue.top();
        if (event_duplicate->is_valid()) {
          if (is_same_event(event_scratch, event_duplicate)) {
            found = true;
          }

          if (found) {
            if (is_event_obsolete(event_duplicate)) {
                CGAL_SS3_CORE_TRACE("Warning: Found event in duplicate queue but it's marked as obsolete");
                CGAL_SS3_CORE_TRACE("Event: " << event_duplicate->to_string());
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
          CGAL_SS3_CORE_TRACE("Event: " << event_scratch->to_string());
          return false;
      }
    }

    CGAL_SS3_CORE_TRACE("OK: found all 'from-scratch' events!");

    return true;
  }

  /**
    * Determines the next event.
    */
  Abstract_event_sptr nextEvent(PQ& queue,
                                const PolyhedronSPtr& polyhedron,
                                const FT& current_time)
  {
    // We could "break" directly in run(), but like this, the save events
    // that are farther than the last event are still processed (arbitrarily decision).
    if (polyhedron->empty()) {
      queue = { };
      // do not return here as to allow save events to exist after the last 'real' event
    }

    Abstract_event_sptr event;
    FT time = 0;

    // If a save event is closest, delay building it in case a const event is even closer
    bool is_save_event = false;

    while (!queue.empty() || !save_times_.empty()) {
      // purge the queue lazily as to avoid wasting time if we stop on the last save event
      if (save_times_.empty()) {
        event = queue.top();
        time = event->time();
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
          if (next_save_time > queue.top()->time()) { // save is strictly earlier
            is_save_event = true;
            time = next_save_time;
          } else {
            // save times exist, but are farther in the future than the next event
            event = queue.top();
            time = event->time();
          }
        }
      }

      // Tentative next event is not a save event, test its sanity.
      // Do this here because we don't want a const event to get created
      // because it's before a zombie event.
      if (!is_save_event) {
        if (!event->is_valid()) {
          CGAL_SS3_CORE_TRACE_V(16, "Skipping invalid event E" << event->get_ID());
          event = {};
          queue.pop();
          continue;
        }

        if (is_event_in_the_past(event, current_time)) {
          CGAL_SS3_CORE_TRACE_V(16, "Skipping event-in-the-past E" << event->get_ID());
          event = {};
          queue.pop();
          continue;
        }

        // This "is_obsolete()" function is only a sufficient condition: if the neighborhoods
        // have changed, then the event should be discarded.
        // "is_actual_..._event()" takes care of the other checks.
        if (is_event_obsolete(event)) {
          CGAL_SS3_CORE_TRACE_V(16, "Skipping obsolete event E" << event->get_ID());
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
    FT pulse = Configuration::get_instance()->get_double("Algorithm", "const_offset");
    if (pulse != 0) {
      FT next_const_time = floor(CGAL::to_double(current_time / pulse) + 1.0) * pulse;
      if (current_time > next_const_time && next_const_time > time) {
        is_save_event = false;
        Const_time_event_sptr const_event = Const_time_event::create();
        const_event->set_time(next_const_time);
        return const_event;
      }
    }

    if (is_save_event) {
      save_times_.erase(save_times_.begin());
      Save_event_sptr save_event = Save_event::create();
      save_event->set_time(time);
      return save_event;
    }

    CGAL_SS3_DEBUG_SPTR(event);

    // if here, the topmost is neither a const nor save event
    queue.pop();

    return event;
  }

  void add_event(const Abstract_event_sptr& event)
  {
    typename std::list<Abstract_event_sptr>::iterator it = events_.insert(events_.end(), event);
    event->setListIt(it);
  }

  bool remove_event(const Abstract_event_sptr& event)
  {
    bool result = false;
    events_.erase(event->getListIt());
    event->setListIt(typename std::list<Abstract_event_sptr>::iterator());
    return result;
  }

  void shit_to_event_time(const PolyhedronSPtr& polyhedron,
                          const FT& current_time,
                          const FT& target_time)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    const FT shift = target_time - current_time;
    CGAL_precondition(!is_zero(shift));

    Transformation::shift_facets(polyhedron, shift);

#ifdef CGAL_SS3_DUMP_FILES
    // below will have degeneracies since we have not yet treated the event
    static int shift_id = -1;
    IO::OBJFile::save("results/shift_" + std::to_string(++shift_id) + ".obj",
                      polyhedron,
                      false /*do not triangulate*/);
#endif
  }

  bool save_polyhedron(PolyhedronSPtr polyhedron,
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

  bool save_skeleton(const StraightSkeletonSPtr& skeleton,
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

  Event_status handle_save_event(const Save_event_sptr& event,
                                 const FT& current_time,
                                 const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#########  Handle Save Event  ##########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->time();
    shit_to_event_time(polyhedron, current_time, event_time);

    bool res = true;

#ifdef CGAL_SS3_DUMP_FILES
    res = save_polyhedron(polyhedron, event_time);
    // res = save_skeleton(skeleton_, event_time) && res;
#endif

    if (res) {
      add_event(event);
    }

    return (res ? Event_status::EVENT_HANDLED : Event_status::EVENT_NOT_HANDLED);
  }

  // This 'handle' is in fact more akin to a collect, but the interesting point
  // is that it happens after pop time
  Event_status handle_const_time_event(Const_time_event_sptr event,
                                       const FT& current_time,
                                       const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#########  Handle Const Event  #########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    const FT& event_time = event->time();
    shit_to_event_time(polyhedron, current_time, event_time);

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  // This 'handle' is in fact more akin to a collect, but the interesting point
  // is that it happens after pop time.
  // This function does NOT shift the polyhedron, it is for the "real" handler
  // to deal with it.
  // This function might not do anything, for example if the vanish event is in fact
  // escalated as a contact event.
  Event_status handle_vanish_event(const Vanish_event_sptr& event,
                                   const FT& current_time,
                                   const std::optional<FT>& time_future_bound,
                                   const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "####### Tentative Vanish Event  ########");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    CGAL_SS3_DEBUG_SPTR(event);

    EdgeSPtr edge = event->get_edge();
    const FT& event_time = event->time();

    event->set_point(vanish_point(edge));
    const Point_3& point = event->point();

    // @todo would nice:
    // - to avoid redundant checks (e.g. is_triangle multiple times)
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
      FacetSPtr facet_l = edge->get_facet_L();
      FacetSPtr facet_r = edge->get_facet_R();
      if (Hds_utils::is_triangle(facet_l, edge) ||
          Hds_utils::is_triangle(facet_r, edge)) {
        // triangle event
        CGAL_SS3_CORE_TRACE_V(8, "Not an Edge Event (one incident facet is a triangle)");
        break;
      }

      FacetSPtr facet_src = edge->get_facet_src();
      FacetSPtr facet_dst = edge->get_facet_dst();

      // This does not work when there is more than one edge between both facets.
      // EdgeSPtr edge_2 = facet_src->find_edge(facet_dst);
      std::list<EdgeSPtr> edges_2 = facet_src->find_edges(facet_dst); // @todo shouldn't this check also happen in other events?...

      bool split_event = false;
      for (const EdgeSPtr& edge_2 : edges_2) {
        FacetSPtr facet_l2 = edge_2->get_facet_L();
        FacetSPtr facet_r2 = edge_2->get_facet_R();
        FacetSPtr facet_2_src = edge_2->get_facet_src();
        FacetSPtr facet_2_dst = edge_2->get_facet_dst();

        const Plane_3& plane_l2 = facet_l2->get_plane();
        CGAL_assertion_code(const FT& l2a = plane_l2.a();)
        CGAL_assertion_code(const FT& l2b = plane_l2.b();)
        CGAL_assertion_code(const FT& l2c = plane_l2.c();)
        CGAL_assertion_code(const FT& l2d = plane_l2.d();)
        CGAL_assertion_code(const FT& speed_l2 = Hds_utils::get_speed(facet_l2);)
        CGAL_assertion_code(FT lt2 = (l2a * point.x() + l2b * point.y() + l2c * point.z() + l2d) / speed_l2);

        CGAL_assertion_code(const Plane_3& plane_r2 = facet_r2->get_plane();)
        CGAL_assertion_code(const FT& r2a = plane_r2.a();)
        CGAL_assertion_code(const FT& r2b = plane_r2.b();)
        CGAL_assertion_code(const FT& r2c = plane_r2.c();)
        CGAL_assertion_code(const FT& r2d = plane_r2.d();)
        CGAL_assertion_code(const FT& speed_r2 = Hds_utils::get_speed(facet_r2);)
        CGAL_assertion_code(FT rt2 = (r2a * point.x() + r2b * point.y() + r2c * point.z() + r2d) / speed_r2);

        CGAL_assertion(lt2 == rt2);
        CGAL_assertion(!is_positive(lt2));

        if (!check_bisectors_V2(edge_2, point, event_time)) {
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
      if (edge_prev->has_same_facets(edge_next)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not an Edge Event (Edge Merge #1)");
        break;
      }
      edge_prev = edge->prev(facet_l)->prev(facet_l);
      edge_next = edge->next(facet_l);
      if (edge_prev->has_same_facets(edge_next)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not an Edge Event (Edge Merge #2)");
        break;
      }
      edge_prev = edge->prev(facet_r);
      edge_next = edge->next(facet_r)->next(facet_r);
      if (edge_prev->has_same_facets(edge_next)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not an Edge Event (Edge Merge #3)");
        break;
      }
      edge_prev = edge->prev(facet_r)->prev(facet_r);
      edge_next = edge->next(facet_r);
      if (edge_prev->has_same_facets(edge_next)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not an Edge Event (Edge Merge #4)");
        break;
      }

      // if here, it's an edge event
      Edge_event_sptr edge_event = Edge_event::create();
      edge_event->set_time(event_time);
      edge_event->set_point(point);
      edge_event->set_edge(edge);

      return handle_edge_event(edge_event, current_time, time_future_bound, polyhedron);
    }

    // Edge_merge_event
    for (;;)
    {
      FacetSPtr facet_l = edge->get_facet_L();
      FacetSPtr facet_r = edge->get_facet_R();
      if (Hds_utils::is_triangle(facet_l, edge) ||
          Hds_utils::is_triangle(facet_r, edge)) {
        // triangle event
        CGAL_SS3_CORE_TRACE_V(8, "Not an EdgeMerge Event (one incident facet is a triangle)");
        break;
      }

      FacetSPtr facet_other = edge->get_facet_L();
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

      facet_other = edge->get_facet_R();
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
      if (edge_prev->has_same_facets(edge_next) && edge_prev != edge_next) {
        facet = facet_l;
        edge_1 = edge_prev;
        edge_2 = edge_next;
      }

      edge_prev = edge->prev(facet_l)->prev(facet_l);
      edge_next = edge->next(facet_l);
      if (edge_prev->has_same_facets(edge_next) && edge_prev != edge_next) {
        facet = facet_l;
        edge_1 = edge_prev;
        edge_2 = edge_next;
      }

      edge_prev = edge->prev(facet_r);
      edge_next = edge->next(facet_r)->next(facet_r);
      if (edge_prev->has_same_facets(edge_next) && edge_prev != edge_next) {
        facet = facet_r;
        edge_1 = edge_prev;
        edge_2 = edge_next;
      }

      edge_prev = edge->prev(facet_r)->prev(facet_r);
      edge_next = edge->next(facet_r);
      if (edge_prev->has_same_facets(edge_next) && edge_prev != edge_next) {
        facet = facet_r;
        edge_1 = edge_prev;
        edge_2 = edge_next;
      }

      if (!(facet && edge_1 && edge_2)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not an EdgeMerge Event (Neighborhood)");
        break;
      }

      // if here, it's an edge merge event
      Edge_merge_event_sptr edge_merge_event = Edge_merge_event::create();
      edge_merge_event->set_time(event_time);
      edge_merge_event->set_point(point);
      edge_merge_event->set_facet(facet);
      edge_merge_event->set_edge_1(edge_1);
      edge_merge_event->set_edge_2(edge_2);

      return handle_edge_merge_event(edge_merge_event, current_time, time_future_bound, polyhedron);
    }

    // Triangle_event
    for (;;)
    {
      if (Hds_utils::is_tetrahedron(edge)) {
        // tetrahedron event
        CGAL_SS3_CORE_TRACE_V(8, "Not a Triangle Event (Tetrahedron)");
        break;
      }

      FacetSPtr facet;
      if (Hds_utils::is_triangle(edge->get_facet_L(), edge)) {
        facet = edge->get_facet_L();
      } else if (Hds_utils::is_triangle(edge->get_facet_R(), edge)) {
        facet = edge->get_facet_R();
      } else {
        CGAL_SS3_CORE_TRACE_V(8, "Not a Triangle Event (not triangle)");
        break;
      }

      bool dbl_triangle_event = false;
      EdgeSPtr edge_tmp = edge;
      for (unsigned int i = 0; i < 3; ++i) {
        FacetSPtr facet_tmp_l = edge_tmp->get_facet_L();
        FacetSPtr facet_tmp_r = edge_tmp->get_facet_R();
        if (facet_tmp_l && facet_tmp_r) {
          if (Hds_utils::is_triangle(facet_tmp_l, edge_tmp) &&
              Hds_utils::is_triangle(facet_tmp_r, edge_tmp)) {
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
      Triangle_event_sptr triangle_event = Triangle_event::create();
      triangle_event->set_time(event_time);
      triangle_event->set_point(point);
      triangle_event->set_facet(facet);
      triangle_event->set_edge_begin(edge);

      return handle_triangle_event(triangle_event, current_time, time_future_bound, polyhedron);
    }

    // Dbl_edge_merge_event
    for (;;)
    {
      // At pop time, edge is degenerate, but by default is_reflex() does a symbolic
      // shift into the future. However, here we want to know if the edge was reflex
      // BEFORE it became degenerate.
      if (!Hds_utils::is_reflex(edge, false /*shift into the past*/)) {
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
      FacetSPtr facet_other = edge->get_facet_L();
      EdgeSPtr edge_next = edge->next(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->prev(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->next(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->prev(facet_other);
      if (edge_next == edge) {
        is_dbl_edge_merge_event = true;
        facet_1 = edge->get_facet_L();
        edge_11 = edge->prev(facet_1);
        edge_12 = edge->next(facet_1)->next(facet_1);
        facet_2 = edge->get_facet_R();
        edge_21 = edge->prev(facet_2);
        edge_22 = edge->next(facet_2)->next(facet_2);
      }

      facet_other = edge->get_facet_R();
      edge_next = edge->prev(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->next(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->prev(facet_other);
      facet_other = edge_next->other(facet_other);
      edge_next = edge_next->next(facet_other);
      if (edge_next == edge) {
        is_dbl_edge_merge_event = true;
        facet_1 = edge->get_facet_R();
        edge_11 = edge->prev(facet_1)->prev(facet_1);
        edge_12 = edge->next(facet_1);
        facet_2 = edge->get_facet_L();
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
      Dbl_edge_merge_event_sptr dbl_edge_merge_event = Dbl_edge_merge_event::create();
      dbl_edge_merge_event->set_time(event_time);
      dbl_edge_merge_event->set_point(point);
      dbl_edge_merge_event->set_facet_1(facet_1);
      dbl_edge_merge_event->set_edge_11(edge_11);
      dbl_edge_merge_event->set_edge_12(edge_12);
      dbl_edge_merge_event->set_facet_2(facet_2);
      dbl_edge_merge_event->set_edge_21(edge_21);
      dbl_edge_merge_event->set_edge_22(edge_22);

      return handle_dbl_edge_merge_event(dbl_edge_merge_event, current_time, time_future_bound, polyhedron);
    }

    // Dbl_triangle_event
    for (;;)
    {
      if (Hds_utils::is_tetrahedron(edge)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not a DblTriangleMerge Event (Tetrahedron)");
        break;
      }
      FacetSPtr facet_l = edge->get_facet_L();
      FacetSPtr facet_r = edge->get_facet_R();
      if (!facet_l || !facet_r) {
        CGAL_SS3_CORE_TRACE_V(8, "Not a DblTriangleMerge Event (neighborhood)");
        break;
      }
      if (!(Hds_utils::is_triangle(facet_l, edge) && Hds_utils::is_triangle(facet_r, edge))) {
        CGAL_SS3_CORE_TRACE_V(8, "Not a DblTriangleMerge Event (not triangles)");
        break;
      }

      // if here, it's a double triangle event
      Dbl_triangle_event_sptr dbl_triangle_event = Dbl_triangle_event::create();
      dbl_triangle_event->set_time(event_time);
      dbl_triangle_event->set_point(point);
      dbl_triangle_event->set_edge(edge);

      return handle_dbl_triangle_event(dbl_triangle_event, current_time, time_future_bound, polyhedron);
    }

    // Tetrahedron_event
    for (;;)
    {
      if (!Hds_utils::is_tetrahedron(edge)) {
        CGAL_SS3_CORE_TRACE_V(8, "Not a Tetrahedron Event");
        break;
      }

      // if here, it's a tetrahedron event
      Tetrahedron_event_sptr tetrahedron_event = Tetrahedron_event::create();
      tetrahedron_event->set_time(event_time);
      tetrahedron_event->set_point(point);
      tetrahedron_event->set_edge_begin(edge);

      return handle_tetrahedron_event(tetrahedron_event, current_time, time_future_bound, polyhedron);
    }

    return Event_status::NON_EVENT;
  }

  Event_status handle_edge_event(const Edge_event_sptr& event,
                                 const FT& current_time,
                                 const std::optional<FT>& time_future_bound,
                                 const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#########  Handle Edge Event  ##########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->time();

    shit_to_event_time(polyhedron, current_time, event_time);

    EdgeSPtr edge = event->get_edge();
    VertexSPtr vertex_src = edge->get_vertex_src();
    VertexSPtr vertex_dst = edge->get_vertex_dst();

    CGAL_SS3_CORE_TRACE_V(4,"Edge:\n" << edge->to_string());

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(event->point());
    skeleton_->add_node(node);

    Hds_utils::get_arc(vertex_src)->close_arc(node);
    Hds_utils::get_arc(vertex_dst)->close_arc(node);

    Hds_utils::get_sheet(edge)->add_node(node);
    Hds_utils::clear_sheet(edge);
#endif

    // grab vertices and facets counter clockwise around node
    std::array<VertexSPtr, 4> vertices;
    vertices[0] = vertex_src->prev(edge->get_facet_L());
    vertices[1] = vertex_src->next(edge->get_facet_R());
    vertices[2] = vertex_dst->prev(edge->get_facet_R());
    vertices[3] = vertex_dst->next(edge->get_facet_L());

    std::array<FacetSPtr, 4> facets;
    facets[1] = edge->get_facet_R();
    facets[3] = edge->get_facet_L();
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
        facets_clone[i]->set_plane(facets[i]->get_plane());
        CGAL_assertion(facets[i]->has_data());
        facets_clone[i]->set_data(facets[i]->get_data()->clone(facets[i]));
      }

      edge_no_flip->set_facet_R(facets_clone[1]);
      edge_no_flip->set_facet_L(facets_clone[3]);

      // In edge events, we have unbounded faces (degree 1 vertices)
      facets_clone[3]->add_edge(edge_no_flip);
      facets_clone[1]->add_edge(edge_no_flip);
      for (unsigned int i = 0; i < 4; ++i) {
        edges[i]->set_facet_R(facets_clone[(i+3)%4]);
        edges[i]->set_facet_L(facets_clone[i]);
        facets_clone[(i+3)%4]->add_edge(edges[i]);
        facets_clone[i]->add_edge(edges[i]);
      }

      PolyhedronSPtr polyhedron_no_flip = Polyhedron::create(facets_clone);
      not_flipped_valid =
          (Transformation::shift_facets_deg1(polyhedron_no_flip, -1.0) &&
           !Self_intersection::has_self_intersecting_surface(polyhedron_no_flip));
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
        facets_clone[i]->set_plane(facets[i]->get_plane());
        CGAL_assertion(facets[i]->has_data());
        facets_clone[i]->set_data(facets[i]->get_data()->clone(facets[i]));
      }

      edge_flipped->set_facet_R(facets_clone[2]);
      edge_flipped->set_facet_L(facets_clone[0]);

      facets_clone[0]->add_edge(edge_flipped);
      facets_clone[2]->add_edge(edge_flipped);
      for (unsigned int i = 0; i < 4; ++i) {
        edges[i]->set_facet_R(facets_clone[(i+3)%4]);
        edges[i]->set_facet_L(facets_clone[i]);
        facets_clone[(i+3)%4]->add_edge(edges[i]);
        facets_clone[i]->add_edge(edges[i]);
      }

      PolyhedronSPtr polyhedron_flipped = Polyhedron::create(facets_clone);
      flipped_valid =
          (Transformation::shift_facets_deg1(polyhedron_flipped, -1.0) &&
           !Self_intersection::has_self_intersecting_surface(polyhedron_flipped));
    }

    CGAL_SS3_CORE_TRACE_V(16, "flipped_valid = " << flipped_valid);

    if (flipped_valid && !not_flipped_valid) {
      flip_edge = true;
    } else if (not_flipped_valid && !flipped_valid) {
      flip_edge = false;
    } else if (flipped_valid && not_flipped_valid) {
      const Vector_3 n0 = facets[0]->get_plane().orthogonal_vector();
      const Vector_3 n2 = facets[2]->get_plane().orthogonal_vector();
      const Vector_3 n1 = facets[1]->get_plane().orthogonal_vector();
      const Vector_3 n3 = facets[3]->get_plane().orthogonal_vector();
      CGAL::Comparison_result dac = CGAL::compare_angle(n0, n2, n1, n3);

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
      throw std::runtime_error("Error: not able to handle Edge_event (2).");
    }

    std::array<EdgeSPtr, 4> edges;
    edges[0] = edge->next(vertex_src);
    edges[1] = edges[0]->next(vertex_src);
    edges[2] = edge->next(vertex_dst);
    edges[3] = edges[2]->next(vertex_dst);

    CGAL_SS3_CORE_TRACE_V(16, "flip_edge = " << flip_edge);

    if (flip_edge) {
      facets[3]->remove_vertex(vertex_src);
      facets[2]->add_vertex(vertex_src);
      facets[1]->remove_vertex(vertex_dst);
      facets[0]->add_vertex(vertex_dst);

      facets[1]->remove_edge(edge);
      facets[3]->remove_edge(edge);
      edge->set_facet_L(facets[0]);
      edge->set_facet_R(facets[2]);
      facets[0]->add_edge(edge);
      facets[2]->add_edge(edge);

      if (edges[0]->get_vertex_src() == vertex_src) {
        edges[0]->replace_vertex_src(vertex_dst);
      } else if (edges[0]->get_vertex_dst() == vertex_src) {
        edges[0]->replace_vertex_dst(vertex_dst);
      }
      if (edges[2]->get_vertex_src() == vertex_dst) {
        edges[2]->replace_vertex_src(vertex_src);
      } else if (edges[2]->get_vertex_dst() == vertex_dst) {
        edges[2]->replace_vertex_dst(vertex_src);
      }

      if (time_future_bound.has_value()) {
        Hds_utils::set_final_point(vertex_src, std::nullopt);
        Hds_utils::set_final_point(vertex_dst, std::nullopt);
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
      Hds_utils::get_sheet(edges[i])->add_node(node);
    }

    // update arcs and sheets
    for (VertexSPtr vertex_ext : { vertex_src, vertex_dst }) {
      Hds_utils::set_node(vertex_ext, node);

      ArcSPtr arc = create_arc(vertex_ext);
      skeleton_->add_arc(arc);

      // 'edge' is common to both extremities, and is lacking a new sheet at this point
      for (EdgeWPtr inc_edge_w : vertex_ext->edges()) {
        if (EdgeSPtr inc_edge = inc_edge_w.lock()) {
          if (inc_edge != edge) {
            Hds_utils::get_sheet(inc_edge)->add_arc(arc);
          }
        }
      }
    }

    // this sets up the two incident arcs, the nodes, etc.
    SheetSPtr sheet = create_sheet(edge);
    skeleton_->add_sheet(sheet);
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
        post_op_vertices_pierce_.insert(poe->get_vertex_src());
        post_op_vertices_pierce_.insert(poe->get_vertex_dst());
    }
    CGAL_postcondition(post_op_vertices_pierce_.size() == 6);

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_edge_merge_event(const Edge_merge_event_sptr& event,
                                       const FT& current_time,
                                       const std::optional<FT>& time_future_bound,
                                       const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Edge Merge Event  #######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->time();

    shit_to_event_time(polyhedron, current_time, event_time);

    FacetSPtr facet = event->get_facet();
    EdgeSPtr edge_1 = event->get_edge_1();
    EdgeSPtr edge_2 = event->get_edge_2();
    FacetSPtr facet_l = edge_1->get_facet_L();
    FacetSPtr facet_r = edge_1->get_facet_R();
    EdgeSPtr edge_toremove_1 = edge_1->next(facet);
    EdgeSPtr edge_toremove_2 = edge_toremove_1->next(facet);

    CGAL_assertion(edge_2 == edge_toremove_2->next(facet));

    CGAL_SS3_CORE_TRACE_V(4, "Edge #1: " << edge_1->get_ID());
    CGAL_SS3_CORE_TRACE_V(4, "Edge #2: " << edge_2->get_ID());

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(event->point());
    skeleton_->add_node(node);

    Hds_utils::get_arc(edge_toremove_1->src(facet))->close_arc(node);
    Hds_utils::get_arc(edge_toremove_1->dst(facet))->close_arc(node);
    Hds_utils::get_arc(edge_toremove_2->dst(facet))->close_arc(node);

    Hds_utils::get_sheet(edge_toremove_1)->add_node(node);
    Hds_utils::get_sheet(edge_toremove_2)->add_node(node);
    Hds_utils::get_sheet(edge_1)->add_node(node);
    // no need for edge_2, since edge_2's and edge_1's sheets are merged

    merge_sheets(edge_1, edge_2);
#endif

    VertexSPtr vertex = edge_toremove_1->dst(facet);
    VertexSPtr vertex_1 = edge_1->dst(facet);
    VertexSPtr vertex_2 = edge_2->src(facet);
    EdgeSPtr edge_b = edge_toremove_1->prev(edge_toremove_1->other(facet));
    EdgeSPtr edge_b1 = edge_1->prev(edge_1->other(facet));
    EdgeSPtr edge_b2 = edge_2->next(edge_2->other(facet));
    facet->remove_vertex(vertex);
    edge_1->other(facet)->add_vertex(vertex);
    if (edge_b1->get_vertex_src() == vertex_1) {
      edge_b1->replace_vertex_src(vertex);
    } else {
      edge_b1->replace_vertex_dst(vertex);
    }
    if (edge_b2->get_vertex_src() == vertex_2) {
      edge_b2->replace_vertex_src(vertex);
    } else {
      edge_b2->replace_vertex_dst(vertex);
    }
    if (edge_1->get_vertex_dst() == vertex_1) {
      if (edge_2->get_vertex_src() == vertex_2) {
        edge_1->replace_vertex_dst(edge_2->get_vertex_dst());
      } else {
        edge_1->replace_vertex_dst(edge_2->get_vertex_src());
      }
    } else {
      if (edge_2->get_vertex_src() == vertex_2) {
        edge_1->replace_vertex_src(edge_2->get_vertex_dst());
      } else {
        edge_1->replace_vertex_src(edge_2->get_vertex_src());
      }
    }
    edge_toremove_1->get_facet_L()->remove_edge(edge_toremove_1);
    edge_toremove_1->get_facet_R()->remove_edge(edge_toremove_1);
    polyhedron->remove_edge(edge_toremove_1);
    edge_toremove_2->get_facet_L()->remove_edge(edge_toremove_2);
    edge_toremove_2->get_facet_R()->remove_edge(edge_toremove_2);
    polyhedron->remove_edge(edge_toremove_2);
    edge_2->get_facet_L()->remove_edge(edge_2);
    edge_2->get_facet_R()->remove_edge(edge_2);
    polyhedron->remove_edge(edge_2);
    for (auto it_f = vertex_1->facets().begin(); it_f != vertex_1->facets().end(); ) { // no C++11
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        facet->remove_vertex(vertex_1);
      }
    }
    polyhedron->remove_vertex(vertex_1);
    for (auto it_f = vertex_2->facets().begin(); it_f != vertex_2->facets().end(); ) { // no C++11
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        facet->remove_vertex(vertex_2);
      }
    }
    polyhedron->remove_vertex(vertex_2);

    if (time_future_bound.has_value()) {
      Hds_utils::set_final_point(vertex, std::nullopt);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    Hds_utils::set_node(vertex, node);

    ArcSPtr arc = create_arc(vertex);
    skeleton_->add_arc(arc);

    for (EdgeWPtr edge_w : vertex->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        Hds_utils::get_sheet(edge)->add_node(node);
        Hds_utils::get_sheet(edge)->add_arc(arc);
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

    CGAL_assertion(!Hds_utils::is_reflex(vertex)); // just to see the configurations where this could not be the case
    post_op_vertices_pierce_.clear();

    // since all faces are getting smaller, we don't need to check unmodified edges
    post_op_edges_edgesplit_ = {{ edge_1, edge_b1, edge_b, edge_b2 }};
    CGAL_postcondition(post_op_edges_edgesplit_.size() == 4);

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_triangle_event(const Triangle_event_sptr& event,
                                     const FT& current_time,
                                     const std::optional<FT>& time_future_bound,
                                     const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#######  Handle Triangle Event  ########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->time();
    const Point_3& point = event->point();

    shit_to_event_time(polyhedron, current_time, event_time);

    std::array<VertexSPtr, 3> vertices = event->get_vertices();
    FacetSPtr facet = event->get_facet();

    CGAL_SS3_CORE_TRACE_V(4, "Facet: " << facet->get_ID());
    CGAL_SS3_CORE_TRACE_V(4, "VS:\n" << vertices[0]->to_string() << "\n"
                                     << vertices[1]->to_string() << "\n"
                                     << vertices[2]->to_string());

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(point);
    skeleton_->add_node(node);

    for (unsigned int i = 0; i < 3; ++i) {
      Hds_utils::get_arc(vertices[i])->close_arc(node);
    }
    std::array<EdgeSPtr, 3> edges = event->get_edges();
    for (unsigned int i = 0; i < 3; ++i) {
      Hds_utils::get_sheet(edges[i])->add_node(node);
    }
#endif

    if (facet->vertices().size() == 3) {
      polyhedron->remove_facet(facet);
    }
    for (unsigned int i = 0; i < 3; ++i) {
      EdgeSPtr edge = vertices[i]->find_edge(vertices[(i+1)%3]);
      if (edge->get_facet_L()) {
        edge->get_facet_L()->remove_edge(edge);
      }
      if (edge->get_facet_R()) {
        edge->get_facet_R()->remove_edge(edge);
      }
      polyhedron->remove_edge(edge);
    }
    VertexSPtr new_vertex = Vertex::create(point);
    Skeleton_vertex_data::create(new_vertex);
    for (unsigned int i = 0; i < 3; ++i) {
      EdgeSPtr edge = vertices[i]->edges().front().lock();
      if (edge->get_vertex_src() == vertices[i]) {
        edge->replace_vertex_src(new_vertex);
      } else if (edge->get_vertex_dst() == vertices[i]) {
        edge->replace_vertex_dst(new_vertex);
      }
      edge->get_facet_L()->remove_vertex(vertices[i]);
      edge->get_facet_R()->remove_vertex(vertices[i]);
      facet->remove_vertex(vertices[i]);
      polyhedron->remove_vertex(vertices[i]);
      if (!edge->get_facet_L()->contains_vertex(new_vertex)) {
        edge->get_facet_L()->add_vertex(new_vertex);
      }
      if (!edge->get_facet_R()->contains_vertex(new_vertex)) {
        edge->get_facet_R()->add_vertex(new_vertex);
      }
    }
    polyhedron->add_vertex(new_vertex);

    if (time_future_bound.has_value()) {
      Hds_utils::set_final_point(new_vertex, std::nullopt);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    Hds_utils::set_node(new_vertex, node);

    ArcSPtr arc = create_arc(new_vertex);
    skeleton_->add_arc(arc);

    for (EdgeWPtr edge_w : new_vertex->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        Hds_utils::get_sheet(edge)->add_node(node);
        Hds_utils::get_sheet(edge)->add_arc(arc);
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
    CGAL_assertion(!Hds_utils::is_reflex(new_vertex));
    post_op_vertices_pierce_.clear();

    // faces are getting smaller so no need to check unmodified edges
    post_op_edges_edgesplit_ = post_op_edges_;
    CGAL_postcondition(post_op_edges_edgesplit_.size() == 3);

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_dbl_edge_merge_event(const Dbl_edge_merge_event_sptr& event,
                                           const FT& current_time,
                                           const std::optional<FT>& /*time_future_bound*/,
                                           const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#######  Handle Dbl Edge Event  ########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->time();

    shit_to_event_time(polyhedron, current_time, event_time);

    EdgeSPtr edge_11 = event->get_edge_11();
    EdgeSPtr edge_12 = event->get_edge_12();
    EdgeSPtr edge_21 = event->get_edge_21();
    EdgeSPtr edge_22 = event->get_edge_22();

    std::array<VertexSPtr, 4> vertices = event->get_vertices();
    std::array<EdgeSPtr, 4> edges = event->get_edges();

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(event->point());
    skeleton_->add_node(node);

    for (unsigned int i = 0; i < 4; ++i) {
      Hds_utils::get_arc(vertices[i])->close_arc(node);
    }
    for (unsigned int i = 0; i < 4; ++i) {
      Hds_utils::get_sheet(edges[i])->add_node(node);
    }

    Hds_utils::get_sheet(edge_11)->add_node(node);
    Hds_utils::get_sheet(edge_21)->add_node(node);
    // no need for _2 because we merge sheets

    merge_sheets(edge_11, edge_12);
    merge_sheets(edge_21, edge_22);
#endif

    for (unsigned int i = 0; i < 4; ++i) {
      EdgeSPtr edge = edges[i];
      edge->get_facet_L()->remove_edge(edge);
      edge->get_facet_R()->remove_edge(edge);
      polyhedron->remove_edge(edge);
    }

    if (edge_11->get_vertex_dst() == vertices[0]) {
      if (edge_12->get_vertex_src() == vertices[2]) {
        edge_11->replace_vertex_dst(edge_12->get_vertex_dst());
      } else {
        edge_11->replace_vertex_dst(edge_12->get_vertex_src());
      }
    } else {
      if (edge_12->get_vertex_src() == vertices[2]) {
        edge_11->replace_vertex_src(edge_12->get_vertex_dst());
      } else {
        edge_11->replace_vertex_src(edge_12->get_vertex_src());
      }
    }
    edge_12->get_facet_L()->remove_edge(edge_12);
    edge_12->get_facet_R()->remove_edge(edge_12);
    polyhedron->remove_edge(edge_12);
    if (edge_21->get_vertex_dst() == vertices[1]) {
      if (edge_22->get_vertex_src() == vertices[3]) {
        edge_21->replace_vertex_dst(edge_22->get_vertex_dst());
      } else {
        edge_21->replace_vertex_dst(edge_22->get_vertex_src());
      }
    } else {
      if (edge_22->get_vertex_src() == vertices[3]) {
        edge_21->replace_vertex_src(edge_22->get_vertex_dst());
      } else {
        edge_21->replace_vertex_src(edge_22->get_vertex_src());
      }
    }
    edge_22->get_facet_L()->remove_edge(edge_22);
    edge_22->get_facet_R()->remove_edge(edge_22);
    polyhedron->remove_edge(edge_22);
    for (unsigned int i = 0; i < 4; ++i) {
      VertexSPtr vertex = vertices[i];
      for (auto it_f = vertex->facets().begin(); it_f != vertex->facets().end(); ) { // no C++11
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          facet->remove_vertex(vertex);
        }
      }
      polyhedron->remove_vertex(vertex);
    }

    // Gather relevant elements for local queue updates
    post_op_vertices_.clear();
    post_op_edges_ = {{ edge_11, edge_21 }};
    post_op_facets_ = {{ edge_11->get_facet_L(),
                         edge_11->get_facet_R(),
                         edge_21->get_facet_L(),
                         edge_21->get_facet_R() }};
    CGAL_postcondition(post_op_vertices_.empty() && post_op_edges_.size() == 2 && post_op_facets_.size() == 4);

    // faces are smaller so nothing from unmodified vertices
    post_op_vertices_VV_.clear();

    // no new vertices & only reducing the size of facets so no edge disconnection
    post_op_vertices_pierce_.clear();

    CGAL_assertion(!Hds_utils::is_reflex(edge_11));
    CGAL_assertion(!Hds_utils::is_reflex(edge_21));
    post_op_edges_edgesplit_.clear();

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_dbl_triangle_event(const Dbl_triangle_event_sptr& event,
                                         const FT& current_time,
                                         const std::optional<FT>& /*time_future_bound*/,
                                         const PolyhedronSPtr& polyhedron) {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#####  Handle Dbl Triangle Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->time();

    shit_to_event_time(polyhedron, current_time, event_time);

    EdgeSPtr edge = event->get_edge();
    CGAL_SS3_CORE_TRACE_V(4,"Edge:\n" << edge->to_string());

    std::array<VertexSPtr, 4> vertices = event->get_vertices();
    std::array<EdgeSPtr, 5> edges = event->get_edges();

    FacetSPtr facet_l = edge->get_facet_L();
    FacetSPtr facet_r = edge->get_facet_R();
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
    node->set_time(event_time);
    node->set_point(event->point());
    skeleton_->add_node(node);

    for (unsigned int i = 0; i < 4; ++i) {
      Hds_utils::get_arc(vertices[i])->close_arc(node);
    }
    for (unsigned int i = 0; i < 5; ++i) {
      Hds_utils::get_sheet(edges[i])->add_node(node);
    }

    Hds_utils::get_sheet(edge_l)->add_node(node);

    merge_sheets(edge_l, edge_r);
#endif

    if (facet_l->edges().size() == 3) {
      polyhedron->remove_facet(facet_l);
    }
    if (facet_r->edges().size() == 3) {
      polyhedron->remove_facet(facet_r);
    }
    for (unsigned int i = 0; i < 5; ++i) {
      EdgeSPtr edge = edges[i];
      if (edge->get_facet_L()) {
        edge->get_facet_L()->remove_edge(edge);
      }
      if (edge->get_facet_R()) {
        edge->get_facet_R()->remove_edge(edge);
      }
      polyhedron->remove_edge(edge);
    }

    if (edge_l->get_vertex_src() == vertex_l) {
      if (edge_r->get_vertex_dst() == vertex_r) {
        edge_l->replace_vertex_src(edge_r->get_vertex_src());
      } else {
        edge_l->replace_vertex_src(edge_r->get_vertex_dst());
      }
    } else {
      if (edge_r->get_vertex_dst() == vertex_r) {
        edge_l->replace_vertex_dst(edge_r->get_vertex_src());
      } else {
        edge_l->replace_vertex_dst(edge_r->get_vertex_dst());
      }
    }
    edge_r->get_facet_L()->remove_edge(edge_r);
    edge_r->get_facet_R()->remove_edge(edge_r);
    polyhedron->remove_edge(edge_r);

    for (unsigned int i = 0; i < 4; ++i) {
      VertexSPtr vertex = vertices[i];
      for (auto it_f = vertex->facets().begin(); it_f != vertex->facets().end(); ) { // no C++11
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          facet->remove_vertex(vertex);
        }
      }
      polyhedron->remove_vertex(vertex);
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

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_tetrahedron_event(const Tetrahedron_event_sptr& event,
                                        const FT& current_time,
                                        const std::optional<FT>& /*time_future_bound*/,
                                        const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Tetrahedron Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->time();

    shit_to_event_time(polyhedron, current_time, event_time);

    std::array<VertexSPtr, 4> vertices = event->get_vertices();
    std::array<EdgeSPtr, 6> edges = event->get_edges();
    std::array<FacetSPtr, 4> facets = event->get_facets();

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(event->point());
    skeleton_->add_node(node);

    for (unsigned int i = 0; i < 4; ++i) {
      Hds_utils::get_arc(vertices[i])->close_arc(node);
    }
    for (unsigned int i = 0; i < 6; ++i) {
      Hds_utils::get_sheet(edges[i])->add_node(node);
    }
#endif

    for (unsigned int i = 0; i < 4; ++i) {
      if (facets[i]->vertices().size() == 3) {
        polyhedron->remove_facet(facets[i]);
      }
    }
    for (unsigned int i = 0; i < 6; ++i) {
      EdgeSPtr edge = edges[i];
      if (edge->get_facet_L()) {
        edge->get_facet_L()->remove_edge(edge);
      }
      if (edge->get_facet_R()) {
        edge->get_facet_R()->remove_edge(edge);
      }
      polyhedron->remove_edge(edge);
    }
    for (unsigned int i = 0; i < 4; ++i) {
      VertexSPtr vertex = vertices[i];
      for (auto it_f = vertex->facets().begin(); it_f != vertex->facets().end(); ) { // no C++11
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          facet->remove_vertex(vertex);
        }
      }
      polyhedron->remove_vertex(vertex);
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

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_vertex_event(const Vertex_event_sptr& event,
                                   const FT& current_time,
                                   const std::optional<FT>& time_future_bound,
                                   const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "########  Handle Vertex Event  #########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->time();

    shit_to_event_time(polyhedron, current_time, event_time);

    VertexSPtr vertex_1 = event->get_vertex_1();
    VertexSPtr vertex_2 = event->get_vertex_2();
    FacetSPtr facet_1 = event->get_facet_1();
    FacetSPtr facet_2 = event->get_facet_2();

    EdgeSPtr edge_tomerge_1 = EdgeSPtr();
    EdgeSPtr edge_11 = EdgeSPtr();
    EdgeSPtr edge_12 = EdgeSPtr();
    EdgeSPtr edge_tomerge_2 = EdgeSPtr();
    EdgeSPtr edge_21 = EdgeSPtr();
    EdgeSPtr edge_22 = EdgeSPtr();
    for (EdgeWPtr edge_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->get_facet_L() == facet_1 && edge->get_facet_R() == facet_2) ||
            (edge->get_facet_L() == facet_2 && edge->get_facet_R() == facet_1)) {
          edge_tomerge_1 = edge;
          continue;
        }
        if (edge->get_facet_L() == facet_1 || edge->get_facet_R() == facet_1) {
          edge_11 = edge;
        }
        if (edge->get_facet_L() == facet_2 || edge->get_facet_R() == facet_2) {
          edge_12 = edge;
        }
      }
    }
    for (EdgeWPtr edge_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->get_facet_L() == facet_1 && edge->get_facet_R() == facet_2) ||
            (edge->get_facet_L() == facet_2 && edge->get_facet_R() == facet_1)) {
          edge_tomerge_2 = edge;
          continue;
        }
        if (edge->get_facet_L() == facet_1 || edge->get_facet_R() == facet_1) {
          edge_21 = edge;
        }
        if (edge->get_facet_L() == facet_2 || edge->get_facet_R() == facet_2) {
          edge_22 = edge;
        }
      }
    }
    FacetSPtr facet_1b = edge_11->get_facet_L();
    if (facet_1b == facet_1 || facet_1b == facet_2) {
      facet_1b = edge_11->get_facet_R();
    }
    FacetSPtr facet_2b = edge_21->get_facet_L();
    if (facet_2b == facet_1 || facet_2b == facet_2) {
      facet_2b = edge_21->get_facet_R();
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(event->point());
    skeleton_->add_node(node);

    Hds_utils::get_arc(vertex_1)->close_arc(node);
    Hds_utils::get_arc(vertex_2)->close_arc(node);

    Hds_utils::get_sheet(edge_11)->add_node(node);
    Hds_utils::get_sheet(edge_12)->add_node(node);
    Hds_utils::get_sheet(edge_21)->add_node(node);
    Hds_utils::get_sheet(edge_22)->add_node(node);

    Hds_utils::get_sheet(edge_tomerge_1)->add_node(node);
    // no need it to edge_tomerge_2's sheet because that sheet gets merged with edge_tomerge_1's

    merge_sheets(edge_tomerge_1, edge_tomerge_2);
#endif

    if (edge_tomerge_1->get_vertex_src() == vertex_1) {
      if (edge_tomerge_2->get_vertex_src() == vertex_2) {
        edge_tomerge_1->replace_vertex_src(edge_tomerge_2->get_vertex_dst());
      } else {
        edge_tomerge_1->replace_vertex_src(edge_tomerge_2->get_vertex_src());
      }
    } else {
      if (edge_tomerge_2->get_vertex_src() == vertex_2) {
        edge_tomerge_1->replace_vertex_dst(edge_tomerge_2->get_vertex_dst());
      } else {
        edge_tomerge_1->replace_vertex_dst(edge_tomerge_2->get_vertex_src());
      }
    }
    facet_1->remove_vertex(vertex_2);
    facet_2->remove_vertex(vertex_1);
    facet_1b->add_vertex(vertex_2);
    facet_2b->add_vertex(vertex_1);
    edge_tomerge_2->replace_vertex_src(vertex_1);
    edge_tomerge_2->replace_vertex_dst(vertex_2);
    edge_tomerge_2->replace_facet_L(facet_1b);
    edge_tomerge_2->replace_facet_R(facet_2b);
    if (edge_12->get_vertex_src() == vertex_1) {
      edge_12->replace_vertex_src(vertex_2);
    } else {
      edge_12->replace_vertex_dst(vertex_2);
    }
    if (edge_21->get_vertex_src() == vertex_2) {
      edge_21->replace_vertex_src(vertex_1);
    } else {
      edge_21->replace_vertex_dst(vertex_1);
    }

    if (time_future_bound.has_value()) {
      Hds_utils::set_final_point(vertex_1, std::nullopt);
      Hds_utils::set_final_point(vertex_2, std::nullopt);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    Hds_utils::clear_sheet(edge_tomerge_2);

    Hds_utils::set_node(vertex_1, node);
    Hds_utils::set_node(vertex_2, node);

    ArcSPtr arc_1 = create_arc(vertex_1);
    skeleton_->add_arc(arc_1);
    ArcSPtr arc_2 = create_arc(vertex_2);
    skeleton_->add_arc(arc_2);

    Hds_utils::get_sheet(edge_11)->add_arc(arc_1);
    Hds_utils::get_sheet(edge_21)->add_arc(arc_1);
    Hds_utils::get_sheet(edge_12)->add_arc(arc_2);
    Hds_utils::get_sheet(edge_22)->add_arc(arc_2);

    SheetSPtr sheet = create_sheet(edge_tomerge_2);
    skeleton_->add_sheet(sheet);

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
      post_op_vertices_pierce_.insert(poe->get_vertex_src());
      post_op_vertices_pierce_.insert(poe->get_vertex_dst());
    }

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_flip_vertex_event(const Flip_vertex_event_sptr& event,
                                        const FT& current_time,
                                        const std::optional<FT>& time_future_bound,
                                        const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Flip Vertex Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->time();

    shit_to_event_time(polyhedron, current_time, event_time);

    VertexSPtr vertex_1 = event->get_vertex_1();
    VertexSPtr vertex_2 = event->get_vertex_2();
    FacetSPtr facet_1 = event->get_facet_1();
    FacetSPtr facet_2 = event->get_facet_2();

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(event->point());
    skeleton_->add_node(node);

    Hds_utils::get_arc(vertex_1)->close_arc(node);
    Hds_utils::get_arc(vertex_2)->close_arc(node);

    for (VertexSPtr vertex : { vertex_1, vertex_2 }) {
      for (EdgeWPtr edge_w : vertex->edges()) {
        if (EdgeSPtr edge = edge_w.lock()) {
          Hds_utils::get_sheet(edge)->add_node(node);
        }
      }
    }
#endif

    EdgeSPtr edge_1;
    for (EdgeWPtr edge_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->get_facet_L() == facet_1 && edge->get_facet_R() == facet_2) ||
            (edge->get_facet_L() == facet_2 && edge->get_facet_R() == facet_1)) {
          edge_1 = edge;
          break;
        }
      }
    }
    EdgeSPtr edge_2;
    for (EdgeWPtr edge_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->get_facet_L() == facet_1 && edge->get_facet_R() == facet_2) ||
            (edge->get_facet_L() == facet_2 && edge->get_facet_R() == facet_1)) {
          edge_2 = edge;
          break;
        }
      }
    }

    if (edge_1->get_vertex_src() == vertex_1) {
      edge_1->replace_vertex_src(vertex_2);
    } else if (edge_1->get_vertex_dst() == vertex_1) {
      edge_1->replace_vertex_dst(vertex_2);
    }
    if (edge_2->get_vertex_src() == vertex_2) {
      edge_2->replace_vertex_src(vertex_1);
    } else if (edge_2->get_vertex_dst() == vertex_2) {
      edge_2->replace_vertex_dst(vertex_1);
    }

    if (time_future_bound.has_value()) {
      Hds_utils::set_final_point(vertex_1, std::nullopt);
      Hds_utils::set_final_point(vertex_2, std::nullopt);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    Hds_utils::set_node(vertex_1, node);
    Hds_utils::set_node(vertex_2, node);

    ArcSPtr arc_1 = create_arc(vertex_1);
    skeleton_->add_arc(arc_1);
    ArcSPtr arc_2 = create_arc(vertex_2);
    skeleton_->add_arc(arc_2);

    for (VertexSPtr vertex : { vertex_1, vertex_2 }) {
      for (EdgeWPtr edge_w : vertex->edges()) {
        if (EdgeSPtr edge = edge_w.lock()) {
          Hds_utils::get_sheet(edge)->add_arc(Hds_utils::get_arc(vertex));
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
      post_op_vertices_pierce_.insert(poe->get_vertex_src());
      post_op_vertices_pierce_.insert(poe->get_vertex_dst());
    }

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_surface_event(const Surface_event_sptr& event,
                                    const FT& current_time,
                                    const std::optional<FT>& time_future_bound,
                                    const PolyhedronSPtr& polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#######  Handle Surface Event  #########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    CGAL_SS3_CORE_TRACE_V(4, "Edge A = " << event->get_edge_1()->to_string());
    CGAL_SS3_CORE_TRACE_V(4, "Edge B = " << event->get_edge_2()->to_string());

    const FT& event_time = event->time();
    const Point_3& point = event->point();

    shit_to_event_time(polyhedron, current_time, event_time);

    EdgeSPtr edge_1 = event->get_edge_1();
    EdgeSPtr edge_2 = event->get_edge_2();
    FacetSPtr facet_1_src = edge_1->get_facet_src();
    FacetSPtr facet_1_dst = edge_1->get_facet_dst();

    VertexSPtr vertex = VertexSPtr();
    EdgeSPtr edge_b1 = EdgeSPtr();
    EdgeSPtr edge_b2 = EdgeSPtr();
    if (edge_2->get_facet_L() == facet_1_src) {
      vertex = edge_1->get_vertex_src();
      edge_b1 = edge_1->prev(edge_1->get_facet_L());
      edge_b2 = edge_1->next(edge_1->get_facet_R());
    } else if (edge_2->get_facet_R() == facet_1_src) {
      vertex = edge_1->get_vertex_src();
      edge_b1 = edge_1->next(edge_1->get_facet_R());
      edge_b2 = edge_1->prev(edge_1->get_facet_L());
    } else if (edge_2->get_facet_L() == facet_1_dst) {
      vertex = edge_1->get_vertex_dst();
      edge_b1 = edge_1->prev(edge_1->get_facet_R());
      edge_b2 = edge_1->next(edge_1->get_facet_L());
    } else if (edge_2->get_facet_R() == facet_1_dst) {
      vertex = edge_1->get_vertex_dst();
      edge_b1 = edge_1->next(edge_1->get_facet_L());
      edge_b2 = edge_1->prev(edge_1->get_facet_R());
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(point);
    skeleton_->add_node(node);

    if (facet_1_src == edge_2->get_facet_L() || facet_1_src == edge_2->get_facet_R()) {
      Hds_utils::get_arc(edge_1->get_vertex_src())->close_arc(node);
    }
    if (facet_1_dst == edge_2->get_facet_L() || facet_1_dst == edge_2->get_facet_R()) {
      Hds_utils::get_arc(edge_1->get_vertex_dst())->close_arc(node);
    }

    for (EdgeWPtr edge_w : vertex->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        Hds_utils::get_sheet(edge)->add_node(node);
      }
    }

    Hds_utils::get_sheet(edge_2)->add_node(node);
#endif

    VertexSPtr vertex_21 = Vertex::create(point);
    VertexSPtr vertex_22 = Vertex::create(point);
    Skeleton_vertex_data::create(vertex_21);
    Skeleton_vertex_data::create(vertex_22);
    polyhedron->add_vertex(vertex_21);
    polyhedron->add_vertex(vertex_22);
    if (edge_b1->get_vertex_src() == vertex) {
      edge_b1->replace_vertex_src(vertex_21);
    } else if (edge_b1->get_vertex_dst() == vertex) {
      edge_b1->replace_vertex_dst(vertex_21);
    }
    if (edge_b2->get_vertex_src() == vertex) {
      edge_b2->replace_vertex_src(vertex_22);
    } else if (edge_b2->get_vertex_dst() == vertex) {
      edge_b2->replace_vertex_dst(vertex_22);
    }
    edge_b1->get_facet_L()->add_vertex(vertex_21);
    edge_b1->get_facet_R()->add_vertex(vertex_21);
    edge_b2->get_facet_L()->add_vertex(vertex_22);
    edge_b2->get_facet_R()->add_vertex(vertex_22);

    EdgeSPtr edge_tmp = edge_2->split(vertex);
    EdgeSPtr edge_21 = edge_2->split(vertex_21);
    EdgeSPtr edge_22 = edge_tmp;
    edge_tmp = edge_22->split(vertex_22);

    Skeleton_edge_data::create(edge_tmp);
    Skeleton_edge_data::create(edge_21);
    Skeleton_edge_data::create(edge_22);

    if (edge_2->get_facet_L() == facet_1_src ||
        edge_2->get_facet_L() == facet_1_dst) {
      edge_2->get_facet_L()->remove_vertex(vertex);
      if (vertex == edge_1->get_vertex_src()) {
        edge_21->replace_facet_L(edge_1->get_facet_L());
        edge_22->replace_facet_L(edge_1->get_facet_R());
      } else if (vertex == edge_1->get_vertex_dst()) {
        edge_21->replace_facet_L(edge_1->get_facet_R());
        edge_22->replace_facet_L(edge_1->get_facet_L());
      }
    } else {
      edge_2->get_facet_R()->remove_vertex(vertex);
      if (vertex == edge_1->get_vertex_src()) {
        edge_21->replace_facet_R(edge_1->get_facet_R());
        edge_22->replace_facet_R(edge_1->get_facet_L());
      } else if (vertex == edge_1->get_vertex_dst()) {
        edge_21->replace_facet_R(edge_1->get_facet_L());
        edge_22->replace_facet_R(edge_1->get_facet_R());
      }
    }

    if (time_future_bound.has_value()) {
      Hds_utils::set_final_point(vertex, std::nullopt);
      Hds_utils::set_final_point(vertex_21, std::nullopt);
      Hds_utils::set_final_point(vertex_22, std::nullopt);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    set_sheet(edge_tmp, Hds_utils::get_sheet(edge_2));

    Hds_utils::set_node(vertex, node);
    ArcSPtr arc = create_arc(vertex);
    skeleton_->add_arc(arc);
    Hds_utils::get_sheet(edge_1)->add_arc(arc); // other 2 sheets are edge_21/edge_22's, created below

    Hds_utils::set_node(vertex_21, node);
    ArcSPtr arc_21 = create_arc(vertex_21);
    skeleton_->add_arc(arc_21);
    Hds_utils::get_sheet(edge_b1)->add_arc(arc_21);
    Hds_utils::get_sheet(edge_2)->add_arc(arc_21); // third sheet is edge_21's, created below

    Hds_utils::set_node(vertex_22, node);
    ArcSPtr arc_22 = create_arc(vertex_22);
    skeleton_->add_arc(arc_22);
    Hds_utils::get_sheet(edge_b2)->add_arc(arc_22);
    Hds_utils::get_sheet(edge_tmp)->add_arc(arc_22); // third sheet is edge_22's, created below

    SheetSPtr sheet_21 = create_sheet(edge_21);
    skeleton_->add_sheet(sheet_21);

    SheetSPtr sheet_22 = create_sheet(edge_22);
    skeleton_->add_sheet(sheet_22);
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

    post_op_facets_ = {{ edge_1->get_facet_L(), edge_1->get_facet_R(),
                         edge_2->get_facet_L(), edge_2->get_facet_R() }};

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
    post_op_vertices_pierce_ = {{ edge_1->get_vertex_src(), edge_1->get_vertex_dst(), vertex_21, vertex_22 }};
    CGAL_postcondition(post_op_vertices_pierce_.size() == 4);

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_polyhedron_split_event(const Polyhedron_split_event_sptr& event,
                                             const FT& current_time,
                                             const std::optional<FT>& time_future_bound,
                                             const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "####  Handle Polyhedron Split Event  ###");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->time();

    shit_to_event_time(polyhedron, current_time, event_time);

    EdgeSPtr edge_1 = event->get_edge_1();
    EdgeSPtr edge_2 = event->get_edge_2();
    FacetSPtr facet_1_src = edge_1->get_facet_src();
    FacetSPtr facet_1_dst = edge_1->get_facet_dst();

    CGAL_SS3_CORE_TRACE_V(4, "Edge 1 = " << edge_1->to_string());
    CGAL_SS3_CORE_TRACE_V(4, "Edge 2 = " << edge_2->to_string());

    VertexSPtr vertex_l = VertexSPtr();
    VertexSPtr vertex_r = VertexSPtr();
    if (edge_2->get_facet_L() == facet_1_src &&
        edge_2->get_facet_R() == facet_1_dst) {
      vertex_l = edge_1->get_vertex_src();
      vertex_r = edge_1->get_vertex_dst();
    }
    if (edge_2->get_facet_L() == facet_1_dst &&
        edge_2->get_facet_R() == facet_1_src) {
      vertex_l = edge_1->get_vertex_dst();
      vertex_r = edge_1->get_vertex_src();
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(event->point());
    skeleton_->add_node(node);

    if (facet_1_src == edge_2->get_facet_L() || facet_1_src == edge_2->get_facet_R()) {
      Hds_utils::get_arc(edge_1->get_vertex_src())->close_arc(node);
    }
    if (facet_1_dst == edge_2->get_facet_L() || facet_1_dst == edge_2->get_facet_R()) {
      Hds_utils::get_arc(edge_1->get_vertex_dst())->close_arc(node);
    }

    std::array<EdgeSPtr, 4> edges;
    edges[0] = edge_1->next(vertex_l);
    edges[1] = edges[0]->next(vertex_l);
    edges[2] = edge_1->next(vertex_r);
    edges[3] = edges[2]->next(vertex_r);

    Hds_utils::get_sheet(edge_1)->add_node(node);
    for (int i = 0; i < 4; ++i) {
      Hds_utils::get_sheet(edges[i])->add_node(node);
    }
    Hds_utils::get_sheet(edge_2)->add_node(node);
#endif

    EdgeSPtr edge_22;
    if (edge_2->get_facet_L() == facet_1_src &&
        edge_2->get_facet_R() == facet_1_dst) {
      EdgeSPtr edge_l = edge_1->prev(edge_1->get_vertex_dst());
      if (edge_l->get_vertex_src() == vertex_r) {
        edge_l->replace_vertex_src(vertex_l);
      } else if (edge_l->get_vertex_dst() == vertex_r) {
        edge_l->replace_vertex_dst(vertex_l);
      }
      EdgeSPtr edge_r = edge_1->prev(edge_1->get_vertex_src());
      if (edge_r->get_vertex_src() == vertex_l) {
        edge_r->replace_vertex_src(vertex_r);
      } else if (edge_r->get_vertex_dst() == vertex_l) {
        edge_r->replace_vertex_dst(vertex_r);
      }

      edge_1->get_facet_R()->remove_vertex(vertex_l);
      edge_2->get_facet_R()->add_vertex(vertex_l);
      edge_1->get_facet_L()->remove_vertex(vertex_r);
      edge_2->get_facet_L()->add_vertex(vertex_r);

      edge_22 = edge_1;
      edge_22->replace_vertex_dst(edge_2->get_vertex_dst());
      edge_2->replace_vertex_dst(vertex_l);
      edge_22->replace_vertex_src(vertex_r);

      edge_22->replace_facet_L(edge_2->get_facet_L());
      edge_22->replace_facet_R(edge_2->get_facet_R());
    }
    // @todo just else it...
    if (edge_2->get_facet_L() == facet_1_dst &&
        edge_2->get_facet_R() == facet_1_src) {
      vertex_l = edge_1->get_vertex_dst();
      vertex_r = edge_1->get_vertex_src();

      EdgeSPtr edge_l = edge_1->next(edge_1->get_vertex_src());
      if (edge_l->get_vertex_src() == vertex_r) {
        edge_l->replace_vertex_src(vertex_l);
      } else if (edge_l->get_vertex_dst() == vertex_r) {
        edge_l->replace_vertex_dst(vertex_l);
      }
      EdgeSPtr edge_r = edge_1->next(edge_1->get_vertex_dst());
      if (edge_r->get_vertex_src() == vertex_l) {
        edge_r->replace_vertex_src(vertex_r);
      } else if (edge_r->get_vertex_dst() == vertex_l) {
        edge_r->replace_vertex_dst(vertex_r);
      }

      edge_1->get_facet_R()->remove_vertex(vertex_l);
      edge_2->get_facet_R()->add_vertex(vertex_l);
      edge_1->get_facet_L()->remove_vertex(vertex_r);
      edge_2->get_facet_L()->add_vertex(vertex_r);

      edge_22 = edge_1;
      edge_22->replace_vertex_dst(edge_2->get_vertex_dst());
      edge_2->replace_vertex_dst(vertex_r);
      edge_22->replace_vertex_src(vertex_l);

      edge_22->replace_facet_L(edge_2->get_facet_L());
      edge_22->replace_facet_R(edge_2->get_facet_R());
    }

    if (time_future_bound.has_value()) {
      Hds_utils::set_final_point(vertex_l, std::nullopt);
      Hds_utils::set_final_point(vertex_r, std::nullopt);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    Hds_utils::set_node(vertex_l, node);
    Hds_utils::set_node(vertex_r, node);

    ArcSPtr arc_l = create_arc(vertex_l);
    skeleton_->add_arc(arc_l);
    ArcSPtr arc_r = create_arc(vertex_r);
    skeleton_->add_arc(arc_r);

    set_sheet(edge_22, Hds_utils::get_sheet(edge_2));

    Hds_utils::get_sheet(edge_2)->add_arc(arc_l);
    Hds_utils::get_sheet(edge_2)->add_arc(arc_r);

    for (EdgeWPtr edge_w : vertex_l->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        if (Hds_utils::get_sheet(edge) != Hds_utils::get_sheet(edge_2)) {
          Hds_utils::get_sheet(edge)->add_arc(arc_l);
        }
      }
    }

    for (EdgeWPtr edge_w : vertex_r->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        if (Hds_utils::get_sheet(edge) != Hds_utils::get_sheet(edge_2)) {
          Hds_utils::get_sheet(edge)->add_arc(arc_r);
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

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_split_merge_event(const Split_merge_event_sptr& event,
                                        const FT& current_time,
                                        const std::optional<FT>& time_future_bound,
                                        const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Split Merge Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->time();

    shit_to_event_time(polyhedron, current_time, event_time);

    VertexSPtr vertex_1 = event->get_vertex_1();
    VertexSPtr vertex_2 = event->get_vertex_2();
    FacetSPtr facet_1 = event->get_facet_1();
    FacetSPtr facet_2 = event->get_facet_2();

    EdgeSPtr edge_tomerge_1 = EdgeSPtr();
    EdgeSPtr edge_11 = EdgeSPtr();
    EdgeSPtr edge_12 = EdgeSPtr();
    EdgeSPtr edge_tomerge_2 = EdgeSPtr();
    EdgeSPtr edge_21 = EdgeSPtr();
    EdgeSPtr edge_22 = EdgeSPtr();
    for (EdgeWPtr edge_wptr : vertex_1->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->get_facet_L() == facet_1 && edge->get_facet_R() == facet_2) ||
            (edge->get_facet_L() == facet_2 && edge->get_facet_R() == facet_1)) {
          edge_tomerge_1 = edge;
          continue;
        }
        if (edge->get_facet_L() == facet_1 || edge->get_facet_R() == facet_1) {
          edge_11 = edge;
        }
        if (edge->get_facet_L() == facet_2 || edge->get_facet_R() == facet_2) {
          edge_12 = edge;
        }
      }
    }
    for (EdgeWPtr edge_wptr : vertex_2->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if ((edge->get_facet_L() == facet_1 && edge->get_facet_R() == facet_2) ||
            (edge->get_facet_L() == facet_2 && edge->get_facet_R() == facet_1)) {
          edge_tomerge_2 = edge;
          continue;
        }
        if (edge->get_facet_L() == facet_1 || edge->get_facet_R() == facet_1) {
          edge_21 = edge;
        }
        if (edge->get_facet_L() == facet_2 || edge->get_facet_R() == facet_2) {
          edge_22 = edge;
        }
      }
    }
    FacetSPtr facet_1b = edge_11->get_facet_L();
    if (facet_1b == facet_1 || facet_1b == facet_2) {
      facet_1b = edge_11->get_facet_R();
    }
    FacetSPtr facet_2b = edge_21->get_facet_L();
    if (facet_2b == facet_1 || facet_2b == facet_2) {
      facet_2b = edge_21->get_facet_R();
    }
    EdgeSPtr edge_tosplit = EdgeSPtr();
    // edge_tosplit = facet_1b->find_edge(facet_2b);
    EdgeSPtr edge_cur = edge_11->next(facet_1b);
    while (edge_cur != edge_11) {
      if ((edge_cur->get_facet_L() == facet_1b && edge_cur->get_facet_R() == facet_2b) ||
          (edge_cur->get_facet_R() == facet_1b && edge_cur->get_facet_L() == facet_2b)) {
        edge_tosplit = edge_cur;
        break;
      }
      edge_cur = edge_cur->next(facet_1b);
    }
    CGAL_SS3_DEBUG_SPTR(edge_tosplit);

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(event->point());
    skeleton_->add_node(node);

    Hds_utils::get_arc(vertex_1)->close_arc(node);
    Hds_utils::get_arc(vertex_2)->close_arc(node);

    for (EdgeSPtr e : {edge_11, edge_12, edge_21, edge_22}) {
      Hds_utils::get_sheet(e)->add_node(node);
    }
    Hds_utils::get_sheet(edge_tosplit)->add_node(node);
    Hds_utils::get_sheet(edge_tomerge_1)->add_node(node);

    merge_sheets(edge_tomerge_1, edge_tomerge_2);
#endif

    if (edge_tomerge_1->get_vertex_src() == vertex_1) {
      if (edge_tomerge_2->get_vertex_src() == vertex_2) {
        edge_tomerge_1->replace_vertex_src(edge_tomerge_2->get_vertex_dst());
      } else {
        edge_tomerge_1->replace_vertex_src(edge_tomerge_2->get_vertex_src());
      }
    } else {
      if (edge_tomerge_2->get_vertex_src() == vertex_2) {
        edge_tomerge_1->replace_vertex_dst(edge_tomerge_2->get_vertex_dst());
      } else {
        edge_tomerge_1->replace_vertex_dst(edge_tomerge_2->get_vertex_src());
      }
    }
    if (edge_12->get_vertex_dst() == vertex_1) {
      edge_12->replace_vertex_dst(vertex_2);
    } else {
      edge_12->replace_vertex_src(vertex_2);
    }
    if (edge_21->get_vertex_dst() == vertex_2) {
      edge_21->replace_vertex_dst(vertex_1);
    } else {
      edge_21->replace_vertex_src(vertex_1);
    }
    facet_1->remove_vertex(vertex_2);
    facet_2->remove_vertex(vertex_1);
    facet_1b->add_vertex(vertex_2);
    facet_2b->add_vertex(vertex_1);
    if (edge_tosplit->get_facet_L() == facet_1b &&
        edge_tosplit->get_facet_R() == facet_2b) {
      edge_tomerge_2->replace_vertex_src(edge_tosplit->get_vertex_src());
      edge_tomerge_2->replace_vertex_dst(vertex_2);
      edge_tosplit->replace_vertex_src(vertex_1);
    } else if (edge_tosplit->get_facet_L() == facet_2b &&
               edge_tosplit->get_facet_R() == facet_1b) {
      edge_tomerge_2->replace_vertex_dst(edge_tosplit->get_vertex_dst());
      edge_tomerge_2->replace_vertex_src(vertex_2);
      edge_tosplit->replace_vertex_dst(vertex_1);
    }
    edge_tomerge_2->replace_facet_L(edge_tosplit->get_facet_L());
    edge_tomerge_2->replace_facet_R(edge_tosplit->get_facet_R());

    if (time_future_bound.has_value()) {
      Hds_utils::set_final_point(vertex_1, std::nullopt);
      Hds_utils::set_final_point(vertex_2, std::nullopt);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    Hds_utils::set_node(vertex_1, node);
    Hds_utils::set_node(vertex_2, node);

    ArcSPtr arc_1 = create_arc(vertex_1);
    skeleton_->add_arc(arc_1);
    ArcSPtr arc_2 = create_arc(vertex_2);
    skeleton_->add_arc(arc_2);

    set_sheet(edge_tomerge_2, Hds_utils::get_sheet(edge_tosplit));

    for (VertexSPtr vertex : { vertex_1, vertex_2 }) {
      for (EdgeWPtr edge_w : vertex->edges()) {
        if (EdgeSPtr edge = edge_w.lock()) {
          Hds_utils::get_sheet(edge)->add_arc(Hds_utils::get_arc(vertex));
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
    CGAL_assertion(!Hds_utils::is_reflex(vertex_1));
    CGAL_assertion(!Hds_utils::is_reflex(vertex_2));
    post_op_vertices_pierce_.clear();

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_edge_split_event(const Edge_split_event_sptr& event,
                                       const FT& current_time,
                                       const std::optional<FT>& time_future_bound,
                                       const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Edge Split Event  #######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    // this must be done **BEFORE** the shift because after the shift there is
    // an intersection between the two edges, and so the orientation becomes coplanar
    int orientation = Kernel_wrapper::orientation(Hds_utils::line(event->get_edge_1()),
                                                 Hds_utils::line(event->get_edge_2()));
    CGAL_assertion(orientation != 0);

    const FT& event_time = event->time();
    const Point_3& point = event->point();

    shit_to_event_time(polyhedron, current_time, event_time);

    EdgeSPtr edge_1 = event->get_edge_1();
    EdgeSPtr edge_2 = event->get_edge_2();

    CGAL_SS3_CORE_TRACE("edge_1 = " << event->get_edge_1()->to_string());
    CGAL_SS3_CORE_TRACE("edge_2 = " << event->get_edge_2()->to_string());

    FacetSPtr facet_l1 = edge_1->get_facet_L();
    FacetSPtr facet_r1 = edge_1->get_facet_R();
    FacetSPtr facet_l2 = edge_2->get_facet_L();
    FacetSPtr facet_r2 = edge_2->get_facet_R();

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(point);
    skeleton_->add_node(node);

    FacetSPtr facet_1_src = edge_1->get_facet_src();
    FacetSPtr facet_1_dst = edge_1->get_facet_dst();

    if (facet_1_src == facet_l2 || facet_1_src == facet_r2) {
      Hds_utils::get_arc(edge_1->get_vertex_src())->close_arc(node);
    }
    if (facet_1_dst == facet_l2 || facet_1_dst == facet_r2) {
      Hds_utils::get_arc(edge_1->get_vertex_dst())->close_arc(node);
    }

    Hds_utils::get_sheet(edge_1)->add_node(node);
    Hds_utils::get_sheet(edge_2)->add_node(node);
#endif

    std::array<VertexSPtr, 4> vertices;
    for (unsigned int i = 0; i < 4; ++i) {
      vertices[i] = Vertex::create(point);
      Skeleton_vertex_data::create(vertices[i]);
      polyhedron->add_vertex(vertices[i]);
    }
    std::array<EdgeSPtr, 4> edges;
    for (unsigned int i = 0; i < 4; ++i) {
      edges[i] = Edge::create(vertices[i], vertices[(i+1)%4]);
      Skeleton_edge_data::create(edges[i]);
      polyhedron->add_edge(edges[i]);
    }

    if (orientation > 0) {
      edges[0]->set_facet_L(edge_2->get_facet_R());
      edges[0]->set_facet_R(edge_1->get_facet_R());
      edges[1]->set_facet_L(edge_2->get_facet_L());
      edges[1]->set_facet_R(edge_1->get_facet_R());
      edges[2]->set_facet_L(edge_2->get_facet_L());
      edges[2]->set_facet_R(edge_1->get_facet_L());
      edges[3]->set_facet_L(edge_2->get_facet_R());
      edges[3]->set_facet_R(edge_1->get_facet_L());
    } else {
      edges[0]->set_facet_L(edge_1->get_facet_L());
      edges[0]->set_facet_R(edge_2->get_facet_L());
      edges[1]->set_facet_L(edge_1->get_facet_L());
      edges[1]->set_facet_R(edge_2->get_facet_R());
      edges[2]->set_facet_L(edge_1->get_facet_R());
      edges[2]->set_facet_R(edge_2->get_facet_R());
      edges[3]->set_facet_L(edge_1->get_facet_R());
      edges[3]->set_facet_R(edge_2->get_facet_L());
    }
    EdgeSPtr edge_12 = edge_1->split(vertices[2]);
    EdgeSPtr edge_22 = edge_2->split(vertices[3]);
    Skeleton_edge_data::create(edge_12);
    Skeleton_edge_data::create(edge_22);
    edge_1->replace_vertex_dst(vertices[0]);
    edge_2->replace_vertex_dst(vertices[1]);
    for (unsigned int i = 0; i < 4; ++i) {
      // adds the vertices also
      edges[i]->get_facet_L()->add_edge(edges[i]);
      edges[i]->get_facet_R()->add_edge(edges[i]);
    }

    if (time_future_bound.has_value()) {
      for (std::size_t i=0; i<4; ++i) {
        Hds_utils::set_final_point(vertices[i], std::nullopt);
      }
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    set_sheet(edge_12, Hds_utils::get_sheet(edge_1));
    set_sheet(edge_22, Hds_utils::get_sheet(edge_2));

    for (unsigned int i = 0; i < 4; ++i) {
      SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertices[i]->get_data());
      vertex_data->set_node(node);
      ArcSPtr arc = create_arc(vertices[i]);
      skeleton_->add_arc(arc);

      for (EdgeWPtr ew : vertices[i]->edges()) {
        if (EdgeSPtr e = ew.lock()) {
          // could abuse the fact that new edges have no data or sheet yet, but this is clearer
          if (std::find(edges.begin(), edges.end(), e) == edges.end()) {
            Hds_utils::get_sheet(e)->add_arc(arc); // other 2 are added when sheets are created below
          }
        }
      }
    }
    for (unsigned int i = 0; i < 4; ++i) {
      SheetSPtr sheet = create_sheet(edges[i]);
      skeleton_->add_sheet(sheet);
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
      post_op_vertices_pierce_.insert(poe->get_vertex_src());
      post_op_vertices_pierce_.insert(poe->get_vertex_dst());
    }
    CGAL_postcondition(post_op_vertices_pierce_.size() == 8);

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_pierce_event(const Pierce_event_sptr& event,
                                   const FT& current_time,
                                   const std::optional<FT>& time_future_bound,
                                   const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "########  Handle Pierce Event  #########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const FT& event_time = event->time();
    const Point_3& point = event->point();

    shit_to_event_time(polyhedron, current_time, event_time);

    VertexSPtr vertex = event->get_vertex();
    FacetSPtr facet = event->get_facet();

    CGAL_SS3_CORE_TRACE("V: " << event->get_vertex()->to_string());
    CGAL_SS3_CORE_TRACE("F: " << event->get_facet()->to_string());

#ifndef CGAL_SS3_NO_SKELETON_DS
    NodeSPtr node = Node::create();
    node->set_time(event_time);
    node->set_point(point);
    skeleton_->add_node(node);

    Hds_utils::get_arc(vertex)->close_arc(node);

    for (EdgeWPtr edge_w : vertex->edges()) {
      if (EdgeSPtr edge = edge_w.lock()) {
        Hds_utils::get_sheet(edge)->add_node(node);
      }
    }
#endif

    // the 3 new vertices cannot be reflex, but since we grow faces,
    // we need to check the other extremities of the edges
    for (EdgeWPtr ew : vertex->edges()) {
      EdgeSPtr edge = ew.lock();
      post_op_vertices_pierce_.insert(edge->other(vertex));
    }

    CGAL_postcondition(post_op_vertices_pierce_.size() == 3);


    std::array<FacetSPtr, 3> facets;
    std::array<EdgeSPtr, 3> edges;
    EdgeSPtr edge = vertex->first_edge();
    for (unsigned int i = 0; i < 3; ++i) {
      edges[i] = edge;
      if (edge->get_vertex_src() == vertex) {
        facets[i] = edge->get_facet_L();
      } else if (edge->get_vertex_dst() == vertex) {
        facets[i] = edge->get_facet_R();
      }
      edge = edge->next(vertex);
    }

    std::array<VertexSPtr, 3> vertices;
    for (unsigned int i = 0; i < 3; ++i) {
      vertices[i] = Vertex::create(point);
      Skeleton_vertex_data::create(vertices[i]);
      facet->add_vertex(vertices[i]);
      polyhedron->add_vertex(vertices[i]);
    }
    for (unsigned int i = 0; i < 3; ++i) {
      EdgeSPtr edge = edges[i];
      if (edge->get_vertex_src() == vertex) {
        edge->replace_vertex_src(vertices[i]);
      } else if (edge->get_vertex_dst() == vertex) {
        edge->replace_vertex_dst(vertices[i]);
      }
      facets[i]->remove_vertex(vertex);
      facets[i]->add_vertex(vertices[i]);
      facets[(i+2)%3]->add_vertex(vertices[i]);
    }
    vertex->facets().clear();
    vertex->edges().clear();
    polyhedron->remove_vertex(vertex);
    std::array<EdgeSPtr, 3> new_edges;
    for (unsigned int i = 0; i < 3; ++i) {
      new_edges[i] = Edge::create(vertices[i], vertices[(i+1)%3]);
      Skeleton_edge_data::create(new_edges[i]);
      new_edges[i]->set_facet_L(facet);
      new_edges[i]->set_facet_R(facets[i]);
      facet->add_edge(new_edges[i]);
      facets[i]->add_edge(new_edges[i]);
      polyhedron->add_edge(new_edges[i]);
    }

    if (time_future_bound.has_value()) {
      for (std::size_t i=0; i<3; ++i) {
        Hds_utils::set_final_point(vertices[i], std::nullopt);
      }
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    for (unsigned int i = 0; i < 3; ++i) {
      Hds_utils::set_node(vertices[i], node);
      ArcSPtr arc = create_arc(vertices[i]);
      skeleton_->add_arc(arc);
      Hds_utils::get_sheet(edges[i])->add_arc(arc); // other 2 sheets are the new edges', created below
    }
    for (unsigned int i = 0; i < 3; ++i) {
      SheetSPtr sheet = create_sheet(new_edges[i]);
      skeleton_->add_sheet(sheet);
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
      post_op_facets_.insert(new_edges[i]->get_facet_L());
      post_op_facets_.insert(new_edges[i]->get_facet_R());
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
    CGAL_postcondition(!Hds_utils::is_reflex(vertices[0]));
    CGAL_postcondition(!Hds_utils::is_reflex(vertices[1]));
    CGAL_postcondition(!Hds_utils::is_reflex(vertices[2]));

    add_event(event);

    return Event_status::EVENT_HANDLED;
  }

  Event_status handle_event(const Abstract_event_sptr& event,
                            const FT& current_time,
                            const std::optional<FT>& time_future_bound,
                            const PolyhedronSPtr& polyhedron)
  {
    Event_status result = Event_status::NON_EVENT;

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
    if (!is_actual_event(event, current_time, time_future_bound)) {
      return result;
    }

    if (event->getType() == Abstract_event::SAVE_EVENT) {
      result = handle_save_event(std::dynamic_pointer_cast<Save_event>(event),
                               current_time, polyhedron);
    } else if (event->getType() == Abstract_event::CONST_TIME_EVENT) {
      result = handle_const_time_event(std::dynamic_pointer_cast<Const_time_event>(event),
                                    current_time, polyhedron);
    } else if (event->getType() == Abstract_event::VANISH_EVENT) {
      result = handle_vanish_event(std::dynamic_pointer_cast<Vanish_event>(event),
                                 current_time, time_future_bound, polyhedron);
    } else if (event->getType() == Abstract_event::VERTEX_EVENT) {
      result = handle_vertex_event(std::dynamic_pointer_cast<Vertex_event>(event),
                                 current_time, time_future_bound, polyhedron);
    } else if (event->getType() == Abstract_event::FLIP_VERTEX_EVENT) {
      result = handle_flip_vertex_event(std::dynamic_pointer_cast<Flip_vertex_event>(event),
                                     current_time, time_future_bound, polyhedron);
    } else if (event->getType() == Abstract_event::SURFACE_EVENT) {
      result = handle_surface_event(std::dynamic_pointer_cast<Surface_event>(event),
                                  current_time, time_future_bound, polyhedron);
    } else if (event->getType() == Abstract_event::POLYHEDRON_SPLIT_EVENT) {
      result = handle_polyhedron_split_event(std::dynamic_pointer_cast<Polyhedron_split_event>(event),
                                          current_time, time_future_bound, polyhedron);
    } else if (event->getType() == Abstract_event::SPLIT_MERGE_EVENT) {
      result = handle_split_merge_event(std::dynamic_pointer_cast<Split_merge_event>(event),
                                     current_time, time_future_bound, polyhedron);
    } else if (event->getType() == Abstract_event::EDGE_SPLIT_EVENT) {
      result = handle_edge_split_event(std::dynamic_pointer_cast<Edge_split_event>(event),
                                    current_time, time_future_bound, polyhedron);
    } else if (event->getType() == Abstract_event::PIERCE_EVENT) {
      result = handle_pierce_event(std::dynamic_pointer_cast<Pierce_event>(event),
                                 current_time, time_future_bound, polyhedron);
    } else {
      CGAL_SS3_CORE_TRACE("Error: Cannot handle event of type " << event->getType());
      CGAL_assertion(false);
      result = Event_status::EVENT_NOT_HANDLED;
    }

    CGAL_SS3_CORE_TRACE_V(4, "-- Finished handling Event --");

    CGAL_postcondition_code(for (const VertexSPtr& v : polyhedron->vertices()))
    CGAL_postcondition(v->get_ID() != -1);
    CGAL_postcondition_code(for (const EdgeSPtr& e : polyhedron->edges()))
    CGAL_postcondition(e->get_ID() != -1);
    CGAL_postcondition_code(for (const FacetSPtr& f : polyhedron->facets()))
    CGAL_postcondition(f->get_ID() != -1);

#ifndef CGAL_SS3_NO_SKELETON_DS
    CGAL_postcondition_code(for (const VertexSPtr& v : polyhedron->vertices()) {)
    CGAL_postcondition(Hds_utils::get_node(v) != NodeSPtr());
    CGAL_postcondition(Hds_utils::get_arc(v) != ArcSPtr());
    CGAL_postcondition_code(})
    CGAL_postcondition_code(for (const EdgeSPtr& e : polyhedron->edges()))
    CGAL_postcondition(Hds_utils::get_sheet(e) != SheetSPtr());
#endif

    return result;
  }

  StraightSkeletonSPtr get_skeleton() const
  {
#ifndef CGAL_SS3_NO_SKELETON_DS
    CGAL_SS3_CORE_TRACE("Warning: no skeleton to return as it was not built");
#endif
    return this->skeleton_;
  }

  std::list<Abstract_event_sptr>& events()
  {
    return this->events_;
  }

  int count_events(int type) const
  {
    int result = 0;
    typename std::list<Abstract_event_sptr>::const_iterator it_e = events_.begin();
    while (it_e != events_.end()) {
      Abstract_event_sptr event = *it_e++;
      if (event->getType() == type) {
        result += 1;
      }
    }
    return result;
  }

  std::string events_summary() const
  {
    std::stringstream sstr;
    sstr << "Events: " << events_.size() << std::endl;
    sstr << "    ConstTimeEvents:       " << count_events(Abstract_event::CONST_TIME_EVENT) << std::endl;
    sstr << "    SaveEvents:            " << count_events(Abstract_event::SAVE_EVENT) << std::endl;
    sstr << "  VanishEvents:" << std::endl;
    sstr << "    Generic VanishEvents:  " << count_events(Abstract_event::VANISH_EVENT) << std::endl;
    sstr << "    EdgeEvents:            " << count_events(Abstract_event::EDGE_EVENT) << std::endl;
    sstr << "    EdgeMergeEvents:       " << count_events(Abstract_event::EDGE_MERGE_EVENT) << std::endl;
    sstr << "    TriangleEvents:        " << count_events(Abstract_event::TRIANGLE_EVENT) << std::endl;
    sstr << "    Dbl_edge_merge_events:    " << count_events(Abstract_event::DBL_EDGE_MERGE_EVENT) << std::endl;
    sstr << "    DblTriangleEvents:     " << count_events(Abstract_event::DBL_TRIANGLE_EVENT) << std::endl;
    sstr << "    TetrahedronEvents:     " << count_events(Abstract_event::TETRAHEDRON_EVENT) << std::endl;
    sstr << "  ContactEvents:" << std::endl;
    sstr << "    VertexEvents:          " << count_events(Abstract_event::VERTEX_EVENT) << std::endl;
    sstr << "    FlipVertexEvents:      " << count_events(Abstract_event::FLIP_VERTEX_EVENT) << std::endl;
    sstr << "    SurfaceEvents:         " << count_events(Abstract_event::SURFACE_EVENT) << std::endl;
    sstr << "    PolyhedronSplitEvents: " << count_events(Abstract_event::POLYHEDRON_SPLIT_EVENT) << std::endl;
    sstr << "    SplitMergeEvents:      " << count_events(Abstract_event::SPLIT_MERGE_EVENT) << std::endl;
    sstr << "    EdgeSplitEvents:       " << count_events(Abstract_event::EDGE_SPLIT_EVENT) << std::endl;
    sstr << "    PierceEvents:          " << count_events(Abstract_event::PIERCE_EVENT) << std::endl;
    sstr << ")" << std::endl;
    return sstr.str();

  }

private:
  PolyhedronSPtr polyhedron_;
  Abstract_vertex_splitter_sptr vertex_splitter_;
  int edge_event_;

  Base_mesh_offset_visitor<Traits>* visitor_ = nullptr;

  std::vector<FT> save_times_;
  std::filesystem::path save_path_;

  std::list<Abstract_event_sptr> events_;
  StraightSkeletonSPtr skeleton_;

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

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_STRAIGHT_SKELETON_BUILDER_3_H */
