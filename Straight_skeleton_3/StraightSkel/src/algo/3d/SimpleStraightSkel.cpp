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
 * @file   algo/3d/SimpleStraightSkel.cpp
 * @author Gernot Walzl
 * @date   2012-03-08
 */

// @fixme yesterday:

// @speed
// - avoid duplicate computations in check_bisectors
// - check_bisector with base time
// - delay crashAt bisector checks till pop time
// - re-enable fast vertex splitter for uniform weights around a vertex

// @fixme:
// - Enable caching + performance model ==> shared pointer is invalid: edge_tosplit
// - customer data w/ no facet merging had bugs

// @fixme later:
// - Fix simultaneous events still happening sometimes (likely the same event multiple times
//   since we don't check the queue before pushing)
// - Fix expired shared ptr errors in prints

// @fixme latest:
// - (a) Random perturbations could still create degenerate configurations if unlucky
// - (b) Random perturbations must be small enough as to not create self-intersections
// - EPECK -> EPICK could create self-intersections

// @todo:
// - must delay more combinatorial checks to pop time?

// @todo: cleaning
// - use CGAL kernel everywhere, get rid of the other kernel
// - get rid of the KernelWrapper
// - assert / CGAL assertion
// - indentation 4 --> 2 spaces
// - merge .h/.cpp into a single file, header only
// - fix all weak --> shared pointer conversions
// - harmonize face/facet
// - nullptr init?
// - check for overly shared objects, redundant function calls (plane normalization, for example)

// @todo later:
// - determine if outward or inward offset. Reject outward if there is no save-offset provided
//   or stop_on_last_save is false
// - Fix straight skeleton construction
// - Lighter data structures
// - implement lazy perturbations where we only perturb if we encounter an issue?
//   * Would cost more (detection of simultaneous events, etc.)
//   * How to detect "hidden" simultaneous events?
// - zero speed (check divisions)
// - debug ptr everywhere

// @todo latest:
// - factorize the three VV events but be very careful with the tiny differences
// - there should be a function "Facet::copyPropertiesAndData(otherFacet)"

// ----

/*
  As to not waste energy building the skeleton if we do not care about it
  The construction is also likely broken since the commit that made it so the polyhedron
  is not rebuilt from scratch at each and every iteration
*/
#define CGAL_SS3_NO_SKELETON_DS

// ----

/*
  Some events can be detected from multime elements. Reduce that to a single element.
  Should not be used when we are not refreshing the queue because when we add local elements
  we might have some representatives in the subset, but not the canonical one.
*/
#define CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS

/*
  Collect a single vanish event type per edge, and determine at queue pop time
  which type of vanish event it is.

  @fixme The point is that for e.g. DblTriangleEvent, the vanish event might change type
  without the event's representative edge being impacted... It'd probably be cheaper to handle the
  DblTriangleEvent particularity, and not have to compute and push a ton of useless vanish events...
*/
#define CGAL_SS3_USE_GENERIC_VANISH_EVENT

/*
  If we are not purging the queue at each iteration, we need to check if the event still exists
  at pop time AND we need not to purge at pop time because the event might become valid
  later in the future but the local updates to pierce events will not re-add it.
*/
#define CGAL_SS3_CHECK_PIERCE_AT_POP_TIME

/*
  Surface events have a filter to determine when an event is actually a split event.
  If we filter this, then it's a lot more work to go back after events and check
  if actually the surface event now is valid.
  Instead, filter at pop time.
*/
#define CGAL_SS3_CHECK_CONV_SPLIT_EVENT_AT_POP_TIME

/*
  Delay the degeneracy check until the event is relevant
*/
#define CGAL_SS3_CHECK_POLYHEDRON_SPLIT_EVENT_AT_POP_TIME

/*
  Delay the side of bisector till pop time
*/
// #define CGAL_SS3_CHECK_BISECTORS_AT_POP_TIME

/*
  Limit the visilibity of other vertices to connected component of the facet for VV events
*/
// #define CGAL_SS3_VV_VERTEX_2_WALK_FACES_FOR_DETECTION

// ----

// Debug

// #define CGAL_SS3_DEBUG_PRINT_QUEUE

// ----

#include "algo/3d/SimpleStraightSkel.h"

#include "debug.h"
#include "algo/Controller.h"
#include "algo/3d/KernelWrapper.h"
#include "algo/3d/LineInFacet.h"
#include "algo/3d/SelfIntersection.h"
#include "algo/3d/PolyhedronTransformation.h"
#include "algo/3d/AbstractVertexSplitter.h"
#include "algo/3d/AngleVertexSplitter.h"
#include "algo/3d/CombiVertexSplitter.h"
#include "algo/3d/ConvexVertexSplitter.h"
#include "algo/3d/VolumeVertexSplitter.h"
#include "algo/3d/WeightVertexSplitter.h"
#include "algo/3d/SphereVertexSplitter.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Triangle.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/skel/StraightSkeleton.h"
#include "data/3d/skel/AbstractEvent.h"
#include "data/3d/skel/ConstOffsetEvent.h"
#include "data/3d/skel/SaveOffsetEvent.h"
#include "data/3d/skel/VanishEvent.h"
#include "data/3d/skel/EdgeEvent.h"
#include "data/3d/skel/EdgeMergeEvent.h"
#include "data/3d/skel/TriangleEvent.h"
#include "data/3d/skel/DblEdgeMergeEvent.h"
#include "data/3d/skel/DblTriangleEvent.h"
#include "data/3d/skel/TetrahedronEvent.h"
#include "data/3d/skel/VertexEvent.h"
#include "data/3d/skel/FlipVertexEvent.h"
#include "data/3d/skel/SurfaceEvent.h"
#include "data/3d/skel/PolyhedronSplitEvent.h"
#include "data/3d/skel/SplitMergeEvent.h"
#include "data/3d/skel/EdgeSplitEvent.h"
#include "data/3d/skel/PierceEvent.h"
#include "data/3d/skel/Sheet.h"
#include "data/3d/skel/Arc.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelVertexData.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "data/3d/skel/SkelFacetData.h"
#include "db/3d/OBJFile.h"
#include "util/Configuration.h"
#include "util/Timer.h"
#include "util/StringFactory.h"

#include <CGAL/assertions.h>
#include <CGAL/Real_timer.h>

#include <limits>
#include <list>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace algo { namespace _3d {

SimpleStraightSkel::SimpleStraightSkel(PolyhedronSPtr polyhedron) {
    polyhedron_ = polyhedron;
    save_path_ = std::filesystem::current_path();
    skel_result_ = StraightSkeleton::create();
#ifndef CGAL_SS3_NO_SKELETON_DS
    skel_result_->setPolyhedron(polyhedron);
#endif
    initVertexSplitter();
    initEdgeEvent();
}

SimpleStraightSkel::SimpleStraightSkel(PolyhedronSPtr polyhedron, ControllerSPtr controller) {
    polyhedron_ = polyhedron;
    controller_ = controller;
    save_path_ = std::filesystem::current_path();
    skel_result_ = StraightSkeleton::create();
#ifndef CGAL_SS3_NO_SKELETON_DS
    skel_result_->setPolyhedron(polyhedron);
#endif
    initVertexSplitter();
    initEdgeEvent();
}

SimpleStraightSkel::SimpleStraightSkel(PolyhedronSPtr polyhedron,
                                       ControllerSPtr controller,
                                       const std::vector<CGAL::FT>& save_offsets,
                                       const std::filesystem::path& save_path)
    : polyhedron_(polyhedron),
      controller_(controller),
      save_offsets_(save_offsets), // intentional copy
      save_path_(save_path)
{
    std::sort(save_offsets_.begin(), save_offsets_.end(), [](const CGAL::FT& a, const CGAL::FT& b) { return CGAL::abs(a) < CGAL::abs(b); });

    skel_result_ = StraightSkeleton::create();
#ifndef CGAL_SS3_NO_SKELETON_DS
    skel_result_->setPolyhedron(polyhedron);
#endif
    initVertexSplitter();
    initEdgeEvent();
}

SimpleStraightSkel::~SimpleStraightSkel() {
    polyhedron_.reset();
    controller_.reset();
    vertex_splitter_.reset();
    skel_result_.reset();
}

SimpleStraightSkelSPtr SimpleStraightSkel::create(PolyhedronSPtr polyhedron) {
    return SimpleStraightSkelSPtr(new SimpleStraightSkel(polyhedron));
}

SimpleStraightSkelSPtr SimpleStraightSkel::create(PolyhedronSPtr polyhedron, ControllerSPtr controller) {
    return SimpleStraightSkelSPtr(new SimpleStraightSkel(polyhedron, controller));
}

SimpleStraightSkelSPtr SimpleStraightSkel::create(PolyhedronSPtr polyhedron,
                                                  ControllerSPtr controller,
                                                  const std::vector<CGAL::FT>& save_offsets,
                                                  const std::filesystem::path& save_path) {
    return SimpleStraightSkelSPtr(new SimpleStraightSkel(polyhedron, controller, save_offsets, save_path));
}

void SimpleStraightSkel::initVertexSplitter() {
    use_fast_vertex_splitter_ = false;

    util::ConfigurationSPtr config = util::Configuration::getInstance();
    std::string s_vertex_splitter;
    if (config->isLoaded()) {
        s_vertex_splitter = config->getString(
                "algo_3d_SimpleStraightSkel", "vertex_splitter");
        if (s_vertex_splitter.compare("AngleVertexSplitter") == 0) {
            vertex_splitter_ = AngleVertexSplitter::create();
        } else if (s_vertex_splitter.compare("CombiVertexSplitter") == 0) {
            vertex_splitter_ = CombiVertexSplitter::create();
        } else if (s_vertex_splitter.compare("ConvexVertexSplitter") == 0) {
            vertex_splitter_ = ConvexVertexSplitter::create();
        } else if (s_vertex_splitter.compare("VolumeVertexSplitter") == 0) {
            vertex_splitter_ = VolumeVertexSplitter::create();
        } else if (s_vertex_splitter.compare("WeightVertexSplitter") == 0) {
            vertex_splitter_ = WeightVertexSplitter::create(controller_);
        } else if (s_vertex_splitter.compare("SphereVertexSplitter") == 0) {
            vertex_splitter_ = SphereVertexSplitter::create(controller_);
        } else {
            CGAL_SS3_SPLITTER_TRACE("Warning: '" << s_vertex_splitter << "' not found.");
            CGAL_SS3_SPLITTER_TRACE("Using 'ConvexVertexSplitter'.");
            vertex_splitter_ = ConvexVertexSplitter::create();
        }
    } else {
        vertex_splitter_ = ConvexVertexSplitter::create();
    }
    skel_result_->appendConfig("vertex_splitter="+vertex_splitter_->toString()+"; ");
}

void SimpleStraightSkel::initEdgeEvent() {
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    std::string s_edge_event;
    if (config->isLoaded()) {
        s_edge_event = util::Configuration::getInstance()->getString(
                "algo_3d_SimpleStraightSkel", "edge_event");
        if (s_edge_event.compare("convex") == 0) {
            edge_event_ = 0;
        } else if (s_edge_event.compare("reflex") == 0) {
            edge_event_ = 1;
        } else if (s_edge_event.compare("flip") == 0) {
            edge_event_ = 2;
        } else if (s_edge_event.compare("sphere") == 0) {
            edge_event_ = 3;
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
    skel_result_->appendConfig("edge_event="+s_edge_event+"; ");
}

bool SimpleStraightSkel::isReflex(EdgeSPtr edge,
                                  const bool future_facing) {
    if (edge->getReflexStatus()) {
        // std::cout << "yes cache (high)" << std::endl;
        return *(edge->getReflexStatus());
    }

    bool result = false;

    VertexSPtr vertex_src = edge->getVertexSrc();
    VertexSPtr vertex_dst = edge->getVertexDst();

    // @todo is degenerate edge? pointer comparison or actual position comparison?
    if (*(vertex_src->getPoint()) == *(vertex_dst->getPoint())) {

        // std::cout << "no cache [degen]" << std::endl;

        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        FacetSPtr facet_src = edge->getFacetSrc();
        FacetSPtr facet_dst = edge->getFacetDst();

        const CGAL::FT& speed_l = std::dynamic_pointer_cast<SkelFacetData>(facet_l->getData())->getSpeed();
        const CGAL::FT& speed_r = std::dynamic_pointer_cast<SkelFacetData>(facet_r->getData())->getSpeed();
        const CGAL::FT& speed_src = std::dynamic_pointer_cast<SkelFacetData>(facet_src->getData())->getSpeed();
        const CGAL::FT& speed_dst = std::dynamic_pointer_cast<SkelFacetData>(facet_dst->getData())->getSpeed();

        CGAL::FT od = future_facing ? -1 : 1;
        Plane3SPtr offset_plane_l = KernelWrapper::offsetPlane(facet_l->getPlane(), od * speed_l);
        Plane3SPtr offset_plane_r = KernelWrapper::offsetPlane(facet_r->getPlane(), od * speed_r);
        Plane3SPtr offset_plane_src = KernelWrapper::offsetPlane(facet_src->getPlane(), od * speed_src);
        Plane3SPtr offset_plane_dst = KernelWrapper::offsetPlane(facet_dst->getPlane(), od * speed_dst);

        Point3SPtr p_src = KernelWrapper::intersection(offset_plane_src, offset_plane_l, offset_plane_r);
        Point3SPtr p_dst = KernelWrapper::intersection(offset_plane_dst, offset_plane_l, offset_plane_r);
        CGAL_assertion(p_src && p_dst);

        Vector3SPtr v_dir = KernelFactory::createVector3((*p_dst) - (*p_src));
        CGAL_assertion(*v_dir != CGAL::NULL_VECTOR);

        Vector3SPtr normal_l = KernelFactory::createVector3(offset_plane_l);
        CGAL_assertion(*normal_l != CGAL::NULL_VECTOR);

        Vector3SPtr v_cross = KernelWrapper::cross(normal_l, v_dir);
        Point3SPtr p = KernelFactory::createPoint3((*p_src) + (*v_cross));
        if (KernelWrapper::side(offset_plane_r, p) > 0) {
            result = true;
        }
    } else {
        result = edge->isReflex();
    }
    return result;
}

bool SimpleStraightSkel::isReflex(VertexSPtr vertex) {
    if (vertex->degree() == 0) {
        return false;
    }
    bool result = true;
    for (EdgeWPtr edge_wptr : vertex->edges()) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
            if (!isReflex(edge)) {
                result = false;
                break;
            }
        }
    }
    return result;
}

bool SimpleStraightSkel::isConvex(VertexSPtr vertex) {
    if (vertex->degree() == 0) {
        return false;
    }
    bool result = true;
    for (EdgeWPtr edge_wptr : vertex->edges()) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
            if (isReflex(edge)) {
                result = false;
                break;
            }
        }
    }
    return result;
}


Line3SPtr SimpleStraightSkel::line(EdgeSPtr edge) {
    Line3SPtr result = Line3SPtr();
    VertexSPtr vertex_src = edge->getVertexSrc();
    VertexSPtr vertex_dst = edge->getVertexDst();
    if (*(vertex_src->getPoint()) == *(vertex_dst->getPoint())) {
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        FacetSPtr facet_src = edge->getFacetSrc();
        FacetSPtr facet_dst = edge->getFacetDst();

        const CGAL::FT& speed_l = std::dynamic_pointer_cast<SkelFacetData>(facet_l->getData())->getSpeed();
        const CGAL::FT& speed_r = std::dynamic_pointer_cast<SkelFacetData>(facet_r->getData())->getSpeed();
        const CGAL::FT& speed_src = std::dynamic_pointer_cast<SkelFacetData>(facet_src->getData())->getSpeed();
        const CGAL::FT& speed_dst = std::dynamic_pointer_cast<SkelFacetData>(facet_dst->getData())->getSpeed();

        Plane3SPtr plane_l = facet_l->getPlane();
        Plane3SPtr plane_r = facet_r->getPlane();

        Plane3SPtr offset_plane_l = KernelWrapper::offsetPlane(plane_l, -speed_l);
        Plane3SPtr offset_plane_r = KernelWrapper::offsetPlane(plane_r, -speed_r);
        Plane3SPtr offset_plane_src = KernelWrapper::offsetPlane(facet_src->getPlane(), -speed_src);
        Plane3SPtr offset_plane_dst = KernelWrapper::offsetPlane(facet_dst->getPlane(), -speed_dst);

        Point3SPtr p_src = KernelWrapper::intersection(offset_plane_src,
                offset_plane_l, offset_plane_r);
        Point3SPtr p_dst = KernelWrapper::intersection(offset_plane_dst,
                offset_plane_l, offset_plane_r);
        Vector3SPtr v_dir = KernelFactory::createVector3((*p_dst) - (*p_src));
        result = KernelFactory::createLine3(vertex_src->getPoint(), v_dir);
    } else {
        result = edge->line();
    }
    return result;
}

Point3SPtr SimpleStraightSkel::getFinalPoint(VertexSPtr vertex,
                                             const CGAL::FT& offset_future_bound) {
    if (!vertex->hasFinalPoint()) {
        PolyhedronTransformation::resetFinalPoint(vertex, offset_future_bound);
    }

    return vertex->getFinalPoint();
}

Plane3SPtr SimpleStraightSkel::getFinalPlane(FacetSPtr facet,
                                             const CGAL::FT& offset_future_bound) {
    if (!facet->hasFinalPlane()) {
        PolyhedronTransformation::resetFinalPlane(facet, offset_future_bound);
    }

    return facet->getFinalPlane();
}

bool SimpleStraightSkel::savePolyhedron(PolyhedronSPtr polyhedron,
                                        const CGAL::FT& current_offset,
                                        const bool attempt_untilting)
{
    bool result = true;

    // attempt naive un-tilting
    if (attempt_untilting) {
        CGAL_SS3_CORE_TRACE_V(4, "Attempting to un-tilt polyhedron...");

#ifdef CGAL_SS3_DUMP_FILES
        db::_3d::OBJFile::save("results/pre-untilt_attempt.obj", polyhedron,
                               false /*do_triangulate*/,
                               true /*convert_to_double*/);
#endif

        PolyhedronSPtr polyhedron_cpy = polyhedron->clone();

        for (FacetSPtr facet : polyhedron_cpy->facets()) {
            // this assumes we have perturbed at the start ('0')
            facet->restorePlaneCoefficients(0, current_offset);
        }

        // As to avoid having a vertex not be defined by e.g. 3 planes with 2 being equal post un-tilt
        // That vertex is then useless, so just remove it
        // Do it here, before we recompute point positions
        db::_3d::AbstractFile::mergeCoplanarFacets(polyhedron_cpy, 0.0);

#ifdef CGAL_SS3_DUMP_FILES
        db::_3d::OBJFile::save("results/restored_merged.obj", polyhedron_cpy,
                               false /*do_triangulate*/,
                               true /*convert_to_double*/);
#endif

        CGAL_assertion(polyhedron_cpy && polyhedron_cpy->isConsistent());

        db::_3d::AbstractFile::sanitize(polyhedron_cpy);

#ifdef CGAL_SS3_DUMP_FILES
        db::_3d::OBJFile::save("results/restored_final.obj", polyhedron_cpy,
                               false /*do_triangulate*/,
                               true /*convert_to_double*/);
#endif

        PolyhedronTransformation::resetPoints(polyhedron_cpy);
        if(!polyhedron_cpy || !polyhedron_cpy->isConsistent()) {
            CGAL_SS3_CORE_TRACE("Warning: failed to un-tilt polyhedron");
            result = false;
        } else {
            polyhedron = polyhedron_cpy;
            if (SelfIntersection::hasSelfIntersectingSurface(polyhedron)) {
                CGAL_SS3_CORE_TRACE("Warning: self-intersections after un-tilting");
            } else {
                CGAL_SS3_CORE_TRACE("Successfully un-tilted polyhedron");
            }
            result = true;
        }
    }

    std::stringstream ss_filename, ss_filename_triangulated, ss_filename_exact;
    ss_filename << save_path_.string() << "/offset_" << current_offset << ".obj";
    ss_filename_triangulated << save_path_.string() << "/offset_" << current_offset << "_triangulated.obj";
    ss_filename_exact << save_path_.string() << "/offset_" << current_offset << "_exact.obj";

    result = (db::_3d::OBJFile::save(ss_filename.str(), polyhedron,
                                     false /*do_triangulate*/,
                                     true /*convert_to_double*/) && result);
    result = (db::_3d::OBJFile::save(ss_filename_triangulated.str(), polyhedron,
                                     true /*do_triangulate*/,
                                     true /*convert_to_double*/) && result);
    result = (db::_3d::OBJFile::save(ss_filename_exact.str(), polyhedron,
                                     do_triangulate,
                                     false /*convert_to_double*/) && result);

    return result;
}

bool SimpleStraightSkel::run() {
    if (controller_) {
        controller_->wait();
        controller_->setDispPolyhedron(polyhedron_);
        controller_->setDispSkel3d(skel_result_);
    }

    CGAL_SS3_CORE_TRACE_V(1, "== Straight Skeleton 3D started ==");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    PolyhedronSPtr polyhedron = polyhedron_->clone();
    CGAL_assertion(polyhedron && polyhedron->isConsistent());

#ifdef CGAL_SS3_DUMP_FILES
    db::_3d::OBJFile::save("results/input.obj", polyhedron, false /*do not triangulate*/);
#endif

    CGAL_SS3_CORE_TRACE_V(1, polyhedron->vertices().size() << " NV " << polyhedron->facets().size() << " NF");

    CGAL_assertion(algo::_3d::PolyhedronTransformation::doAll2PlanesIntersect(polyhedron_));
    CGAL_assertion(algo::_3d::PolyhedronTransformation::doAll3PlanesIntersect(polyhedron_));
    CGAL_assertion(!algo::_3d::SelfIntersection::hasSelfIntersectingSurface(polyhedron_));

    // store base plane coefficients
    cacheBasePlanes(polyhedron);

// @tmp some hardcoded weights for specific inputs
// #define CGAL_SS3_ACUTE_WEIGHTS
// #define CGAL_SS3_MERGING_WEIGHTS
// #define CGAL_SS3_PERFORMANCE_WEIGHTS
#if defined(CGAL_SS3_ACUTE_WEIGHTS) || defined(CGAL_SS3_MERGING_WEIGHTS) || defined(CGAL_SS3_PERFORMANCE_WEIGHTS)
# ifdef CGAL_SS3_ACUTE_WEIGHTS
    const CGAL::FT x_speed = 20;
    const CGAL::FT y_speed = 20;
    const CGAL::FT z_speed = 20;
    const CGAL::FT other_speed = 18.7939;
# elif defined(CGAL_SS3_MERGING_WEIGHTS)
    const CGAL::FT x_speed = 20;
    const CGAL::FT y_speed = 20;
    const CGAL::FT z_speed = 20;
    const CGAL::FT other_speed = 19.8777;
# elif defined(CGAL_SS3_PERFORMANCE_WEIGHTS)
    const CGAL::FT x_speed = 5;
    const CGAL::FT y_speed = 5;
    const CGAL::FT z_speed = 2;
    const CGAL::FT other_speed = 5;
# else
#  error
# endif

    for (FacetSPtr facet : polyhedron->facets()) {
        CGAL::FT speed = other_speed;
        const auto pl = facet->getPlane();
        const auto normal = KernelFactory::createVector3(pl);
        // DEBUG_PRINT("SP X " << CGAL::scalar_product(*normal, Vector3(1,0,0)));
        // DEBUG_PRINT("SP Y " << CGAL::scalar_product(*normal, Vector3(0,1,0)));
        // DEBUG_PRINT("SP Z " << CGAL::scalar_product(*normal, Vector3(0,0,1)));
        if(CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector3(1,0,0))) - 1) < 1e-3)
          speed = x_speed;
        if(CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector3(0,1,0))) - 1) < 1e-3)
          speed = y_speed;
        if(CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector3(0,0,1))) - 1) < 1e-3)
          speed = z_speed;

        SkelFacetDataSPtr data = SkelFacetData::create(facet);
        data->setSpeed(speed);
        // DEBUG_PRINT("speed to " << speed);
    }
#endif

    if (init(polyhedron, vertex_splitter_, use_fast_vertex_splitter_, controller_)) {
        if (controller_) {
            controller_->wait();
        }

        // if we stop immediately after the last save event, there is no point investigating events
        // that are farther away
        std::optional<CGAL::FT> offset_future_bound;
        if (!save_offsets_.empty()) {
            util::ConfigurationSPtr config = util::Configuration::getInstance();
            if (config->isLoaded()) {
                if ((config->contains("main", "stop_after_last_save_event") &&
                    config->getBool("main", "stop_after_last_save_event"))) {
                    offset_future_bound = save_offsets_.back();
                }
            }
        }

        step_id_ = -1;
        CGAL::FT current_offset = 0;
        CGAL::FT upcoming_offset;

        CGAL_assertion_code(const bool is_emptiness_expected = save_offsets_.empty();)

        PQ queue;
        collectEvents(polyhedron, current_offset, offset_future_bound, queue);

        for(;;) {
            ++step_id_;

            CGAL_SS3_CORE_TRACE_V(2, "\n=========== ITERATION #" << step_id_ << " AT OFFSET " << current_offset);
            CGAL_SS3_CORE_TRACE_V(2, polyhedron->vertices().size() << " NV " << polyhedron->facets().size() << " NF");

            if (visitor_) {
                if (!visitor_->go_further(step_id_, polyhedron, current_offset)) {
                    CGAL_SS3_CORE_TRACE_V(2, "Stopping on visitor request");
                    break;
                }
            }

            CGAL_assertion_code(for (FacetSPtr facet : polyhedron->facets()) {)
            CGAL_assertion(facet->getPlane()->a() == facet->getBasePlane()->a());
            CGAL_assertion(facet->getPlane()->b() == facet->getBasePlane()->b());
            CGAL_assertion(facet->getPlane()->c() == facet->getBasePlane()->c());
            CGAL_assertion_code(CGAL::FT speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();)
            CGAL_assertion(facet->getPlane()->d() == facet->getBasePlane()->d() - speed * current_offset);
            CGAL_assertion_code(})

            AbstractEventSPtr event = nextEvent(queue, current_offset);
            if (!event) {
                CGAL_SS3_CORE_TRACE_V(2, "No more events to treat");
                break;
            }

            CGAL_SS3_CORE_TRACE_V(2, "popped E" << event->getID() << " Type [" << event->getType() << "]");

            static int event_id = -1;
            CGAL_SS3_CORE_TRACE_V(2, "--> Accepted event #" << ++event_id << " " << event->toString() << " --");

            upcoming_offset = event->getOffset();

            // the next event should be at a time that is further away than the current one
            bool simultaneousEvents = (current_offset == upcoming_offset);
            if (simultaneousEvents) {
                CGAL_SS3_CORE_TRACE("Warning: there should not be any simultaneous events");
                // @fixme there are still some configurations where we have non canonicality
                // return false;
            }

#ifdef CGAL_SS3_RUN_TIMERS
            CGAL_SS3_CORE_TRACE_V(2, "current elapsed time: " << timer.time());
#endif

            if (visitor_) {
                visitor_->before_offset_event(polyhedron, current_offset, event);
            }

            if (controller_) {
                controller_->wait();
            }

            CGAL_assertion_code(Point3SPtr p_box_min = PolyhedronTransformation::boundingBoxMin(polyhedron);)
            CGAL_assertion_code(Point3SPtr p_box_max = PolyhedronTransformation::boundingBoxMax(polyhedron);)

            // Here is the event treatment: if the event is valid,
            // shift the polyhedron + apply the combinatorial changes
            EventStatus es = handleEvent(event, current_offset, offset_future_bound, polyhedron);
            CGAL_assertion(es != EventStatus::EVENT_NOT_HANDLED);
            if (es == EventStatus::NON_EVENT) {
                continue;
            }

            current_offset = upcoming_offset;

#ifdef CGAL_SS3_DUMP_FILES
            db::_3d::OBJFile::save("results/event_" + std::to_string(event_id) + ".obj", polyhedron, false /*do triangulate*/);
            db::_3d::OBJFile::save("results/event_" + std::to_string(event_id) + "_triangulated.obj", polyhedron);
#endif

            if(visitor_ && event->getType() == AbstractEvent::SAVE_OFFSET_EVENT) {
                visitor_->on_save_offset_event(polyhedron, current_offset);
            }

            CGAL_assertion(polyhedron->isConsistent());
#ifndef CGAL_SS3_NO_SKELETON_DS
            CGAL_assertion(skel_result_->isConsistent());
#endif
            CGAL_assertion(p_box_min && p_box_max);

            // @fixme for outward offsets without an enclosing bbox, this is wrong
            CGAL_warning(PolyhedronTransformation::isInsideBox(polyhedron, p_box_min, p_box_max));

            // this is tempting, but the mesh is usually not in a nice state here
            // CGAL_assertion(!SelfIntersection::hasSelfIntersectingSurface(polyhedron));

            if (controller_) {
              controller_->wait();
            }

            if (visitor_) {
                visitor_->after_offset_event(polyhedron, current_offset);
            }

            // If we are interested in a specific offset, there is no point going further
            if (event->getType() == AbstractEvent::SAVE_OFFSET_EVENT && save_offsets_.empty()) {
                util::ConfigurationSPtr config = util::Configuration::getInstance();
                if (config->isLoaded() &&
                    config->contains("main", "stop_after_last_save_event") &&
                    config->getBool("main", "stop_after_last_save_event")) {
                    break;
                }
            }

            // Update the event priority queue
            collectLocalEvents(polyhedron, current_offset, offset_future_bound, queue);

            post_op_vertices_.clear();
            post_op_edges_.clear();
            post_op_facets_.clear();

            post_op_vertices_VV_.clear();
            post_op_vertices_pierce_.clear();
            post_op_edges_edgesplit_.clear();
        }

        CGAL_SS3_CORE_TRACE_V(1, "== Straight Skeleton 3D finished ==");

        CGAL_warning(!is_emptiness_expected || polyhedron->empty());

#ifdef CGAL_SS3_RUN_TIMERS
        timer.stop();
        skel_result_->appendDescription("time=" + std::to_string(timer.time()) + "; ");
#endif

#ifndef CGAL_SS3_NO_SKELETON_DS
        skel_result_->appendDescription("controller=" +
                util::StringFactory::fromBoolean(controller_) + "; ");
#endif
        CGAL_SS3_CORE_TRACE_V(2, skel_result_->toString());
    } else {
        CGAL_SS3_CORE_TRACE_V(1, "Error: failed to initialize polyhedron");
        return false;
    }

    return true;
}

ThreadSPtr SimpleStraightSkel::startThread() {
    return ThreadSPtr(new std::thread(
            std::bind(&SimpleStraightSkel::run, this)));
}


NodeSPtr SimpleStraightSkel::createNode(VertexSPtr vertex) {
    NodeSPtr result = NodeSPtr();
    SkelVertexDataSPtr data;
    if (vertex->hasData()) {
        data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
    } else {
        data = SkelVertexData::create(vertex);
    }
    result = Node::create(vertex->getPoint());
    data->setNode(result);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
}

ArcSPtr SimpleStraightSkel::createArc(VertexSPtr vertex) {
    ArcSPtr result = ArcSPtr();
    if (vertex->degree() == 3) {
        SkelVertexDataSPtr data;
        if (vertex->hasData()) {
            data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
        } else {
            data = SkelVertexData::create(vertex);
        }
        FacetSPtr facets[3];
        for (unsigned int i = 0; i < 3; i++) {
            facets[i] = FacetSPtr();
        }
        unsigned int i = 0;
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (i < 3 && it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (FacetSPtr facet = facet_wptr.lock()) {
                facets[i] = facet;
                ++i;
            }
        }

        if (i >= 3) {
            Vector3SPtr direction;
            Plane3SPtr plane_1 = facets[0]->getPlane();
            Plane3SPtr plane_2 = facets[1]->getPlane();
            Plane3SPtr plane_3 = facets[2]->getPlane();

            CGAL::FT speed_1 = 1.0;
            if (facets[0]->hasData()) {
                speed_1 = std::dynamic_pointer_cast<SkelFacetData>(
                        facets[0]->getData())->getSpeed();
            }
            CGAL::FT speed_2 = 1.0;
            if (facets[1]->hasData()) {
                speed_2 = std::dynamic_pointer_cast<SkelFacetData>(
                        facets[1]->getData())->getSpeed();
            }
            CGAL::FT speed_3 = 1.0;
            if (facets[2]->hasData()) {
                speed_3 = std::dynamic_pointer_cast<SkelFacetData>(
                        facets[2]->getData())->getSpeed();
            }

#if 0 // supporting planes might not intersect in a point
            Point3SPtr src = vertex->getPoint();

            Plane3SPtr off_1 = KernelWrapper::offsetPlane(plane_1, -speed_1);
            Plane3SPtr off_2 = KernelWrapper::offsetPlane(plane_2, -speed_2);
            Plane3SPtr off_3 = KernelWrapper::offsetPlane(plane_3, -speed_3);
            Point3SPtr dst = KernelWrapper::intersection(off_1, off_2, off_3);
            if (src && dst) {
                direction = KernelFactory::createVector3(*dst - *src);
            }
#else
            Vector3SPtr n_1 = KernelFactory::createVector3(plane_1);
            Vector3SPtr n_2 = KernelFactory::createVector3(plane_2);
            Vector3SPtr n_3 = KernelFactory::createVector3(plane_3);

            // @fixme is division by "speed" correct?
            direction = KernelFactory::createVector3(speed_1 * (*n_1) + speed_2 * (*n_2) + speed_3 * (*n_3));
#endif

            if (direction) {
                result = Arc::create(data->getNode(), direction);
                data->setArc(result);

                for (EdgeWPtr edge_wptr : vertex->edges()) {
                    if (EdgeSPtr edge = edge_wptr.lock()) {
                        if (edge->hasData()) {
                            SkelEdgeDataSPtr edge_data =
                                    std::dynamic_pointer_cast<SkelEdgeData>(edge->getData());
                            SheetSPtr sheet = edge_data->getSheet();
                            if (sheet) {
                                sheet->addArc(result);
                            }
                        }
                    }
                }
            }
        }
    } else {
        CGAL_SS3_SKEL_DS_TRACE("Vertex has degree: " << vertex->degree());
    }
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
}

SheetSPtr SimpleStraightSkel::createSheet(EdgeSPtr edge) {
    SheetSPtr result = SheetSPtr();
    SkelEdgeDataSPtr data;
    if (edge->hasData()) {
        data = std::dynamic_pointer_cast<SkelEdgeData>(edge->getData());
    } else {
        data = SkelEdgeData::create(edge);
    }
    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    if (facet_l && facet_r) {
        Plane3SPtr plane_l = facet_l->getPlane();
        Plane3SPtr plane_r = facet_r->getPlane();
        FacetSPtr facet_b = facet_l;
        CGAL::FT speed_l = 1.0;
        if (facet_l->hasData()) {
            SkelFacetDataSPtr data_l = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_l->getData());
            facet_b = data_l->getFacetOrigin();
            speed_l = data_l->getSpeed();
        }
        FacetSPtr facet_f = facet_r;
        CGAL::FT speed_r = 1.0;
        if (facet_r->hasData()) {
            SkelFacetDataSPtr data_r = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_r->getData());
            facet_f = data_r->getFacetOrigin();
            speed_r = data_r->getSpeed();
        }
        Plane3SPtr plane_sheet;
        if (speed_l == speed_r) {
            if (isReflex(edge)) {
                plane_sheet = KernelWrapper::bisector(
                    KernelWrapper::opposite(plane_l), plane_r);
            } else {
                plane_sheet = KernelWrapper::bisector(
                    plane_l, KernelWrapper::opposite(plane_r));
            }
        } else {
            Line3SPtr line = KernelWrapper::intersection(plane_l, plane_r);
            Plane3SPtr offset_l = KernelWrapper::offsetPlane(plane_l, -speed_l);
            Plane3SPtr offset_r = KernelWrapper::offsetPlane(plane_r, -speed_r);
            Line3SPtr line_offset = KernelWrapper::intersection(offset_l, offset_r);
            Point3SPtr point_1 = KernelFactory::createPoint3(line->point());
            Vector3SPtr direction = KernelFactory::createVector3(line);
            Point3SPtr point_2 = KernelFactory::createPoint3(*point_1 + *direction);
            Point3SPtr point_3 = KernelFactory::createPoint3(line_offset->point());
            if (isReflex(edge)) {
                plane_sheet = KernelFactory::createPlane3(point_3, point_2, point_1);
            } else {
                plane_sheet = KernelFactory::createPlane3(point_1, point_2, point_3);
            }
        }

        result = Sheet::create();
        result->setPlane(plane_sheet);
        result->setFacetB(facet_b);
        result->setFacetF(facet_f);
        data->setSheet(result);
        SkelVertexDataSPtr data_src = std::dynamic_pointer_cast<SkelVertexData>(edge->getVertexSrc()->getData());
        SkelVertexDataSPtr data_dst = std::dynamic_pointer_cast<SkelVertexData>(edge->getVertexDst()->getData());
        if (data_src) {
            result->addNode(data_src->getNode());
            result->addArc(data_src->getArc());
        }
        if (data_dst) {
            result->addNode(data_dst->getNode());
            result->addArc(data_dst->getArc());
        }
    }
    return result;
}

void SimpleStraightSkel::cacheBasePlanes(PolyhedronSPtr polyhedron) {
    for (FacetSPtr facet : polyhedron->facets()) {
        if (!facet->hasData()) {
            SkelFacetData::create(facet);
        }

        Plane3SPtr base_plane = KernelFactory::createPlane3(*(facet->getPlane()));
        facet->setBasePlane(base_plane);
    }
}

bool SimpleStraightSkel::init(PolyhedronSPtr polyhedron,
                              AbstractVertexSplitterSPtr vertex_splitter,
                              const bool use_fast_vertex_splitter,
                              ControllerSPtr controller) {
    bool result = true;

    CGAL_SS3_CORE_TRACE("Input: " << polyhedron->vertices().size() << " NV " << polyhedron->facets().size() << " NF");

    for (VertexSPtr vertex : polyhedron->vertices()) {
        if (vertex->degree() < 3) {
            continue;
        }
        if (!vertex->hasData()) {
            SkelVertexData::create(vertex);
        }
#ifndef CGAL_SS3_NO_SKELETON_DS
        NodeSPtr node = createNode(vertex);
        if (node) {
            skel_result_->addNode(node);
        } else {
            result = false;
        }
#endif
    }

    CGAL_SS3_CORE_TRACE_V(2, "Using " << vertex_splitter->toString() << " to split vertices.");

    std::list<VertexSPtr> vertices_tosplit;
    for (VertexSPtr vertex : polyhedron->vertices()) {
        if (vertex->degree() > 3) {
            vertices_tosplit.push_back(vertex);
            SkelVertexDataSPtr data;
            if (vertex->hasData()) {
                data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
            } else {
                data = SkelVertexData::create(vertex);
            }
            data->setHighlight(true);
        }
    }

    CGAL_SS3_CORE_TRACE_V(2, vertices_tosplit.size() << " vertices to split");

    if (controller) {
        controller->setDispPolyhedron(polyhedron);
        controller->wait();
    }
    for (VertexSPtr vertex : vertices_tosplit) {
        vertex->getData()->setHighlight(false);
    }

    for (VertexSPtr vertex : vertices_tosplit) {
        bool equal_speeds = false;

        if (use_fast_vertex_splitter) {
            equal_speeds = true;
            CGAL::FT first_speed = 1.0;
            bool first_speed_set = false;

            for (FacetWPtr facet_wptr : vertex->facets()) {
                if (FacetSPtr facet = facet_wptr.lock()) {
                    const CGAL::FT& speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();

                    if (!first_speed_set) {
                        first_speed = speed;
                        first_speed_set = true;
                    } else if (speed != first_speed) {
                        equal_speeds = false;
                        break;
                    }
                }
            }
        }

        if (equal_speeds && vertex->isConvex()) {
            AbstractVertexSplitter::splitConvexVertex(vertex);
        } else if (equal_speeds && vertex->isReflex()) {
            AbstractVertexSplitter::splitReflexVertex(vertex);
        } else {
            CGAL_SS3_SPLITTER_TRACE("Generic split vertex:\n" << vertex->toString());

#ifdef CGAL_SS3_USE_COMBINATORIAL_SPLITTER_FOR_HIGH_DEGREE_VERTICES
            if (vertex->degree() > 15) {
                CGAL_SS3_SPLITTER_TRACE("Warning: degree is so high that even a combinatorial split will take forever.");
            }

            if (vertex->degree() > 10) {
#ifdef CGAL_SS3_RUN_TIMERS
                CGAL::Real_timer timer;
                timer.start();
#endif

                CGAL_SS3_SPLITTER_TRACE("High degree, use a combinatorial splitter");
                AbstractVertexSplitterSPtr combi_splitter = CombiVertexSplitter::create();
                combi_splitter->splitVertex(vertex);

#ifdef CGAL_SS3_RUN_TIMERS
                CGAL_SS3_SPLITTER_TRACE("Time taken to split vertex #" << vertex->getID() << " = " << timer.time());
#endif
            } else
#endif // CGAL_SS3_USE_COMBINATORIAL_SPLITTER_FOR_HIGH_DEGREE_VERTICES

                vertex_splitter->splitVertex(vertex);

            if (controller) {
                controller->setDispPolyhedron(polyhedron);
#ifndef CGAL_SS3_NO_SKELETON_DS
                controller->setDispSkel3d(skel_result_);
#endif
            }
        }
    }

#ifdef CGAL_SS3_DUMP_FILES
    db::_3d::OBJFile::save("results/split.obj", polyhedron, false /*do not triangulate*/);
#endif

    CGAL_postcondition_code(for(auto v : polyhedron->vertices()))
    CGAL_postcondition(v->getID() != -1);
    CGAL_postcondition_code(for(auto e : polyhedron->edges()))
    CGAL_postcondition(e->getID() != -1);
    CGAL_postcondition_code(for(auto f : polyhedron->facets()))
    CGAL_postcondition(f->getID() != -1);

#ifndef CGAL_SS3_NO_SKELETON_DS
    for (VertexSPtr vertex : polyhedron->vertices()) {
        if (vertex->degree() < 3) {
            continue;
        }
        ArcSPtr arc = createArc(vertex);
        if (arc) {
            skel_result_->addArc(arc);
        } else {
            result = false;
        }
    }

    for (EdgeSPtr edge : polyhedron->edges()) {
        if (!edge->hasData()) {
            SkelEdgeData::create(edge);
        }
        SheetSPtr sheet = createSheet(edge);
        if (sheet) {
            skel_result_->addSheet(sheet);
        } else {
            result = false;
        }
    }

    for (FacetSPtr facet : polyhedron->facets()) {
        if (!facet->hasData()) {
            SkelFacetData::create(facet);
        }
    }

    CGAL_postcondition(skel_result_->isConsistent());
#endif

    return result;
}

bool SimpleStraightSkel::isTriangle(FacetSPtr facet, EdgeSPtr edge_begin) {
    bool result = false;
    if (!facet->containsEdge(edge_begin)) {
        return false;
    }
    unsigned int i = 0;
//    VertexSPtr vertices[3];
    EdgeSPtr edge = edge_begin;
    for (i = 0; i < 3; i++) {
//        vertices[i] = edge->src(facet);
        edge = edge->next(facet);
        if (!edge) {
            break;
        }
    }
    if (i == 3 && edge_begin == edge) {
        result = true;

        // careful about ptr equality if below is uncommented
//        if (vertices[0]->getPoint() != vertices[1]->getPoint() &&
//                vertices[1]->getPoint() != vertices[2]->getPoint() &&
//                vertices[2]->getPoint() != vertices[0]->getPoint()) {
//            // check if triangle is a hole
//            Plane3SPtr plane = KernelFactory::createPlane3(
//                    vertices[0]->getPoint(),
//                    vertices[1]->getPoint(),
//                    vertices[2]->getPoint());
//            Vector3SPtr v1 = KernelFactory::createVector3(plane);
//            Vector3SPtr v2 = KernelFactory::createVector3(facet->getPlane());
//            if (((*v1) * (*v2)) > 0.0) {
//                // angle between orthogonal vectors of planes < CGAL_PI/2.0
//                // not a hole
//                result = true;
//            }
//        }
    }
    return result;
}

bool SimpleStraightSkel::isTetrahedron(EdgeSPtr edge_begin) {
    bool result = true;
    VertexSPtr vertices[4];
    for (unsigned int i = 0; i < 4; i++) {
        vertices[i] = VertexSPtr();
    }
    vertices[0] = edge_begin->getVertexSrc();
    vertices[1] = edge_begin->getVertexDst();
    EdgeSPtr edge_l = edge_begin->next(edge_begin->getFacetL());
    vertices[2] = edge_l->dst(edge_begin->getFacetL());
    EdgeSPtr edge_r = edge_begin->next(edge_begin->getFacetR());
    vertices[3] = edge_r->dst(edge_begin->getFacetR());
    for (unsigned int i = 0; i < 4; i++) {
        if (!vertices[i]) {
            result = false;
            return result;
        }
    }
    for (unsigned int i = 0; i < 4; i++) {
        for (unsigned int j = 1; j < 4; j++) {
            EdgeSPtr edge = vertices[i]->findEdge(vertices[(i+j)%4]);
            if (!edge) {
                result = false;
            }
        }
    }
    return result;
}

std::pair<Point3SPtr, CGAL::FT>
SimpleStraightSkel::intersectionPointAndTimeOffsetPlanes(FacetSPtr facet_0,
                                                         FacetSPtr facet_1,
                                                         FacetSPtr facet_2,
                                                         FacetSPtr facet_3,
                                                         const std::optional<CGAL::FT>& offset_past_bound,
                                                         const std::optional<CGAL::FT>& offset_future_bound)
{
    CGAL_SS3_DEBUG_SPTR(facet_0);
    CGAL_SS3_DEBUG_SPTR(facet_1);
    CGAL_SS3_DEBUG_SPTR(facet_2);
    CGAL_SS3_DEBUG_SPTR(facet_3);

    Plane3SPtr plane_0 = facet_0->getBasePlane();
    Plane3SPtr plane_1 = facet_1->getBasePlane();
    Plane3SPtr plane_2 = facet_2->getBasePlane();
    Plane3SPtr plane_3 = facet_3->getBasePlane();

    const CGAL::FT& speed_0 = std::dynamic_pointer_cast<SkelFacetData>(facet_0->getData())->getSpeed();
    const CGAL::FT& speed_1 = std::dynamic_pointer_cast<SkelFacetData>(facet_1->getData())->getSpeed();
    const CGAL::FT& speed_2 = std::dynamic_pointer_cast<SkelFacetData>(facet_2->getData())->getSpeed();
    const CGAL::FT& speed_3 = std::dynamic_pointer_cast<SkelFacetData>(facet_3->getData())->getSpeed();

    return KernelWrapper::intersectionPointAndTimeOffsetPlanes(plane_0, speed_0, plane_1, speed_1,
                                                               plane_2, speed_2, plane_3, speed_3,
                                                               offset_past_bound, offset_future_bound);
}

// @speed how about a first filter before Point&Time computation: IsEdgeGrowing
// if not, then there is definitely no intersection
std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::vanishesAt(EdgeSPtr edge,
                                                               const std::optional<CGAL::FT>& offset_past_bound,
                                                               const std::optional<CGAL::FT>& offset_future_bound)
{
    Point3SPtr point = Point3SPtr();
    CGAL::FT offset_event;

    CGAL_SS3_CORE_TRACE_V(16, "vanishesAt " << edge->toString());

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

    return intersectionPointAndTimeOffsetPlanes(facetL, facetP, facetR, facetN,
                                                offset_past_bound, offset_future_bound);
}

// returns 'true' if the point is on f's side of the bisector between f and f third
// edge: edge that is vanishing or crashing into another edge
// f: one of the faces incident to the edge
// t: the time of vanishing or the time of crash
// f_third: edge shared between 'f' and either the dst of the edge seen in f
//
// @todo this could be a predicate (oriented_side_of_event_point_wrt_bisectorC2)
bool SimpleStraightSkel::check_bisector(EdgeSPtr edge,
                                        FacetSPtr f,
                                        const CGAL::FT& t,
                                        FacetSPtr f_third, // @todo superfluous parameter
                                        Point3SPtr point)
{
    // @todo for speeds 0, when this function is called from crashAt, the t
    // is the event time, but maybe if f or f_third have speed 0, we could quickly
    // check like below?

    Plane3SPtr plane_third = f_third->getPlane();
    const CGAL::FT& speed_third = std::dynamic_pointer_cast<SkelFacetData>(f_third->getData())->getSpeed();
    if (is_zero(speed_third)) {
        return (KernelWrapper::side(plane_third, point) <= 0);
    }

    EdgeSPtr edge_f_f_third = edge->next(f);
    CGAL_precondition(edge_f_f_third->other(f) == f_third);

    // Determine which side of f-third is legal; that is determined by the angle
    // that the facet makes at the common vertex.
    //
    // We can't use the actual geometry of the edges because they might be degenerate,
    // so everything must be done with planes.
    //
    // @todo predicates

    FacetSPtr f_other = edge->other(f);
    Plane3SPtr plane_f = f->getPlane();
    Plane3SPtr plane_f_other = f_other->getPlane();
    Vector3SPtr n_f = KernelFactory::createVector3(plane_f);
    Vector3SPtr n_f_other = KernelFactory::createVector3(plane_f_other);
    Vector3SPtr n_f_third = KernelFactory::createVector3(plane_third);

    Vector3 v_f_f_other = isReflex(edge) ? CGAL::cross_product(*n_f_other, *n_f)
                                         : CGAL::cross_product(*n_f, *n_f_other);
    Vector3 v_f_f_third = isReflex(edge_f_f_third) ? CGAL::cross_product(*n_f_third, *n_f)
                                                   : CGAL::cross_product(*n_f, *n_f_third);

    Point3 pvp = *(edge->getVertexSrc()->getPoint()); // the position doesn't matter
    Point3 vp = pvp + v_f_f_other;
    Point3 nvp = vp + v_f_f_third;
    CGAL_assertion(pvp != vp && vp != nvp);

    auto n = CGAL::cross_product(nvp - vp, pvp - vp);
    CGAL::FT sp = CGAL::scalar_product(*n_f, n);

    CGAL_SS3_CORE_TRACE_V(128, "f = " << f->getID());
    CGAL_SS3_CORE_TRACE_V(128, "f_other = " << f_other->getID());
    CGAL_SS3_CORE_TRACE_V(128, "f_third = " << f_third->getID());
    CGAL_SS3_CORE_TRACE_V(128, "pv = " << pvp);
    CGAL_SS3_CORE_TRACE_V(128, "v = " << vp);
    CGAL_SS3_CORE_TRACE_V(128, "nv = " << nvp);
    CGAL_SS3_CORE_TRACE_V(128, "n_f = " << *n_f);
    CGAL_SS3_CORE_TRACE_V(128, "n_f_other = " << *n_f_other);
    CGAL_SS3_CORE_TRACE_V(128, "n_f_third = " << *n_f_third);
    CGAL_SS3_CORE_TRACE_V(128, "n = " << n);
    CGAL_SS3_CORE_TRACE_V(128, "sp = " << sp);

    // @todo below means that the edge is aligned, so the two bisectors are coplanar
    // and it'll be annoying. In that case, we can simply check on which side the
    // other vertex is and use the third bisector... but currently it's not possible
    // because no 2 faces are coplanar, so below should always be true!
    CGAL_assertion(sp != 0);

    // now, check on which side of the bisector we are
    const CGAL::FT& a_third = plane_third->a();
    const CGAL::FT& b_third = plane_third->b();
    const CGAL::FT& c_third = plane_third->c();
    const CGAL::FT& d_third = plane_third->d();
    CGAL::FT t_third = (a_third * point->x() + b_third * point->y() + c_third * point->z() + d_third) / speed_third;

    CGAL_SS3_CORE_TRACE_V(128, "edge = " << edge->toString());
    CGAL_SS3_CORE_TRACE_V(128, "edge_f_f_third = " << edge_f_f_third->toString());
    CGAL_SS3_CORE_TRACE_V(128, "time from third = " << t_third);
    CGAL_assertion(edge_f_f_third != edge);

    // we want SMALLER to mean "on f's side of the bisector"
    // and     LARGER  to mean "on f_third's side of the bisector"
    CGAL::Comparison_result f_third_point_side;
    if (isReflex(edge_f_f_third)) {
        CGAL_SS3_CORE_TRACE_V(128, "F1O: reflex edge ");
        f_third_point_side = CGAL::compare(t, t_third);
    } else {
        CGAL_SS3_CORE_TRACE_V(128, "F1O: convex edge ");
        f_third_point_side = CGAL::compare(t_third, t);
    }

    CGAL_SS3_CORE_TRACE_V(128, "f_third_point_side = " << f_third_point_side);

    // Finally, combine both info to know on which side we are
    if (sp < 0) {
        CGAL_SS3_CORE_TRACE_V(128, "SP: right turn");
        if (f_third_point_side == CGAL::SMALLER) {
            CGAL_SS3_CORE_TRACE_V(128, "reject!");
            return false;
        }
    } else {
        CGAL_SS3_CORE_TRACE_V(128, "SP: left turn");
        if (f_third_point_side == CGAL::LARGER) {
            CGAL_SS3_CORE_TRACE_V(128, "reject!");
            return false;
        }
    }

    return true;
}

// Check that the point is inside bounds
bool
SimpleStraightSkel::check_bisectors(EdgeSPtr edge_1,
                                    EdgeSPtr edge_2,
                                    Point3SPtr point,
                                    const CGAL::FT& t)
{
    FacetSPtr facet_l1 = edge_1->getFacetL();
    FacetSPtr facet_r1 = edge_1->getFacetR();
    FacetSPtr facet_l2 = edge_2->getFacetL();
    FacetSPtr facet_r2 = edge_2->getFacetR();

    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();
    FacetSPtr facet_2_src = edge_2->getFacetSrc();
    FacetSPtr facet_2_dst = edge_2->getFacetDst();

    CGAL_SS3_CORE_TRACE_V(16, "Facet 1 SRC = " << facet_1_src->getID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet 1 DST = " << facet_1_dst->getID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet 2 SRC = " << facet_2_src->getID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet 2 DST = " << facet_2_dst->getID());

    // We want the point to be left of the right arc, and right of the left arc
    //
    // We can just check the position of the query 'point' with respect
    // to one of the other bisector. A few tricky parts:
    // - The side of bisector is easy to know: it's comparing times, but that
    //   depends on whether the edge is reflex or convex
    // - Once we know on which side of the bisector it is, it's not done yet:
    //   we need to know which side of the bisector is the correct one and since
    //   faces can be non-convex polygons, we need to check if we are in a concave
    //   (within the facet) vertex to know if the clipping bisector is inverted

#ifndef CGAL_SS3_EXIT_ASAP
    bool reject_2b = false;
    bool reject_3b = false;
    bool reject_5b = false;
    bool reject_6b = false;
#endif

    if (!(facet_1_src == facet_l2 || facet_1_src == facet_r2)) {
        // src is the target of the edge when in the right facet
        if (!check_bisector(edge_1, facet_r1, t /*rt1*/, facet_1_src, point)) {
#ifdef CGAL_SS3_EXIT_ASAP
            return false;
#else
            reject_2b = true;
#endif
        }
    }

    if (!(facet_1_dst == facet_l2 || facet_1_dst == facet_r2)) {
        if (!check_bisector(edge_1, facet_l1, t /*lt1*/, facet_1_dst, point)) {
#ifdef CGAL_SS3_EXIT_ASAP
            return false;
#else
            reject_3b = true;
#endif
        }
    }

    if (!(facet_2_src == facet_l1 || facet_2_src == facet_r1)) {
        if (!check_bisector(edge_2, facet_r2, t /*rt2*/, facet_2_src, point)) {
#ifdef CGAL_SS3_EXIT_ASAP
            return false;
#else
            reject_5b = true;
#endif
        }
    }

    if (!(facet_2_dst == facet_l1 || facet_2_dst == facet_r1)) {
        if (!check_bisector(edge_2, facet_l2, t /*lt2*/, facet_2_dst, point)) {
#ifdef CGAL_SS3_EXIT_ASAP
            return false;
#else
            reject_6b = true;
#endif
        }
    }

#ifndef CGAL_SS3_EXIT_ASAP
    const bool reject = (reject_2b || reject_3b || reject_5b || reject_6b);
    if(reject) {
        return false;
    }
#endif

    return true;
}

// this function does an early exit if the result is irrelevant (in the past or too far in the future)
//
// @speed, should be able to not solve the system but just exit early if the 4 planes are clearly
// not intersecting (diametral spheres around the edges of size something?)
std::pair<Point3SPtr, CGAL::FT>
SimpleStraightSkel::crashAt(EdgeSPtr edge_1, EdgeSPtr edge_2,
                            const std::optional<CGAL::FT>& offset_past_bound,
                            const std::optional<CGAL::FT>& offset_future_bound)
{
    CGAL_SS3_CORE_TRACE_V(16, "-- Crash At\n    " << edge_1->toString() << "\n    " << edge_2->toString());

    FacetSPtr facet_l1 = edge_1->getFacetL();
    FacetSPtr facet_r1 = edge_1->getFacetR();
    FacetSPtr facet_l2 = edge_2->getFacetL();
    FacetSPtr facet_r2 = edge_2->getFacetR();

    CGAL_SS3_CORE_TRACE_V(16, "Facet L1 = " << facet_l1->getID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet R1 = " << facet_r1->getID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet L2 = " << facet_l2->getID());
    CGAL_SS3_CORE_TRACE_V(16, "Facet R2 = " << facet_r2->getID());

    // @speed we could delay point computation and orientation checks below to queue pop time,
    // but it probably does not gain much: 99.99% of the time, we compute an intersection time
    // and it is filtered, so the times where we compute a (lazy) point and use it is somewhat
    // negligible (currently).
    Point3SPtr point;
    CGAL::FT event_offset;
    std::tie(point, event_offset) = intersectionPointAndTimeOffsetPlanes(facet_l1, facet_r1, facet_l2, facet_r2,
                                                                         offset_past_bound, offset_future_bound);

    if (!point) {
        return { };
    }

    CGAL_SS3_CORE_TRACE_V(16, "Intersection: " << *point << " @ " << event_offset);

    const CGAL::FT& speed_l1 = std::dynamic_pointer_cast<SkelFacetData>(facet_l1->getData())->getSpeed();
    const CGAL::FT& speed_r1 = std::dynamic_pointer_cast<SkelFacetData>(facet_r1->getData())->getSpeed();
    const CGAL::FT& speed_l2 = std::dynamic_pointer_cast<SkelFacetData>(facet_l2->getData())->getSpeed();
    const CGAL::FT& speed_r2 = std::dynamic_pointer_cast<SkelFacetData>(facet_r2->getData())->getSpeed();
    CGAL_USE(speed_l1);
    CGAL_USE(speed_r1);
    CGAL_USE(speed_l2);
    CGAL_USE(speed_r2);

    CGAL_SS3_CORE_TRACE_CODE(CGAL::FT current_offset = (facet_l1->getBasePlane()->d()
                                                        - facet_l1->getPlane()->d()) / speed_l1);
    CGAL_SS3_CORE_TRACE_CODE(CGAL::FT shift_offset = event_offset - current_offset);
    CGAL_SS3_CORE_TRACE_V(16, "current offset " << current_offset);
    CGAL_SS3_CORE_TRACE_V(16, "shift offset " << shift_offset);
    CGAL_SS3_CORE_TRACE_CODE(Segment3SPtr offset_e1 = PolyhedronTransformation::shiftEdge(edge_1, shift_offset);)
    CGAL_SS3_CORE_TRACE_CODE(Segment3SPtr offset_e2 = PolyhedronTransformation::shiftEdge(edge_2, shift_offset);)
    CGAL_SS3_CORE_TRACE_V(16, "Offset edge 1: " << *offset_e1);
    CGAL_SS3_CORE_TRACE_V(16, "Offset edge 2: " << *offset_e2);

    Plane3SPtr plane_l1 = facet_l1->getPlane();
    Plane3SPtr plane_r1 = facet_r1->getPlane();
    Plane3SPtr plane_l2 = facet_l2->getPlane();
    Plane3SPtr plane_r2 = facet_r2->getPlane();

    // Since the algorithm shrinks the polyhedron, "not to be in the past" is equivalent to
    // "being on the negative side of the planes of the facets incident to the edge" (assuming positive weights)
    //
    // Since we have filtered positive times, this is already checked
    CGAL_assertion(!(KernelWrapper::side(plane_l1, point) > 0 ||
                     KernelWrapper::side(plane_r1, point) > 0));
    CGAL_assertion(!(KernelWrapper::side(plane_l2, point) > 0 ||
                     KernelWrapper::side(plane_r2, point) > 0));

    const CGAL::FT& l1a = plane_l1->a();
    const CGAL::FT& l1b = plane_l1->b();
    const CGAL::FT& l1c = plane_l1->c();
    const CGAL::FT& l1d = plane_l1->d();
    const CGAL::FT& r1a = plane_r1->a();
    const CGAL::FT& r1b = plane_r1->b();
    const CGAL::FT& r1c = plane_r1->c();
    const CGAL::FT& r1d = plane_r1->d();
    const CGAL::FT& l2a = plane_l2->a();
    const CGAL::FT& l2b = plane_l2->b();
    const CGAL::FT& l2c = plane_l2->c();
    const CGAL::FT& l2d = plane_l2->d();
    const CGAL::FT& r2a = plane_r2->a();
    const CGAL::FT& r2b = plane_r2->b();
    const CGAL::FT& r2c = plane_r2->c();
    const CGAL::FT& r2d = plane_r2->d();

    CGAL_SS3_CORE_TRACE_V(16, "time: " << event_offset);

    CGAL_assertion_code(CGAL::FT lt1 = (l1a * point->x() + l1b * point->y() + l1c * point->z() + l1d) / speed_l1;)
    CGAL_assertion_code(CGAL::FT rt1 = (r1a * point->x() + r1b * point->y() + r1c * point->z() + r1d) / speed_r1;)
    CGAL_assertion_code(CGAL::FT lt2 = (l2a * point->x() + l2b * point->y() + l2c * point->z() + l2d) / speed_l2;)
    CGAL_assertion_code(CGAL::FT rt2 = (r2a * point->x() + r2b * point->y() + r2c * point->z() + r2d) / speed_r2;)
    CGAL_assertion(lt1 == rt1 && lt1 == lt2 && lt1 == rt2);

#ifndef CGAL_SS3_CHECK_BISECTORS_AT_POP_TIME
    // @speed use base planes and event_offset
    CGAL::FT t = (l1a * point->x() + l1b * point->y() + l1c * point->z() + l1d) / speed_l1;
    if (!check_bisectors(edge_1, edge_2, point, t)) {
        return { };
    }
#endif // CGAL_SS3_CHECK_BISECTORS_AT_POP_TIME

    return { point, event_offset };
}

CGAL::FT SimpleStraightSkel::offsetDist(FacetSPtr facet, Point3SPtr point) {
    Plane3SPtr plane = facet->getPlane();
    CGAL::FT result = KernelWrapper::distance(plane, point);
    if (KernelWrapper::side(plane, point) < 0) {
        result *= -1.0;
    }
    if (facet->hasData()) {
        const CGAL::FT& speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
        result /= speed;
    }
    return result;
}

bool SimpleStraightSkel::isEventInThePast(const CGAL::FT& current_offset,
                                          AbstractEventSPtr event)
{
    CGAL_precondition(event->isValid());
    return event->getOffset() >= current_offset;
}

bool SimpleStraightSkel::isEventObsolete(AbstractEventSPtr event)
{
    CGAL_precondition(event->isValid());
    return event->isObsolete();
}

bool SimpleStraightSkel::isActualEvent(const CGAL::FT& current_offset,
                                       AbstractEventSPtr event,
                                       PolyhedronSPtr polyhedron)
{
    CGAL_precondition(event->isValid());

    bool result = true;

    if (event->getType() == AbstractEvent::VERTEX_EVENT) {
        result = isActualVertexEvent(std::dynamic_pointer_cast<VertexEvent>(event), polyhedron);
    } else if (event->getType() == AbstractEvent::FLIP_VERTEX_EVENT) {
        result = isActualFlipVertexEvent(std::dynamic_pointer_cast<FlipVertexEvent>(event), polyhedron);
    } else if (event->getType() == AbstractEvent::SURFACE_EVENT) {
        result = isActualSurfaceEvent(std::dynamic_pointer_cast<SurfaceEvent>(event), polyhedron);
    } else if (event->getType() == AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
        result = isActualPolyhedronSplitEvent(std::dynamic_pointer_cast<PolyhedronSplitEvent>(event), current_offset, polyhedron);
    } else if (event->getType() == AbstractEvent::SPLIT_MERGE_EVENT) {
        result = isActualSplitMergeEvent(std::dynamic_pointer_cast<SplitMergeEvent>(event), polyhedron);
    } else if (event->getType() == AbstractEvent::PIERCE_EVENT) {
        result = isActualPierceEvent(std::dynamic_pointer_cast<PierceEvent>(event), current_offset, polyhedron);
    }

    return result;
}

// @speed we could mark edges to say "vanish event already in the queue" (or for other events)
void SimpleStraightSkel::collectVanishEvents(const std::list<EdgeSPtr>& edges,
                                             PolyhedronSPtr polyhedron,
                                             const CGAL::FT& current_offset,
                                             const std::optional<CGAL::FT>& offset_future_bound,
                                             PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Vanish Events [" << current_offset << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    for (EdgeSPtr edge : edges) {
        CGAL_assertion(edge->getID() != -1);

        VertexSPtr vertex_src = edge->getVertexSrc();
        VertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->getPoint() == vertex_dst->getPoint()) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge, current_offset, offset_future_bound);
        if (!point) {
            continue;
        }

        CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

        NodeSPtr node = Node::create();
        VanishEventSPtr event = VanishEvent::create();
        event->setStepID(step_id_);
        event->setNode(node);
        node->clear();
        node->setOffset(offset_event);
        node->setPoint(point);
        event->setEdge(edge);

#ifndef CGAL_SS3_NO_SKELETON_DS
# error "todo"
#endif

        // CGAL_SS3_CORE_TRACE("Edge " << edge->getID() << " vanishes at " << *point << " @ " << offset_event);

        queue.push(event);
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Vanish Events in: " << timer.time());
#endif
}


void SimpleStraightSkel::collectVanishEvents(PolyhedronSPtr polyhedron,
                                             const CGAL::FT& current_offset,
                                             const std::optional<CGAL::FT>& offset_future_bound,
                                             PQ& queue)
{
    return collectVanishEvents(polyhedron->edges(), polyhedron, current_offset, offset_future_bound, queue);
}

void SimpleStraightSkel::collectEdgeEvents(const std::list<EdgeSPtr>& edges,
                                           PolyhedronSPtr polyhedron,
                                           const CGAL::FT& current_offset,
                                           const std::optional<CGAL::FT>& offset_future_bound,
                                           PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Edge Events [" << current_offset << "]");

    for (EdgeSPtr edge : edges) {
        CGAL_assertion(edge->getID() != -1);

        VertexSPtr vertex_src = edge->getVertexSrc();
        VertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->getPoint() == vertex_dst->getPoint()) {
            continue;
        }

        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (isTriangle(facet_l, edge) || isTriangle(facet_r, edge)) {
            // triangle event
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge, current_offset, offset_future_bound);
        if (!point) {
            continue;
        }

        CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

        FacetSPtr facet_src = edge->getFacetSrc();
        FacetSPtr facet_dst = edge->getFacetDst();

        // This does not work when there is more than one edge between both facets.
        // EdgeSPtr edge_2 = facet_src->findEdge(facet_dst);
        std::list<EdgeSPtr> edges_2 = facet_src->findEdges(facet_dst); // @todo shouldn't this check also happen in other events?...

        bool split_event = false;
        for (EdgeSPtr edge_2 : edges_2) {
#if defined(CGAL_SS3_OLD_CODE_BOUND_CHECKS) || defined(CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS)
            bool split_event_current_1 = true;
            bool split_event_current_2 = true;
            bool split_event_current_3 = true;

            SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(edge_2->getData());
            Vector3SPtr normal_2 = KernelFactory::createVector3(data_2->getSheet()->getPlane());
            Line3SPtr line_normal_2 = KernelFactory::createLine3(point, normal_2);
            if (KernelWrapper::orientation(line(edge_2), line_normal_2) < 0) {
                // out of bounded area
                split_event_current_1 = false;
            }
            SkelVertexDataSPtr data_2_src = std::dynamic_pointer_cast<SkelVertexData>(
                edge_2->getVertexSrc()->getData());
            ArcSPtr arc_2_src = data_2_src->getArc();
            if (KernelWrapper::orientation(arc_2_src->line(), line_normal_2) > 0) {
                // out of bounded area
                split_event_current_2 = false;
            }
            SkelVertexDataSPtr data_2_dst = std::dynamic_pointer_cast<SkelVertexData>(
                edge_2->getVertexDst()->getData());
            ArcSPtr arc_2_dst = data_2_dst->getArc();
            if (KernelWrapper::orientation(arc_2_dst->line(), line_normal_2) < 0) {
                // out of bounded area
                split_event_current_3 = false;
            }

            const bool split_event_current = (split_event_current_1 &&
                                              split_event_current_2 &&
                                              split_event_current_3);
#endif

#if !defined(CGAL_SS3_OLD_CODE_BOUND_CHECKS) || defined(CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS)
            bool split_event_current_1_b = true;
            bool split_event_current_2_b = true;
            bool split_event_current_3_b = true;

            FacetSPtr facet_l2 = edge_2->getFacetL();
            FacetSPtr facet_r2 = edge_2->getFacetR();
            FacetSPtr facet_2_src = edge_2->getFacetSrc();
            FacetSPtr facet_2_dst = edge_2->getFacetDst();

            Plane3SPtr plane_l2 = facet_l2->getPlane();

            const CGAL::FT& l2a = plane_l2->a();
            const CGAL::FT& l2b = plane_l2->b();
            const CGAL::FT& l2c = plane_l2->c();
            const CGAL::FT& l2d = plane_l2->d();
            const CGAL::FT& speed_l2 = std::dynamic_pointer_cast<SkelFacetData>(facet_l2->getData())->getSpeed();
            CGAL::FT t = (l2a * point->x() + l2b * point->y() + l2c * point->z() + l2d) / speed_l2;

            CGAL_assertion_code(Plane3SPtr plane_r2 = facet_r2->getPlane();)
            CGAL_assertion_code(const CGAL::FT& r2a = plane_r2->a();)
            CGAL_assertion_code(const CGAL::FT& r2b = plane_r2->b();)
            CGAL_assertion_code(const CGAL::FT& r2c = plane_r2->c();)
            CGAL_assertion_code(const CGAL::FT& r2d = plane_r2->d();)
            CGAL_assertion_code(const CGAL::FT& speed_r2 = std::dynamic_pointer_cast<SkelFacetData>(facet_r2->getData())->getSpeed();)
            CGAL_assertion_code(CGAL::FT lt2 = (l2a * point->x() + l2b * point->y() + l2c * point->z() + l2d) / speed_l2);
            CGAL_assertion_code(CGAL::FT rt2 = (r2a * point->x() + r2b * point->y() + r2c * point->z() + r2d) / speed_r2);
            CGAL_assertion(lt2 == rt2);

            if (is_positive(t)) {
#ifdef CGAL_SS3_EXIT_ASAP
                // can 'break' directly because it's the same value for all 'edge_2's
                break;
#else
                split_event_current_1_b = false;
#endif
            }

            if (!check_bisector(edge_2, facet_r2, t /*rt2*/, facet_2_src, point)) {
#ifdef CGAL_SS3_EXIT_ASAP
                continue;
#else
                split_event_current_2_b = false;
#endif
            }

            if (!check_bisector(edge_2, facet_l2, t /*lt2*/, facet_2_dst, point)) {
#ifdef CGAL_SS3_EXIT_ASAP
                continue;
#else
                split_event_current_3_b = false;
#endif
            }

            const bool split_event_current_b = (split_event_current_1_b &&
                                                split_event_current_2_b &&
                                                split_event_current_3_b);
#endif

#ifdef CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS
# ifdef CGAL_SS3_EXIT_ASAP
#  error "Cannot compare if doing early exits"
# endif
            CGAL_assertion(split_event_current_1 == split_event_current_1_b);
            CGAL_assertion(split_event_current_2 == split_event_current_2_b);
            CGAL_assertion(split_event_current_3 == split_event_current_3_b);

            CGAL_assertion(split_event_current == split_event_current_b);
#endif

            if (split_event_current_b) {
                split_event = true;
                break;
            }
        }
        if (split_event) {
            continue;
        }
        // edge merge event
        EdgeSPtr edge_prev = edge->prev(facet_l);
        EdgeSPtr edge_next = edge->next(facet_l)->next(facet_l);
        if (edge_prev->hasSameFacets(edge_next)) {
            continue;
        }
        edge_prev = edge->prev(facet_l)->prev(facet_l);
        edge_next = edge->next(facet_l);
        if (edge_prev->hasSameFacets(edge_next)) {
            continue;
        }
        edge_prev = edge->prev(facet_r);
        edge_next = edge->next(facet_r)->next(facet_r);
        if (edge_prev->hasSameFacets(edge_next)) {
            continue;
        }
        edge_prev = edge->prev(facet_r)->prev(facet_r);
        edge_next = edge->next(facet_r);
        if (edge_prev->hasSameFacets(edge_next)) {
            continue;
        }

        NodeSPtr node = Node::create();
        EdgeEventSPtr event = EdgeEvent::create();
        event->setStepID(step_id_);
        event->setNode(node);
        node->clear();
        node->setOffset(offset_event);
        node->setPoint(point);
        event->setEdge(edge);

#ifndef CGAL_SS3_NO_SKELETON_DS
        SkelVertexDataSPtr data_src = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexSrc()->getData());
        SkelVertexDataSPtr data_dst = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexDst()->getData());
        node->addArc(data_src->getArc());
        node->addArc(data_dst->getArc());
        SkelEdgeDataSPtr data_edge = std::dynamic_pointer_cast<SkelEdgeData>(
                edge->getData());
        node->addSheet(data_edge->getSheet());
#endif

        queue.push(event);
    }
}

void SimpleStraightSkel::collectEdgeEvents(PolyhedronSPtr polyhedron,
                                           const CGAL::FT& current_offset,
                                           const std::optional<CGAL::FT>& offset_future_bound,
                                           PQ& queue)
{
    return collectEdgeEvents(polyhedron->edges(), polyhedron,
                             current_offset, offset_future_bound, queue);
}

void SimpleStraightSkel::collectEdgeMergeEvents(const std::list<EdgeSPtr>& edges,
                                                PolyhedronSPtr polyhedron,
                                                const bool use_canonical_event_reps,
                                                const CGAL::FT& current_offset,
                                                const std::optional<CGAL::FT>& offset_future_bound,
                                                PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Edge Merge Events [" << current_offset << "]");

    for (EdgeSPtr edge : edges) {
        CGAL_assertion(edge->getID() != -1);

        VertexSPtr vertex_src = edge->getVertexSrc();
        VertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->getPoint() == vertex_dst->getPoint()) {
            continue;
        }

        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (isTriangle(facet_l, edge) || isTriangle(facet_r, edge)) {
            // triangle event
            continue;
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
            continue;
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
            continue;
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

#ifndef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        if (use_canonical_event_reps) {
            // @fixme not sure about all of this...
            // this is essentially the same as above but if 'edge' were 'edge->prev(facet_l)'
            // so when we pass by again with this edge, then we will meet it in the check just above.
            // BUT: the dbl edge merge event filter above applies only to 'edge' so we might filter
            // for edge->prev(facet_l) which we wouldn't have filtered when seeing it from 'edge'?

            edge_prev = edge->prev(facet_l)->prev(facet_l);
            edge_next = edge->next(facet_l);
            if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
                facet = facet_l;
                edge_1 = edge_prev;
                edge_2 = edge_next;
            }
        }
#endif
        edge_prev = edge->prev(facet_r);
        edge_next = edge->next(facet_r)->next(facet_r);
        if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
            facet = facet_r;
            edge_1 = edge_prev;
            edge_2 = edge_next;
        }
#ifndef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        // @fixme not sure about this - see above
        if (use_canonical_event_reps) {
            edge_prev = edge->prev(facet_r)->prev(facet_r);
            edge_next = edge->next(facet_r);
            if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
                facet = facet_r;
                edge_1 = edge_prev;
                edge_2 = edge_next;
            }
        }
#endif

        if (!(facet && edge_1 && edge_2)) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge, current_offset, offset_future_bound);
        if (!point) {
            continue;
        }

        CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

        NodeSPtr node = Node::create();
        EdgeMergeEventSPtr event = EdgeMergeEvent::create();
        event->setStepID(step_id_);
        event->setNode(node);
        node->clear();
        node->setOffset(offset_event);
        node->setPoint(point);
        event->setFacet(facet);
        event->setEdge1(edge_1);
        event->setEdge2(edge_2);

#ifndef CGAL_SS3_NO_SKELETON_DS
        EdgeSPtr edge_toremove_1 = edge_1->next(facet);
        EdgeSPtr edge_toremove_2 = edge_toremove_1->next(facet);
        SkelVertexDataSPtr data_vertex = std::dynamic_pointer_cast<SkelVertexData>(
                edge_toremove_1->src(facet)->getData());
        node->addArc(data_vertex->getArc());
        data_vertex = std::dynamic_pointer_cast<SkelVertexData>(
                edge_toremove_1->dst(facet)->getData());
        node->addArc(data_vertex->getArc());
        data_vertex = std::dynamic_pointer_cast<SkelVertexData>(
                edge_toremove_2->dst(facet)->getData());
        node->addArc(data_vertex->getArc());
        SkelEdgeDataSPtr data_edge = std::dynamic_pointer_cast<SkelEdgeData>(
                edge_toremove_1->getData());
        node->addSheet(data_edge->getSheet());
        data_edge = std::dynamic_pointer_cast<SkelEdgeData>(
                edge_toremove_2->getData());
        node->addSheet(data_edge->getSheet());
#endif

        queue.push(event);
    }
}

void SimpleStraightSkel::collectEdgeMergeEvents(PolyhedronSPtr polyhedron,
                                                const CGAL::FT& current_offset,
                                                const std::optional<CGAL::FT>& offset_future_bound,
                                                PQ& queue)
{
    return collectEdgeMergeEvents(polyhedron->edges(), polyhedron, true /*canonical reps*/,
                                  current_offset, offset_future_bound, queue);
}

void SimpleStraightSkel::collectTriangleEvents(const std::list<EdgeSPtr>& edges,
                                               PolyhedronSPtr polyhedron,
                                               const bool use_canonical_event_reps,
                                               const CGAL::FT& current_offset,
                                               const std::optional<CGAL::FT>& offset_future_bound,
                                               PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Triangle Event [" << current_offset << "]");

    for (EdgeSPtr edge : edges) {
        CGAL_assertion(edge->getID() != -1);

        if (edge->getVertexSrc()->getPoint() == edge->getVertexDst()->getPoint()) {
            continue;
        }

        if (isTetrahedron(edge)) {
            // tetrahedron event
            continue;
        }

        FacetSPtr facet;
        if (isTriangle(edge->getFacetL(), edge)) {
            facet = edge->getFacetL();
        } else if (isTriangle(edge->getFacetR(), edge)) {
            facet = edge->getFacetR();
        } else {
            continue;
        }

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        if (use_canonical_event_reps) {
            // rep canonicity: only check event existence when at the smallest edge of the facet
            if (edge->getID() > edge->prev(facet)->getID() ||
                edge->getID() > edge->next(facet)->getID()) {
                continue;
            }
        }
#endif

        bool dbl_triangle_event = false;
        EdgeSPtr edge_tmp = edge;
        for (unsigned int i = 0; i < 3; i++) {
            FacetSPtr facet_tmp_l = edge_tmp->getFacetL();
            FacetSPtr facet_tmp_r = edge_tmp->getFacetR();
            if (facet_tmp_l && facet_tmp_r) {
                if (isTriangle(facet_tmp_l, edge_tmp) &&
                        isTriangle(facet_tmp_r, edge_tmp)) {
                    dbl_triangle_event = true;
                    break;
                }
            }
            edge_tmp = edge_tmp->next(facet);
        }
        if (dbl_triangle_event) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge, current_offset, offset_future_bound);
        if (!point) {
            continue;
        }

        CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

        CGAL_assertion(KernelWrapper::side(edge->getFacetL()->getPlane(), point) < 0);
        CGAL_assertion(KernelWrapper::side(edge->getFacetR()->getPlane(), point) < 0);

        NodeSPtr node = Node::create();
        TriangleEventSPtr event = TriangleEvent::create();
        event->setStepID(step_id_);
        event->setNode(node);
        node->clear();
        node->setOffset(offset_event);
        node->setPoint(point);
        event->setFacet(facet);
        event->setEdgeBegin(edge);

#ifndef CGAL_SS3_NO_SKELETON_DS
        VertexSPtr vertices[3];
        event->getVertices(vertices);
        for (unsigned int i = 0; i < 3; i++) {
            SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(
                    vertices[i]->getData());
            ArcSPtr arc = data->getArc();
            node->addArc(arc);
        }
        EdgeSPtr edges[3];
        event->getEdges(edges);
        for (unsigned int i = 0; i < 3; i++) {
            SkelEdgeDataSPtr data = std::dynamic_pointer_cast<SkelEdgeData>(
                    edges[i]->getData());
            SheetSPtr sheet = data->getSheet();
            node->addSheet(sheet);
        }
#endif

        queue.push(event);
    }
}

void SimpleStraightSkel::collectTriangleEvents(PolyhedronSPtr polyhedron,
                                               const CGAL::FT& current_offset,
                                               const std::optional<CGAL::FT>& offset_future_bound,
                                               PQ& queue)
{
    return collectTriangleEvents(polyhedron->edges(), polyhedron, true /*use canonical reps*/,
                                 current_offset, offset_future_bound, queue);
}

void SimpleStraightSkel::collectDblEdgeMergeEvents(const std::list<EdgeSPtr>& edges,
                                                   PolyhedronSPtr polyhedron,
                                                   const bool use_canonical_event_reps,
                                                   const CGAL::FT& current_offset,
                                                   const std::optional<CGAL::FT>& offset_future_bound,
                                                   PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Dbl Edge Merge Events [" << current_offset << "]");

    for (EdgeSPtr edge : edges) {
        CGAL_assertion(edge->getID() != -1);

        if (!isReflex(edge)) {
            continue;
        }

        bool is_dbl_edge_merge_event = false;
        FacetSPtr facet_1;
        EdgeSPtr edge_11;
        EdgeSPtr edge_12;
        FacetSPtr facet_2;
        EdgeSPtr edge_21;
        EdgeSPtr edge_22;
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
            continue;
        }

        if (!is_dbl_edge_merge_event) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge, current_offset, offset_future_bound);
        if (!point) {
            continue;
        }

        CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

        NodeSPtr node = Node::create();
        DblEdgeMergeEventSPtr event = DblEdgeMergeEvent::create();
        event->setStepID(step_id_);
        event->setNode(node);
        node->clear();
        node->setOffset(offset_event);
        node->setPoint(point);
        event->setFacet1(facet_1);
        event->setEdge11(edge_11);
        event->setEdge12(edge_12);
        event->setFacet2(facet_2);
        event->setEdge21(edge_21);
        event->setEdge22(edge_22);

#ifndef CGAL_SS3_NO_SKELETON_DS
        VertexSPtr vertices[4];
        event->getVertices(vertices);
        for (unsigned int i = 0; i < 4; i++) {
            SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                    vertices[i]->getData());
            ArcSPtr arc = vertex_data->getArc();
            node->addArc(arc);
        }
        EdgeSPtr edges[4];
        event->getEdges(edges);
        for (unsigned int i = 0; i < 4; i++) {
            SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                    edges[i]->getData());
            SheetSPtr sheet = edge_data->getSheet();
            node->addSheet(sheet);
        }
#endif

        queue.push(event);
    }
}

void SimpleStraightSkel::collectDblEdgeMergeEvents(PolyhedronSPtr polyhedron,
                                                   const CGAL::FT& current_offset,
                                                   const std::optional<CGAL::FT>& offset_future_bound,
                                                   PQ& queue)
{
    return collectDblEdgeMergeEvents(polyhedron->edges(), polyhedron, true /*use canonical reps*/,
                                     current_offset, offset_future_bound, queue);
}

// contrary to all the other vanish events, this event is not detected from all
// of the edges it involves, but only from the main, "ridge", edge...
void SimpleStraightSkel::collectDblTriangleEvents(const std::list<EdgeSPtr>& edges,
                                                  PolyhedronSPtr polyhedron,
                                                  const bool /*use_canonical_event_reps*/,
                                                  const CGAL::FT& current_offset,
                                                  const std::optional<CGAL::FT>& offset_future_bound,
                                                  PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Dbl Triangle Events [" << current_offset << "]");

    for (EdgeSPtr edge : edges) {
        CGAL_assertion(edge->getID() != -1);

        if (isTetrahedron(edge)) {
            continue;
        }
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (!facet_l || !facet_r) {
            continue;
        }

        // use_canonical_event_reps
        // rep canonicity: this event is already only seen from its ridge edge
        // because of this check
        // @fixme which is something that is in fact bad in fact if we are updating the queue
        if (!(isTriangle(facet_l, edge) &&
                isTriangle(facet_r, edge))) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge, current_offset, offset_future_bound);
        if (!point) {
            continue;
        }

        CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

        NodeSPtr node = Node::create();
        DblTriangleEventSPtr event = DblTriangleEvent::create();
        event->setStepID(step_id_);
        event->setNode(node);
        node->clear();
        node->setOffset(offset_event);
        node->setPoint(point);
        event->setEdge(edge);

#ifndef CGAL_SS3_NO_SKELETON_DS
        VertexSPtr vertices[4];
        event->getVertices(vertices);
        for (unsigned int i = 0; i < 4; i++) {
            SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(
                    vertices[i]->getData());
            ArcSPtr arc = data->getArc();
            node->addArc(arc);
        }
        EdgeSPtr edges[5];
        event->getEdges(edges);
        for (unsigned int i = 0; i < 5; i++) {
            SkelEdgeDataSPtr data = std::dynamic_pointer_cast<SkelEdgeData>(
                    edges[i]->getData());
            SheetSPtr sheet = data->getSheet();
            node->addSheet(sheet);
        }
#endif

        queue.push(event);
    }
}

void SimpleStraightSkel::collectDblTriangleEvents(PolyhedronSPtr polyhedron,
                                                  const CGAL::FT& current_offset,
                                                  const std::optional<CGAL::FT>& offset_future_bound,
                                                  PQ& queue)
{
    return collectDblTriangleEvents(polyhedron->edges(), polyhedron, true /*use canonical reps*/,
                                    current_offset, offset_future_bound, queue);
}

void SimpleStraightSkel::collectTetrahedronEvents(const std::list<EdgeSPtr>& edges,
                                                  PolyhedronSPtr polyhedron,
                                                  const bool use_canonical_event_reps,
                                                  const CGAL::FT& current_offset,
                                                  const std::optional<CGAL::FT>& offset_future_bound,
                                                  PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Tetrahedron Events [" << current_offset << "]");

    for (EdgeSPtr edge : edges) {
        CGAL_assertion(edge->getID() != -1);

        if (isTetrahedron(edge)) {

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
            if (use_canonical_event_reps) {
                  // rep canonicity: only investigate this event if we are in the smallest edge
                  // of the smallest facet of the tetrahedron
                  FacetSPtr facet = (edge->getFacetL()->getID() < edge->getFacetR()->getID()) ? edge->getFacetL()
                                                                                              : edge->getFacetR();

                  // smallest incident facet is the smallest of the tetrahedron
                  if (facet->getID() > edge->prev(facet)->other(facet)->getID() ||
                      facet->getID() > edge->other(facet)->getID() ||
                      facet->getID() > edge->next(facet)->other(facet)->getID()) {
                      continue;
                  }

                  // edge is the smallest within the facet
                  if (edge->getID() > edge->next(facet)->getID() ||
                      edge->getID() > edge->prev(facet)->getID()) {
                      continue;
                  }
            }
#endif

            Point3SPtr point = Point3SPtr();
            CGAL::FT offset_event;
            std::tie(point, offset_event) = vanishesAt(edge, current_offset, offset_future_bound);
            if (!point) {
                continue;
            }

            CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

            NodeSPtr node = Node::create();
            TetrahedronEventSPtr event = TetrahedronEvent::create();
            event->setStepID(step_id_);
            event->setNode(node);
            node->clear();
            node->setOffset(offset_event);
            node->setPoint(point);
            event->setEdgeBegin(edge);

#ifndef CGAL_SS3_NO_SKELETON_DS
            VertexSPtr vertices[4];
            event->getVertices(vertices);
            for (unsigned int i = 0; i < 4; i++) {
                SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                        vertices[i]->getData());
                ArcSPtr arc = vertex_data->getArc();
                node->addArc(arc);
            }
            EdgeSPtr edges[6];
            event->getEdges(edges);
            for (unsigned int i = 0; i < 6; i++) {
                SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                        edges[i]->getData());
                SheetSPtr sheet = edge_data->getSheet();
                node->addSheet(sheet);
            }
#endif

            queue.push(event);
        }
    }
}

void SimpleStraightSkel::collectTetrahedronEvents(PolyhedronSPtr polyhedron,
                                                  const CGAL::FT& current_offset,
                                                  const std::optional<CGAL::FT>& offset_future_bound,
                                                  PQ& queue)
{
    return collectTetrahedronEvents(polyhedron->edges(), polyhedron, true /*use canonical reps*/,
                                    current_offset, offset_future_bound, queue);
}

void SimpleStraightSkel::collectVertexEvents(const std::list<VertexSPtr>& vertices,
                                             PolyhedronSPtr polyhedron,
                                             const bool use_canonical_event_reps,
                                             const CGAL::FT& current_offset,
                                             const std::optional<CGAL::FT>& offset_future_bound,
                                             PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Vertex Events [" << current_offset << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    for (VertexSPtr vertex_1 : vertices) {
        CGAL_assertion(vertex_1->getID() != -1);

        if (isConvex(vertex_1)) {
            continue;
        }

// this only seeks the first vertex_2 that share the same incident facets as vertex_1 and
// is in the cycle (instead of all vertex_2 from all cycles of all incident facets...)
// #define CGAL_SS3_NEWER_VV_VERTEX_2_DETECTION
#ifdef CGAL_SS3_NEWER_VV_VERTEX_2_DETECTION
        I do not think whatever is below is enough in case of some weird tangled faces.
        Indeed, we want the first vertices left and right _on the supporting line_ of the edge
        but with some weird enough facet, we could have a vertex that is not the first one
        but that is still the closest to the edge on the supporting line.

        also it is pointless to walk both left and right, just walk one direction for each facet

        also the other_facet looks kinda weird, why edge_r?

        // @todo rewrite all of that stuff not to waste the facet_1/facet_2 information
        std::set<VertexSPtr> vertices_2;
        for (FacetWPtr facet_wptr : vertex_1->facets()) {
            if (FacetSPtr facet = facet_wptr.lock()) {
                CGAL_SS3_CORE_TRACE("walking facet: " << facet->getID());

                EdgeSPtr edge_l = vertex_1->findEdge(facet);
                CGAL_SS3_DEBUG_SPTR(edge_l);
                CGAL_assertion(edge_l->src(facet) == vertex_1);
                EdgeSPtr edge_r = edge_l->prev(facet);
                CGAL_assertion(edge_r->dst(facet) == vertex_1);

                // seek to the left
                CGAL_SS3_CORE_TRACE("------> LEFT WALK");
                FacetSPtr other_facet = start_e->other(facet);
                CGAL_SS3_CORE_TRACE("other facet: " << other_facet->getID());

                CGAL_assertion(other_facet != edge_l->other(facet));

                // skip edge_l directly, and ignore one more because we want v-v events, not vanish.
                // could skip two, but then it becomes annoying with triangular faces
                CGAL_assertion(edge_l->next(facet) != edge_l);
                CGAL_assertion(edge_l->next(facet)->next(facet) != edge_l);
                EdgeSPtr edge_l1 = edge_l->next(facet);
                EdgeSPtr edge_l2 = edge_l1->next(facet);

                // skip edge_r directly, and ignore one more because we want v-v events, not vanish.
                // could skip two, but then it becomes annoying with triangular faces
                CGAL_assertion(edge_r->prev(facet) != edge_r);
                CGAL_assertion(edge_r->prev(facet)->prev(facet) != edge_r);
                EdgeSPtr edge_r1 = edge_r->prev(facet);
                EdgeSPtr edge_r2 = edge_r1->prev(facet);

                EdgeSPtr current_edge = edge_l2;
                for (;;) {
                    CGAL_SS3_CORE_TRACE("current_edge: " << current_edge->getVertexSrc()->getID()
                                                         << " " << current_edge->getVertexDst()->getID()
                                                         << " faces: " << current_edge->getFacetL()->getID()
                                                         << " " << current_edge->getFacetR()->getID());

                    if (current_edge == edge_r1 || current_edge == edge_r) {
                        break;
                    }

                    if (current_edge->other(facet) == other_facet) {

                      // we have found another vertex that is along an edge that has the same
                      // incident facets. But to be of interest, this needs to be on the correct side
                      // and not a full loop around
                      if (CGAL::collinear_are_ordered_along_line(*(edge_l->src(facet)->getPoint()),
                                                                 *(edge_l->dst(facet)->getPoint()),
                                                                 *(current_edge->src(facet)->getPoint()))) {
                            CGAL_SS3_CORE_TRACE("left neighbor: " << current_edge->src(facet)->getID());
                            vertices_2.insert(current_edge->src(facet));
                            break;
                        }
                      }

                    current_edge = current_edge->next(facet);
                }

                // seek to the right
                CGAL_SS3_CORE_TRACE("------> RIGHT WALK");
                other_facet = edge_l->other(facet);
                CGAL_assertion(other_facet != edge_r->other(facet));

                current_edge = edge_r2;
                for (;;) {
                    CGAL_SS3_CORE_TRACE("current_edge: " << current_edge->getVertexSrc()->getID() << " " << current_edge->getVertexDst()->getID() << " faces: " << current_edge->getFacetL()->getID() << " " << current_edge->getFacetR()->getID());

                    if (current_edge == edge_l1 || current_edge == edge_l) {
                        break;
                    }

                    if (current_edge->other(facet) == other_facet) {
                        // we have found another vertex that is along an edge that has the same
                        // incident facets. But to be of interest, this needs to be on the correct side
                        // and not a full loop around
                        if (CGAL::collinear_are_ordered_along_line(*(edge_r->dst(facet)->getPoint()),
                                                                   *(edge_r->src(facet)->getPoint()),
                                                                   *(current_edge->dst(facet)->getPoint()))) {
                            CGAL_SS3_CORE_TRACE("right neighbor: " << current_edge->dst(facet)->getID());
                            vertices_2.insert(current_edge->dst(facet));
                            break;
                        }
                    }

                    current_edge = current_edge->prev(facet);
                }
            }
        }
#else // CGAL_SS3_NEWER_VV_VERTEX_2_DETECTION
        std::set<VertexSPtr> vertices_2;
        for (FacetWPtr facet_wptr : vertex_1->facets()) {
            if (FacetSPtr facet = facet_wptr.lock()) {
# ifndef CGAL_SS3_VV_VERTEX_2_WALK_FACES_FOR_DETECTION
                vertices_2.insert(facet->vertices().begin(), facet->vertices().end());
# else
                // Should be enough to check within the same border because different CCs
                // of the same facet cannot merge again

                // Subtlety: within a facet we are only checking the edge whose source is vertex_1
                // That's enough for two reasons:
                // - VV events are symmetrical
                // - The other edge will be checked when seen from the other facet
                EdgeSPtr start_e = vertex_1->findEdge(facet);
                FacetSPtr other_facet = start_e->other(facet);
                EdgeSPtr current_edge = start_e->next(facet);
                while (current_edge != start_e) {
                    if (current_edge->other(facet) == other_facet) {
                        vertices_2.insert(current_edge->src(facet));
                        vertices_2.insert(current_edge->dst(facet));
                    }

                    current_edge = current_edge->next(facet);
                }
# endif
                // @speed sort vertices_2 as to get the first vertex after start_e->dst(facet)
            }
        }
#endif // CGAL_SS3_NEWER_VV_VERTEX_2_DETECTION

        for (VertexSPtr vertex_2 : vertices_2) {
            CGAL_assertion(vertex_2->getID() != -1);
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
            if (isConvex(vertex_2)) {
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
            std::list<FacetWPtr>::iterator it_f1 = v1->facets().begin();
            while (it_f1 != v1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = v2->facets().begin();
                    while (it_f2 != v2->facets().end()) {
                        FacetWPtr facet_2_wptr = *it_f2++;
                        if (facet_1_wptr == facet_2_wptr) {
                            if (num_equal_facets == 0) {
                                facet_1 = FacetSPtr(facet_1_wptr);
                            } else {
                                facet_2 = FacetSPtr(facet_2_wptr);
                            }
                            num_equal_facets++;
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

#ifndef CGAL_SS3_CHECK_CONV_SPLIT_EVENT_AT_POP_TIME
            bool conv_split_event = false;
            FacetSPtr facet_1b = facet_2->next(v1);
            FacetSPtr facet_2b = facet_1->next(v2);
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
                continue;
            }
#endif

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset, offset_future_bound);
            if (!point) {
                continue;
            }

            CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

            NodeSPtr node = Node::create();
            VertexEventSPtr event = VertexEvent::create();
            event->setStepID(step_id_);
            event->setNode(node);
            node->clear();

#ifndef CGAL_SS3_NO_SKELETON_DS
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(v1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(->getData());
            node->addArc(data_2->getArc());
#endif
            node->setOffset(offset_event);
            node->setPoint(point);
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

void SimpleStraightSkel::collectVertexEvents(PolyhedronSPtr polyhedron,
                                             const CGAL::FT& current_offset,
                                             const std::optional<CGAL::FT>& offset_future_bound,
                                             PQ& queue)
{
    return collectVertexEvents(polyhedron->vertices(), polyhedron, true /*use canonical reps*/,
                               current_offset, offset_future_bound, queue);
}

void SimpleStraightSkel::collectFlipVertexEvents(const std::list<VertexSPtr>& vertices,
                                                 PolyhedronSPtr polyhedron,
                                                 const bool use_canonical_event_reps,
                                                 const CGAL::FT& current_offset,
                                                 const std::optional<CGAL::FT>& offset_future_bound,
                                                 PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Flip Vertex Events [" << current_offset << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    // @speed seems like this could be rewritten in a much cheaper manner (and storing
    // the convexity property). But currently, it's a cheap collect().
    for (VertexSPtr vertex_1 : vertices) {
        CGAL_assertion(vertex_1->getID() != -1);

        if (isConvex(vertex_1)) {
            continue;
        }

#ifdef CGAL_SS3_NEWER_VV_VERTEX_2_DETECTION
        take code from collectVertexEvents
#else // CGAL_SS3_NEWER_VV_VERTEX_2_DETECTION
        std::set<VertexSPtr> vertices_2;
        for (FacetWPtr facet_wptr : vertex_1->facets()) {
            if (FacetSPtr facet = facet_wptr.lock()) {
# ifndef CGAL_SS3_VV_VERTEX_2_WALK_FACES_FOR_DETECTION
                vertices_2.insert(facet->vertices().begin(), facet->vertices().end());
# else
                // Should be enough to check within the same border because different CCs
                // of the same facet cannot merge again
                EdgeSPtr start_e = vertex_1->findEdge(facet);
                FacetSPtr other_facet = start_e->other(facet);
                EdgeSPtr current_edge = start_e->next(facet);
                while (current_edge != start_e) {
                    if (current_edge->other(facet) == other_facet) {
                        vertices_2.insert(current_edge->src(facet));
                        vertices_2.insert(current_edge->dst(facet));
                    }

                    current_edge = current_edge->next(facet);
                }
# endif
            }
        }
#endif // CGAL_SS3_NEWER_VV_VERTEX_2_DETECTION

        for (VertexSPtr vertex_2 : vertices_2) {
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
            if (isConvex(vertex_2)) {
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
            std::list<FacetWPtr>::iterator it_f1 = v1->facets().begin();
            while (it_f1 != v1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = v2->facets().begin();
                    while (it_f2 != v2->facets().end()) {
                        FacetWPtr facet_2_wptr = *it_f2++;
                        if (facet_1_wptr == facet_2_wptr) {
                            if (num_equal_facets == 0) {
                                facet_1 = FacetSPtr(facet_1_wptr);
                            } else {
                                facet_2 = FacetSPtr(facet_2_wptr);
                            }
                            num_equal_facets++;
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

#ifndef CGAL_SS3_CHECK_CONV_SPLIT_EVENT_AT_POP_TIME
            bool conv_split_event = false;
            FacetSPtr facet_1b = facet_2->next(v1);
            FacetSPtr facet_2b = facet_2->next(v2);
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
                continue;
            }
#endif

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset, offset_future_bound);
            if (!point) {
                continue;
            }

            CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

            NodeSPtr node = Node::create();
            FlipVertexEventSPtr event = FlipVertexEvent::create();
            event->setStepID(step_id_);
            event->setNode(node);
            node->clear();

#ifndef CGAL_SS3_NO_SKELETON_DS
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(v1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(v2->getData());
            node->addArc(data_2->getArc());
#endif

            node->setOffset(offset_event);
            node->setPoint(point);
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

void SimpleStraightSkel::collectFlipVertexEvents(PolyhedronSPtr polyhedron,
                                                 const CGAL::FT& current_offset,
                                                 const std::optional<CGAL::FT>& offset_future_bound,
                                                 PQ& queue)
{
    return collectFlipVertexEvents(polyhedron->vertices(), polyhedron, true /*use canonical reps*/,
                                   current_offset, offset_future_bound, queue);
}

void SimpleStraightSkel::collectSurfaceEvent(EdgeSPtr edge_1,
                                             EdgeSPtr edge_2,
                                             PolyhedronSPtr polyhedron,
                                             const CGAL::FT& current_offset,
                                             const std::optional<CGAL::FT>& offset_future_bound,
                                             PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(8, ">>> Collect Surface Event [\n  " << edge_1->toString() << "\n  " << edge_2->toString() << "]");

    CGAL_assertion(edge_1->getID() != -1);
    CGAL_assertion(edge_2->getID() != -1);

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

#ifndef CGAL_SS3_CHECK_CONV_SPLIT_EVENT_AT_POP_TIME
    bool is_conv_split_event = false;
    std::list<EdgeSPtr> common_edges = edge_1->getFacetL()->findEdges(edge_1->getFacetR());
    for (EdgeSPtr edge : common_edges) {
        if (edge == edge_1) {
            continue;
        }
        FacetSPtr facet_src = edge->getFacetSrc();
        FacetSPtr facet_dst = edge->getFacetDst();
        if (facet_1_src == edge_2->getFacetL() ||
                facet_1_dst == edge_2->getFacetL()) {
            if (facet_src == edge_2->getFacetR() ||
                    facet_dst == edge_2->getFacetR()) {
                is_conv_split_event = true;
                break;
            }
        } else if (facet_1_src == edge_2->getFacetR() ||
                facet_1_dst == edge_2->getFacetR()) {
            if (facet_src == edge_2->getFacetL() ||
                    facet_dst == edge_2->getFacetL()) {
                is_conv_split_event = true;
                break;
            }
        }
    }

    if (is_conv_split_event) {
        return;
    }
#endif

    // not outside of the loop just because maybe one day this will be called
    // as the first collect function with an initial bound that gets updated...

    // let's just check if bboxes overlap first
    if (offset_future_bound) {
        CGAL::Bbox_3 b1;
        b1 += edge_1->getVertexSrc()->getPoint()->bbox();
        b1 += edge_1->getVertexDst()->getPoint()->bbox();
        b1 += getFinalPoint(edge_1->getVertexSrc(), *offset_future_bound)->bbox();
        b1 += getFinalPoint(edge_1->getVertexDst(), *offset_future_bound)->bbox();

        CGAL::Bbox_3 b2;
        b2 += edge_2->getVertexSrc()->getPoint()->bbox();
        b2 += edge_2->getVertexDst()->getPoint()->bbox();
        b2 += getFinalPoint(edge_2->getVertexSrc(), *offset_future_bound)->bbox();
        b2 += getFinalPoint(edge_2->getVertexDst(), *offset_future_bound)->bbox();

        if (!CGAL::do_overlap(b1, b2)) {
            CGAL_SS3_CORE_TRACE_V(32, "Filtered possible surface event candidates\n\t" << edge_1->toString() << "\n\t"
                                                                                       << edge_2->toString());
            return;
        } else {
            CGAL_SS3_CORE_TRACE_V(32, "Checking possible surface event event\n\t" << edge_1->toString() << "\n\t"
                                                                                  << edge_2->toString());
        }
    }

    // calculate intersection point
    Point3SPtr point;
    CGAL::FT offset_event;
    std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset, offset_future_bound);
    if (!point) {
        return;
    }

    CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

    NodeSPtr node = Node::create();
    SurfaceEventSPtr event = SurfaceEvent::create();
    event->setStepID(step_id_);
    event->setNode(node);
    node->clear();
    node->setOffset(offset_event);
    node->setPoint(point);
    event->setEdge1(edge_1);
    event->setEdge2(edge_2);

#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
            edge_1->getData());
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
            edge_2->getData());
    node->addSheet(data_1->getSheet());
    node->addSheet(data_2->getSheet());

    if (facet_1_src == edge_2->getFacetL() ||
            facet_1_src == edge_2->getFacetR()) {
        SkelVertexDataSPtr data_1_src = std::dynamic_pointer_cast<SkelVertexData>(
            edge_1->getVertexSrc()->getData());
        node->addArc(data_1_src->getArc());
    }
    if (facet_1_dst == edge_2->getFacetL() ||
            facet_1_dst == edge_2->getFacetR()) {
        SkelVertexDataSPtr data_1_dst = std::dynamic_pointer_cast<SkelVertexData>(
            edge_1->getVertexDst()->getData());
        node->addArc(data_1_dst->getArc());
    }
#endif

    queue.push(event);
}

void SimpleStraightSkel::collectSurfaceEvents(const std::list<EdgeSPtr>& edges,
                                              PolyhedronSPtr polyhedron,
                                              const CGAL::FT& current_offset,
                                              const std::optional<CGAL::FT>& offset_future_bound,
                                              PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Surface Events [" << current_offset << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    unsigned int filtered_candidates = 0;
    unsigned int tested_candidates = 0;
#endif

    for (EdgeSPtr edge_1 : edges) {
        CGAL_assertion(edge_1->getID() != -1);

        FacetSPtr facet_1_src = edge_1->getFacetSrc();
        FacetSPtr facet_1_dst = edge_1->getFacetDst();
        std::list<EdgeSPtr> edges_2;
        edges_2.insert(edges_2.end(), facet_1_src->edges().begin(), facet_1_src->edges().end());
        edges_2.insert(edges_2.end(), facet_1_dst->edges().begin(), facet_1_dst->edges().end());

        // not outside of the loop just because maybe one day this will be called
        // as the first collect function with an initial bound that gets updated...
        CGAL::Bbox_3 b1;
        if (offset_future_bound) {
            b1 += edge_1->getVertexSrc()->getPoint()->bbox();
            b1 += edge_1->getVertexDst()->getPoint()->bbox();
            b1 += getFinalPoint(edge_1->getVertexSrc(), *offset_future_bound)->bbox();
            b1 += getFinalPoint(edge_1->getVertexDst(), *offset_future_bound)->bbox();
        }

        for (EdgeSPtr edge_2 : edges_2) {
            CGAL_assertion(edge_2->getID() != -1);

            CGAL_SS3_CORE_TRACE_V(32, "Possible surface event:\n\t" << edge_1->toString() << "\n\t"
                                                                    << edge_2->toString());

            if (edge_1 == edge_2) {
                continue;
            }

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

#ifndef CGAL_SS3_CHECK_CONV_SPLIT_EVENT_AT_POP_TIME
            bool is_conv_split_event = false;
            std::list<EdgeSPtr> common_edges = edge_1->getFacetL()->findEdges(edge_1->getFacetR());
            for (EdgeSPtr edge : common_edges) {
                if (edge == edge_1) {
                    continue;
                }
                FacetSPtr facet_src = edge->getFacetSrc();
                FacetSPtr facet_dst = edge->getFacetDst();
                if (facet_1_src == edge_2->getFacetL() ||
                        facet_1_dst == edge_2->getFacetL()) {
                    if (facet_src == edge_2->getFacetR() ||
                            facet_dst == edge_2->getFacetR()) {
                        is_conv_split_event = true;
                        break;
                    }
                } else if (facet_1_src == edge_2->getFacetR() ||
                        facet_1_dst == edge_2->getFacetR()) {
                    if (facet_src == edge_2->getFacetL() ||
                            facet_dst == edge_2->getFacetL()) {
                        is_conv_split_event = true;
                        break;
                    }
                }
            }

            if (is_conv_split_event) {
                continue;
            }
#endif

            // let's just check if bboxes overlap first
            if (offset_future_bound) {
                CGAL::Bbox_3 b2;
                b2 += edge_2->getVertexSrc()->getPoint()->bbox();
                b2 += edge_2->getVertexDst()->getPoint()->bbox();
                b2 += getFinalPoint(edge_2->getVertexSrc(), *offset_future_bound)->bbox();
                b2 += getFinalPoint(edge_2->getVertexDst(), *offset_future_bound)->bbox();

                if (!CGAL::do_overlap(b1, b2)) {
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
                    ++filtered_candidates;
#endif
                    CGAL_SS3_CORE_TRACE_V(32, "Filtered possible surface event candidates\n\t" << edge_1->toString() << "\n\t"
                                                                                               << edge_2->toString());
                    continue;
                } else {
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
                    ++tested_candidates;
#endif
                }
            }

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset, offset_future_bound);
            if (!point) {
                continue;
            }

            CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

            NodeSPtr node = Node::create();
            SurfaceEventSPtr event = SurfaceEvent::create();
            event->setStepID(step_id_);
            event->setNode(node);
            node->clear();
            node->setOffset(offset_event);
            node->setPoint(point);
            event->setEdge1(edge_1);
            event->setEdge2(edge_2);

#ifndef CGAL_SS3_NO_SKELETON_DS
            SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_1->getData());
            SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_2->getData());
            node->addSheet(data_1->getSheet());
            node->addSheet(data_2->getSheet());

            if (facet_1_src == edge_2->getFacetL() ||
                    facet_1_src == edge_2->getFacetR()) {
                SkelVertexDataSPtr data_1_src = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_1->getVertexSrc()->getData());
                node->addArc(data_1_src->getArc());
            }
            if (facet_1_dst == edge_2->getFacetL() ||
                    facet_1_dst == edge_2->getFacetR()) {
                SkelVertexDataSPtr data_1_dst = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_1->getVertexDst()->getData());
                node->addArc(data_1_dst->getArc());
            }
#endif

            queue.push(event);
        }
    }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    CGAL_SS3_CORE_TRACE_V(16, "  " << filtered_candidates << " filtered, " << tested_candidates << " tests");
#endif

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Surface Events in: " << timer.time());
#endif
}

void SimpleStraightSkel::collectSurfaceEvents(PolyhedronSPtr polyhedron,
                                              const CGAL::FT& current_offset,
                                              const std::optional<CGAL::FT>& offset_future_bound,
                                              PQ& queue)
{
    return collectSurfaceEvents(polyhedron->edges(), polyhedron,
                                current_offset, offset_future_bound, queue);
}

void SimpleStraightSkel::collectPolyhedronSplitEvent(EdgeSPtr edge_1,
                                                     EdgeSPtr edge_2,
                                                     PolyhedronSPtr polyhedron,
                                                     const CGAL::FT& current_offset,
                                                     const std::optional<CGAL::FT>& offset_future_bound,
                                                     PQ& queue)
{
  CGAL_SS3_DEBUG_SPTR(edge_1);
  CGAL_SS3_DEBUG_SPTR(edge_2);

  CGAL_assertion(edge_1->getID() != -1);
  CGAL_assertion(edge_2->getID() != -1);

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

    // calculate intersection point
    Point3SPtr point;
    CGAL::FT offset_event;
    std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset, offset_future_bound);
    if (!point) {
        return;
    }

    CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

#ifndef CGAL_SS3_CHECK_POLYHEDRON_SPLIT_EVENT_AT_POP_TIME
    CGAL::FT shift = offset_event - current_offset;
    Segment3SPtr e1o = PolyhedronTransformation::shiftEdge(edge_1, shift);
    if (!e1o->is_degenerate()) {
        // not a polyhedron split (edge split?)
        return;
    }
#endif

    NodeSPtr node = Node::create();
    PolyhedronSplitEventSPtr event = PolyhedronSplitEvent::create();
    event->setStepID(step_id_);
    event->setNode(node);
    node->clear();
    node->setOffset(offset_event);
    node->setPoint(point);
    event->setEdge1(edge_1);
    event->setEdge2(edge_2);

#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(edge_1->getData());
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(edge_2->getData());
    node->addSheet(data_1->getSheet());
    node->addSheet(data_2->getSheet());

    if (facet_1_src == edge_2->getFacetL() ||
        facet_1_src == edge_2->getFacetR()) {
        SkelVertexDataSPtr data_1_src = std::dynamic_pointer_cast<SkelVertexData>(
            edge_1->getVertexSrc()->getData());
        node->addArc(data_1_src->getArc());
    }
    if (facet_1_dst == edge_2->getFacetL() ||
        facet_1_dst == edge_2->getFacetR()) {
        SkelVertexDataSPtr data_1_dst = std::dynamic_pointer_cast<SkelVertexData>(
          edge_1->getVertexDst()->getData());
        node->addArc(data_1_dst->getArc());
    }
#endif

    queue.push(event);
}

void SimpleStraightSkel::collectPolyhedronSplitEvents(const std::list<EdgeSPtr>& edges,
                                                      PolyhedronSPtr polyhedron,
                                                      const CGAL::FT& current_offset,
                                                      const std::optional<CGAL::FT>& offset_future_bound,
                                                      PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Polyhedron Split Events [" << current_offset << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    for (EdgeSPtr edge_1 : edges) {
        CGAL_assertion(edge_1->getID() != -1);

        if (!isReflex(edge_1)) {
            continue;
        }

        FacetSPtr facet_1_src = edge_1->getFacetSrc();
        CGAL_assertion(bool(facet_1_src));
        for (EdgeSPtr edge_2 : facet_1_src->edges()) {
            collectPolyhedronSplitEvent(edge_1, edge_2, polyhedron,
                                        current_offset, offset_future_bound,
                                        queue);
        }
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Polyhedron Split Events in: " << timer.time());
#endif
}

void SimpleStraightSkel::collectPolyhedronSplitEvents(PolyhedronSPtr polyhedron,
                                                      const CGAL::FT& current_offset,
                                                      const std::optional<CGAL::FT>& offset_future_bound,
                                                      PQ& queue)
{
    return collectPolyhedronSplitEvents(polyhedron->edges(), polyhedron,
                                        current_offset, offset_future_bound, queue);
}

// @todo large amount of code duplicated between VertexEvent, FlipVertexEvent and SplitMergeEvent
void SimpleStraightSkel::collectSplitMergeEvents(const std::list<VertexSPtr>& vertices,
                                                 PolyhedronSPtr polyhedron,
                                                 const bool use_canonical_event_reps,
                                                 const CGAL::FT& current_offset,
                                                 const std::optional<CGAL::FT>& offset_future_bound,
                                                 PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Split Merge Events [" << current_offset << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    for (VertexSPtr vertex_1 : vertices) {
        CGAL_assertion(vertex_1->getID() != -1);

        if (isConvex(vertex_1)) {
            continue;
        }

#ifdef CGAL_SS3_NEWER_VV_VERTEX_2_DETECTION
        take code from collectVertexEvents
#else // CGAL_SS3_NEWER_VV_VERTEX_2_DETECTION
        std::set<VertexSPtr> vertices_2;
        for (FacetWPtr facet_wptr : vertex_1->facets()) {
            if (FacetSPtr facet = facet_wptr.lock()) {
# ifndef CGAL_SS3_VV_VERTEX_2_WALK_FACES_FOR_DETECTION
                vertices_2.insert(facet->vertices().begin(), facet->vertices().end());
# else
                // Should be enough to check within the same border because different CCs
                // of the same facet cannot merge again

                currently this improvement creates an annoying fake bug in the queue check:
                by walking incident faces in a single direction, this introduces
                some asymetry in the detection of the event: e.g. if the event were 5-6,
                we could detect it when looking from 5 but not from 6 because the facet walked
                from 5 has 6, but the facet walked from 6 does not have 5 (5 is in a difference_type CC of the facet) and then we could be at 6 in local, but 5
                when filling the queue globally

                EdgeSPtr start_e = vertex_1->findEdge(facet);
                FacetSPtr other_facet = start_e->other(facet);
                EdgeSPtr current_edge = start_e->next(facet);
                while (current_edge != start_e) {
                    if (current_edge->other(facet) == other_facet) {
                        vertices_2.insert(current_edge->src(facet));
                        vertices_2.insert(current_edge->dst(facet));
                    }

                    current_edge = current_edge->next(facet);
                }
# endif
            }
        }

# ifdef CGAL_SS3_RUN_TIMERS
        // CGAL_SS3_CORE_TRACE("filled vertices_2 after " << timer.time() << "s");
# endif
#endif // CGAL_SS3_NEWER_VV_VERTEX_2_DETECTION

        for (VertexSPtr vertex_2 : vertices_2) {
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
            if (isConvex(vertex_2)) {
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
            std::list<FacetWPtr>::iterator it_f1 = v1->facets().begin();
            while (it_f1 != v1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = v2->facets().begin();
                    while (it_f2 != v2->facets().end()) {
                        FacetWPtr facet_2_wptr = *it_f2++;
                        if (facet_1_wptr == facet_2_wptr) {
                            if (num_equal_facets == 0) {
                                facet_1 = FacetSPtr(facet_1_wptr);
                            } else {
                                facet_2 = FacetSPtr(facet_2_wptr);
                            }
                            num_equal_facets++;
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

#ifndef CGAL_SS3_CHECK_CONV_SPLIT_EVENT_AT_POP_TIME
            bool conv_split_event = false;
            FacetSPtr facet_1b = facet_2->next(v1);
            FacetSPtr facet_2b = facet_1->next(v2);
            if (facet_2b == facet_2) {
                // flip vertex event
                facet_2b = facet_2b->next(v2);
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
                continue;
            }
#endif

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset, offset_future_bound);
            if (!point) {
                continue;
            }

            CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

            NodeSPtr node = Node::create();
            SplitMergeEventSPtr event = SplitMergeEvent::create();
            event->setStepID(step_id_);
            event->setNode(node);
            node->clear();

#ifndef CGAL_SS3_NO_SKELETON_DS
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(v1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(v2->getData());
            node->addArc(data_2->getArc());
#endif

            node->setOffset(offset_event);
            node->setPoint(point);
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

void SimpleStraightSkel::collectSplitMergeEvents(PolyhedronSPtr polyhedron,
                                                 const CGAL::FT& current_offset,
                                                 const std::optional<CGAL::FT>& offset_future_bound,
                                                 PQ& queue)
{
    return collectSplitMergeEvents(polyhedron->vertices(), polyhedron, true /*use canonical reps*/,
                                   current_offset, offset_future_bound, queue);
}

void SimpleStraightSkel::collectEdgeSplitEvents(const std::list<EdgeSPtr>& edges_1,
                                                const std::list<EdgeSPtr>& edges_2,
                                                PolyhedronSPtr polyhedron,
                                                const bool use_canonical_event_reps,
                                                const CGAL::FT& current_offset,
                                                const std::optional<CGAL::FT>& offset_future_bound,
                                                PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Edge Split Events [" << current_offset << "]");

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    unsigned int filtered_candidates = 0;
    unsigned int tested_candidates = 0;
#endif

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    std::list<EdgeSPtr> edges_reflex_1, edges_reflex_2;
    auto fill_reflex_edges = [&](const std::list<EdgeSPtr>& edges,
                                 std::list<EdgeSPtr>& edges_reflex) {
        for (EdgeSPtr edge : edges) {
            if (isReflex(edge)) {
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
    CGAL_SS3_CORE_TRACE_V(8, "  Collect reflex edges: " << timer.time());
#endif

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    CGAL_SS3_CORE_TRACE_V(16, "  " << edges_reflex_1.size() << " reflex edges [1]");
    CGAL_SS3_CORE_TRACE_V(16, "  " << edges_reflex_2.size() << " reflex edges [2]");
#endif

    // maybe we will use canonical reps with non symmetrical ranges one day, but not for now
    CGAL_warning(!use_canonical_event_reps || edges_reflex_1 == edges_reflex_2);

    for (EdgeSPtr edge_1 : edges_reflex_1) {
        CGAL_assertion(edge_1->getID() != -1);

        FacetSPtr facet_1_src = edge_1->getFacetSrc();
        FacetSPtr facet_1_dst = edge_1->getFacetDst();

        // not outside of the loop just because maybe one day this will be called
        // as the first collect function with an initial bound that gets updated...

        CGAL::Bbox_3 b1;
        if (offset_future_bound) {
            b1 += edge_1->getVertexSrc()->getPoint()->bbox();
            b1 += edge_1->getVertexDst()->getPoint()->bbox();
            b1 += getFinalPoint(edge_1->getVertexSrc(), *offset_future_bound)->bbox();
            b1 += getFinalPoint(edge_1->getVertexDst(), *offset_future_bound)->bbox();
        }

        for (EdgeSPtr edge_2 : edges_reflex_2) {
            CGAL_assertion(edge_2->getID() != -1);

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
            if (offset_future_bound) {
                CGAL::Bbox_3 b2;
                b2 += edge_2->getVertexSrc()->getPoint()->bbox();
                b2 += edge_2->getVertexDst()->getPoint()->bbox();
                b2 += getFinalPoint(edge_2->getVertexSrc(), *offset_future_bound)->bbox();
                b2 += getFinalPoint(edge_2->getVertexDst(), *offset_future_bound)->bbox();

                if (!CGAL::do_overlap(b1, b2)) {
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
                    ++filtered_candidates;
#endif
                    CGAL_SS3_CORE_TRACE_V(64, "Filtered edge split candidates\n\t" << edge_1->toString() << "\n\t"
                                                                                   << edge_2->toString() << "\n\t"
                                                                                   << b1 << "\n\t" << b2);
                    continue;
                } else {
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
                    ++tested_candidates;
#endif
                    CGAL_SS3_CORE_TRACE_V(64, "Checking possible edge split event\n\t" << edge_1->toString() << "\n\t"
                                                                                       << edge_2->toString() << "\n\t"
                                                                                       << b1 << "\n\t" << b2);
                }
            }

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset, offset_future_bound);
            if (!point) {
                continue;
            }

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
            // @fixme below needs to be double checked
            if (use_canonical_event_reps) {
                CGAL::FT shift = offset_event - current_offset;
                Segment3SPtr e1o = PolyhedronTransformation::shiftEdge(edge_1, shift);
                if (e1o->is_degenerate()) {
                    // polyhedron split
                    continue;
                }
                Segment3SPtr e2o = PolyhedronTransformation::shiftEdge(edge_2, shift);
                if (e2o->is_degenerate()) {
                    // polyhedron split
                    continue;
                }
            }
#endif

            NodeSPtr node = Node::create();
            EdgeSplitEventSPtr event = EdgeSplitEvent::create();
            event->setStepID(step_id_);
            event->setNode(node);
            node->clear();
            node->setOffset(offset_event);
            node->setPoint(point);
            event->setEdge1(edge_1);
            event->setEdge2(edge_2);
            event->setEdgeOrientation(KernelWrapper::orientation(line(edge_1), line(edge_2)));

#ifndef CGAL_SS3_NO_SKELETON_DS
            SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_1->getData());
            SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_2->getData());
            node->addSheet(data_1->getSheet());
            node->addSheet(data_2->getSheet());

            if (facet_1_src == edge_2->getFacetL() ||
                    facet_1_src == edge_2->getFacetR()) {
                SkelVertexDataSPtr data_1_src = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_1->getVertexSrc()->getData());
                node->addArc(data_1_src->getArc());
            }
            if (facet_1_dst == edge_2->getFacetL() ||
                    facet_1_dst == edge_2->getFacetR()) {
                SkelVertexDataSPtr data_1_dst = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_1->getVertexDst()->getData());
                node->addArc(data_1_dst->getArc());
            }
#endif

            queue.push(event);
        }
    }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    CGAL_SS3_CORE_TRACE_V(16, "  " << filtered_candidates << " filtered, " << tested_candidates << " tests");
#endif

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Edge Split Events in: " << timer.time());
#endif
}

void SimpleStraightSkel::collectEdgeSplitEvents(PolyhedronSPtr polyhedron,
                                                const CGAL::FT& current_offset,
                                                const std::optional<CGAL::FT>& offset_future_bound,
                                                PQ& queue)
{
    return collectEdgeSplitEvents(polyhedron->edges(), polyhedron->edges(), polyhedron, true /*use canonical reps*/,
                                  current_offset, offset_future_bound, queue);
}

// @speed should keep a (sorted) list of vertex ---> facet of known contact events
// to filter how far we need to look and also avoid checking again multiple times
// something like SLS2
void SimpleStraightSkel::collectPierceEvents(const std::list<VertexSPtr>& vertices,
                                             const std::list<FacetSPtr>& facets,
                                             PolyhedronSPtr polyhedron,
                                             const CGAL::FT& current_offset,
                                             const std::optional<CGAL::FT>& offset_future_bound,
                                             PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(4, ">>> Collect Pierce Events [" << current_offset << "]");

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    unsigned int pierce_vertex_counter = 0;
    unsigned int filtered_candidates = 0;
    unsigned int tested_candidates = 0;
#endif

    for (VertexSPtr vertex : vertices) {
        CGAL_assertion(vertex->getID() != -1);

        // actual check
        if (isReflex(vertex)) {
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
            ++pierce_vertex_counter;
#endif

            for (FacetSPtr facet : facets) {
                // @todo contains_vertex is redundant with has_edge_to_facet?
                bool contains_vertex = false;
                for (VertexSPtr vertex_2 : facet->vertices()) {
                    if (vertex_2->getPoint() == vertex->getPoint()) {
                        contains_vertex = true;
                        break;
                    }
                }
                if (contains_vertex) {
                    continue;
                }

#ifndef CGAL_SS3_CHECK_PIERCE_AT_POP_TIME
                bool has_edge_to_facet = false;
                for (EdgeWPtr edge_wptr : vertex->edges()) {
                    if (EdgeSPtr edge = edge_wptr.lock()) {
                        FacetSPtr facet_src = edge->getFacetSrc();
                        FacetSPtr facet_dst = edge->getFacetDst();
                        if (facet == facet_src || facet == facet_dst) {
                            has_edge_to_facet = true;
                            break;
                        }
                    }
                }
                if (has_edge_to_facet) {
                    continue;
                }
#endif

                // shrinking, so the vertex must be on the backside of the plane
                if (KernelWrapper::side(facet->getPlane(), vertex->getPoint()) > 0) {
                    continue;
                }

                // If the facet is so far that even when shifting point and plane by the current
                // best lower bound on offset delta, the vertex has not crossed it yet, then we are done

                // not outside of the loop just because maybe one day this might get called
                // as the first collect function with an initial bound that gets updated...
                if (offset_future_bound) {
                    Point3SPtr shifted_pt = getFinalPoint(vertex, *offset_future_bound);
                    Plane3SPtr shifted_plane = getFinalPlane(facet, *offset_future_bound);

                    if (KernelWrapper::side(shifted_plane, shifted_pt) < 0) {
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
                        ++filtered_candidates;
#endif
                        continue;
                    } else {
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
                        ++tested_candidates;
#endif
                    }

                    // It would be tempting to use a bbox around the vertex and facet's trajectories,
                    // but if we did, then we would have to detect potential piercing events involving
                    // that facet after each modification of the facet, and this we obviously do not want.
                }

                CGAL_assertion(vertex->facets().size() >= 3);
                FacetSPtr fs[3];
                for (int i = 0; i < 3; ++i) {
                    FacetWPtr wf = *(std::next(vertex->facets().begin(), i));
                    CGAL_assertion(!wf.expired());
                    fs[i] = wf.lock();
                }

                Point3SPtr point;
                CGAL::FT offset_event;
                std::tie(point, offset_event) = intersectionPointAndTimeOffsetPlanes(facet, fs[0], fs[1], fs[2], current_offset, offset_future_bound);
                if (!point) {
                    continue;
                }

                CGAL_assertion(offset_event < current_offset && offset_event > offset_future_bound);

#ifndef CGAL_SS3_CHECK_PIERCE_AT_POP_TIME
                // Filter if the event point is on an edge (and a fortiori on a vertex)
                // as it will be a different kind of event
                FacetSPtr facet_offset = facet->clone();

                CGAL::FT shift = offset_event - current_offset;
                const CGAL::FT& speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
                Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->getPlane(), shift*speed);
                facet_offset->setPlane(offset_plane);

                // abusing the fact that vertices will have the same order in both facets
                std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
                std::list<VertexSPtr>::iterator it_v_offset = facet_offset->vertices().begin();
                while (it_v != facet->vertices().end()) {
                    VertexSPtr vertex = *it_v++;
                    VertexSPtr offset_vertex = *it_v_offset++;
                    Point3SPtr point_offset = PolyhedronTransformation::shiftPoint(vertex, shift);
                    offset_vertex->setPoint(point_offset);
                }

                // Note that this result could be meaningless if the offset facet
                // is not a simple polygon. However, if it's not simple, then some event
                // has happened before the pierce, and the pierce event - if whitelisted -
                // would be checked again later, thus it's safe to call.
                if (!SelfIntersection::isInsideWithRayShootingV2(point, facet_offset)) {
                    continue;
                }

                // @todo would be good if it could be merged with the function above...
                bool boundary_rejection = false;
                for (EdgeSPtr edge : facet_offset->edges()) {
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
                    continue;
                }
#endif

                NodeSPtr node = Node::create();
                PierceEventSPtr event = PierceEvent::create();
                event->setStepID(step_id_);
                event->setNode(node);
                node->clear();
#ifndef CGAL_SS3_NO_SKELETON_DS
                node->addArc(arc);
#endif
                node->setOffset(offset_event);
                node->setPoint(point);
                event->setFacet(facet);
                event->setVertex(vertex);

                queue.push(event);
            }
        }
    }

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    CGAL_SS3_CORE_TRACE_V(16, "  " << pierce_vertex_counter << " reflex vertices");
    CGAL_SS3_CORE_TRACE_V(16, "  " << filtered_candidates << " filtered, " << tested_candidates << " tests");
#endif

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(4, "  Sought Pierce Events in: " << timer.time());
#endif
}

void SimpleStraightSkel::collectPierceEvents(PolyhedronSPtr polyhedron,
                                             const CGAL::FT& current_offset,
                                             const std::optional<CGAL::FT>& offset_future_bound,
                                             PQ& queue)
{
    return collectPierceEvents(polyhedron->vertices(), polyhedron->facets(), polyhedron,
                               current_offset, offset_future_bound, queue);
}

void SimpleStraightSkel::printQueue(const PQ& queue) {
    CGAL_SS3_CORE_TRACE("-------------------------------------------------");
    CGAL_SS3_CORE_TRACE("--- Event queue (size = " << queue.size() << "; iter = " << step_id_ << ") ---");
    CGAL_SS3_CORE_TRACE("-------------------------------------------------");

    PQ duplicate_queue = queue;
    while (!duplicate_queue.empty()) {
        AbstractEventSPtr event = duplicate_queue.top();
        CGAL_SS3_CORE_TRACE("Event E" << event->getID()
                            << " T" << event->getType()
                            << " @ " << event->getOffset());
        if(event->isValid()) {
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
    CGAL_SS3_CORE_TRACE_CODE(for (CGAL::FT save_offset : save_offsets_))
    CGAL_SS3_CORE_TRACE_CODE(ss << " " << save_offset;)
    CGAL_SS3_CORE_TRACE(ss.str());
    CGAL_SS3_CORE_TRACE("-------------------------------------------------");
    CGAL_SS3_CORE_TRACE("-------------------------------------------------");
}

// this function checks the correctness of the local queue by computing the queue from scratch
// and checking that the first valid event is the same for both queues
bool SimpleStraightSkel::checkQueueCorrectness(const PQ& queue,
                                               PolyhedronSPtr polyhedron,
                                               const CGAL::FT& current_offset,
                                               const std::optional<CGAL::FT>& offset_future_bound) {
    CGAL_SS3_CORE_TRACE("checkQueueCorrectness()");

    // Compute a queue from scratch using collectEvents()
    PQ queue_from_scratch;
    collectEvents(polyhedron, current_offset, offset_future_bound, queue_from_scratch);

    // Duplicate the queue since we need to pop events
    PQ duplicate_queue = queue;

    // Tops should be the same
    auto purge_top = [&](PQ& q)
    {
        while (!q.empty()) {
            AbstractEventSPtr event = q.top();
            if (!event->isValid() ||
                isEventInThePast(current_offset, event) ||
                isEventObsolete(event) ||
                !isActualEvent(current_offset, event, polyhedron)) {
                q.pop();
            } else {
                break;
            }
        }
    };

    auto is_same_event = [](const AbstractEventSPtr& event_1, const AbstractEventSPtr& event_2) {
        if (event_1->getType() != event_2->getType()) {
            return false;
        }

        switch (event_1->getType()) {
            case AbstractEvent::SAVE_OFFSET_EVENT: {
                auto save_event_1 = std::dynamic_pointer_cast<SaveOffsetEvent>(event_1);
                auto save_event_2 = std::dynamic_pointer_cast<SaveOffsetEvent>(event_2);
                return *save_event_1 == *save_event_2;
            }
            case AbstractEvent::CONST_OFFSET_EVENT: {
                auto const_event_1 = std::dynamic_pointer_cast<ConstOffsetEvent>(event_1);
                auto const_event_2 = std::dynamic_pointer_cast<ConstOffsetEvent>(event_2);
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

        CGAL_SS3_CORE_TRACE("Seek event @ " << event_scratch->getOffset() << " Type " << event_scratch->getType());

        CGAL_SS3_CORE_TRACE("Event E" << event_scratch->getID()
                            << " T" << event_scratch->getType()
                            << " @ " << event_scratch->getOffset());
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

void SimpleStraightSkel::collectLocalEvents(PolyhedronSPtr polyhedron,
                                            const CGAL::FT& current_offset,
                                            const std::optional<CGAL::FT>& offset_future_bound,
                                            PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(2, "collectLocalEvents(" << current_offset << ")");

    // return collectEvents(polyhedron, current_offset, offset_future_bound, queue);

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    CGAL_SS3_CORE_TRACE_V(16, "Past bound = " << current_offset);
    CGAL_SS3_CORE_TRACE_IF(offset_future_bound, 4, "Initial future bound = " << *offset_future_bound);

    {
        std::list<EdgeSPtr> local_edges(post_op_edges_.begin(), post_op_edges_.end());

        CGAL_SS3_CORE_TRACE_V(8, "Local Edges for Vanish Events (" << local_edges.size() << ")");
        CGAL_SS3_CORE_TRACE_CODE(for(EdgeSPtr e : local_edges))
        CGAL_SS3_CORE_TRACE_V(8, "\t" << e->toString());

#ifdef CGAL_SS3_USE_GENERIC_VANISH_EVENT
        collectVanishEvents(local_edges, polyhedron, current_offset, offset_future_bound, queue);
#else
        collectEdgeEvents(local_edges, polyhedron,
                          current_offset, offset_future_bound, queue);
        collectEdgeMergeEvents(local_edges, polyhedron, false /*do not use canonical reps*/,
                              current_offset, offset_future_bound, queue);
        collectTriangleEvents(local_edges, polyhedron, false /*do not use canonical reps*/,
                              current_offset, offset_future_bound, queue);
        collectDblEdgeMergeEvents(local_edges, polyhedron,
                                  current_offset, offset_future_bound, queue);
        collectDblTriangleEvents(local_edges, polyhedron,
                                current_offset, offset_future_bound, queue);
        collectTetrahedronEvents(local_edges, polyhedron, false /*do not use canonical reps*/,
                                current_offset, offset_future_bound, queue);
#endif
    }

    // V-V events
    {
#if 1
        // for these three events, we need to look farther than the vertices that were involved
        // in the event because an event might add new edges to an edge cycle and so what was
        // before maybe a simple edge merge event now becomes a vertex event...
        std::list<VertexSPtr> local_vertices_VV(post_op_vertices_VV_.begin(),
                                                post_op_vertices_VV_.end());

        CGAL_SS3_CORE_TRACE_V(8, "Local Vertices for Vertex-Vertex Events (" << local_vertices_VV.size() << "):")
        CGAL_SS3_CORE_TRACE_CODE(std::stringstream ss;)
        CGAL_SS3_CORE_TRACE_CODE(for(VertexSPtr v : local_vertices_VV))
        CGAL_SS3_CORE_TRACE_CODE(ss << " " << v->getID();)
        CGAL_SS3_CORE_TRACE_V(8, ss.str());

        const bool use_canonical_reps = false;
#else
        const std::list<VertexSPtr>& local_vertices_VV = polyhedron->vertices();
        const bool use_canonical_reps = true;
#endif

        collectVertexEvents(local_vertices_VV, polyhedron, use_canonical_reps,
                            current_offset, offset_future_bound, queue);
        collectFlipVertexEvents(local_vertices_VV, polyhedron, use_canonical_reps,
                                current_offset, offset_future_bound, queue);
        collectSplitMergeEvents(local_vertices_VV, polyhedron, use_canonical_reps,
                                current_offset, offset_future_bound, queue);
    }

    // == POLYHEDRON SPLIT EVENTS ==
    {
#if 1
        // Polyhedron Split Events are not symmetrical, so we need to check both modified edges
        // as the first edge parameter, but also grab all the cases where a modified
        // edge is the second edge.

        std::list<EdgeSPtr> local_edges_EE(post_op_edges_.begin(), post_op_edges_.end());

        CGAL_SS3_CORE_TRACE_V(8, "Local Edges for Polyhedron Events (" << local_edges_EE.size() << ")");
        CGAL_SS3_CORE_TRACE_CODE(for(EdgeSPtr e : local_edges_EE))
        CGAL_SS3_CORE_TRACE_V(8, "\t" << e->toString());

        // this is the modified edges as 'edge_1'
        collectPolyhedronSplitEvents(local_edges_EE, polyhedron,
                                     current_offset, offset_future_bound, queue);

        // this is the modified edges as 'edge_2'
        // since we know that for a polyhedron split event, edge_2's LR facets
        // are the src and dst of edge_1, we can build a subset of all edges for edge_1s
        for (EdgeSPtr edge_2 : local_edges_EE) {
            std::set<EdgeSPtr> edges_1;
            for (FacetSPtr facet_2 : {edge_2->getFacetL(), edge_2->getFacetR()}) {
                for (VertexSPtr vertex_1 : facet_2->vertices()) {
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

            for (EdgeSPtr edge_1 : edges_1) {
              collectPolyhedronSplitEvent(edge_1, edge_2, polyhedron,
                                          current_offset, offset_future_bound, queue);
            }
        }


#else
        collectPolyhedronSplitEvents(polyhedron->edges(), polyhedron,
                                     current_offset, offset_future_bound, queue);
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
        CGAL_SS3_CORE_TRACE_CODE(for(VertexSPtr v : local_vertices_VF))
        CGAL_SS3_CORE_TRACE_CODE(ss << " " << v->getID());
        CGAL_SS3_CORE_TRACE_V(8, ss.str());
#else
        std::list<VertexSPtr> local_vertices_VF = polyhedron->vertices();
#endif

        // Pierce events must check with all faces, no choice about this
        collectPierceEvents(local_vertices_VF, polyhedron->facets(), polyhedron,
                            current_offset, offset_future_bound, queue);
    }

    // == SURFACE EVENTS ==
    {
        // Surface events are not symmetrical, so we need to check both modified edges
        // as the first edge parameter, but also grab all the cases where a modified
        // edge is the second edge.
#if 1
        std::list<EdgeSPtr> local_edges_EE(post_op_edges_.begin(), post_op_edges_.end());

        CGAL_SS3_CORE_TRACE_V(8, "Local Edges for Surface Events (" << local_edges_EE.size() << ")");
        CGAL_SS3_CORE_TRACE_CODE(for(EdgeSPtr e : local_edges_EE))
        CGAL_SS3_CORE_TRACE_V(8, "\t" << e->toString());

        // this is the modified edges as 'edge_1'
        collectSurfaceEvents(local_edges_EE, polyhedron,
                             current_offset, offset_future_bound, queue);

        // this is the modified edges as 'edge_2'
        for (EdgeSPtr edge_2 : local_edges_EE) {
            std::set<EdgeSPtr> edges_1;
            for (FacetSPtr facet_2 : {edge_2->getFacetL(), edge_2->getFacetR()}) {
                for (VertexSPtr vertex_1 : facet_2->vertices()) {
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

            for (EdgeSPtr edge_1 : edges_1) {
                collectSurfaceEvent(edge_1, edge_2, polyhedron,
                                    current_offset, offset_future_bound, queue);
            }
        }
#else
        collectSurfaceEvents(polyhedron->edges(), polyhedron,
                             current_offset, offset_future_bound, queue);
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
        CGAL_SS3_CORE_TRACE_CODE(for(EdgeSPtr e : local_edges_EE))
        CGAL_SS3_CORE_TRACE_V(8, "\t" << e->toString());

        const bool use_canonical_reps = false;
#else
        std::list<EdgeSPtr> local_edges_EE = polyhedron->edges();
        const bool use_canonical_reps = true;
#endif
        collectEdgeSplitEvents(local_edges_EE, polyhedron->edges(), polyhedron, use_canonical_reps,
                               current_offset, offset_future_bound, queue);
    }

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(2, "Sought All Local Events in: " << timer.time());
#endif

#ifdef CGAL_SS3_DEBUG_PRINT_QUEUE
    printQueue(queue);
#endif

    CGAL_postcondition(checkQueueCorrectness(queue, polyhedron, current_offset, offset_future_bound));
}

// two types of useless events:
// - events that are in the past:
//     offset > current_offset <--- values are negative and decreasing!
// - events that are stricly later than the current next tentative offset:
//     offset < curr_earliest_next_offset
void SimpleStraightSkel::collectEvents(PolyhedronSPtr polyhedron,
                                       const CGAL::FT& current_offset,
                                       const std::optional<CGAL::FT>& offset_future_bound,
                                       PQ& queue)
{
    CGAL_SS3_CORE_TRACE_V(2, "collectEvents(offset = " << current_offset << ")");

    AbstractEventSPtr result = AbstractEventSPtr();
    if (!polyhedron || polyhedron->facets().size() == 0) {
        return;
    }

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    CGAL_SS3_CORE_TRACE_V(8, "Past bound = " << current_offset);
    CGAL_SS3_CORE_TRACE_IF(offset_future_bound, 8, "Initial future bound = " << *offset_future_bound);

    // --- Vanish Events
#ifdef CGAL_SS3_USE_GENERIC_VANISH_EVENT
    collectVanishEvents(polyhedron, current_offset, offset_future_bound, queue);
#else
    collectEdgeEvents(polyhedron, current_offset, offset_future_bound, queue);
    collectEdgeMergeEvents(polyhedron, current_offset, offset_future_bound, queue);
    collectTriangleEvents(polyhedron, current_offset, offset_future_bound, queue);
    collectDblEdgeMergeEvents(polyhedron, current_offset, offset_future_bound, queue);
    collectDblTriangleEvents(polyhedron, current_offset, offset_future_bound, queue);
    collectTetrahedronEvents(polyhedron, current_offset, offset_future_bound, queue);
#endif

    // --- Contact Event
    collectVertexEvents(polyhedron, current_offset, offset_future_bound, queue);
    collectFlipVertexEvents(polyhedron, current_offset, offset_future_bound, queue);
    collectPolyhedronSplitEvents(polyhedron, current_offset, offset_future_bound, queue);
    collectSplitMergeEvents(polyhedron, current_offset, offset_future_bound, queue);

    // the next event types are particularly slow, so reduce the bound by doing them last
    // so other events lower the bound
    collectPierceEvents(polyhedron, current_offset, offset_future_bound, queue);
    collectSurfaceEvents(polyhedron, current_offset, offset_future_bound, queue);
    collectEdgeSplitEvents(polyhedron, current_offset, offset_future_bound, queue);

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_CORE_TRACE_V(2, "Sought All Events in: " << timer.time());
#endif

#ifdef CGAL_SS3_DEBUG_PRINT_QUEUE
    printQueue(queue);
#endif
}

AbstractEventSPtr SimpleStraightSkel::nextEvent(PQ& queue,
                                                const CGAL::FT& current_offset)
{
    if (queue.empty() && save_offsets_.empty())
        return { };

    AbstractEventSPtr event;
    CGAL::FT offset_next = 0;

    // If a save event is closest, delay building it in case a const event is even closer
    bool do_save = false;

    // purge the queue lazily as to avoid wasting time if we stop on the last save offset
    if (save_offsets_.empty()) {
        event = queue.top();
        offset_next = event->getOffset();
    } else {
        // If we have upcoming save events, compare
        CGAL::FT next_save_offset = save_offsets_.front();

        if (queue.empty()) {
            // queue is empty but save_offsets cannot be empty as well
            do_save = true;
            offset_next = save_offsets_.front();
        } else {
            // neither queue nor save_offsets are empty, take the earliest
            if (next_save_offset > queue.top()->getOffset()) { // save is strictly earlier
                do_save = true;
                offset_next = next_save_offset;
            } else {
                // save offsets exist, but are farther in the future than the next event
                event = queue.top();
                offset_next = event->getOffset();
            }
        }
    }

    // Tentative next event is not a save event, test its sanity.
    // Do this here because we don't want a const event to get created
    // because it's before a zombie event.
    if (!do_save) {
        if (!event->isValid()) {
            CGAL_SS3_CORE_TRACE_V(8, "Skipping invalid event E" << event->getID());
            queue.pop();
            return nextEvent(queue, current_offset);
        }

        if (isEventInThePast(current_offset, event)) {
            CGAL_SS3_CORE_TRACE_V(8, "Skipping event-in-the-past E" << event->getID());
            queue.pop();
            return nextEvent(queue, current_offset);
        }

        // @fixme the current "isObsolete" is only a sufficient condition: if the neighborhoods
        // have changed, then the event should be discarded.
        // But we likely also need the necessary condition, i.e. "if it is not obsolete, then
        // we know it is valid". Is it just doing _again_ all combinatorial checks to pop time?
        if (isEventObsolete(event)) {
            CGAL_SS3_CORE_TRACE_V(8, "Skipping obsolete event E" << event->getID());
            queue.pop();
            return nextEvent(queue, current_offset);
        }
    }

    CGAL_assertion(do_save || bool(event));
    CGAL_assertion(offset_next != 0);

    // Check if the next const event would be (strictly) closer
    CGAL::FT const_offset = util::Configuration::getInstance()->getDouble(
            "algo_3d_SimpleStraightSkel", "const_offset");
    if (const_offset != 0) {
        CGAL::FT next_const_offset = floor(CGAL::to_double(current_offset / const_offset) + 1.0) * const_offset;
        if (current_offset > next_const_offset && next_const_offset > offset_next) {
            do_save = false;
            return ConstOffsetEvent::create(next_const_offset);
        }
    }

    if (do_save) { // save event is the topmost
        save_offsets_.erase(save_offsets_.begin()); // @todo awkward
        return SaveOffsetEvent::create(offset_next);
    }

    CGAL_assertion(bool(event));

    // neither a ConstOffsetEvent nor a SaveOffsetEvent
    queue.pop();

    return event;
}

void SimpleStraightSkel::appendEventNode(NodeSPtr node) {
    std::list<ArcWPtr>::iterator it_a = node->arcs().begin();
    while (it_a != node->arcs().end()) {
        ArcWPtr arc_wptr = *it_a++;
        if (ArcSPtr arc = arc_wptr.lock()) {
            arc->setNodeDst(node);
            arc->setNodeDstListIt(util::weak_find(node->arcs().begin(), node->arcs().end(), arc_wptr));
        }
    }
    for (SheetWPtr sheet_wptr : node->sheets()) {
        if (SheetSPtr sheet = sheet_wptr.lock()) {
            sheet->addNode(node);
        }
    }
    skel_result_->addNode(node);
}

PolyhedronSPtr SimpleStraightSkel::shiftToEventOffset(PolyhedronSPtr polyhedron,
                                                      const CGAL::FT& start_offset,
                                                      const CGAL::FT& target_offset)
{
  const CGAL::FT shift = target_offset - start_offset;
  CGAL_precondition(!is_zero(shift));

  // @speed don't actually shift at all intermediate steps because we can do everything with base planes
  PolyhedronTransformation::shiftFacetsInPlace(polyhedron, shift);

#ifdef CGAL_SS3_DUMP_FILES
  // below will have degeneracies since we have not yet treated the event
  static int shift_id = -1;
  db::_3d::OBJFile::save("results/shift_" + std::to_string(++shift_id) + ".obj",
                         polyhedron,
                         false /*do not triangulate*/);
#endif

  return polyhedron;
}

// This 'handle' is in fact more akin to a collect, but the interesting point
// is that it happens after pop time
SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleSaveOffsetEvent(SaveOffsetEventSPtr event,
                                          const CGAL::FT& current_offset,
                                          PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#########  Handle Save Event  ##########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    skel_result_->addEvent(event);

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    bool res = savePolyhedron(polyhedron, event_offset, false /*attempt untilting*/);

#ifdef CGAL_SS3_DUMP_FILES
    db::_3d::OBJFile::save("results/last_save.obj", polyhedron);
#endif

    return (res ? EventStatus::EVENT_HANDLED : EventStatus::EVENT_NOT_HANDLED);
}

// This 'handle' is in fact more akin to a collect, but the interesting point
// is that it happens after pop time
SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleConstOffsetEvent(ConstOffsetEventSPtr event,
                                           const CGAL::FT& current_offset,
                                           PolyhedronSPtr polyhedron)
{
  CGAL_SS3_CORE_TRACE_V(4, "########################################");
  CGAL_SS3_CORE_TRACE_V(4, "#########  Handle Const Event  #########");
  CGAL_SS3_CORE_TRACE_V(4, "########################################");

  const CGAL::FT& event_offset = event->getOffset();
  polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

#ifndef CGAL_SS3_NO_SKELETON_DS
    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    bool screenshot_on_const_offset_event =
            util::Configuration::getInstance()->getBool(
            "algo_3d_SimpleStraightSkel", "screenshot_on_const_offset_event");
    if (controller_ && screenshot_on_const_offset_event) {
        controller_->screenshot();
    }

    return EventStatus::EVENT_HANDLED;
}

// This 'handle' is in fact more akin to a collect, but the interesting point
// is that it happens after pop time.
// This function does NOT shift the polyhedron, it is for the "real" handler
// to deal with it.
// This function might not do anything, for example if the vanish event is in fact
// escalated as a contact event.
SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleVanishEvent(VanishEventSPtr event,
                                      const CGAL::FT& current_offset,
                                      const std::optional<CGAL::FT>& offset_future_bound,
                                      PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "####### Tentative Vanish Event  ########");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

    NodeSPtr node = event->getNode();
    EdgeSPtr edge = event->getEdge();
    Point3SPtr point = node->getPoint();

    // @todo would nice:
    // - to avoid redundant checks (e.g. isTriangle multiple times)
    // - to avoid code duplication with collectXYZEvents()
    //
    // To avoid redundant code, it could be re-ordered by descending amount of collapsing edges (6 --> 1)?
    // Might not speed things up if majority of events are a small number of edge collapses.
    // Anyway, vanish events are cheap to collect, and cheap to process.

    // The fake infinite loops are only there to enable 'break's such that the code
    // is identical to that of the respective collectXYZEvent() functions

    // Edge Event
    for(;;)
    {
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (isTriangle(facet_l, edge) || isTriangle(facet_r, edge)) {
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
        for (EdgeSPtr edge_2 : edges_2) {
            bool split_event_current_1_b = true;
            bool split_event_current_2_b = true;
            bool split_event_current_3_b = true;

            FacetSPtr facet_l2 = edge_2->getFacetL();
            FacetSPtr facet_r2 = edge_2->getFacetR();
            FacetSPtr facet_2_src = edge_2->getFacetSrc();
            FacetSPtr facet_2_dst = edge_2->getFacetDst();

            Plane3SPtr plane_l2 = facet_l2->getPlane();

            const CGAL::FT& l2a = plane_l2->a();
            const CGAL::FT& l2b = plane_l2->b();
            const CGAL::FT& l2c = plane_l2->c();
            const CGAL::FT& l2d = plane_l2->d();
            const CGAL::FT& speed_l2 = std::dynamic_pointer_cast<SkelFacetData>(facet_l2->getData())->getSpeed();
            CGAL::FT t = (l2a * point->x() + l2b * point->y() + l2c * point->z() + l2d) / speed_l2;

            CGAL_assertion_code(Plane3SPtr plane_r2 = facet_r2->getPlane();)
            CGAL_assertion_code(const CGAL::FT& r2a = plane_r2->a();)
            CGAL_assertion_code(const CGAL::FT& r2b = plane_r2->b();)
            CGAL_assertion_code(const CGAL::FT& r2c = plane_r2->c();)
            CGAL_assertion_code(const CGAL::FT& r2d = plane_r2->d();)
            CGAL_assertion_code(const CGAL::FT& speed_r2 = std::dynamic_pointer_cast<SkelFacetData>(facet_r2->getData())->getSpeed();)
            CGAL_assertion_code(CGAL::FT lt2 = (l2a * point->x() + l2b * point->y() + l2c * point->z() + l2d) / speed_l2);
            CGAL_assertion_code(CGAL::FT rt2 = (r2a * point->x() + r2b * point->y() + r2c * point->z() + r2d) / speed_r2);
            CGAL_assertion(lt2 == rt2);

            if (is_positive(t)) {
                // can 'break' directly because it's the same value for all 'edge_2's
                break;
            }

            if (!check_bisector(edge_2, facet_r2, t/*rt2*/, facet_2_src, point)) {
                continue;
            }

            if (!check_bisector(edge_2, facet_l2, t/*lt2*/, facet_2_dst, point)) {
                continue;
            }

            const bool split_event_current_b = (split_event_current_1_b &&
                                                split_event_current_2_b &&
                                                split_event_current_3_b);
            if (split_event_current_b) {
                split_event = true;
                break;
            }
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
        edge_event->setStepID(step_id_);
        edge_event->setNode(node);
        edge_event->setEdge(edge);

        return handleEdgeEvent(edge_event, current_offset, offset_future_bound, polyhedron);
    }

    // EdgeMergeEvent
    for(;;)
    {
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (isTriangle(facet_l, edge) || isTriangle(facet_r, edge)) {
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
        edge_merge_event->setStepID(step_id_);
        edge_merge_event->setNode(node);
        edge_merge_event->setFacet(facet);
        edge_merge_event->setEdge1(edge_1);
        edge_merge_event->setEdge2(edge_2);

        return handleEdgeMergeEvent(edge_merge_event, current_offset, offset_future_bound, polyhedron);
    }

    // TriangleEvent
    for(;;)
    {
        if (isTetrahedron(edge)) {
            // tetrahedron event
            CGAL_SS3_CORE_TRACE_V(8, "Not a Triangle Event (Tetrahedron)");
            break;
        }

        FacetSPtr facet;
        if (isTriangle(edge->getFacetL(), edge)) {
            facet = edge->getFacetL();
        } else if (isTriangle(edge->getFacetR(), edge)) {
            facet = edge->getFacetR();
        } else {
            CGAL_SS3_CORE_TRACE_V(8, "Not a Triangle Event (not triangle)");
            break;
        }

        bool dbl_triangle_event = false;
        EdgeSPtr edge_tmp = edge;
        for (unsigned int i = 0; i < 3; i++) {
            FacetSPtr facet_tmp_l = edge_tmp->getFacetL();
            FacetSPtr facet_tmp_r = edge_tmp->getFacetR();
            if (facet_tmp_l && facet_tmp_r) {
                if (isTriangle(facet_tmp_l, edge_tmp) &&
                        isTriangle(facet_tmp_r, edge_tmp)) {
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
        triangle_event->setStepID(step_id_);
        triangle_event->setNode(node);
        triangle_event->setFacet(facet);
        triangle_event->setEdgeBegin(edge);

        return handleTriangleEvent(triangle_event, current_offset, offset_future_bound, polyhedron);
    }

    // DblEdgeMergeEvent
    for(;;)
    {
        // At pop time, edge is degenerate, but by default isReflex() does a symbolic
        // shift into the future. However, here we want to know if the edge was reflex
        // BEFORE it became degenerate.
        if (!isReflex(edge, false /*shift into the past*/)) {
            CGAL_SS3_CORE_TRACE_V(8, "Not a DblEdgeMerge Event (not reflex)");
            break;
        }

        bool is_dbl_edge_merge_event = false;
        FacetSPtr facet_1;
        EdgeSPtr edge_11;
        EdgeSPtr edge_12;
        FacetSPtr facet_2;
        EdgeSPtr edge_21;
        EdgeSPtr edge_22;
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
        dbl_edge_merge_event->setStepID(step_id_);
        dbl_edge_merge_event->setNode(node);
        dbl_edge_merge_event->setFacet1(facet_1);
        dbl_edge_merge_event->setEdge11(edge_11);
        dbl_edge_merge_event->setEdge12(edge_12);
        dbl_edge_merge_event->setFacet2(facet_2);
        dbl_edge_merge_event->setEdge21(edge_21);
        dbl_edge_merge_event->setEdge22(edge_22);

        return handleDblEdgeMergeEvent(dbl_edge_merge_event, current_offset, offset_future_bound, polyhedron);
    }

    // DblTriangleEvent
    for(;;)
    {
        if (isTetrahedron(edge)) {
            CGAL_SS3_CORE_TRACE_V(8, "Not a DblTriangleMerge Event (Tetrahedron)");
            break;
        }
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (!facet_l || !facet_r) {
            CGAL_SS3_CORE_TRACE_V(8, "Not a DblTriangleMerge Event (neighborhood)");
            break;
        }
        if (!(isTriangle(facet_l, edge) &&
                isTriangle(facet_r, edge))) {
            CGAL_SS3_CORE_TRACE_V(8, "Not a DblTriangleMerge Event (not triangles)");
            break;
        }

        // if here, it's a double triangle event
        DblTriangleEventSPtr dbl_triangle_event = DblTriangleEvent::create();
        dbl_triangle_event->setStepID(step_id_);
        dbl_triangle_event->setNode(node);
        dbl_triangle_event->setEdge(edge);

        return handleDblTriangleEvent(dbl_triangle_event, current_offset, offset_future_bound, polyhedron);
    }

    // TetrahedronEvent
    for(;;)
    {
        if (!isTetrahedron(edge)) {
            CGAL_SS3_CORE_TRACE_V(8, "Not a Tetrahedron Event");
            break;
        }

        // if here, it's a tetrahedron event
        TetrahedronEventSPtr tetrahedron_event = TetrahedronEvent::create();
        tetrahedron_event->setStepID(step_id_);
        tetrahedron_event->setNode(node);
        tetrahedron_event->setEdgeBegin(edge);

        return handleTetrahedronEvent(tetrahedron_event, current_offset, offset_future_bound, polyhedron);
    }

    return EventStatus::NON_EVENT;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleEdgeEvent(EdgeEventSPtr event,
                                    const CGAL::FT& current_offset,
                                    const std::optional<CGAL::FT>& offset_future_bound,
                                    PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#########  Handle Edge Event  ##########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    EdgeSPtr edge = event->getEdge();
    SkelVertexDataSPtr data_src = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexSrc()->getData());
    VertexSPtr vertex_src_offset = data_src->getOffsetVertex();
    SkelVertexDataSPtr data_dst = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexDst()->getData());
    VertexSPtr vertex_dst_offset = data_dst->getOffsetVertex();
    EdgeSPtr edge_offset = vertex_src_offset->findEdge(vertex_dst_offset);
    vertex_src_offset->setPoint(node->getPoint());
    vertex_dst_offset->setPoint(node->getPoint());

    CGAL_SS3_CORE_TRACE_V(4,"Node:\n" << node->toString());
    CGAL_SS3_CORE_TRACE_V(4,"Edge:\n" << edge->toString());

    // grab vertices and facets counter clockwise around node
    VertexSPtr vertices[4];
    for (unsigned int i = 0; i < 4; i++) {
        vertices[i] = VertexSPtr();
    }

    vertices[0] = vertex_src_offset->prev(edge_offset->getFacetL());
    vertices[1] = vertex_src_offset->next(edge_offset->getFacetR());
    vertices[2] = vertex_dst_offset->prev(edge_offset->getFacetR());
    vertices[3] = vertex_dst_offset->next(edge_offset->getFacetL());

    FacetSPtr facets[4];
    for (unsigned int i = 0; i < 4; i++) {
        facets[i] = FacetSPtr();
    }
    facets[1] = edge_offset->getFacetR();
    facets[3] = edge_offset->getFacetL();
    facets[0] = facets[3]->next(vertex_src_offset);
    facets[2] = facets[1]->next(vertex_dst_offset);

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

        VertexSPtr vertex_src_clone = vertex_src_offset->clone();
        VertexSPtr vertex_dst_clone = vertex_dst_offset->clone();
        EdgeSPtr edge_no_flip = Edge::create(vertex_src_clone, vertex_dst_clone);
        EdgeSPtr edges[4];
        for (unsigned int i = 0; i < 4; i++) {
            edges[i] = EdgeSPtr();
        }
        edges[0] = Edge::create(vertex_src_clone, vertices[0]->clone());
        edges[1] = Edge::create(vertex_src_clone, vertices[1]->clone());
        edges[2] = Edge::create(vertex_dst_clone, vertices[2]->clone());
        edges[3] = Edge::create(vertex_dst_clone, vertices[3]->clone());
        std::vector<FacetSPtr> facets_clone(4);
        for (unsigned int i = 0; i < 4; i++) {
            facets_clone[i] = Facet::create();
            facets_clone[i]->setPlane(facets[i]->getPlane());
            facets_clone[i]->setBasePlane(facets[i]->getBasePlane()); // @todo useless but just in case for now...
            if (facets[i]->hasData()) {
                SkelFacetDataSPtr data_clone = SkelFacetData::create(facets_clone[i]);
                data_clone->setSpeed(std::dynamic_pointer_cast<SkelFacetData>(
                        facets[i]->getData())->getSpeed());
            }
        }

        edge_no_flip->setFacetR(facets_clone[1]);
        edge_no_flip->setFacetL(facets_clone[3]);

        // In edge events, we have unbounded faces (degree 1 vertices)
        facets_clone[3]->addEdge(edge_no_flip);
        facets_clone[1]->addEdge(edge_no_flip);
        for (unsigned int i = 0; i < 4; i++) {
            edges[i]->setFacetR(facets_clone[(i+3)%4]);
            edges[i]->setFacetL(facets_clone[i]);
            facets_clone[(i+3)%4]->addEdge(edges[i]);
            facets_clone[i]->addEdge(edges[i]);
        }

        PolyhedronSPtr polyhedron_no_flip = Polyhedron::create(facets_clone);
        PolyhedronSPtr polyhedron_no_flip_offset = PolyhedronTransformation::shiftFacets(polyhedron_no_flip, -1.0);

        if (polyhedron_no_flip_offset && !SelfIntersection::hasSelfIntersectingSurface(polyhedron_no_flip_offset)) {
            not_flipped_valid = true;
        }
    }

    CGAL_SS3_CORE_TRACE_V(16, "not_flipped_valid = " << not_flipped_valid);

    bool flipped_valid = false;
#if 1
    // if one fails, the other one must be valid and we can avoid self-intersections
    if (!not_flipped_valid) { // redundant but simply for clarity
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

        VertexSPtr vertex_src_clone = vertex_src_offset->clone();
        VertexSPtr vertex_dst_clone = vertex_dst_offset->clone();
        EdgeSPtr edge_flipped = Edge::create(vertex_src_clone, vertex_dst_clone);
        EdgeSPtr edges[4];
        for (unsigned int i = 0; i < 4; i++) {
            edges[i] = EdgeSPtr();
        }
        edges[0] = Edge::create(vertex_dst_clone, vertices[0]->clone());
        edges[1] = Edge::create(vertex_src_clone, vertices[1]->clone());
        edges[2] = Edge::create(vertex_src_clone, vertices[2]->clone());
        edges[3] = Edge::create(vertex_dst_clone, vertices[3]->clone());
        std::vector<FacetSPtr> facets_clone(4);
        for (unsigned int i = 0; i < 4; i++) {
            facets_clone[i] = Facet::create();
            facets_clone[i]->setPlane(facets[i]->getPlane());
            facets_clone[i]->setBasePlane(facets[i]->getBasePlane());
            if (facets[i]->hasData()) {
                SkelFacetDataSPtr data_clone = SkelFacetData::create(facets_clone[i]);
                data_clone->setSpeed(std::dynamic_pointer_cast<SkelFacetData>(
                        facets[i]->getData())->getSpeed());
            }
        }

        edge_flipped->setFacetR(facets_clone[2]);
        edge_flipped->setFacetL(facets_clone[0]);

        facets_clone[0]->addEdge(edge_flipped);
        facets_clone[2]->addEdge(edge_flipped);
        for (unsigned int i = 0; i < 4; i++) {
            edges[i]->setFacetR(facets_clone[(i+3)%4]);
            edges[i]->setFacetL(facets_clone[i]);
            facets_clone[(i+3)%4]->addEdge(edges[i]);
            facets_clone[i]->addEdge(edges[i]);
        }

        PolyhedronSPtr polyhedron_flipped = Polyhedron::create(facets_clone);
        PolyhedronSPtr polyhedron_flipped_offset = PolyhedronTransformation::shiftFacets(polyhedron_flipped, -1.0);

        if (polyhedron_flipped_offset && !SelfIntersection::hasSelfIntersectingSurface(polyhedron_flipped_offset)) {
            flipped_valid = true;
        }
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
        } else if (edge_event_ == 2) {
            // flip_when_possible
            flip_edge = true;
        } else if (edge_event_ == 3) {
            // sphere
            VertexSPtr vertex_c = Vertex::create(node->getPoint());
            EdgeSPtr edges[4];
            for (unsigned int i = 0; i < 4; i++) {
                edges[i] = Edge::create(vertex_c, vertices[i]->clone());
            }
            std::vector<FacetSPtr> facets_c(4);
            for (unsigned int i = 0; i < 4; i++) {
                facets_c[i] = Facet::create();
                facets_c[i]->setPlane(facets[i]->getPlane());
                facets_c[i]->setBasePlane(facets[i]->getBasePlane());
            }
            for (unsigned int i = 0; i < 4; i++) {
                edges[i]->setFacetR(facets_c[(i+3)%4]);
                edges[i]->setFacetL(facets_c[i]);
                facets_c[(i+3)%4]->addEdge(edges[i]);
                facets_c[i]->addEdge(edges[i]);
            }
            PolyhedronSPtr polyhedron_sphere = Polyhedron::create(facets_c);
            SphereVertexSplitterSPtr splitter = SphereVertexSplitter::create();
            splitter->splitVertex(vertex_c);
            if (facets_c[0]->findEdge(facets_c[2])) {
                flip_edge = true;
            } else if (facets_c[1]->findEdge(facets_c[3])) {
                flip_edge = false;
            } else {
                throw std::runtime_error("Error: Not able to handle EdgeEvent (1).");
            }
        }
    } else {
        throw std::runtime_error("Error: not able to handle EdgeEvent (2).");
    }

    CGAL_SS3_CORE_TRACE_V(16, "flip_edge = " << flip_edge);

    if (flip_edge) {
        EdgeSPtr edges[4];
        for (unsigned int i = 0; i < 4; i++) {
            edges[i] = EdgeSPtr();
        }
        edges[0] = edge_offset->next(vertex_src_offset);
        edges[1] = edges[0]->next(vertex_src_offset);
        edges[2] = edge_offset->next(vertex_dst_offset);
        edges[3] = edges[2]->next(vertex_dst_offset);

        facets[3]->removeVertex(vertex_src_offset);
        facets[2]->addVertex(vertex_src_offset);
        facets[1]->removeVertex(vertex_dst_offset);
        facets[0]->addVertex(vertex_dst_offset);

        facets[1]->removeEdge(edge_offset);
        facets[3]->removeEdge(edge_offset);
        edge_offset->setFacetL(facets[0]);
        edge_offset->setFacetR(facets[2]);
        facets[0]->addEdge(edge_offset);
        facets[2]->addEdge(edge_offset);

        if (edges[0]->getVertexSrc() == vertex_src_offset) {
            edges[0]->replaceVertexSrc(vertex_dst_offset);
        } else if (edges[0]->getVertexDst() == vertex_src_offset) {
            edges[0]->replaceVertexDst(vertex_dst_offset);
        }
        if (edges[2]->getVertexSrc() == vertex_dst_offset) {
            edges[2]->replaceVertexSrc(vertex_src_offset);
        } else if (edges[2]->getVertexDst() == vertex_dst_offset) {
            edges[2]->replaceVertexDst(vertex_src_offset);
        }

        if (offset_future_bound) {
            vertex_src_offset->final_point_ = nullptr;
            vertex_dst_offset->final_point_ = nullptr;
        }

        // @fixme duplicates with below
        post_op_vertices_VV_ = {{ vertex_src_offset, vertex_dst_offset }};

        // now, in the facets that have grown in size, we also need to collect
        // vertices that might create new events because they are now separated-enough
        // from other vertices in the cycle
        post_op_vertices_VV_.insert(vertex_src_offset->prev(facets[0]));
        post_op_vertices_VV_.insert(vertex_dst_offset->next(facets[0]));
        post_op_vertices_VV_.insert(vertex_src_offset->next(facets[2]));
        post_op_vertices_VV_.insert(vertex_dst_offset->prev(facets[2]));
        CGAL_assertion(post_op_vertices_VV_.size() == 6);

        // if there was no flip, then unmodified vertices should not have new events
    }

    // update arcs and sheets
#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(edge_offset->getData());
    edge_data->setSheet(SheetSPtr());

    data_src = std::dynamic_pointer_cast<SkelVertexData>(vertex_src_offset->getData());
    data_dst = std::dynamic_pointer_cast<SkelVertexData>(vertex_dst_offset->getData());
    data_src->setNode(node);
    data_dst->setNode(node);

    ArcSPtr arc_src = createArc(vertex_src_offset);
    skel_result_->addArc(arc_src);
    ArcSPtr arc_dst = createArc(vertex_dst_offset);
    skel_result_->addArc(arc_dst);
    SheetSPtr sheet = createSheet(edge_offset);
    skel_result_->addSheet(sheet);

    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_ = {{ vertex_src_offset, vertex_dst_offset }};
    post_op_edges_ = {{ edge_offset,
                        edge_offset->next(vertex_src_offset),
                        edge_offset->prev(vertex_src_offset),
                        edge_offset->next(vertex_dst_offset),
                        edge_offset->prev(vertex_dst_offset) }};
    post_op_facets_ = {{ facets[0], facets[1], facets[2], facets[3] }};
    CGAL_postcondition(post_op_vertices_.size() == 2 &&
                       post_op_edges_.size() == 5 &&
                       post_op_facets_.size() == 4);

    for (EdgeSPtr poe : post_op_edges_) {
        post_op_vertices_pierce_.insert(poe->getVertexSrc());
        post_op_vertices_pierce_.insert(poe->getVertexDst());
    }
    CGAL_postcondition(post_op_vertices_pierce_.size() == 6);

    return EventStatus::EVENT_HANDLED;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleEdgeMergeEvent(EdgeMergeEventSPtr event,
                                         const CGAL::FT& current_offset,
                                         const std::optional<CGAL::FT>& offset_future_bound,
                                         PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Edge Merge Event  #######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    appendEventNode(event->getNode());

    SkelFacetDataSPtr facet_data = std::dynamic_pointer_cast<SkelFacetData>(
                event->getFacet()->getData());
    FacetSPtr facet = facet_data->getOffsetFacet();
    SkelEdgeDataSPtr edge_data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
                event->getEdge1()->getData());
    SkelEdgeDataSPtr edge_data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
                event->getEdge2()->getData());
    EdgeSPtr edge_1 = edge_data_1->getOffsetEdge();
    EdgeSPtr edge_2 = edge_data_2->getOffsetEdge();

    FacetSPtr facet_l = edge_1->getFacetL();
    FacetSPtr facet_r = edge_1->getFacetR();
    EdgeSPtr edge_toremove_1 = edge_1->next(facet);
    EdgeSPtr edge_toremove_2 = edge_toremove_1->next(facet);

    assert(edge_2 == edge_toremove_2->next(facet));

    VertexSPtr vertex = edge_toremove_1->dst(facet);
    vertex->setPoint(event->getNode()->getPoint());
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

    if (offset_future_bound) {
        vertex->final_point_ = nullptr;
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
    vertex_data->setNode(event->getNode());
    ArcSPtr arc = createArc(vertex);
    skel_result_->addArc(arc);
    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_ = {{ vertex }};
    post_op_edges_ = {{ edge_1, edge_b1, edge_b, edge_b2 }};
    post_op_facets_ = {{ facet_l, facet_r }};
    for (FacetWPtr wf : vertex->facets()) { post_op_facets_.insert(wf.lock()); }
    CGAL_postcondition(post_op_vertices_.size() == 1 && post_op_edges_.size() == 4 && post_op_facets_.size() == 4);

    // faces are smaller so nothing from unmodified vertices
    post_op_vertices_VV_ = {{ vertex }};

    CGAL_assertion(!isReflex(vertex)); // just to see the configurations where this could not be the case
    post_op_vertices_pierce_.clear();

    // since all faces are getting smaller, we don't need to check unmodified edges
    post_op_edges_edgesplit_ = {{ edge_1, edge_b1, edge_b, edge_b2 }};
    CGAL_postcondition(post_op_edges_edgesplit_.size() == 4);

    return EventStatus::EVENT_HANDLED;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleTriangleEvent(TriangleEventSPtr event,
                                        const CGAL::FT& current_offset,
                                        const std::optional<CGAL::FT>& offset_future_bound,
                                        PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#######  Handle Triangle Event  ########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    appendEventNode(event->getNode());

    VertexSPtr vertices[3];
    event->getVertices(vertices);

    VertexSPtr vertices_offset[3];
    for (unsigned int i = 0; i < 3; i++) {
        SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(
                vertices[i]->getData());
        vertices_offset[i] = data->getOffsetVertex();
    }
    SkelFacetDataSPtr facet_data = std::dynamic_pointer_cast<SkelFacetData>(
            event->getFacet()->getData());
    FacetSPtr facet_offset = facet_data->getOffsetFacet();

    CGAL_SS3_CORE_TRACE_V(4, "Facet: " << facet_offset->getID());
    CGAL_SS3_CORE_TRACE_V(4, "VS:\n" << vertices[0]->toString() << "\n"
                                     << vertices[1]->toString() << "\n"
                                     << vertices[2]->toString());
    CGAL_SS3_CORE_TRACE_V(4, "VSO:\n" << vertices_offset[0]->toString() << "\n"
                                      << vertices_offset[1]->toString() << "\n"
                                      << vertices_offset[2]->toString());

    if (facet_offset->vertices().size() == 3) {
        polyhedron->removeFacet(facet_offset);
    }
    for (unsigned int i = 0; i < 3; i++) {
        EdgeSPtr edge = vertices_offset[i]->findEdge(vertices_offset[(i+1)%3]);
        if (edge->getFacetL()) {
            edge->getFacetL()->removeEdge(edge);
        }
        if (edge->getFacetR()) {
            edge->getFacetR()->removeEdge(edge);
        }
        polyhedron->removeEdge(edge);
    }
    VertexSPtr vertex_offset = Vertex::create(event->getNode()->getPoint());
    for (unsigned int i = 0; i < 3; i++) {
        EdgeSPtr edge_offset = vertices_offset[i]->edges().front().lock();
        if (edge_offset->getVertexSrc() == vertices_offset[i]) {
            edge_offset->replaceVertexSrc(vertex_offset);
        } else if (edge_offset->getVertexDst() == vertices_offset[i]) {
            edge_offset->replaceVertexDst(vertex_offset);
        }
        edge_offset->getFacetL()->removeVertex(vertices_offset[i]);
        edge_offset->getFacetR()->removeVertex(vertices_offset[i]);
        facet_offset->removeVertex(vertices_offset[i]);
        polyhedron->removeVertex(vertices_offset[i]);
        if (!edge_offset->getFacetL()->containsVertex(vertex_offset)) {
            edge_offset->getFacetL()->addVertex(vertex_offset);
        }
        if (!edge_offset->getFacetR()->containsVertex(vertex_offset)) {
            edge_offset->getFacetR()->addVertex(vertex_offset);
        }
    }
    polyhedron->addVertex(vertex_offset);

    if (offset_future_bound) {
        vertex_offset->final_point_ = nullptr;
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelVertexDataSPtr data_offset = SkelVertexData::create(vertex_offset);
    data_offset->setNode(event->getNode());
    ArcSPtr arc = createArc(vertex_offset);
    skel_result_->addArc(arc);
    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_ = {{ vertex_offset }};
    for (EdgeWPtr we : vertex_offset->edges()) {
        post_op_edges_.insert(EdgeSPtr(we.lock()));
    }
    for (FacetWPtr wf : vertex_offset->facets()) {
        post_op_facets_.insert(wf.lock());
    }
    CGAL_postcondition(post_op_vertices_.size() == 1 && post_op_edges_.size() == 3 && post_op_facets_.size() == 3);

    // faces are smaller so nothing from unmodified vertices
    post_op_vertices_VV_ = {{ vertex_offset }};

    // Assume a triangle event with a facet whose incident edges are all reflex.
    // At a reflex edge, an epsilon offset is a growth of the facet.
    // If all edges are reflex, the facet will grow and there could not have been
    // a triangle event. So at least one edge is convex and the resulting vertex
    // is not reflex. Since it is not reflex, we don't need to look at its pierce events.
    //
    // @todo checking for reflexness is (relatively) so cheap that we should insert it
    // anyway, in case the reasoning above is wrong!
    CGAL_assertion(!isReflex(vertex_offset));
    post_op_vertices_pierce_.clear();

    // faces are getting smaller so no need to check unmodified edges
    post_op_edges_edgesplit_ = post_op_edges_;
    CGAL_postcondition(post_op_edges_edgesplit_.size() == 3);

    return EventStatus::EVENT_HANDLED;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleDblEdgeMergeEvent(DblEdgeMergeEventSPtr event,
                                            const CGAL::FT& current_offset,
                                            const std::optional<CGAL::FT>& offset_future_bound,
                                            PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#######  Handle Dbl Edge Event  ########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    appendEventNode(event->getNode());

    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge11()->getData());
    EdgeSPtr edge_offset_11 = edge_data->getOffsetEdge();
    edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge12()->getData());
    EdgeSPtr edge_offset_12 = edge_data->getOffsetEdge();
    edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge21()->getData());
    EdgeSPtr edge_offset_21 = edge_data->getOffsetEdge();
    edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge22()->getData());
    EdgeSPtr edge_offset_22 = edge_data->getOffsetEdge();

    VertexSPtr vertices[4];
    event->getVertices(vertices);
    VertexSPtr vertices_offset[4];
    for (unsigned int i = 0; i < 4; i++) {
        SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                vertices[i]->getData());
        vertices_offset[i] = vertex_data->getOffsetVertex();
    }

    EdgeSPtr edges[4];
    event->getEdges(edges);
    EdgeSPtr edges_offset[4];
    for (unsigned int i = 0; i < 4; i++) {
        SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                edges[i]->getData());
        edges_offset[i] = edge_data->getOffsetEdge();
    }

    // @fixme what's going on here? Aren't we deleting the offset edges multiple times?
    // What's the point of changing vertices of edge_offset_xy if we are deleting them anyway
    for (unsigned int i = 0; i < 4; i++) {
        EdgeSPtr edge = edges_offset[i];
        edge->getFacetL()->removeEdge(edge);
        edge->getFacetR()->removeEdge(edge);
        polyhedron->removeEdge(edge);
    }

    if (edge_offset_11->getVertexDst() == vertices_offset[0]) {
        if (edge_offset_12->getVertexSrc() == vertices_offset[2]) {
            edge_offset_11->replaceVertexDst(edge_offset_12->getVertexDst());
        } else {
            edge_offset_11->replaceVertexDst(edge_offset_12->getVertexSrc());
        }
    } else {
        if (edge_offset_12->getVertexSrc() == vertices_offset[2]) {
            edge_offset_11->replaceVertexSrc(edge_offset_12->getVertexDst());
        } else {
            edge_offset_11->replaceVertexSrc(edge_offset_12->getVertexSrc());
        }
    }
    edge_offset_12->getFacetL()->removeEdge(edge_offset_12);
    edge_offset_12->getFacetR()->removeEdge(edge_offset_12);
    polyhedron->removeEdge(edge_offset_12);
    if (edge_offset_21->getVertexDst() == vertices_offset[1]) {
        if (edge_offset_22->getVertexSrc() == vertices_offset[3]) {
            edge_offset_21->replaceVertexDst(edge_offset_22->getVertexDst());
        } else {
            edge_offset_21->replaceVertexDst(edge_offset_22->getVertexSrc());
        }
    } else {
        if (edge_offset_22->getVertexSrc() == vertices_offset[3]) {
            edge_offset_21->replaceVertexSrc(edge_offset_22->getVertexDst());
        } else {
            edge_offset_21->replaceVertexSrc(edge_offset_22->getVertexSrc());
        }
    }
    edge_offset_22->getFacetL()->removeEdge(edge_offset_22);
    edge_offset_22->getFacetR()->removeEdge(edge_offset_22);
    polyhedron->removeEdge(edge_offset_22);
    for (unsigned int i = 0; i < 4; i++) {
        VertexSPtr vertex = vertices_offset[i];
        for (auto it_f = vertex->facets().begin(); it_f != vertex->facets().end(); ) { // no C++11
            FacetWPtr facet_wptr = *it_f++;
            if (FacetSPtr facet = facet_wptr.lock()) {
                facet->removeVertex(vertex);
            }
        }
        polyhedron->removeVertex(vertex);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_.clear();
    post_op_edges_ = {{ edge_offset_11, edge_offset_21 }};
    post_op_facets_ = {{ edge_offset_11->getFacetL(),
                         edge_offset_11->getFacetR(),
                         edge_offset_21->getFacetL(),
                         edge_offset_21->getFacetR() }};
    CGAL_postcondition(post_op_vertices_.empty() && post_op_edges_.size() == 2 && post_op_facets_.size() == 4);

    // faces are smaller so nothing from unmodified vertices
    post_op_vertices_VV_.clear();

    // no new vertices & only reducing the size of facets so no edge disconnection
    post_op_vertices_pierce_.clear();

    CGAL_assertion(!isReflex(edge_offset_11));
    CGAL_assertion(!isReflex(edge_offset_21));
    post_op_edges_edgesplit_.clear();

    return EventStatus::EVENT_HANDLED;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleDblTriangleEvent(DblTriangleEventSPtr event,
                                           const CGAL::FT& current_offset,
                                           const std::optional<CGAL::FT>& offset_future_bound,
                                           PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#####  Handle Dbl Triangle Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    appendEventNode(event->getNode());

    EdgeSPtr edge = event->getEdge();

    VertexSPtr vertices[4];
    event->getVertices(vertices);
    VertexSPtr vertices_offset[4];
    for (unsigned int i = 0; i < 4; i++) {
        SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                vertices[i]->getData());
        vertices_offset[i] = vertex_data->getOffsetVertex();
    }
    EdgeSPtr edges[5];
    event->getEdges(edges);
    EdgeSPtr edges_offset[5];
    for (unsigned int i = 0; i < 5; i++) {
        SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                edges[i]->getData());
        edges_offset[i] = edge_data->getOffsetEdge();
    }

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

    SkelFacetDataSPtr facet_data_l = std::dynamic_pointer_cast<SkelFacetData>(
            facet_l->getData());
    FacetSPtr facet_offset_l = facet_data_l->getOffsetFacet();
    SkelFacetDataSPtr facet_data_r = std::dynamic_pointer_cast<SkelFacetData>(
            facet_r->getData());
    FacetSPtr facet_offset_r = facet_data_r->getOffsetFacet();
    SkelVertexDataSPtr vertex_data_l = std::dynamic_pointer_cast<SkelVertexData>(
            vertex_l->getData());
    VertexSPtr vertex_offset_l = vertex_data_l->getOffsetVertex();
    SkelVertexDataSPtr vertex_data_r = std::dynamic_pointer_cast<SkelVertexData>(
            vertex_r->getData());
    VertexSPtr vertex_offset_r = vertex_data_r->getOffsetVertex();
    SkelEdgeDataSPtr edge_data_l = std::dynamic_pointer_cast<SkelEdgeData>(
            edge_l->getData());
    EdgeSPtr edge_offset_l = edge_data_l->getOffsetEdge();
    SkelEdgeDataSPtr edge_data_r = std::dynamic_pointer_cast<SkelEdgeData>(
            edge_r->getData());
    EdgeSPtr edge_offset_r = edge_data_r->getOffsetEdge();

    if (facet_offset_l->edges().size() == 3) {
        polyhedron->removeFacet(facet_offset_l);
    }
    if (facet_offset_r->edges().size() == 3) {
        polyhedron->removeFacet(facet_offset_r);
    }
    for (unsigned int i = 0; i < 5; i++) {
        EdgeSPtr edge = edges_offset[i];
        if (edge->getFacetL()) {
            edge->getFacetL()->removeEdge(edge);
        }
        if (edge->getFacetR()) {
            edge->getFacetR()->removeEdge(edge);
        }
        polyhedron->removeEdge(edge);
    }

    if (edge_offset_l->getVertexSrc() == vertex_offset_l) {
        if (edge_offset_r->getVertexDst() == vertex_offset_r) {
            edge_offset_l->replaceVertexSrc(edge_offset_r->getVertexSrc());
        } else {
            edge_offset_l->replaceVertexSrc(edge_offset_r->getVertexDst());
        }
    } else {
        if (edge_offset_r->getVertexDst() == vertex_offset_r) {
            edge_offset_l->replaceVertexDst(edge_offset_r->getVertexSrc());
        } else {
            edge_offset_l->replaceVertexDst(edge_offset_r->getVertexDst());
        }
    }
    edge_offset_r->getFacetL()->removeEdge(edge_offset_r);
    edge_offset_r->getFacetR()->removeEdge(edge_offset_r);
    polyhedron->removeEdge(edge_offset_r);

    for (unsigned int i = 0; i < 4; i++) {
        VertexSPtr vertex = vertices_offset[i];
        for (auto it_f = vertex->facets().begin(); it_f != vertex->facets().end(); ) { // no C++11
            FacetWPtr facet_wptr = *it_f++;
            if (FacetSPtr facet = facet_wptr.lock()) {
                facet->removeVertex(vertex);
            }
        }
        polyhedron->removeVertex(vertex);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_.clear();
    post_op_edges_ = {{ edge_offset_l }};
    post_op_facets_ = {{ facet_ll, facet_rr }};
    CGAL_postcondition(post_op_vertices_.empty() && post_op_edges_.size() == 1 && post_op_facets_.size() == 2);

    // faces are smaller so nothing from unmodified vertices
    post_op_vertices_VV_.clear();

    // no new vertices & only reducing the size of facets so no edge disconnection
    post_op_vertices_pierce_.clear();

    // faces are getting smaller so no need to check unmodified edges
    post_op_edges_edgesplit_ = {{ edge_offset_l }};

    return EventStatus::EVENT_HANDLED;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleTetrahedronEvent(TetrahedronEventSPtr event,
                                           const CGAL::FT& current_offset,
                                           const std::optional<CGAL::FT>& offset_future_bound,
                                           PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Tetrahedron Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    appendEventNode(event->getNode());

    VertexSPtr vertices[4];
    event->getVertices(vertices);
    VertexSPtr vertices_offset[4];
    for (unsigned int i = 0; i < 4; i++) {
        SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                vertices[i]->getData());
        vertices_offset[i] = vertex_data->getOffsetVertex();
    }
    EdgeSPtr edges[6];
    event->getEdges(edges);
    EdgeSPtr edges_offset[6];
    for (unsigned int i = 0; i < 6; i++) {
        SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                edges[i]->getData());
        edges_offset[i] = edge_data->getOffsetEdge();
    }
    FacetSPtr facets[4];
    event->getFacets(facets);
    FacetSPtr facets_offset[4];
    for (unsigned int i = 0; i < 4; i++) {
        SkelFacetDataSPtr facet_data = std::dynamic_pointer_cast<SkelFacetData>(
                facets[i]->getData());
        facets_offset[i] = facet_data->getOffsetFacet();
    }

    for (unsigned int i = 0; i < 4; i++) {
        if (facets_offset[i]->vertices().size() == 3) {
            polyhedron->removeFacet(facets_offset[i]);
        }
    }
    for (unsigned int i = 0; i < 6; i++) {
        EdgeSPtr edge = edges_offset[i];
        if (edge->getFacetL()) {
            edge->getFacetL()->removeEdge(edge);
        }
        if (edge->getFacetR()) {
            edge->getFacetR()->removeEdge(edge);
        }
        polyhedron->removeEdge(edge);
    }
    for (unsigned int i = 0; i < 4; i++) {
        VertexSPtr vertex = vertices_offset[i];
        for (auto it_f = vertex->facets().begin(); it_f != vertex->facets().end(); ) { // no C++11
            FacetWPtr facet_wptr = *it_f++;
            if (FacetSPtr facet = facet_wptr.lock()) {
                facet->removeVertex(vertex);
            }
        }
        polyhedron->removeVertex(vertex);
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_.clear();
    post_op_edges_.clear();
    post_op_facets_.clear();
    CGAL_postcondition(post_op_vertices_.empty() && post_op_edges_.empty() && post_op_facets_.empty());

    // no new vertices & only reducing the size of facets so no edge disconnection
    post_op_vertices_VV_.clear();
    post_op_vertices_pierce_.clear();
    post_op_edges_edgesplit_.clear();

    return EventStatus::EVENT_HANDLED;
}

bool
SimpleStraightSkel::isActualVertexEvent(VertexEventSPtr event,
                                        PolyhedronSPtr polyhedron)
{
    VertexSPtr vertex_1 = event->getVertex1();
    VertexSPtr vertex_2 = event->getVertex2();
    FacetSPtr facet_1 = event->getFacet1();
    FacetSPtr facet_2 = event->getFacet2();

    // @todo avoid all this duplication...
    EdgeSPtr edge_11 = EdgeSPtr();
    EdgeSPtr edge_12 = EdgeSPtr();
    for (EdgeWPtr edge_1_wptr : vertex_1->edges()) {
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
        CGAL_SS3_CORE_TRACE_V(8, "VE Convex split event detected");
        return false;
    }

    CGAL_SS3_CORE_TRACE_V(8, "Vertex event accepted");
    return true;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleVertexEvent(VertexEventSPtr event,
                                      const CGAL::FT& current_offset,
                                      const std::optional<CGAL::FT>& offset_future_bound,
                                      PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "#######  Tentative Vertex Event  #######");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

#ifdef CGAL_SS3_CHECK_CONV_SPLIT_EVENT_AT_POP_TIME
    CGAL_SS3_CORE_TRACE_V(8, "Checking vertex event at pop time...");

    if (!isActualVertexEvent(event, polyhedron)) {
        return EventStatus::NON_EVENT;
    }
#endif

    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "########  Handle Vertex Event  #########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    appendEventNode(event->getNode());

    SkelVertexDataSPtr vertex_data_1 = std::dynamic_pointer_cast<SkelVertexData>(
            event->getVertex1()->getData());
    SkelVertexDataSPtr vertex_data_2 = std::dynamic_pointer_cast<SkelVertexData>(
            event->getVertex2()->getData());
    VertexSPtr vertex_1 = vertex_data_1->getOffsetVertex();
    VertexSPtr vertex_2 = vertex_data_2->getOffsetVertex();
    vertex_1->setPoint(event->getNode()->getPoint());
    vertex_2->setPoint(event->getNode()->getPoint());
    SkelFacetDataSPtr facet_data_1 = std::dynamic_pointer_cast<SkelFacetData>(
            event->getFacet1()->getData());
    SkelFacetDataSPtr facet_data_2 = std::dynamic_pointer_cast<SkelFacetData>(
            event->getFacet2()->getData());
    FacetSPtr facet_1 = facet_data_1->getOffsetFacet();
    FacetSPtr facet_2 = facet_data_2->getOffsetFacet();

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
                edge_22 = edge; // @fixme this is unused...
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

    if (offset_future_bound) {
        vertex_1->final_point_ = nullptr;
        vertex_2->final_point_ = nullptr;
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(edge_tomerge_2->getData());
    edge_data->setSheet(SheetSPtr());

    vertex_data_1 = std::dynamic_pointer_cast<SkelVertexData>(vertex_1->getData());
    vertex_data_1->setNode(event->getNode());
    vertex_data_2 = std::dynamic_pointer_cast<SkelVertexData>(vertex_2->getData());
    vertex_data_2->setNode(event->getNode());

    ArcSPtr arc_1 = createArc(vertex_1);
    skel_result_->addArc(arc_1);
    ArcSPtr arc_2 = createArc(vertex_2);
    skel_result_->addArc(arc_2);
    SheetSPtr sheet = createSheet(edge_tomerge_2);
    skel_result_->addSheet(sheet);

    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

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
    for (EdgeSPtr poe : post_op_edges_) {
        post_op_vertices_pierce_.insert(poe->getVertexSrc());
        post_op_vertices_pierce_.insert(poe->getVertexDst());
    }

    return EventStatus::EVENT_HANDLED;
}

bool
SimpleStraightSkel::isActualFlipVertexEvent(FlipVertexEventSPtr event,
                                            PolyhedronSPtr polyhedron)
{
  VertexSPtr vertex_1 = event->getVertex1();
  VertexSPtr vertex_2 = event->getVertex2();
  FacetSPtr facet_1 = event->getFacet1();
  FacetSPtr facet_2 = event->getFacet2();

  EdgeSPtr edge_11 = EdgeSPtr();
  EdgeSPtr edge_12 = EdgeSPtr();
  for (EdgeWPtr edge_1_wptr : vertex_1->edges()) {
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
      CGAL_SS3_CORE_TRACE_V(8, "FV Convex split event detected");
      return false;
  }

  CGAL_SS3_CORE_TRACE_V(8, "Flip vertex event accepted");
  return true;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleFlipVertexEvent(FlipVertexEventSPtr event,
                                          const CGAL::FT& current_offset,
                                          const std::optional<CGAL::FT>& offset_future_bound,
                                          PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "#####  Tentative Flip Vertex Event  ####");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

#ifdef CGAL_SS3_CHECK_CONV_SPLIT_EVENT_AT_POP_TIME
    CGAL_SS3_CORE_TRACE_V(8, "Checking flip vertex event at pop time...");

    if (!isActualFlipVertexEvent(event, polyhedron)) {
        return EventStatus::NON_EVENT;
    }
#endif

    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Flip Vertex Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    appendEventNode(event->getNode());

    SkelVertexDataSPtr vertex_data_1 = std::dynamic_pointer_cast<SkelVertexData>(
            event->getVertex1()->getData());
    SkelVertexDataSPtr vertex_data_2 = std::dynamic_pointer_cast<SkelVertexData>(
            event->getVertex2()->getData());
    VertexSPtr vertex_1 = vertex_data_1->getOffsetVertex();
    VertexSPtr vertex_2 = vertex_data_2->getOffsetVertex();
    SkelFacetDataSPtr facet_data_1 = std::dynamic_pointer_cast<SkelFacetData>(
            event->getFacet1()->getData());
    SkelFacetDataSPtr facet_data_2 = std::dynamic_pointer_cast<SkelFacetData>(
            event->getFacet2()->getData());
    FacetSPtr facet_1 = facet_data_1->getOffsetFacet();
    FacetSPtr facet_2 = facet_data_2->getOffsetFacet();

    vertex_1->setPoint(event->getNode()->getPoint());
    vertex_2->setPoint(event->getNode()->getPoint());

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

    if (offset_future_bound) {
        vertex_1->final_point_ = nullptr;
        vertex_2->final_point_ = nullptr;
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    vertex_data_1 = std::dynamic_pointer_cast<SkelVertexData>(vertex_1->getData());
    vertex_data_1->setNode(event->getNode());
    vertex_data_2 = std::dynamic_pointer_cast<SkelVertexData>(vertex_2->getData());
    vertex_data_2->setNode(event->getNode());
    ArcSPtr arc_1 = createArc(vertex_1);
    skel_result_->addArc(arc_1);
    ArcSPtr arc_2 = createArc(vertex_2);
    skel_result_->addArc(arc_2);

    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_ = {{ vertex_1, vertex_2 }};
    for (EdgeWPtr we : vertex_1->edges()) { post_op_edges_.insert(EdgeSPtr(we.lock())); }
    for (EdgeWPtr we : vertex_2->edges()) { post_op_edges_.insert(EdgeSPtr(we.lock())); }
    for (FacetWPtr wf : vertex_1->facets()) { post_op_facets_.insert(wf.lock()); }
    for (FacetWPtr wf : vertex_2->facets()) { post_op_facets_.insert(wf.lock()); }
    CGAL_postcondition(post_op_vertices_.size() == 2 && post_op_edges_.size() == 6 && post_op_facets_.size() == 4);

    // facets did not grow
    post_op_vertices_VV_ = {{ vertex_1, vertex_2 }};

    // probably could improve this but flip vertex events are rare_, so we are rarely
    // looking for new pierce events after a vertex event so it doesn't matter much
    for (EdgeSPtr poe : post_op_edges_) {
      post_op_vertices_pierce_.insert(poe->getVertexSrc());
      post_op_vertices_pierce_.insert(poe->getVertexDst());
    }

    return EventStatus::EVENT_HANDLED;
}

bool
SimpleStraightSkel::isActualSurfaceEvent(SurfaceEventSPtr event,
                                         PolyhedronSPtr polyhedron)
{
    EdgeSPtr edge_1 = event->getEdge1();
    EdgeSPtr edge_2 = event->getEdge2();
    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();

    bool is_conv_split_event = false;
    std::list<EdgeSPtr> common_edges = edge_1->getFacetL()->findEdges(edge_1->getFacetR());
    for (EdgeSPtr edge : common_edges) {
        if (edge == edge_1) {
            continue;
        }
        FacetSPtr facet_src = edge->getFacetSrc();
        FacetSPtr facet_dst = edge->getFacetDst();
        if (facet_1_src == edge_2->getFacetL() ||
                facet_1_dst == edge_2->getFacetL()) {
            if (facet_src == edge_2->getFacetR() ||
                    facet_dst == edge_2->getFacetR()) {
                is_conv_split_event = true;
                break;
            }
        } else if (facet_1_src == edge_2->getFacetR() ||
                facet_1_dst == edge_2->getFacetR()) {
            if (facet_src == edge_2->getFacetL() ||
                    facet_dst == edge_2->getFacetL()) {
                is_conv_split_event = true;
                break;
            }
        }
    }

    if (is_conv_split_event) {
        CGAL_SS3_CORE_TRACE_V(8, "SE Convex split event detected");
        return false;
    }

    CGAL_SS3_CORE_TRACE_V(8, "Surface event accepted");
    return true;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleSurfaceEvent(SurfaceEventSPtr event,
                                       const CGAL::FT& current_offset,
                                       const std::optional<CGAL::FT>& offset_future_bound,
                                       PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "######  Tentative Surface Event  #######");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

#ifdef CGAL_SS3_CHECK_CONV_SPLIT_EVENT_AT_POP_TIME
    CGAL_SS3_CORE_TRACE_V(8, "Checking surface event at pop time...");

    if (!isActualSurfaceEvent(event, polyhedron)) {
        return EventStatus::NON_EVENT;
    }
#endif

    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "#######  Handle Surface Event  #########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    CGAL_SS3_CORE_TRACE_V(4, "Edge A = " << event->getEdge1()->toString());
    CGAL_SS3_CORE_TRACE_V(4, "Edge B = " << event->getEdge2()->toString());

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
      event->getEdge1()->getData());
    EdgeSPtr edge_1 = data_1->getOffsetEdge();
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
      event->getEdge2()->getData());
    EdgeSPtr edge_2 = data_2->getOffsetEdge();

    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();

    VertexSPtr vertex;
    EdgeSPtr edge_b1;
    EdgeSPtr edge_b2;
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

    vertex->setPoint(node->getPoint());
    VertexSPtr vertex_21 = Vertex::create(node->getPoint());
    VertexSPtr vertex_22 = Vertex::create(node->getPoint());
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

    if (offset_future_bound) {
        vertex->final_point_ = nullptr;
        vertex_21->final_point_ = nullptr;
        vertex_22->final_point_ = nullptr;
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(edge_2->getData());
    SheetSPtr sheet = edge_data->getSheet();
    edge_data = SkelEdgeData::create(edge_tmp);
    edge_data->setSheet(sheet);

    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
    vertex_data->setNode(node);

    ArcSPtr arc = createArc(vertex);
    skel_result_->addArc(arc);
    vertex_data = SkelVertexData::create(vertex_21);
    vertex_data->setNode(node);
    arc = createArc(vertex_21);
    skel_result_->addArc(arc);
    vertex_data = SkelVertexData::create(vertex_22);
    vertex_data->setNode(node);
    arc = createArc(vertex_22);
    skel_result_->addArc(arc);

    SkelEdgeData::create(edge_21);
    sheet = createSheet(edge_21);
    skel_result_->addSheet(sheet);
    SkelEdgeData::create(edge_22);
    sheet = createSheet(edge_22);
    skel_result_->addSheet(sheet);

    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_ = {{ vertex, vertex_21, vertex_22 }};
    for (VertexSPtr v : post_op_vertices_) {
        for (EdgeWPtr we : v->edges()) {
            post_op_edges_.insert(EdgeSPtr(we.lock()));
        }
    }

    post_op_facets_ = {{ edge_1->getFacetL(), edge_1->getFacetR(),
                         edge_2->getFacetL(), edge_2->getFacetR() }};

    // not sure if we need that much stuff
    post_op_vertices_VV_ = {{ vertex, vertex_21, vertex_22 }};
    for (VertexSPtr v : { vertex, vertex_21, vertex_22 }) {
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

    return EventStatus::EVENT_HANDLED;
}

bool
SimpleStraightSkel::isActualPolyhedronSplitEvent(PolyhedronSplitEventSPtr event,
                                                 const CGAL::FT& current_offset,
                                                 PolyhedronSPtr polyhedron)
{
    // @speed is_degenerate(4 planes)? But it would be a nice filter failure, so costly
    const CGAL::FT& offset_event = event->getOffset();
    CGAL::FT shift = offset_event - current_offset;
    Segment3SPtr e1o = PolyhedronTransformation::shiftEdge(event->getEdge1(), shift);
    if (!e1o->is_degenerate()) {
        return false;
    }

    return true;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handlePolyhedronSplitEvent(PolyhedronSplitEventSPtr event,
                                               const CGAL::FT& current_offset,
                                               const std::optional<CGAL::FT>& offset_future_bound,
                                               PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "####  Tentative Split Merge Event  #####");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

#ifdef CGAL_SS3_CHECK_POLYHEDRON_SPLIT_EVENT_AT_POP_TIME
    CGAL_SS3_CORE_TRACE_V(8, "Checking polyhedron split at pop time...");

    if (!isActualPolyhedronSplitEvent(event, current_offset, polyhedron)) {
        return EventStatus::NON_EVENT;
    }
#endif

    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "####  Handle Polyhedron Split Event  ###");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    CGAL_SS3_CORE_TRACE_V(4, "Node = " << *(node->getPoint()));
    CGAL_SS3_CORE_TRACE_V(4, "Edge A = " << event->getEdge1()->toString());
    CGAL_SS3_CORE_TRACE_V(4, "Edge B = " << event->getEdge2()->toString());

    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge1()->getData());
    EdgeSPtr edge_1 = data_1->getOffsetEdge();
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge2()->getData());
    EdgeSPtr edge_2 = data_2->getOffsetEdge();

    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();

    VertexSPtr vertex_l;
    VertexSPtr vertex_r;
    EdgeSPtr edge_22;
    if (edge_2->getFacetL() == facet_1_src &&
            edge_2->getFacetR() == facet_1_dst) {
        vertex_l = edge_1->getVertexSrc();
        vertex_r = edge_1->getVertexDst();
        vertex_l->setPoint(node->getPoint());
        vertex_r->setPoint(node->getPoint());

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
    if (edge_2->getFacetL() == facet_1_dst &&
            edge_2->getFacetR() == facet_1_src) {
        vertex_l = edge_1->getVertexDst();
        vertex_r = edge_1->getVertexSrc();
        vertex_l->setPoint(node->getPoint());
        vertex_r->setPoint(node->getPoint());

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

    if (offset_future_bound) {
        vertex_l->final_point_ = nullptr;
        vertex_r->final_point_ = nullptr;
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelVertexDataSPtr data_l = std::dynamic_pointer_cast<SkelVertexData>(vertex_l->getData());
    data_l->setNode(node);
    SkelVertexDataSPtr data_r = std::dynamic_pointer_cast<SkelVertexData>(vertex_r->getData());
    data_r->setNode(node);
    SkelEdgeDataSPtr data_22 = std::dynamic_pointer_cast<SkelEdgeData>(edge_22->getData());
    data_22->setSheet(data_2->getSheet());
    ArcSPtr arc_l = createArc(vertex_l);
    ArcSPtr arc_r = createArc(vertex_r);
    skel_result_->addArc(arc_l);
    skel_result_->addArc(arc_r);

    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_ = {{ vertex_l, vertex_r }};
    for (VertexSPtr v : post_op_vertices_) {
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

    return EventStatus::EVENT_HANDLED;
}

bool
SimpleStraightSkel::isActualSplitMergeEvent(SplitMergeEventSPtr event,
                                            PolyhedronSPtr polyhedron)
{
    VertexSPtr vertex_1 = event->getVertex1();
    VertexSPtr vertex_2 = event->getVertex2();
    FacetSPtr facet_1 = event->getFacet1();
    FacetSPtr facet_2 = event->getFacet2();

    EdgeSPtr edge_11 = EdgeSPtr();
    EdgeSPtr edge_12 = EdgeSPtr();
    for (EdgeWPtr edge_1_wptr : vertex_1->edges()) {
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
    if (!conv_split_event) {
        CGAL_SS3_CORE_TRACE_V(8, "Convex split event detected");
        return false;
    }

    CGAL_SS3_CORE_TRACE_V(8, "Split merge event accepted");
    return true;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleSplitMergeEvent(SplitMergeEventSPtr event,
                                          const CGAL::FT& current_offset,
                                          const std::optional<CGAL::FT>& offset_future_bound,
                                          PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "####  Tentative Split Merge Event  #####");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

#ifdef CGAL_SS3_CHECK_CONV_SPLIT_EVENT_AT_POP_TIME
    CGAL_SS3_CORE_TRACE_V(8, "Checking split merge at pop time...");

    if (!isActualSplitMergeEvent(event, polyhedron)) {
        return EventStatus::NON_EVENT;
    }
#endif

    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Split Merge Event  ######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    appendEventNode(event->getNode());

    SkelVertexDataSPtr vertex_data_1 = std::dynamic_pointer_cast<SkelVertexData>(
            event->getVertex1()->getData());
    SkelVertexDataSPtr vertex_data_2 = std::dynamic_pointer_cast<SkelVertexData>(
            event->getVertex2()->getData());
    VertexSPtr vertex_1 = vertex_data_1->getOffsetVertex();
    VertexSPtr vertex_2 = vertex_data_2->getOffsetVertex();
    vertex_1->setPoint(event->getNode()->getPoint());
    vertex_2->setPoint(event->getNode()->getPoint());
    SkelFacetDataSPtr facet_data_1 = std::dynamic_pointer_cast<SkelFacetData>(
            event->getFacet1()->getData());
    SkelFacetDataSPtr facet_data_2 = std::dynamic_pointer_cast<SkelFacetData>(
            event->getFacet2()->getData());
    FacetSPtr facet_1 = facet_data_1->getOffsetFacet();
    FacetSPtr facet_2 = facet_data_2->getOffsetFacet();

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
    EdgeSPtr edge_tosplit;
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

    if (offset_future_bound) {
        vertex_1->final_point_ = nullptr;
        vertex_2->final_point_ = nullptr;
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(edge_tosplit->getData());
    SheetSPtr sheet = edge_data->getSheet();
    edge_data = std::dynamic_pointer_cast<SkelEdgeData>(edge_tomerge_2->getData());
    edge_data->setSheet(sheet);

    vertex_data_1 = std::dynamic_pointer_cast<SkelVertexData>(vertex_1->getData());
    vertex_data_1->setNode(event->getNode());
    vertex_data_2 = std::dynamic_pointer_cast<SkelVertexData>(vertex_2->getData());
    vertex_data_2->setNode(event->getNode());
    ArcSPtr arc_1 = createArc(vertex_1);
    skel_result_->addArc(arc_1);
    ArcSPtr arc_2 = createArc(vertex_2);
    skel_result_->addArc(arc_2);

    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_ = {{ vertex_1, vertex_2 }};
    for (VertexSPtr v : post_op_vertices_) {
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
    CGAL_assertion(!isReflex(vertex_1));
    CGAL_assertion(!isReflex(vertex_2));
    post_op_vertices_pierce_.clear();

    return EventStatus::EVENT_HANDLED;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleEdgeSplitEvent(EdgeSplitEventSPtr event,
                                         const CGAL::FT& current_offset,
                                         const std::optional<CGAL::FT>& offset_future_bound,
                                         PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "######  Handle Edge Split Event  #######");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    CGAL_SS3_CORE_TRACE("edge_1 = " << event->getEdge1()->toString());
    CGAL_SS3_CORE_TRACE("edge_2 = " << event->getEdge2()->toString());

    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge1()->getData());
    EdgeSPtr edge_1 = data_1->getOffsetEdge();
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge2()->getData());
    EdgeSPtr edge_2 = data_2->getOffsetEdge();

    FacetSPtr facet_l1 = edge_1->getFacetL();
    FacetSPtr facet_r1 = edge_1->getFacetR();
    FacetSPtr facet_l2 = edge_2->getFacetL();
    FacetSPtr facet_r2 = edge_2->getFacetR();

    VertexSPtr vertices[4];
    for (unsigned int i = 0; i < 4; i++) {
        vertices[i] = Vertex::create(node->getPoint());
        polyhedron->addVertex(vertices[i]);
    }
    EdgeSPtr edges[4];
    for (unsigned int i = 0; i < 4; i++) {
        edges[i] = Edge::create(vertices[i], vertices[(i+1)%4]);
        polyhedron->addEdge(edges[i]);
    }
    if (event->getEdgeOrientation() > 0) {
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
    edge_1->replaceVertexDst(vertices[0]);
    edge_2->replaceVertexDst(vertices[1]);
    for (unsigned int i = 0; i < 4; i++) {
        // adds the vertices also
        edges[i]->getFacetL()->addEdge(edges[i]);
        edges[i]->getFacetR()->addEdge(edges[i]);
    }

    if (offset_future_bound) {
        for (std::size_t i=0; i<4; ++i) {
            vertices[i]->final_point_ = nullptr;
        }
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelEdgeDataSPtr edge_data_1 = std::dynamic_pointer_cast<SkelEdgeData>(edge_1->getData());
    SkelEdgeDataSPtr edge_data_12 = SkelEdgeData::create(edge_12);
    edge_data_12->setSheet(edge_data_1->getSheet());
    SkelEdgeDataSPtr edge_data_2 = std::dynamic_pointer_cast<SkelEdgeData>(edge_2->getData());
    SkelEdgeDataSPtr edge_data_22 = SkelEdgeData::create(edge_22);
    edge_data_22->setSheet(edge_data_2->getSheet());
    for (unsigned int i = 0; i < 4; i++) {
        SkelVertexDataSPtr vertex_data = SkelVertexData::create(vertices[i]);
        vertex_data->setNode(node);
        ArcSPtr arc = createArc(vertices[i]);
        skel_result_->addArc(arc);
    }
    for (unsigned int i = 0; i < 4; i++) {
        SkelEdgeData::create(edges[i]);
        SheetSPtr sheet = createSheet(edges[i]);
        skel_result_->addSheet(sheet);
    }

    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_ = {{ vertices[0], vertices[1], vertices[2], vertices[3] }};
    for (VertexSPtr v : post_op_vertices_) {
        for (EdgeWPtr we : v->edges()) {
            post_op_edges_.insert(EdgeSPtr(we.lock()));
        }
    }
    post_op_facets_ = {{ facet_l1, facet_r1, facet_l2, facet_r2 }};
    CGAL_postcondition(post_op_vertices_.size() == 4 && post_op_edges_.size() == 8 && post_op_facets_.size() == 4);

    // @todo do it a bit more elegantly
    post_op_vertices_VV_ = {{ vertices[0], vertices[1], vertices[2], vertices[3] }};
    for (VertexSPtr v : { vertices[0], vertices[1], vertices[2], vertices[3] }) {
        for (EdgeWPtr we : v->edges()) {
            EdgeSPtr e = we.lock();
            post_op_vertices_VV_.insert(e->other(v));
        }
    }

    // facets grow so we also need to check vertices that are extremities of edges
    // being subdivided
    for (EdgeSPtr poe : post_op_edges_) {
        post_op_vertices_pierce_.insert(poe->getVertexSrc());
        post_op_vertices_pierce_.insert(poe->getVertexDst());
    }
    CGAL_postcondition(post_op_vertices_pierce_.size() == 8);

    return EventStatus::EVENT_HANDLED;
}

bool
SimpleStraightSkel::isActualPierceEvent(PierceEventSPtr event,
                                        const CGAL::FT& current_offset,
                                        PolyhedronSPtr polyhedron)

{
    CGAL_SS3_CORE_TRACE_V(8, "Is actual Pierce Event?");

    VertexSPtr pv = event->getVertex();
    FacetSPtr pf = event->getFacet();

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
        CGAL_SS3_CORE_TRACE_V(8, "Pierce rejected at pop time (a)");
        return false;
    }

    // Filter if the event point is on an edge (and a fortiori on a vertex)
    // as it will be a different kind of event
    Point3SPtr point = event->getNode()->getPoint();
    FacetSPtr facet_clone = pf->clone();

    CGAL::FT shift = event->getOffset() - current_offset;
    const CGAL::FT& speed = std::dynamic_pointer_cast<SkelFacetData>(pf->getData())->getSpeed();
    Plane3SPtr offset_plane = KernelWrapper::offsetPlane(pf->getPlane(), shift*speed);
    facet_clone->setPlane(offset_plane);

    // abusing the fact that vertices will have the same order in both facets
    std::list<VertexSPtr>::iterator it_v = pf->vertices().begin();
    std::list<VertexSPtr>::iterator it_v_offset = facet_clone->vertices().begin();
    while (it_v != pf->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        VertexSPtr offset_vertex = *it_v_offset++;
        Point3SPtr point_offset = PolyhedronTransformation::shiftPoint(vertex, shift);
        offset_vertex->setPoint(point_offset);
    }

#define CGAL_SLS3_NEW_IS_INSIDE
#ifdef CGAL_SLS3_NEW_IS_INSIDE
    if (!SelfIntersection::isInsideWithRayShootingV2(point, facet_clone)) {
        CGAL_SS3_CORE_TRACE_V(8, "Pierce rejected at pop time (b)");
        return false;
    }
#else
    if (!SelfIntersection::isInsideWithRayShootingV2(point, facet_clone)) {
        CGAL_SS3_CORE_TRACE_V(8, "Pierce rejected at pop time (b)");
        return false;
    }

    bool boundary_rejection = false;
    for (EdgeSPtr edge : facet_clone->edges()) {
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
        CGAL_SS3_CORE_TRACE_V(8, "Pierce rejected at pop time (c)");
        return false;
    }
#endif

    CGAL_SS3_CORE_TRACE_V(8, "Pierce event accepted");
    return true;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handlePierceEvent(PierceEventSPtr event,
                                      const CGAL::FT& current_offset,
                                      const std::optional<CGAL::FT>& offset_future_bound,
                                      PolyhedronSPtr polyhedron)
{
    CGAL_SS3_CORE_TRACE_V(8, "########################################");
    CGAL_SS3_CORE_TRACE_V(8, "######  Tentative Pierce Event  ########");
    CGAL_SS3_CORE_TRACE_V(8, "########################################");

#ifdef CGAL_SS3_CHECK_PIERCE_AT_POP_TIME
    CGAL_SS3_CORE_TRACE_V(8, "Checking pierce event at pop time...");

    if (!isActualPierceEvent(event, current_offset, polyhedron)) {
        return EventStatus::NON_EVENT;
    }
#endif

    CGAL_SS3_CORE_TRACE_V(4, "########################################");
    CGAL_SS3_CORE_TRACE_V(4, "########  Handle Pierce Event  #########");
    CGAL_SS3_CORE_TRACE_V(4, "########################################");

    const CGAL::FT& event_offset = event->getOffset();
    polyhedron = shiftToEventOffset(polyhedron, current_offset, event_offset);

    NodeSPtr node = event->getNode();

    CGAL_SS3_CORE_TRACE("Node: " << node->toString());
    CGAL_SS3_CORE_TRACE("V: " << event->getVertex()->toString());
    CGAL_SS3_CORE_TRACE("F: " << event->getFacet()->toString());

    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
            event->getVertex()->getData());
    VertexSPtr vertex_offset = vertex_data->getOffsetVertex();
    SkelFacetDataSPtr facet_data = std::dynamic_pointer_cast<SkelFacetData>(
            event->getFacet()->getData());
    FacetSPtr facet_offset = facet_data->getOffsetFacet();

    // the 3 new vertices cannot be reflex, but since we grow faces,
    // we need to check the other extremities of the edges
    for (EdgeWPtr ew : vertex_offset->edges()) {
        EdgeSPtr edge = ew.lock();
        if (edge->getVertexSrc() != vertex_offset) {
            post_op_vertices_pierce_.insert(edge->getVertexSrc());
        } else {
            post_op_vertices_pierce_.insert(edge->getVertexDst());
        }
    }

    CGAL_postcondition(post_op_vertices_pierce_.size() == 3);

    appendEventNode(node);

    FacetSPtr facets[3];
    EdgeSPtr edges[3];
    EdgeSPtr edge = vertex_offset->firstEdge();
    for (unsigned int i = 0; i < 3; i++) {
        edges[i] = edge;
        if (edge->getVertexSrc() == vertex_offset) {
            facets[i] = edge->getFacetL();
        } else if (edge->getVertexDst() == vertex_offset) {
            facets[i] = edge->getFacetR();
        }
        edge = edge->next(vertex_offset);
    }

    VertexSPtr vertices[3];
    for (unsigned int i = 0; i < 3; i++) {
        vertices[i] = Vertex::create(node->getPoint());
        facet_offset->addVertex(vertices[i]);
        polyhedron->addVertex(vertices[i]);
    }
    for (unsigned int i = 0; i < 3; i++) {
        EdgeSPtr edge = edges[i];
        if (edge->getVertexSrc() == vertex_offset) {
            edge->replaceVertexSrc(vertices[i]);
        } else if (edge->getVertexDst() == vertex_offset) {
            edge->replaceVertexDst(vertices[i]);
        }
        facets[i]->removeVertex(vertex_offset);
        facets[i]->addVertex(vertices[i]);
        facets[(i+2)%3]->addVertex(vertices[i]);
    }
    vertex_offset->facets().clear();
    vertex_offset->edges().clear();
    polyhedron->removeVertex(vertex_offset);
    for (unsigned int i = 0; i < 3; i++) {
        edges[i] = Edge::create(vertices[i], vertices[(i+1)%3]);
        edges[i]->setFacetL(facet_offset);
        edges[i]->setFacetR(facets[i]);
        facet_offset->addEdge(edges[i]);
        facets[i]->addEdge(edges[i]);
        polyhedron->addEdge(edges[i]);
    }

    if (offset_future_bound) {
        for (std::size_t i=0; i<3; ++i) {
            vertices[i]->final_point_ = nullptr;
        }
    }

#ifndef CGAL_SS3_NO_SKELETON_DS
    for (unsigned int i = 0; i < 3; i++) {
        SkelVertexDataSPtr vertex_data = SkelVertexData::create(vertices[i]);
        vertex_data->setNode(event->getNode());
        ArcSPtr arc = createArc(vertices[i]);
        skel_result_->addArc(arc);
    }
    for (unsigned int i = 0; i < 3; i++) {
        SkelEdgeData::create(edges[i]);
        SheetSPtr sheet = createSheet(edges[i]);
        skel_result_->addSheet(sheet);
    }

    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    post_op_vertices_ = {{ vertices[0], vertices[1], vertices[2] }};
    for (VertexSPtr v : post_op_vertices_) {
        for (EdgeWPtr we : v->edges()) {
            post_op_edges_.insert(EdgeSPtr(we.lock()));
        }
    }
    for (unsigned int i = 0; i < 3; i++) {
      post_op_facets_.insert(edges[i]->getFacetL());
      post_op_facets_.insert(edges[i]->getFacetR());
    }
    CGAL_postcondition(post_op_vertices_.size() == 3 && post_op_edges_.size() == 6 && post_op_facets_.size() == 4);

    post_op_vertices_VV_ = {{ vertices[0], vertices[1], vertices[2] }};
    for (VertexSPtr v : { vertices[0], vertices[1], vertices[2] }) {
        for (EdgeWPtr we : v->edges()) {
            EdgeSPtr e = we.lock();
            post_op_vertices_VV_.insert(e->other(v));
        }
    }

    // the edge with 'facet offset' is convex so these vertices are not interesting
    CGAL_assertion(!isReflex(vertices[0]));
    CGAL_assertion(!isReflex(vertices[1]));
    CGAL_assertion(!isReflex(vertices[2]));

    return EventStatus::EVENT_HANDLED;
}

SimpleStraightSkel::EventStatus
SimpleStraightSkel::handleEvent(AbstractEventSPtr event,
                                const CGAL::FT& current_offset,
                                const std::optional<CGAL::FT>& offset_future_bound,
                                PolyhedronSPtr polyhedron)
{
    EventStatus result = EventStatus::NON_EVENT;

    if (event->getType() == AbstractEvent::SAVE_OFFSET_EVENT) {
        result = handleSaveOffsetEvent(std::dynamic_pointer_cast<SaveOffsetEvent>(event),
                                       current_offset, polyhedron);
    } else if (event->getType() == AbstractEvent::CONST_OFFSET_EVENT) {
        result = handleConstOffsetEvent(std::dynamic_pointer_cast<ConstOffsetEvent>(event),
                                        current_offset, polyhedron);
#ifdef CGAL_SS3_USE_GENERIC_VANISH_EVENT
    } else if (event->getType() == AbstractEvent::VANISH_EVENT) {
        result = handleVanishEvent(std::dynamic_pointer_cast<VanishEvent>(event),
                                   current_offset, offset_future_bound, polyhedron);
#else
    } else if (event->getType() == AbstractEvent::EDGE_EVENT) {
        result = handleEdgeEvent(std::dynamic_pointer_cast<EdgeEvent>(event),
                                 current_offset, offset_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::EDGE_MERGE_EVENT) {
        result = handleEdgeMergeEvent(std::dynamic_pointer_cast<EdgeMergeEvent>(event),
                                      current_offset, offset_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
        result = handleTriangleEvent(std::dynamic_pointer_cast<TriangleEvent>(event),
                                     current_offset, offset_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::DBL_EDGE_MERGE_EVENT) {
        result = handleDblEdgeMergeEvent(std::dynamic_pointer_cast<DblEdgeMergeEvent>(event),
                                         current_offset, offset_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::DBL_TRIANGLE_EVENT) {
        result = handleDblTriangleEvent(std::dynamic_pointer_cast<DblTriangleEvent>(event),
                                        current_offset, offset_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::TETRAHEDRON_EVENT) {
        result = handleTetrahedronEvent(std::dynamic_pointer_cast<TetrahedronEvent>(event),
                                        current_offset, offset_future_bound, polyhedron);
#endif
    } else if (event->getType() == AbstractEvent::VERTEX_EVENT) {
        result = handleVertexEvent(std::dynamic_pointer_cast<VertexEvent>(event),
                                   current_offset, offset_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::FLIP_VERTEX_EVENT) {
        result = handleFlipVertexEvent(std::dynamic_pointer_cast<FlipVertexEvent>(event),
                                       current_offset, offset_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::SURFACE_EVENT) {
        result = handleSurfaceEvent(std::dynamic_pointer_cast<SurfaceEvent>(event),
                                    current_offset, offset_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
        result = handlePolyhedronSplitEvent(std::dynamic_pointer_cast<PolyhedronSplitEvent>(event),
                                            current_offset, offset_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::SPLIT_MERGE_EVENT) {
        result = handleSplitMergeEvent(std::dynamic_pointer_cast<SplitMergeEvent>(event),
                                       current_offset, offset_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::EDGE_SPLIT_EVENT) {
        result = handleEdgeSplitEvent(std::dynamic_pointer_cast<EdgeSplitEvent>(event),
                                      current_offset, offset_future_bound, polyhedron);
    } else if (event->getType() == AbstractEvent::PIERCE_EVENT) {
        result = handlePierceEvent(std::dynamic_pointer_cast<PierceEvent>(event),
                                   current_offset, offset_future_bound, polyhedron);
    } else {
        CGAL_SS3_CORE_TRACE("Error: Cannot handle event of type " << event->getType());
        CGAL_assertion(false);
        result = EventStatus::EVENT_NOT_HANDLED;
    }

    // Only two event types are currently allowed not to handle:
    // - Vanish events: since we do not filter during collect time, they can be requalified
    //                  as contact events at pop time, even while still being valid
    // - Surface events: other events can change whether an event is a surface event or a split event
    //                   so we need to consider the possible event until it's popped
    // - Pierce events: the filtering (i.e., is the contact point actually on the facet)
    //                  is done at pop time, so the event could in fact not exist
#if 0 // @todo VV also
    CGAL_postcondition(event->getType() == AbstractEvent::VANISH_EVENT ||
                       event->getType() == AbstractEvent::SURFACE_EVENT ||
                       event->getType() == AbstractEvent::PIERCE_EVENT ||
                       result == EventStatus::EVENT_HANDLED);
#endif

    CGAL_postcondition_code(for(auto v : polyhedron->vertices()))
    CGAL_postcondition(v->getID() != -1);
    CGAL_postcondition_code(for(auto e : polyhedron->edges()))
    CGAL_postcondition(e->getID() != -1);
    CGAL_postcondition_code(for(auto f : polyhedron->facets()))
    CGAL_postcondition(f->getID() != -1);

    CGAL_SS3_CORE_TRACE_V(4, "-- Finished handling Event --");
    return result;
}

StraightSkeletonSPtr SimpleStraightSkel::getResult() const {
#ifndef CGAL_SS3_NO_SKELETON_DS
    CGAL_SS3_CORE_TRACE("Error: no skeleton to return: it was not built");
    return { };
#endif
    return this->skel_result_;
}

} }
