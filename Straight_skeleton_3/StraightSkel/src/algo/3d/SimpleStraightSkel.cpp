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
#include "data/3d/Polyhedron.h"
#include "data/3d/skel/StraightSkeleton.h"
#include "data/3d/skel/AbstractEvent.h"
#include "data/3d/skel/ConstOffsetEvent.h"
#include "data/3d/skel/SaveOffsetEvent.h"
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

namespace algo { namespace _3d {

SimpleStraightSkel::SimpleStraightSkel(PolyhedronSPtr polyhedron) {
    polyhedron_ = polyhedron;
    skel_result_ = StraightSkeleton::create();
    skel_result_->setPolyhedron(polyhedron);
    initVertexSplitter();
    initEdgeEvent();
}

SimpleStraightSkel::SimpleStraightSkel(PolyhedronSPtr polyhedron, ControllerSPtr controller) {
    polyhedron_ = polyhedron;
    controller_ = controller;
    skel_result_ = StraightSkeleton::create();
    skel_result_->setPolyhedron(polyhedron);
    initVertexSplitter();
    initEdgeEvent();
}

SimpleStraightSkel::SimpleStraightSkel(PolyhedronSPtr polyhedron, ControllerSPtr controller, const std::list<CGAL::FT>& save_offsets) {
    polyhedron_ = polyhedron;
    controller_ = controller;
    save_offsets_ = save_offsets;
    skel_result_ = StraightSkeleton::create();
    skel_result_->setPolyhedron(polyhedron);
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

SimpleStraightSkelSPtr SimpleStraightSkel::create(PolyhedronSPtr polyhedron, ControllerSPtr controller, const std::list<CGAL::FT>& save_offsets) {
    return SimpleStraightSkelSPtr(new SimpleStraightSkel(polyhedron, controller, save_offsets));
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
            DEBUG_VAL("Warning: '" << s_vertex_splitter << "' not found.");
            DEBUG_VAL("Using 'ConvexVertexSplitter'.");
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
            DEBUG_VAL("Warning: option '" << s_edge_event << "' not found.");
            DEBUG_VAL("Using 'convex'.");
            edge_event_ = 0;
            s_edge_event = "convex";
        }
    } else {
        edge_event_ = 0;
        s_edge_event = "convex";
    }
    skel_result_->appendConfig("edge_event="+s_edge_event+"; ");
}


// @fixme return three states to handle flat cases?
bool SimpleStraightSkel::isReflex(EdgeSPtr edge) {
    bool result = false;
    VertexSPtr vertex_src = edge->getVertexSrc();
    VertexSPtr vertex_dst = edge->getVertexDst();
    if (vertex_src->getPoint() == vertex_dst->getPoint()) {
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        FacetSPtr facet_src = getFacetSrc(edge);
        FacetSPtr facet_dst = getFacetDst(edge);
        CGAL::FT speed_l = 1.0;
        if (facet_l->hasData()) {
            speed_l = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_l->getData())->getSpeed();
        }
        CGAL::FT speed_r = 1.0;
        if (facet_r->hasData()) {
            speed_r = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_r->getData())->getSpeed();
        }
        CGAL::FT speed_src = 1.0;
        if (facet_src->hasData()) {
            speed_src = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_src->getData())->getSpeed();
        }
        CGAL::FT speed_dst = 1.0;
        if (facet_dst->hasData()) {
            speed_dst = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_dst->getData())->getSpeed();
        }
        Plane3SPtr offset_plane_l = KernelWrapper::offsetPlane(facet_l->plane(), -speed_l);
        Plane3SPtr offset_plane_r = KernelWrapper::offsetPlane(facet_r->plane(), -speed_r);
        Plane3SPtr offset_plane_src = KernelWrapper::offsetPlane(facet_src->plane(), -speed_src);
        Plane3SPtr offset_plane_dst = KernelWrapper::offsetPlane(facet_dst->plane(), -speed_dst);
        Point3SPtr p_src = KernelWrapper::intersection(offset_plane_src,
                offset_plane_l, offset_plane_r);
        Point3SPtr p_dst = KernelWrapper::intersection(offset_plane_dst,
                offset_plane_l, offset_plane_r);
        Vector3SPtr v_dir = KernelFactory::createVector3((*p_dst) - (*p_src));
        Vector3SPtr normal_l = KernelFactory::createVector3(offset_plane_l);
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
    std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
    while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge = EdgeSPtr(edge_wptr);
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
    std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
    while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge = EdgeSPtr(edge_wptr);
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
    if (vertex_src->getPoint() == vertex_dst->getPoint()) {
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        FacetSPtr facet_src = getFacetSrc(edge);
        FacetSPtr facet_dst = getFacetDst(edge);
        CGAL::FT speed_l = 1.0;
        if (facet_l->hasData()) {
            speed_l = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_l->getData())->getSpeed();
        }
        CGAL::FT speed_r = 1.0;
        if (facet_r->hasData()) {
            speed_r = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_r->getData())->getSpeed();
        }
        CGAL::FT speed_src = 1.0;
        if (facet_src->hasData()) {
            speed_src = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_src->getData())->getSpeed();
        }
        CGAL::FT speed_dst = 1.0;
        if (facet_dst->hasData()) {
            speed_dst = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_dst->getData())->getSpeed();
        }
        Plane3SPtr plane_l = facet_l->plane();
        Plane3SPtr plane_r = facet_r->plane();

        Plane3SPtr offset_plane_l = KernelWrapper::offsetPlane(plane_l, -speed_l);
        Plane3SPtr offset_plane_r = KernelWrapper::offsetPlane(plane_r, -speed_r);
        Plane3SPtr offset_plane_src = KernelWrapper::offsetPlane(facet_src->plane(), -speed_src);
        Plane3SPtr offset_plane_dst = KernelWrapper::offsetPlane(facet_dst->plane(), -speed_dst);

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


FacetSPtr SimpleStraightSkel::getFacetSrc(EdgeSPtr edge) {
    FacetSPtr result = FacetSPtr();
    VertexSPtr vertex_src = edge->getVertexSrc();
    if (vertex_src->degree() == 3) {
        result = edge->getFacetL()->next(vertex_src);
    }
    DEBUG_SPTR(result);
    return result;
}

FacetSPtr SimpleStraightSkel::getFacetDst(EdgeSPtr edge) {
    FacetSPtr result = FacetSPtr();
    VertexSPtr vertex_dst = edge->getVertexDst();
    if (vertex_dst->degree() == 3) {
        result = edge->getFacetR()->next(vertex_dst);
    }
    DEBUG_SPTR(result);
    return result;
}

PolyhedronSPtr SimpleStraightSkel::enablePerturbedMode(PolyhedronSPtr polyhedron,
                                                       const CGAL::FT simultaneousOffset) {
    std::cout << "Enabling perturbed mode..." << std::endl;
    CGAL_precondition(!usingPerturbedMode_);

    db::_3d::OBJFile::save("results/pertubation_pre_enable.obj", polyhedron, false /*do not triangulate*/);

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
#define CGAL_SS3_PERTURB_PLANE_COEFFICIENTS
#ifdef CGAL_SS3_PERTURB_PLANE_COEFFICIENTS
        facet->storePlaneCoefficients();
        facet->perturbPlaneCoefficients();
#else
        SkelFacetDataSPtr data;
        if (facet->hasData()) {
            data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
        } else {
            data = SkelFacetData::create(facet);
            data->setSpeed(1.0);
        }

        // store
        facet->cachedSpeed_ = data->getSpeed();
        std::cout << "caching speed: " << facet->cachedSpeed_ << std::endl;

        // perturb
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dist(0.0, 1e-10);

        data->setSpeed(facet->cachedSpeed_ * (1 + dist(gen)));
        std::cout << "speed is now: " << data->getSpeed() << std::endl;
#endif
    }

    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, 0.0, true /*force recomputation of positions*/);
    db::_3d::OBJFile::save("results/pertubation_post_enable.obj", polyhedron, false /*do not triangulate*/);
    CGAL_assertion(bool(polyhedron));
    CGAL_assertion(polyhedron->isConsistent());

    usingPerturbedMode_ = true;
    simultaneousOffset_ = simultaneousOffset;

    return polyhedron;
}

PolyhedronSPtr SimpleStraightSkel::disablePerturbedMode(PolyhedronSPtr polyhedron,
                                                        CGAL::FT currentOffset,
                                                        CGAL::FT nextEventOffset) {
    std::cout << "Disabling perturbed mode..." << std::endl;
    CGAL_precondition(usingPerturbedMode_);

    db::_3d::OBJFile::save("results/pertubation_pre_disable.obj", polyhedron, false /*do not triangulate*/);

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
#ifdef CGAL_SS3_PERTURB_PLANE_COEFFICIENTS
        facet->restorePlaneCoefficients();
#else
        CGAL_assertion(facet->hasData());
        SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
        data->setSpeed(facet->cachedSpeed_);
        std::cout << "speed is back to: " << data->getSpeed() << std::endl;
#endif
    }

    // shift a little more to get some waylay
    CGAL::FT shift = (currentOffset + nextEventOffset) / 2 - currentOffset;
    CGAL_assertion(!is_zero(shift));
    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, shift);
    db::_3d::OBJFile::save("results/pertubation_post_disable.obj", polyhedron, false /*do not triangulate*/);
    CGAL_assertion(bool(polyhedron));

    usingPerturbedMode_ = false;
    simultaneousOffset_ = 0;

    return polyhedron;
}

bool SimpleStraightSkel::run() {
    if (controller_) {
        controller_->wait();
        controller_->setDispPolyhedron(polyhedron_);
        controller_->setDispSkel3d(skel_result_);
    }

    DEBUG_PRINT("== Straight Skeleton 3D started ==");

    double t_start = util::Timer::now();

    PolyhedronSPtr polyhedron = polyhedron_->clone();

    db::_3d::OBJFile::save("results/init_pre.obj", polyhedron, false /*do not triangulate*/);

    DEBUG_VAL("Using " << vertex_splitter_->toString() << " to initialize polyhedron.");
    if (init(polyhedron)) {

        if (controller_) {
            controller_->wait();
        }

        CGAL::FT offset = 0.0;
        CGAL::FT offset_prev = 0.0;

        db::_3d::OBJFile::save("results/init_post.obj", polyhedron, false /*do not triangulate*/);

        for(;;) {
            PQ queue;
            collectEvents(polyhedron, offset, queue);

            // print the queue
            std::cout << "-------------------" << std::endl;
            std::cout << "--- Event queue ---" << std::endl;
            std::cout << "-------------------" << std::endl;
            PQ duplicate_queue = queue;
            while (!duplicate_queue.empty()) {
                AbstractEventSPtr event = duplicate_queue.top();
                std::cout << event->toString() << std::endl;
                duplicate_queue.pop();
            }

            AbstractEventSPtr event;
            bool simultaneousEvents;
            std::tie(event, simultaneousEvents) = nextEvent(queue);

            if (!event) {
                break; // we are done
            }

            std::cout << " current offset: " << offset << "\n"
                      << " next offset: " << event->getOffset() << " (type " << event->getType() << ")\n"
                      << " simultaneous? " << simultaneousEvents << "\n"
                      << " already simultaneous? " << usingPerturbedMode_
                      << " @ " << simultaneousOffset_ << std::endl;

            CGAL::FT delta = 1e-5; // @fixme don't hardcode this value
            if (usingPerturbedMode_) {
                // if the next event is far and we are in perturbed mode, disable perturbation
                if (simultaneousOffset_ - delta > event->getOffset()) { // offsets are negative
                    polyhedron = disablePerturbedMode(polyhedron, offset, event->getOffset());
                    continue; // recompute the polyhedron and its events
                }
            }

            if (simultaneousEvents) {
                if (event->getOffset() == 0) {
                    std::cerr << "NYI: first event simultaneous at T=0" << std::endl;
                    std::exit(1);
                }

                if (usingPerturbedMode_) {
                    std::cerr << "Error: met a simultaneous event in perturbed mode!" << std::endl;
                    // std::exit(1); // @tmp nuke it, for now; do a different perturbation?
                } else {
                    std::cout << "simultaneous events ahead @ offset: " << event->getOffset() << std::endl;
                    polyhedron = enablePerturbedMode(polyhedron, event->getOffset());
                    continue; // recompute events
                }
            }

            static int event_id = 0;
            std::cout << "\n-----------------------------------------------------" << std::endl;
            std::cout << "-- Event #" << event_id << " " << event->toString() << " --" << std::endl;
            std::cout << "-----------------------------------------------------\n" << std::endl;

            if (controller_) {
                controller_->wait();
            }

            offset_prev = offset;
            offset = event->getOffset();

            Point3SPtr p_box_min = PolyhedronTransformation::boundingBoxMin(polyhedron);
            Point3SPtr p_box_max = PolyhedronTransformation::boundingBoxMax(polyhedron);

            const CGAL::FT shift = offset - offset_prev;
            const bool recompute_positions = (shift != 0);
            polyhedron = PolyhedronTransformation::shiftFacets(polyhedron, shift, recompute_positions);

            // will have degeneracies since we haven't treated the event yet
            db::_3d::OBJFile::save("results/shift_" + std::to_string(event_id) + ".obj", polyhedron, false /*do not triangulate*/);

            if (event->getType() == AbstractEvent::SAVE_OFFSET_EVENT) {
                event->setPolyhedronResult(polyhedron);
                skel_result_->addEvent(event);
                std::stringstream ss_filename;
                ss_filename << "offset_" << offset << ".obj";
                db::_3d::OBJFile::save(ss_filename.str(), polyhedron);
                save_offsets_.pop_front();
                if (save_offsets_.empty()) { // @todo ought to be a config flag
                    break;
                }
            } else if (event->getType() == AbstractEvent::CONST_OFFSET_EVENT) {
                event->setPolyhedronResult(polyhedron);
                skel_result_->addEvent(event);
                bool screenshot_on_const_offset_event =
                        util::Configuration::getInstance()->getBool(
                        "algo_3d_SimpleStraightSkel", "screenshot_on_const_offset_event");
                if (controller_ && screenshot_on_const_offset_event) {
                    controller_->screenshot();
                }
            } else if (event->getType() == AbstractEvent::EDGE_EVENT) {
                handleEdgeEvent(std::dynamic_pointer_cast<EdgeEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::EDGE_MERGE_EVENT) {
                handleEdgeMergeEvent(std::dynamic_pointer_cast<EdgeMergeEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
                handleTriangleEvent(std::dynamic_pointer_cast<TriangleEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::DBL_EDGE_MERGE_EVENT) {
                handleDblEdgeMergeEvent(std::dynamic_pointer_cast<DblEdgeMergeEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::DBL_TRIANGLE_EVENT) {
                handleDblTriangleEvent(std::dynamic_pointer_cast<DblTriangleEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::TETRAHEDRON_EVENT) {
                handleTetrahedronEvent(std::dynamic_pointer_cast<TetrahedronEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::VERTEX_EVENT) {
                handleVertexEvent(std::dynamic_pointer_cast<VertexEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::FLIP_VERTEX_EVENT) {
                handleFlipVertexEvent(std::dynamic_pointer_cast<FlipVertexEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::SURFACE_EVENT) {
                handleSurfaceEvent(std::dynamic_pointer_cast<SurfaceEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
                handlePolyhedronSplitEvent(std::dynamic_pointer_cast<PolyhedronSplitEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::SPLIT_MERGE_EVENT) {
                handleSplitMergeEvent(std::dynamic_pointer_cast<SplitMergeEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::EDGE_SPLIT_EVENT) {
                handleEdgeSplitEvent(std::dynamic_pointer_cast<EdgeSplitEvent>(event), polyhedron);
            } else if (event->getType() == AbstractEvent::PIERCE_EVENT) {
                handlePierceEvent(std::dynamic_pointer_cast<PierceEvent>(event), polyhedron);
            }

            DEBUG_PRINT("-- Finished handling Event --");

            db::_3d::OBJFile::save("results/iter_" + std::to_string(event_id) + ".obj", polyhedron, false /*do triangulate*/);
            db::_3d::OBJFile::save("results/iter_" + std::to_string(event_id) + "_triangulated.obj", polyhedron);

            DEBUG_PRINT("-- Check consistency... --");

            assert(polyhedron->isConsistent());
#ifndef CGAL_SS3_NO_ARC_NO_SHEET
            assert(skel_result_->isConsistent());
#endif
            CGAL_assertion(p_box_min && p_box_max);
            CGAL_assertion(PolyhedronTransformation::isInsideBox(polyhedron, p_box_min, p_box_max));

            // this is tempting, but the mesh is usually not in a nice state here
            // assert(!SelfIntersection::hasSelfIntersectingSurface(polyhedron));

            if (controller_) {
                controller_->wait();
            }

            // Can't perturb to treat simultaneous events AND get a meaningful "SAVE" result
            if (simultaneousEvents && event->getType() == AbstractEvent::SAVE_OFFSET_EVENT) {
                std::cout << "WARNING: simultaneous event @ save time, not proceeding further" << std::endl;
                break;
            }

            ++event_id;
        }


        DEBUG_PRINT("== Straight Skeleton 3D finished ==");
        double time = util::Timer::now() - t_start;
        skel_result_->appendDescription("time=" +
                util::StringFactory::fromDouble(time) + "; ");
        //skel_result_->appendDescription("controller=" +
        //        util::StringFactory::fromBoolean(controller_) + "; ");
        DEBUG_VAR(skel_result_->toString());
    } else {
        DEBUG_PRINT("Error: Failed to initialize");
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
    DEBUG_SPTR(result);
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
            if (!facet_wptr.expired()) {
                facets[i] = FacetSPtr(facet_wptr);
                i++;
            }
        }

        // std::cout << "Facet #" << facets[0]->getID() << std::endl;
        // std::cout << "Facet #" << facets[1]->getID() << std::endl;
        // std::cout << "Facet #" << facets[2]->getID() << std::endl;

        if (i >= 3) {
            Vector3SPtr direction;
            Plane3SPtr plane_1 = facets[0]->plane();
            Plane3SPtr plane_2 = facets[1]->plane();
            Plane3SPtr plane_3 = facets[2]->plane();
            // std::cout << "Plane #" << facets[0]->getID() << " --> " << *plane_1 << std::endl;
            // std::cout << "\t" << plane_1->point() << std::endl;
            // std::cout << "\t" << plane_1->point() + plane_1->base1() << std::endl;
            // std::cout << "\t" << plane_1->point() + plane_1->base2() << std::endl;
            // std::cout << "Plane #" << facets[1]->getID() << " --> " << *plane_2 << std::endl;
            // std::cout << "\t" << plane_2->point() << std::endl;
            // std::cout << "\t" << plane_2->point() + plane_2->base1() << std::endl;
            // std::cout << "\t" << plane_2->point() + plane_2->base2() << std::endl;
            // std::cout << "Plane #" << facets[2]->getID() << " --> " << *plane_3 << std::endl;
            // std::cout << "\t" << plane_3->point() << std::endl;
            // std::cout << "\t" << plane_3->point() + plane_3->base1() << std::endl;
            // std::cout << "\t" << plane_3->point() + plane_3->base2() << std::endl;

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
            std::cout << "Offset Plane #:" << facets[0]->getID() << " --> " << *off_1 << std::endl;
            std::cout << "\t" << off_1->point() << std::endl;
            std::cout << "\t" << off_1->point() + off_1->base1() << std::endl;
            std::cout << "\t" << off_1->point() + off_1->base2() << std::endl;
            std::cout << "Offset Plane #:" << facets[1]->getID() << " --> " << *off_2 << std::endl;
            std::cout << "\t" << off_2->point() << std::endl;
            std::cout << "\t" << off_2->point() + off_2->base1() << std::endl;
            std::cout << "\t" << off_2->point() + off_2->base2() << std::endl;
            std::cout << "Offset Plane #:" << facets[2]->getID() << " --> " << *off_3 << std::endl;
            std::cout << "\t" << off_3->point() << std::endl;
            std::cout << "\t" << off_3->point() + off_3->base1() << std::endl;
            std::cout << "\t" << off_3->point() + off_3->base2() << std::endl;

            Point3SPtr dst = KernelWrapper::intersection(off_1, off_2, off_3);
            if (src && dst) {
                direction = KernelFactory::createVector3(*dst - *src);
            }
#else
            Vector3SPtr n_1 = KernelFactory::createVector3(plane_1);
            Vector3SPtr n_2 = KernelFactory::createVector3(plane_2);
            Vector3SPtr n_3 = KernelFactory::createVector3(plane_3);

            // "/ speed"? Is it even correct anyway?
            CGAL_warning_msg(false, "check this");
            direction = KernelFactory::createVector3(speed_1 * (*n_1) + speed_2 * (*n_2) + speed_3 * (*n_3));
#endif

            if (direction) {
                result = Arc::create(data->getNode(), direction);
                data->setArc(result);

                std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
                while (it_e != vertex->edges().end()) {
                    EdgeWPtr edge_wptr = *it_e++;
                    if (!edge_wptr.expired()) {
                        EdgeSPtr edge(edge_wptr);
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
        DEBUG_VAR(vertex->degree());
    }
    DEBUG_SPTR(result);
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
        Plane3SPtr plane_l = facet_l->plane();
        Plane3SPtr plane_r = facet_r->plane();
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
        SkelVertexDataSPtr data_src = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexSrc()->getData());
        SkelVertexDataSPtr data_dst = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexDst()->getData());
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

bool SimpleStraightSkel::init(PolyhedronSPtr polyhedron) {
    std::cout << "Simple Skeleton -- Initialization" << std::endl;

    WriteLock l(polyhedron->mutex());
    bool result = true;

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (vertex->degree() < 3) {
            continue;
        }
        if (!vertex->hasData()) {
            SkelVertexData::create(vertex);
        }
        NodeSPtr node = createNode(vertex);
        if (node) {
            skel_result_->addNode(node);
        } else {
            result = false;
        }
    }

    std::list<VertexSPtr> vertices_tosplit;
    it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
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

    std::cout << vertices_tosplit.size() << " vertices to split" << std::endl;

    if (controller_) {
        l.unlock();
        controller_->setDispPolyhedron(polyhedron);
        controller_->wait();
        l.lock();
    }
    it_v = vertices_tosplit.begin();
    while (it_v != vertices_tosplit.end()) {
        VertexSPtr vertex = *it_v++;
        vertex->getData()->setHighlight(false);
    }

    unsigned int v2s_i = 0;
    it_v = vertices_tosplit.begin();
    while (it_v != vertices_tosplit.end()) {
        VertexSPtr vertex = *it_v++;
        std::cout << "Split... " << v2s_i++ << " pos: " << *(vertex->getPoint()) << std::endl;

        bool equal_speeds = true;
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (facet_wptr.expired()) {
                continue;
            }
            FacetSPtr facet(facet_wptr);
            CGAL::FT speed = 1.0;
            if (facet->hasData()) {
                speed = std::dynamic_pointer_cast<SkelFacetData>(
                        facet->getData())->getSpeed();
            }
            if (speed != 1.0) {
                equal_speeds = false;
                break;
            }
        }

        if (use_fast_vertex_splitter_ && equal_speeds && vertex->isConvex()) {
            std::cout << "Split convex vertex" << std::endl;
            AbstractVertexSplitter::splitConvexVertex(vertex);
        } else if (use_fast_vertex_splitter_ && equal_speeds && vertex->isReflex()) {
            std::cout << "Split reflex vertex" << std::endl;
            AbstractVertexSplitter::splitReflexVertex(vertex);
        } else {
            l.unlock();
            DEBUG_VAR(vertex->toString());
            std::cout << "Split generic vertex" << std::endl;
            vertex_splitter_->splitVertex(vertex);
            if (controller_) {
                controller_->setDispPolyhedron(polyhedron);
                controller_->setDispSkel3d(skel_result_);
            }
            l.lock();
        }
    }

    it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
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
    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
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
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        if (!facet->hasData()) {
            SkelFacetData::create(facet);
        }
    }
    polyhedron->initializeAllIDs();
    CGAL_postcondition(skel_result_->isConsistent());
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
//            Vector3SPtr v2 = KernelFactory::createVector3(facet->plane());
//            if (((*v1) * (*v2)) > 0.0) {
//                // angle between orthogonal vectors of planes < M_PI/2.0
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

/*
  vanishAt:

      \                      /
       \       F0           /
        \                  /
   F1    ------------------   F3
        /                  \
       /       F2           \
      /                      \

  crashAt:

      \                      /                  \                      /
       \       F0           /                    \       F2           /
        \                  /                      \                  /
         ------------------                          ----------------
        /                  \                      /                  \
       /       F1           \                    /       F3           \
      /                      \                  /                      \
*/
enum Quadplane_parallelism {
    PLANES_PARALLELISM_NONE = 0, // No planes are parallel
    PLANES_PARALLEL_01,          // Planes 0 and 1 are parallel (1)
    PLANES_PARALLEL_02,          // Planes 0 and 2 are parallel (2)
    PLANES_PARALLEL_03,          // Planes 0 and 3 are parallel (3)
    PLANES_PARALLEL_12,          // Planes 1 and 2 are parallel (4)
    PLANES_PARALLEL_13,          // Planes 1 and 3 are parallel (5)
    PLANES_PARALLEL_23,          // Planes 2 and 3 are parallel (6)
    PLANES_PARALLEL_012,         // Planes 0, 1, and 2 are parallel (7)
    PLANES_PARALLEL_013,         // Planes 0, 1, and 3 are parallel (8)
    PLANES_PARALLEL_023,         // Planes 0, 2, and 3 are parallel (9)
    PLANES_PARALLEL_123,         // Planes 1, 2, and 3 are parallel (10)
    PLANES_PARALLEL_01_23,       // Planes 0 and 1 are parallel, and planes 2 and 3 are parallel (11)
    PLANES_PARALLEL_02_13,       // Planes 0 and 2 are parallel, and planes 1 and 3 are parallel (12)
    PLANES_PARALLEL_03_12,       // Planes 0 and 3 are parallel, and planes 1 and 2 are parallel (13)
    PLANES_PARALLEL_0123         // All planes are parallel (14)
};

Quadplane_parallelism quadplane_parallelism(FacetSPtr facet_0,
                                            FacetSPtr facet_1,
                                            FacetSPtr facet_2,
                                            FacetSPtr facet_3)
{
    Plane3SPtr plane_0 = facet_0->getPlane();
    Plane3SPtr plane_1 = facet_1->getPlane();
    Plane3SPtr plane_2 = facet_2->getPlane();
    Plane3SPtr plane_3 = facet_3->getPlane();

    // Function to check if two planes are parallel
    auto are_planes_parallel = [](const Plane3SPtr pl1, const Plane3SPtr pl2) {
        // plane coefficients are normalized
        return ((pl1->a() == pl2->a() && pl1->b() == pl2->b() && pl1->c() == pl2->c()) ||
                (pl1->a() == - pl2->a() && pl1->b() == - pl2->b() && pl1->c() == - pl2->c()));
    };

    // @todo could avoid some computations: 1&2 + 2&3 ==> no need to do 1&3,
    // but it would be a lot more code and branches
    int mask = (are_planes_parallel(plane_0, plane_1) << 5) |
               (are_planes_parallel(plane_0, plane_2) << 4) |
               (are_planes_parallel(plane_0, plane_3) << 3) |
               (are_planes_parallel(plane_1, plane_2) << 2) |
               (are_planes_parallel(plane_1, plane_3) << 1) |
               (are_planes_parallel(plane_2, plane_3));

    // Switch based on the bitmask to return the appropriate enum
    switch (mask) {
        case 0b111111: return PLANES_PARALLEL_0123;    // All planes are parallel
        case 0b111000: return PLANES_PARALLEL_012;     // 0, 1, 2 are parallel
        case 0b101100: return PLANES_PARALLEL_013;     // 0, 1, 3 are parallel
        case 0b011010: return PLANES_PARALLEL_023;     // 0, 2, 3 are parallel
        case 0b000111: return PLANES_PARALLEL_123;     // 1, 2, 3 are parallel
        case 0b100001: return PLANES_PARALLEL_01_23;   // 0, 1 are parallel, and 2, 3 are parallel
        case 0b001001: return PLANES_PARALLEL_02_13;   // 0, 2 are parallel, and 1, 3 are parallel
        case 0b010010: return PLANES_PARALLEL_03_12;   // 0, 3 are parallel, and 1, 2 are parallel
        case 0b100000: return PLANES_PARALLEL_01;      // 0, 1 are parallel
        case 0b010000: return PLANES_PARALLEL_02;      // 0, 2 are parallel
        case 0b001000: return PLANES_PARALLEL_03;      // 0, 3 are parallel
        case 0b000100: return PLANES_PARALLEL_12;      // 1, 2 are parallel
        case 0b000010: return PLANES_PARALLEL_13;      // 1, 3 are parallel
        case 0b000001: return PLANES_PARALLEL_23;      // 2, 3 are parallel
        default:       return PLANES_PARALLELISM_NONE; // No planes are parallel
    }
}

std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::vanishesAtGeneric(FacetSPtr facet_0,
                                                                      FacetSPtr facet_1,
                                                                      FacetSPtr facet_2,
                                                                      FacetSPtr facet_3)
{
    Plane3SPtr plane_0 = facet_0->plane();
    Plane3SPtr plane_1 = facet_1->plane();
    Plane3SPtr plane_2 = facet_2->plane();
    Plane3SPtr plane_3 = facet_3->plane();

    // @fixme in theory, one needs to check if the data exists etc. etc.
    CGAL::FT speed_0 = std::dynamic_pointer_cast<SkelFacetData>(facet_0->getData())->getSpeed();
    CGAL::FT speed_1 = std::dynamic_pointer_cast<SkelFacetData>(facet_1->getData())->getSpeed();
    CGAL::FT speed_2 = std::dynamic_pointer_cast<SkelFacetData>(facet_2->getData())->getSpeed();
    CGAL::FT speed_3 = std::dynamic_pointer_cast<SkelFacetData>(facet_3->getData())->getSpeed();

    Point3SPtr point = Point3SPtr();
    CGAL::FT event_offset;
    std::tie(point, event_offset) = KernelWrapper::intersectionAndTimeOffsetPlanes(
      plane_0, speed_0, plane_1, speed_1, plane_2, speed_2, plane_3, speed_3);

    if(point && event_offset <= 0) {
        return { point, event_offset };
    } else {
        return { };
    }
}

std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::vanishesAtOnePairOpposite(FacetSPtr facet_0,
                                                                              FacetSPtr facet_1,
                                                                              FacetSPtr facet_2,
                                                                              FacetSPtr facet_3)
{
    std::cout << "vanishesAtOnePairOpposite()" << std::endl;

    Point3SPtr point = Point3SPtr();

    Plane3SPtr plane_0 = facet_0->plane();
    Plane3SPtr plane_1 = facet_1->plane();
    Plane3SPtr plane_2 = facet_2->plane();
    Plane3SPtr plane_3 = facet_3->plane();

    CGAL_precondition(CGAL::parallel(*plane_1, *plane_3));

    CGAL::FT speed_0 = std::dynamic_pointer_cast<SkelFacetData>(facet_0->getData())->getSpeed();
    CGAL::FT speed_1 = std::dynamic_pointer_cast<SkelFacetData>(facet_1->getData())->getSpeed();
    CGAL::FT speed_2 = std::dynamic_pointer_cast<SkelFacetData>(facet_2->getData())->getSpeed();
    CGAL::FT speed_3 = std::dynamic_pointer_cast<SkelFacetData>(facet_3->getData())->getSpeed();

    // The facet #1 and #3 are parallel
    // - Get the time of intersection of #1 and #3 as the intersection is at that time
    // - Intersec the three non planar planes at that time

    Vector3SPtr n_1 = KernelFactory::createVector3(plane_1);
    Vector3SPtr n_3 = KernelFactory::createVector3(plane_3);
    Point3 seed_1 = plane_1->point();
    Line3 orth_line (seed_1, *n_1);

#ifdef CGAL_SS3_DEBUG_VANISH_AT
    std::cout << "facet1 = " << facet_1->getID() << std::endl;
    std::cout << "facet3 = " << facet_3->getID() << std::endl;
    std::cout << "plane1 = " << *plane_1 << std::endl;
    std::cout << "plane3 = " << *plane_3 << std::endl;
    std::cout << "n1 = " << *n_1 << std::endl;
    std::cout << "n3 = " << *n_3 << std::endl;
    std::cout << "speed1 = " << speed_1 << std::endl;
    std::cout << "speed3 = " << speed_3 << std::endl;
    std::cout << "seed_1 = " << seed_1 << std::endl;
#endif

    CGAL::Object obj = CGAL::intersection(orth_line, *plane_3);
    if (const CGAL::Point3* seed_3_ptr = CGAL::object_cast<CGAL::Point3>(&obj)) {
#ifdef CGAL_SS3_DEBUG_VANISH_AT
        std::cout << "seed_3 = " << *seed_3_ptr << std::endl;
#endif

        // (1): p = seed_1 + t * speed_1 * n_1
        // (2): p = seed_3 + t * speed_3 * n_3
        // t = (seed_3.i - seed_1.i) / (speed_1 * n_1.i - speed_3 * n_3.i)

        CGAL::FT t;
        for (int i=0; i<3; ++i) {
            if(is_zero(n_1->operator[](i)) && is_zero(n_3->operator[](i)))
                continue;

            const CGAL::FT den = speed_1 * n_1->operator[](i) - speed_3 * n_3->operator[](i);
            if (is_zero(den)) {
                std::cerr << "Moving in the same direction, no event" << std::endl;
                return { };
            }

            t = (seed_3_ptr->operator[](i) - seed_1[i]) / den;
            break;
        }

        if (t > 0) {
            std::cerr << "Event in the past" << std::endl;
            return { };
        }

        Plane3SPtr shifted_planed_0 = KernelWrapper::offsetPlane(plane_0, t*speed_0);
        Plane3SPtr shifted_planed_1 = KernelWrapper::offsetPlane(plane_1, t*speed_1);
        Plane3SPtr shifted_planed_2 = KernelWrapper::offsetPlane(plane_2, t*speed_2);

        point = KernelWrapper::intersection(shifted_planed_0, shifted_planed_1, shifted_planed_2);
        if (!point) {
            std::cerr << "No intersection of shifted planes" << std::endl;
            return { };
        } else {
            std::cout << "---> " << *point << " at " << t << std::endl;
            return { point, t };
        }
    } else {
        CGAL_assertion_msg(false, "Unknown second point...?");
    }

    return { };
}

std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::vanishesAtOnePairContiguous(FacetSPtr,
                                                                                FacetSPtr,
                                                                                FacetSPtr,
                                                                                FacetSPtr)
{
    CGAL_assertion(false);
    return {};
}

std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::vanishesAtTwoPairs(FacetSPtr,
                                                                       FacetSPtr,
                                                                       FacetSPtr,
                                                                       FacetSPtr)
{
    // the only possibility is if the middle edge is degenerate?
    // otherwise, the bisecting planes going through edges O2 and 13 are parallel and not intersecting?

    CGAL_assertion(false);
    return {};
}

std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::vanishesAt(EdgeSPtr edge)
{
    Point3SPtr point = Point3SPtr();
    CGAL::FT offset_event;

// #define CGAL_SS3_DEBUG_VANISH_AT
#ifdef CGAL_SS3_DEBUG_VANISH_AT
    std::cout << "vanishesAt " << edge->toString() << std::endl;
#endif

    // @fixme this seemingly assumes that extremities of 'edge' have degree 3?
    // what happens if it's not the case, shouldn't we consider both fL fR fLP fLN,
    // but also fL fR fRP fRN? OR even fL fR & (2 out of all faces incident to any
    // extremity of 'edge')?

// #define CGAL_SS3_OLD_CODE_VANISH_AT
#ifdef CGAL_SS3_OLD_CODE_VANISH_AT
    SheetSPtr sheets[3];
    for (unsigned int i = 0; i < 3; i++) {
        sheets[i] = SheetSPtr();
    }
    FacetSPtr facet = edge->getFacetL();
    if (!facet) {
        facet = edge->getFacetR();
    }
    SkelEdgeDataSPtr data = std::dynamic_pointer_cast<SkelEdgeData>(edge->getData());
    if (data) {
        sheets[0] = data->getSheet();
    }
    EdgeSPtr edge_prev = edge->prev(facet);
    data = std::dynamic_pointer_cast<SkelEdgeData>(edge_prev->getData());
    if (data) {
        sheets[1] = data->getSheet();
    }
    EdgeSPtr edge_next = edge->next(facet);
    data = std::dynamic_pointer_cast<SkelEdgeData>(edge_next->getData());
    if (data) {
        sheets[2] = data->getSheet();
    }
    if (sheets[0] && sheets[1] && sheets[2]) {
        Point3SPtr p_intersect = KernelWrapper::intersection(
                sheets[0]->getPlane(),
                sheets[1]->getPlane(),
                sheets[2]->getPlane());
        if (p_intersect) {
            // @fixme? negative offsets only?
            // @fixme should it check with the other incident facet as well?
            if (KernelWrapper::side(facet->plane(), p_intersect) <= 0) {
                // inside polyhedron
                point = p_intersect;
                offset_event = offsetDist(facet_l, point);
            }
        }
    }
#else // CGAL_SS3_OLD_CODE_VANISH_AT
    {
        FacetSPtr facetL = edge->getFacetL();
        FacetSPtr facetR = edge->getFacetR();
        CGAL_assertion(facetL && facetR && facetL != facetR);
# ifdef CGAL_SS3_DEBUG_VANISH_AT
        std::cout << "facetL: " << facetL->getID() << std::endl;
        std::cout << "facetR: " << facetR->getID() << std::endl;
# endif
        FacetSPtr facetP = edge->prev(facetL)->other(facetL);
        FacetSPtr facetN = edge->next(facetL)->other(facetL);
# ifdef CGAL_SS3_DEBUG_VANISH_AT
        if (facetP) {
            std::cout << "facetP: " << facetP->getID() << std::endl;
        }
        if (facetN) {
            std::cout << "facetN: " << facetN->getID() << std::endl;
        }
# endif
        CGAL_assertion(facetP && facetP != facetL && facetP != facetR);
        CGAL_assertion(facetN && facetN != facetL && facetN != facetR && facetN != facetP);

        Quadplane_parallelism parallelism = quadplane_parallelism(facetL, facetP, facetR, facetN);

        if(parallelism != PLANES_PARALLELISM_NONE)
          std::cout << "parallelism = " << parallelism << std::endl;

        if (parallelism == PLANES_PARALLELISM_NONE) {
            // generic case
            std::tie(point, offset_event) = vanishesAtGeneric(facetL, facetP, facetR, facetN);
        } else if (parallelism == PLANES_PARALLEL_01) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtOnePairContiguous();
        } else if (parallelism == PLANES_PARALLEL_02) {
            // if the two main faces are parallel, everything should be parallel
            CGAL_assertion_msg(false, "This configuration shouldn't be possible with vanish events");
        } else if (parallelism == PLANES_PARALLEL_03) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtOnePairContiguous();
        } else if (parallelism == PLANES_PARALLEL_12) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtOnePairContiguous();
        } else if (parallelism == PLANES_PARALLEL_13) {
            std::tie(point, offset_event) = vanishesAtOnePairOpposite(facetL, facetP, facetR, facetN);
        } else if (parallelism == PLANES_PARALLEL_23) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtOnePairContiguous();
        } else if (parallelism == PLANES_PARALLEL_012 ||
                   parallelism == PLANES_PARALLEL_013 ||
                   parallelism == PLANES_PARALLEL_023 ||
                   parallelism == PLANES_PARALLEL_123) {
         // if 3 are parallel, it should be all 4 since the faces are contiguous
          CGAL_assertion_msg(false, "This configuration shouldn't be possible with vanish events");
        } else if (parallelism == PLANES_PARALLEL_02_13) {
            CGAL_assertion_msg(false, "This configuration shouldn't be possible with vanish events");
        } else if (parallelism == PLANES_PARALLEL_03_12) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtTwoPairs();
        } else if (parallelism == PLANES_PARALLEL_01_23) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtTwoPairs();
        } else if (parallelism == PLANES_PARALLEL_0123) {
            // nothing happens
        } else {
            CGAL_assertion_msg(false, "Unknown parallelism enum value");
        }
    }
#endif // CGAL_SS3_OLD_CODE_VANISH_AT

    return { point, offset_event };
}

// edge: edge that is vanishing or crashing into another edge
// f: one of the faces incident to the edge
// t: the time of vanishing or the time of crash
// f_third:
bool
SimpleStraightSkel::check_bisector(EdgeSPtr edge,
                                   FacetSPtr f,
                                   CGAL::FT t,
                                   FacetSPtr f_third,
                                   Point3SPtr point)
{
    Plane3SPtr plane_third = f_third->plane();
    CGAL::FT a_third = plane_third->a();
    CGAL::FT b_third = plane_third->b();
    CGAL::FT c_third = plane_third->c();
    CGAL::FT d_third = plane_third->d();
    CGAL::FT speed_third = std::dynamic_pointer_cast<SkelFacetData>(f_third->getData())->getSpeed();
    CGAL::FT t_third = (a_third * point->x() + b_third * point->y() + c_third * point->z() + d_third) / speed_third;

    EdgeSPtr edge_f_f_third = edge->next(f);

    // std::cout << "edge = " << edge->toString() << std::endl;
    // std::cout << "edge_f_f_third = " << edge_f_f_third->toString() << std::endl;
    // std::cout << "time from third = " << t_third << std::endl;
    CGAL_assertion(edge_f_f_third != edge);

    // we want SMALLER to mean "on F1 side of the bisector"
    // and     LARGER  to mean "on third side of the bisector"
    CGAL::Comparison_result f_third_point_side;
    if (isReflex(edge_f_f_third)) {
        // std::cout << "F1O: reflex edge " << std::endl;
        f_third_point_side = CGAL::compare(t, t_third);
    } else {
        // std::cout << "F1O: convex edge " << std::endl;
        f_third_point_side = CGAL::compare(t_third, t);
    }

    // std::cout << "f_third_point_side = " << f_third_point_side << std::endl;

    // Now, we have to determine which side of f-third is legal
    // and that is determined by the angle that the face makes at the vertex
    // common
    //
    // We can't use the actual geometry of the edges because they might be degenerate,
    // so everything must be done with planes
    //
    // @todo construction-less

    FacetSPtr f_other = edge->other(f);
    Plane3SPtr plane_f = f->plane();
    Plane3SPtr plane_f_other = f_other->plane();
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

    // std::cout << "pv = " << pvp << std::endl;
    // std::cout << "v = " << vp << std::endl;
    // std::cout << "nv = " << nvp << std::endl;
    // std::cout << "n_f = " << *n_f << std::endl;
    // std::cout << "n_f_other = " << *n_f_other << std::endl;
    // std::cout << "n_f_third = " << *n_f_third << std::endl;
    // std::cout << "n = " << n << std::endl;
    // std::cout << "sp = " << sp << std::endl;

    // @todo below means that the edge is aligned, so the two bisectors are coplanar
    // and it'll be annoying. In that case, we can simply check on which side the
    // other vertex is and use the third bisector... but currently it's not possible
    // because no 2 faces are coplanar, so below should be always true!
    CGAL_assertion(sp != 0);

    if (sp < 0) {
        // std::cout << "SP: right turn" << std::endl;
        if (f_third_point_side == CGAL::SMALLER) {
            // std::cout << "reject!" << std::endl;
            return false;
        }
    } else {
        // std::cout << "SP: left turn" << std::endl;
        if (f_third_point_side == CGAL::LARGER) {
            // std::cout << "reject!" << std::endl;
            return false;
        }
    }

    return true;
}

// this one does an early exit if the result is irrelevant (in the past or too far in the)
std::pair<Point3SPtr, CGAL::FT>
SimpleStraightSkel::crashAt(EdgeSPtr edge_1, EdgeSPtr edge_2,
                            const std::optional<CGAL::FT> offset_max)
{
    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(edge_1->getData());
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(edge_2->getData());

    FacetSPtr facet_l1 = edge_1->getFacetL();
    FacetSPtr facet_r1 = edge_1->getFacetR();
    FacetSPtr facet_l2 = edge_2->getFacetL();
    FacetSPtr facet_r2 = edge_2->getFacetR();
    Plane3SPtr plane_l1 = facet_l1->plane();
    Plane3SPtr plane_r1 = facet_r1->plane();
    Plane3SPtr plane_l2 = facet_l2->plane();
    Plane3SPtr plane_r2 = facet_r2->plane();
    CGAL::FT speed_l1 = std::dynamic_pointer_cast<SkelFacetData>(facet_l1->getData())->getSpeed();
    CGAL::FT speed_r1 = std::dynamic_pointer_cast<SkelFacetData>(facet_r1->getData())->getSpeed();
    CGAL::FT speed_l2 = std::dynamic_pointer_cast<SkelFacetData>(facet_l2->getData())->getSpeed();
    CGAL::FT speed_r2 = std::dynamic_pointer_cast<SkelFacetData>(facet_r2->getData())->getSpeed();

    Quadplane_parallelism parallelism = quadplane_parallelism(facet_l1, facet_r1, facet_l2, facet_r2);

// #define CGAL_SS3_DEBUG_CRASH_AT
#ifdef CGAL_SS3_DEBUG_CRASH_AT
    std::cout << "    -- Crash At\n    " << edge_1->toString() << "\n    " << edge_2->toString() << std::endl;

    CGAL_warning(*(edge_1->getVertexSrc()->getPoint()) == *(edge_1->getVertexDst()->getPoint()));
    CGAL_warning(*(edge_2->getVertexSrc()->getPoint()) == *(edge_2->getVertexDst()->getPoint()));

    std::cout << "Facet L1" << std::endl;
    std::cout << facet_l1->getID() << std::endl;
    std::cout << "Facet R1" << std::endl;
    std::cout << facet_r1->getID() << std::endl;
    std::cout << "Facet L2" << std::endl;
    std::cout << facet_l2->getID() << std::endl;
    std::cout << "Facet R2" << std::endl;
    std::cout << facet_r2->getID() << std::endl;

    std::cout << "Parallelism: " << parallelism << std::endl;
#endif

    if (parallelism != PLANES_PARALLELISM_NONE) {
        std::cerr << "WARNING: QUADFACES WITH PARALLEL PLANES!" << std::endl;
        std::cerr << "parallelism = " << parallelism << std::endl;
        if (parallelism == PLANES_PARALLEL_01 ||
            parallelism == PLANES_PARALLEL_23 ||
            parallelism == PLANES_PARALLEL_012 ||
            parallelism == PLANES_PARALLEL_013 ||
            parallelism == PLANES_PARALLEL_023 ||
            parallelism == PLANES_PARALLEL_123 ||
            parallelism == PLANES_PARALLEL_01_23 ||
            parallelism == PLANES_PARALLEL_0123) {
            std::exit(1);
        }
    }

    Point3SPtr point;
    CGAL::FT offset_event;
    std::tie(point, offset_event) = KernelWrapper::intersectionAndTimeOffsetPlanes(
        plane_l1, speed_l1, plane_r1, speed_r1, plane_l2, speed_l2, plane_r2, speed_r2);

    if (!point) {
        std::cerr << "Warning: no crashing?" << std::endl;
        return { };
    }

    CGAL_assertion(*point == *(KernelWrapper::intersectionOffsetPlanes(
        plane_l1, speed_l1, plane_r1, speed_r1, plane_l2, speed_l2, plane_r2, speed_r2)));

    // Ignore the intersection if it is in the past
    // This is equivalent to the check that the point lives
    // on the inward half plane (split by the edge's supporting line) of the bisector plane
    if (offset_event > 0) {
      return { };
    }

    // Ignore the event if it happens further the future compared to the soonest top event.
    if(offset_max && offset_max.value() > offset_event) {
#ifdef CGAL_SS3_DEBUG_CRASH_AT
        std::cout << "event is strictly in the future" << std::endl;
#endif
        return { };
    }

    // Check that the point is inside bounds
    FacetSPtr facet_1_src = getFacetSrc(edge_1);
    FacetSPtr facet_1_dst = getFacetDst(edge_1);
    FacetSPtr facet_2_src = getFacetSrc(edge_2);
    FacetSPtr facet_2_dst = getFacetDst(edge_2);

#ifdef CGAL_SS3_DEBUG_CRASH_AT
    std::cout << "Facet 1 SRC" << std::endl;
    std::cout << facet_1_src->toString() << std::endl;
    std::cout << "Facet 1 DST" << std::endl;
    std::cout << facet_1_dst->toString() << std::endl;
    std::cout << "Facet 2 SRC" << std::endl;
    std::cout << facet_2_src->toString() << std::endl;
    std::cout << "Facet 2 DST" << std::endl;
    std::cout << facet_2_dst->toString() << std::endl;

    std::cout << "Point: " << *point << std::endl;
#endif

// #define CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS
#if defined(CGAL_SS3_OLD_CODE_BOUND_CHECKS) || defined(CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS)
    bool reject_1 = false;
    bool reject_2 = false;
    bool reject_3 = false;
    bool reject_4 = false;
    bool reject_5 = false;
    bool reject_6 = false;
    Vector3SPtr normal_1 = KernelFactory::createVector3(data_1->getSheet()->getPlane());
    Line3SPtr line_normal_1 = KernelFactory::createLine3(point, normal_1);
    if (KernelWrapper::orientation(line(edge_1), line_normal_1) <= 0) {
        std::cout << "reject #1" << std::endl;
        reject_1 = true;
    }
    if (!(facet_1_src == facet_l2 ||
            facet_1_src == facet_r2)) {
        SkelVertexDataSPtr data_1_src = std::dynamic_pointer_cast<SkelVertexData>(
            edge_1->getVertexSrc()->getData());
        ArcSPtr arc_1_src = data_1_src->getArc();
        if (KernelWrapper::orientation(arc_1_src->line(), line_normal_1) > 0) {
            std::cout << "reject #2" << std::endl;
            reject_2 = true;
        }
    }
    if (!(facet_1_dst == facet_l2 ||
            facet_1_dst == facet_r2)) {
        SkelVertexDataSPtr data_1_dst = std::dynamic_pointer_cast<SkelVertexData>(
            edge_1->getVertexDst()->getData());
        ArcSPtr arc_1_dst = data_1_dst->getArc();
        if (KernelWrapper::orientation(arc_1_dst->line(), line_normal_1) < 0) {
            std::cout << "reject #3" << std::endl;
            reject_3 = true;
        }
    }
    Vector3SPtr normal_2 = KernelFactory::createVector3(data_2->getSheet()->getPlane());
    Line3SPtr line_normal_2 = KernelFactory::createLine3(point, normal_2);
    if (KernelWrapper::orientation(line(edge_2), line_normal_2) <= 0) {
            std::cout << "reject #4" << std::endl;
            reject_4 = true;
    }
    if (!(facet_2_src == facet_l1 ||
            facet_2_src == facet_r1)) {
        SkelVertexDataSPtr data_2_src = std::dynamic_pointer_cast<SkelVertexData>(
            edge_2->getVertexSrc()->getData());
        ArcSPtr arc_2_src = data_2_src->getArc();
        if (KernelWrapper::orientation(arc_2_src->line(), line_normal_2) > 0) {
            std::cout << "reject #5" << std::endl;
            reject_5 = true;
        }
    }
    if (!(facet_2_dst == facet_l1 ||
            facet_2_dst == facet_r1)) {
        SkelVertexDataSPtr data_2_dst = std::dynamic_pointer_cast<SkelVertexData>(
            edge_2->getVertexDst()->getData());
        ArcSPtr arc_2_dst = data_2_dst->getArc();
        if (KernelWrapper::orientation(arc_2_dst->line(), line_normal_2) < 0) {
            std::cout << "reject #6" << std::endl;
            reject_6 = true;
        }
    }

    const bool reject_old = (reject_1 || reject_2 || reject_3 || reject_4 || reject_5 || reject_6);
#endif

#if !defined(CGAL_SS3_OLD_CODE_BOUND_CHECKS) || defined(GAL_SS3_COMPARE_BOTH_BOUND_CHECKS)
    bool reject_2b = false;
    bool reject_3b = false;
    bool reject_5b = false;
    bool reject_6b = false;

    // Not to be in the past is equivalent to be on the negative side (we are shrinking)
    // of the planes of the facets incident to the edge
    //
    // Since we have filtered positive times, this shouldn't be needed (assuming positive weights)

    // @fixme assumes positive weights
    CGAL_assertion(!(KernelWrapper::side(plane_l1, point) > 0 ||
                     KernelWrapper::side(plane_r1, point) > 0));
    CGAL_assertion(!(KernelWrapper::side(plane_l2, point) > 0 ||
                     KernelWrapper::side(plane_r2, point) > 0));

    // @fixme this should be a predicate (oriented_side_of_event_point_wrt_bisectorC2)
    CGAL::FT l1a = plane_l1->a();
    CGAL::FT l1b = plane_l1->b();
    CGAL::FT l1c = plane_l1->c();
    CGAL::FT l1d = plane_l1->d();
    CGAL::FT r1a = plane_r1->a();
    CGAL::FT r1b = plane_r1->b();
    CGAL::FT r1c = plane_r1->c();
    CGAL::FT r1d = plane_r1->d();
    CGAL::FT l2a = plane_l2->a();
    CGAL::FT l2b = plane_l2->b();
    CGAL::FT l2c = plane_l2->c();
    CGAL::FT l2d = plane_l2->d();
    CGAL::FT r2a = plane_r2->a();
    CGAL::FT r2b = plane_r2->b();
    CGAL::FT r2c = plane_r2->c();
    CGAL::FT r2d = plane_r2->d();

    CGAL::FT lt1 = (l1a * point->x() + l1b * point->y() + l1c * point->z() + l1d) / speed_l1;
    CGAL::FT rt1 = (r1a * point->x() + r1b * point->y() + r1c * point->z() + r1d) / speed_r1;
    CGAL::FT lt2 = (l2a * point->x() + l2b * point->y() + l2c * point->z() + l2d) / speed_l2;
    CGAL::FT rt2 = (r2a * point->x() + r2b * point->y() + r2c * point->z() + r2d) / speed_r2;
    // std::cout << "time from l1 " << lt1 << std::endl;
    // std::cout << "time from r1 " << rt1 << std::endl;
    // std::cout << "time from l2 " << lt2 << std::endl;
    // std::cout << "time from r2 " << rt2 << std::endl;
    CGAL_assertion(lt1 == rt1 && lt1 == lt2 && lt1 == rt2);


    // We want the point to be left of the right arc, and right of the left arc
    //
    // We can just check the position of the query 'point' with respect
    // to one of the other bisector. A few tricky parts:
    // - The side of bisector is easy to know: it's comparing times, but that
    //   depends on whether the edge is reflex or convex
    // - Once we know on which side of the bisector it is, it's not done yet:
    //   we need to know which side of the bisector is the correct one and since
    //   faces can be non-convex polygons, we need to check if we are in a concave
    //   (within the face) vertex to know if the clipping bisector is inverted

    // @todo lots of duplicate computations
    // @todo exit as soon as there is a single rejection (to be done once the old code is safe to remove)

    if (!(facet_1_src == facet_l2 || // @fixme are these facet checks required?
            facet_1_src == facet_r2)) {
        // std::cout << "-- check 1_src" << std::endl;
        // src is the target of the edge when in the right face
        if (!check_bisector(edge_1, facet_r1, rt1, facet_1_src, point)) {
            // std::cout << "reject #2b" << std::endl;
            reject_2b = true;
        }
    }

    if (!(facet_1_dst == facet_l2 ||
            facet_1_dst == facet_r2)) {
        // std::cout << "-- check 1_dst" << std::endl;
        if (!check_bisector(edge_1, facet_l1, lt1, facet_1_dst, point)) {
            // std::cout << "reject #3b" << std::endl;
            reject_3b = true;
        }
    }

    if (!(facet_2_src == facet_l1 ||
            facet_2_src == facet_r1)) {
        // std::cout << "-- check 2_src" << std::endl;
        if (!check_bisector(edge_2, facet_r2, rt2, facet_2_src, point)) {
            // std::cout << "reject #5b" << std::endl;
            reject_5b = true;
        }
    }

    if (!(facet_2_dst == facet_l1 ||
            facet_2_dst == facet_r1)) {
          // std::cout << "-- check 2_dst" << std::endl;
          if (!check_bisector(edge_2, facet_l2, lt2, facet_2_dst, point)) {
              // std::cout << "reject #6b" << std::endl;
              reject_6b = true;
        }
    }


    const bool reject_new = (reject_2b || reject_3b || reject_5b || reject_6b);
#endif

#ifdef CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS
    CGAL_assertion(reject_2 == reject_2b);
    CGAL_assertion(reject_3 == reject_3b);
    CGAL_assertion(reject_5 == reject_5b);
    CGAL_assertion(reject_6 == reject_6b);

    CGAL_assertion(reject_old == reject_new);
#endif

    if(reject_new) {
        return { };
    }

    return { point, offset_event };
}

CGAL::FT SimpleStraightSkel::offsetDist(FacetSPtr facet, Point3SPtr point) {
    Plane3SPtr plane = facet->plane();
    CGAL::FT result = KernelWrapper::distance(plane, point);
    if (KernelWrapper::side(plane, point) < 0) {
        result *= -1.0;
    }
    if (facet->hasData()) {
        CGAL::FT speed = std::dynamic_pointer_cast<SkelFacetData>(
                facet->getData())->getSpeed();
        result /= speed;
    }
    return result;
}

#if 0 // individual events of a given type
EdgeEventSPtr SimpleStraightSkel::nextEdgeEvent(PolyhedronSPtr polyhedron,
                                                const CGAL::FT current_offset,
                                                CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Edge Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    EdgeEventSPtr result = EdgeEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
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
        std::tie(point, offset_event) = vanishesAt(edge);
        if (!point) {
            continue;
        }

        std::cout << "Tentative edge event @ " << edge->toString() << std::endl;
        std::cout << "Would give: " << *point << " (at t=" << offset_event << ")" << std::endl;

        FacetSPtr facet_src = getFacetSrc(edge);
        FacetSPtr facet_dst = getFacetDst(edge);

        // This does not work when there is more than one edge between both facets.
        // EdgeSPtr edge_2 = facet_src->findEdge(facet_dst);
        std::list<EdgeSPtr> edges_2 = facet_src->findEdges(facet_dst);

        bool split_event = false;
        std::list<EdgeSPtr>::iterator it_e2 = edges_2.begin();
        while (it_e2 != edges_2.end()) {
            EdgeSPtr edge_2 = *it_e2++;

            // @todo exit as soon as there is a single rejection (to be done once the old code is safe to remove)
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
            FacetSPtr facet_2_src = getFacetSrc(edge_2);
            FacetSPtr facet_2_dst = getFacetDst(edge_2);

            Plane3SPtr plane_l2 = facet_l2->plane();
            Plane3SPtr plane_r2 = facet_r2->plane();
            CGAL::FT speed_l2 = std::dynamic_pointer_cast<SkelFacetData>(facet_l2->getData())->getSpeed();
            CGAL::FT speed_r2 = std::dynamic_pointer_cast<SkelFacetData>(facet_r2->getData())->getSpeed();

            CGAL::FT l2a = plane_l2->a();
            CGAL::FT l2b = plane_l2->b();
            CGAL::FT l2c = plane_l2->c();
            CGAL::FT l2d = plane_l2->d();
            CGAL::FT r2a = plane_r2->a();
            CGAL::FT r2b = plane_r2->b();
            CGAL::FT r2c = plane_r2->c();
            CGAL::FT r2d = plane_r2->d();

            CGAL::FT lt2 = (l2a * point->x() + l2b * point->y() + l2c * point->z() + l2d) / speed_l2;
            CGAL::FT rt2 = (r2a * point->x() + r2b * point->y() + r2c * point->z() + r2d) / speed_r2;
            CGAL_assertion(lt2 == rt2);

            if ((lt2 > 0) || (rt2 > 0)) {
                split_event_current_1_b = false;
            }

            if (!check_bisector(edge_2, facet_r2, rt2, facet_2_src, point)) {
                split_event_current_2_b = false;
            }

            if (!check_bisector(edge_2, facet_l2, lt2, facet_2_dst, point)) {
                split_event_current_3_b = false;
            }

            const bool split_event_current_b = (split_event_current_1_b &&
                                                split_event_current_2_b &&
                                                split_event_current_3_b);
#endif

#ifdef CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS
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
        // @todo could move up this check
        if (offset_event <= 0 && offset_event > current_offset_to_nearest_event)
        {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = EdgeEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setEdge(edge);
            SkelVertexDataSPtr data_src = std::dynamic_pointer_cast<SkelVertexData>(
                    edge->getVertexSrc()->getData());
            SkelVertexDataSPtr data_dst = std::dynamic_pointer_cast<SkelVertexData>(
                    edge->getVertexDst()->getData());
            node->addArc(data_src->getArc());
            node->addArc(data_dst->getArc());
            SkelEdgeDataSPtr data_edge = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge->getData());
            node->addSheet(data_edge->getSheet());
            current_offset_to_nearest_event = offset_event;
        }
    }

    return result;
}

EdgeMergeEventSPtr SimpleStraightSkel::nextEdgeMergeEvent(PolyhedronSPtr polyhedron,
                                                          const CGAL::FT current_offset,
                                                          CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Edge Merge Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    EdgeMergeEventSPtr result = EdgeMergeEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;

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
        EdgeSPtr edge_prev = edge->prev(facet_l);
        edge_next = edge->next(facet_l)->next(facet_l);
        if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
            facet = facet_l;
            edge_1 = edge_prev;
            edge_2 = edge_next;
        }
        // @todo do we need to test these other combinations if above matched?
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
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge);
        if (!point) {
            continue;
        }
        if (offset_event <= 0 && offset_event > current_offset_to_nearest_event) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = EdgeMergeEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setFacet(facet);
            result->setEdge1(edge_1);
            result->setEdge2(edge_2);
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
            current_offset_to_nearest_event = offset_event;
        }
    }

    if(result)
    {
      std::cout << "Best Edge Merge Event with time increment of: " << result->getOffset() << std::endl;
    }

    return result;
}

TriangleEventSPtr SimpleStraightSkel::nextTriangleEvent(PolyhedronSPtr polyhedron,
                                                        const CGAL::FT current_offset,
                                                        CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Triangle Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    TriangleEventSPtr result = TriangleEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;

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
        std::tie(point, offset_event) = vanishesAt(edge);
        if (!point) {
            continue;
        }

        if ((KernelWrapper::side(edge->getFacetL()->plane(), point) > 0) ||
                KernelWrapper::side(edge->getFacetR()->plane(), point) > 0) {
            // triangle may not be a hole
            // after pierce event
            continue;
        }

        if (offset_event <= 0 && offset_event > current_offset_to_nearest_event) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = TriangleEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setFacet(facet);
            result->setEdgeBegin(edge);

            VertexSPtr vertices[3];
            result->getVertices(vertices);
            for (unsigned int i = 0; i < 3; i++) {
                SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(
                        vertices[i]->getData());
                ArcSPtr arc = data->getArc();
                node->addArc(arc);
            }
            EdgeSPtr edges[3];
            result->getEdges(edges);
            for (unsigned int i = 0; i < 3; i++) {
                SkelEdgeDataSPtr data = std::dynamic_pointer_cast<SkelEdgeData>(
                        edges[i]->getData());
                SheetSPtr sheet = data->getSheet();
                node->addSheet(sheet);
            }

            current_offset_to_nearest_event = offset_event;
        }
    }
    return result;
}

DblEdgeMergeEventSPtr SimpleStraightSkel::nextDblEdgeMergeEvent(PolyhedronSPtr polyhedron,
                                                                const CGAL::FT current_offset,
                                                                CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Dbl Edge Merge Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    DblEdgeMergeEventSPtr result = DblEdgeMergeEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;

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
        std::tie(point, offset_event) = vanishesAt(edge);
        if (!point) {
            continue;
        }
        if (offset_event <= 0 && offset_event > current_offset_to_nearest_event) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = DblEdgeMergeEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setFacet1(facet_1);
            result->setEdge11(edge_11);
            result->setEdge12(edge_12);
            result->setFacet2(facet_2);
            result->setEdge21(edge_21);
            result->setEdge22(edge_22);
            VertexSPtr vertices[4];
            result->getVertices(vertices);
            for (unsigned int i = 0; i < 4; i++) {
                SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                        vertices[i]->getData());
                ArcSPtr arc = vertex_data->getArc();
                node->addArc(arc);
            }
            EdgeSPtr edges[4];
            result->getEdges(edges);
            for (unsigned int i = 0; i < 4; i++) {
                SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                        edges[i]->getData());
                SheetSPtr sheet = edge_data->getSheet();
                node->addSheet(sheet);
            }
            current_offset_to_nearest_event = offset_event;
        }
    }
    return result;
}

DblTriangleEventSPtr SimpleStraightSkel::nextDblTriangleEvent(PolyhedronSPtr polyhedron,
                                                              const CGAL::FT current_offset,
                                                              CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Dbl Triangle Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    DblTriangleEventSPtr result = DblTriangleEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (isTetrahedron(edge)) {
            continue;
        }
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (!facet_l || !facet_r) {
            continue;
        }
        if (!(isTriangle(facet_l, edge) &&
                isTriangle(facet_r, edge))) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge);
        if (!point) {
            continue;
        }
        if (offset_event <= 0 && offset_event > current_offset_to_nearest_event) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = DblTriangleEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setEdge(edge);

            VertexSPtr vertices[4];
            result->getVertices(vertices);
            for (unsigned int i = 0; i < 4; i++) {
                SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(
                        vertices[i]->getData());
                ArcSPtr arc = data->getArc();
                node->addArc(arc);
            }
            EdgeSPtr edges[5];
            result->getEdges(edges);
            for (unsigned int i = 0; i < 5; i++) {
                SkelEdgeDataSPtr data = std::dynamic_pointer_cast<SkelEdgeData>(
                        edges[i]->getData());
                SheetSPtr sheet = data->getSheet();
                node->addSheet(sheet);
            }
            current_offset_to_nearest_event = offset_event;
        }
    }
    return result;
}

TetrahedronEventSPtr SimpleStraightSkel::nextTetrahedronEvent(PolyhedronSPtr polyhedron,
                                                              const CGAL::FT current_offset,
                                                              CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Tetrahedron Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    TetrahedronEventSPtr result = TetrahedronEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (isTetrahedron(edge)) {
            FacetSPtr facet = edge->getFacetL();
            Point3SPtr point = Point3SPtr();
            CGAL::FT offset_event;
            std::tie(point, offset_event) = vanishesAt(edge);
            if (!point) {
                continue;
            }
            if (offset_event <= 0 && offset_event > current_offset_to_nearest_event) {
                NodeSPtr node;
                if (!result) {
                    node = Node::create(point);
                    result = TetrahedronEvent::create(polyhedron);
                    result->setNode(node);
                }
                node = result->getNode();
                node->clear();
                node->setOffset(current_offset + offset_event);
                node->setPoint(point);
                result->setEdgeBegin(edge);
                VertexSPtr vertices[4];
                result->getVertices(vertices);
                for (unsigned int i = 0; i < 4; i++) {
                    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                            vertices[i]->getData());
                    ArcSPtr arc = vertex_data->getArc();
                    node->addArc(arc);
                }
                EdgeSPtr edges[6];
                result->getEdges(edges);
                for (unsigned int i = 0; i < 6; i++) {
                    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                            edges[i]->getData());
                    SheetSPtr sheet = edge_data->getSheet();
                    node->addSheet(sheet);
                }

                current_offset_to_nearest_event = offset_event;
            }
        }
    }
    return result;
}

VertexEventSPtr SimpleStraightSkel::nextVertexEvent(PolyhedronSPtr polyhedron,
                                                    const CGAL::FT current_offset,
                                                    CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Vertex Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    VertexEventSPtr result = VertexEventSPtr();

    std::list<VertexSPtr>::iterator it_v1 = polyhedron->vertices().begin();
    while (it_v1 != polyhedron->vertices().end()) {
        VertexSPtr vertex_1 = *it_v1++;
        if (isConvex(vertex_1)) {
            continue;
        }

        std::list<VertexSPtr> vertices_2;
        std::list<FacetWPtr>::iterator it_f = vertex_1->facets().begin();
        while (it_f != vertex_1->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                vertices_2.insert(vertices_2.end(),
                        facet->vertices().begin(), facet->vertices().end());
            }
        }
        std::list<VertexSPtr>::iterator it_v2 = vertices_2.begin();
        while (it_v2 != vertices_2.end()) {
            VertexSPtr vertex_2 = *it_v2++;
            if (vertex_1 == vertex_2) {
                continue;
            }
            if (vertex_1->getPoint() == vertex_2->getPoint()) {
                continue;
            }
            if (isConvex(vertex_2)) {
                continue;
            }

            if (vertex_1->findEdge(vertex_2)) {
                // edge event
                continue;
            }

            FacetSPtr facet_1;
            FacetSPtr facet_2;
            int num_equal_facets = 0;
            std::list<FacetWPtr>::iterator it_f1 = vertex_1->facets().begin();
            while (it_f1 != vertex_1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = vertex_2->facets().begin();
                    while (it_f2 != vertex_2->facets().end()) {
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
            if (facet_1->next(vertex_1) != facet_2) {
                FacetSPtr facet_tmp = facet_1;
                facet_1 = facet_2;
                facet_2 = facet_tmp;
            }
            if (vertex_1->next(facet_1)->next(facet_1) == vertex_2 ||
                    vertex_1->next(facet_2)->next(facet_2) == vertex_2 ||
                    vertex_1->prev(facet_1)->prev(facet_1) == vertex_2 ||
                    vertex_1->prev(facet_2)->prev(facet_2) == vertex_2) {
                // edge merge event
                continue;
            }

            EdgeSPtr edge_11 = EdgeSPtr();
            EdgeSPtr edge_12 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
            while (it_e1 != vertex_1->edges().end()) {
                EdgeWPtr edge_1_wptr = *it_e1++;
                if (!edge_1_wptr.expired()) {
                    EdgeSPtr edge_1 = EdgeSPtr(edge_1_wptr);
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
            std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
            while (it_e2 != vertex_2->edges().end()) {
                EdgeWPtr edge_2_wptr = *it_e2++;
                if (!edge_2_wptr.expired()) {
                    EdgeSPtr edge_2 = EdgeSPtr(edge_2_wptr);
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
            if (!((edge_11->next(vertex_1) == edge_12 && edge_22->next(vertex_2) == edge_21) ||
                    (edge_12->next(vertex_1) == edge_11 && edge_21->next(vertex_2) == edge_22))) {
                // flip vertex event
                continue;
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
                continue;
            }

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion(offset_event > current_offset_to_nearest_event);

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = VertexEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_2->getData());
            node->addArc(data_2->getArc());
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setVertex1(vertex_1);
            result->setVertex2(vertex_2);
            result->setFacet1(facet_1);
            result->setFacet2(facet_2);
            current_offset_to_nearest_event = offset_event;
        }
    }
    return result;
}

FlipVertexEventSPtr SimpleStraightSkel::nextFlipVertexEvent(PolyhedronSPtr polyhedron,
                                                            const CGAL::FT current_offset,
                                                            CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Flip Vertex Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    FlipVertexEventSPtr result = FlipVertexEventSPtr();

    std::list<VertexSPtr>::iterator it_v1 = polyhedron->vertices().begin();
    while (it_v1 != polyhedron->vertices().end()) {
        VertexSPtr vertex_1 = *it_v1++;
        if (isConvex(vertex_1)) {
            continue;
        }

        std::list<VertexSPtr> vertices_2;
        std::list<FacetWPtr>::iterator it_f = vertex_1->facets().begin();
        while (it_f != vertex_1->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                vertices_2.insert(vertices_2.end(),
                        facet->vertices().begin(), facet->vertices().end());
            }
        }
        std::list<VertexSPtr>::iterator it_v2 = vertices_2.begin();
        while (it_v2 != vertices_2.end()) {
            VertexSPtr vertex_2 = *it_v2++;
            if (vertex_1 == vertex_2) {
                continue;
            }
            if (vertex_1->getPoint() == vertex_2->getPoint()) {
                continue;
            }
            if (isConvex(vertex_2)) {
                continue;
            }

            if (vertex_1->findEdge(vertex_2)) {
                // edge event
                continue;
            }

            FacetSPtr facet_1;
            FacetSPtr facet_2;
            int num_equal_facets = 0;
            std::list<FacetWPtr>::iterator it_f1 = vertex_1->facets().begin();
            while (it_f1 != vertex_1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = vertex_2->facets().begin();
                    while (it_f2 != vertex_2->facets().end()) {
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
            if (facet_1->next(vertex_1) != facet_2) {
                FacetSPtr facet_tmp = facet_1;
                facet_1 = facet_2;
                facet_2 = facet_tmp;
            }
            if (vertex_1->next(facet_1)->next(facet_1) == vertex_2 ||
                    vertex_1->next(facet_2)->next(facet_2) == vertex_2 ||
                    vertex_1->prev(facet_1)->prev(facet_1) == vertex_2 ||
                    vertex_1->prev(facet_2)->prev(facet_2) == vertex_2) {
                // edge merge event
                continue;
            }

            EdgeSPtr edge_11 = EdgeSPtr();
            EdgeSPtr edge_12 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
            while (it_e1 != vertex_1->edges().end()) {
                EdgeWPtr edge_1_wptr = *it_e1++;
                if (!edge_1_wptr.expired()) {
                    EdgeSPtr edge_1 = EdgeSPtr(edge_1_wptr);
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
            std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
            while (it_e2 != vertex_2->edges().end()) {
                EdgeWPtr edge_2_wptr = *it_e2++;
                if (!edge_2_wptr.expired()) {
                    EdgeSPtr edge_2 = EdgeSPtr(edge_2_wptr);
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
            if (!(edge_12->next(vertex_1) == edge_11 && edge_22->next(vertex_2) == edge_21)) {
                // vertex event
                continue;
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
                continue;
            }

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion(offset_event > current_offset_to_nearest_event);

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = FlipVertexEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_2->getData());
            node->addArc(data_2->getArc());
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setVertex1(vertex_1);
            result->setVertex2(vertex_2);
            result->setFacet1(facet_1);
            result->setFacet2(facet_2);
            current_offset_to_nearest_event = offset_event;
        }
    }
    return result;
}

SurfaceEventSPtr SimpleStraightSkel::nextSurfaceEvent(PolyhedronSPtr polyhedron,
                                                      const CGAL::FT current_offset,
                                                      CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Surface Event" << std::endl;
    CGAL::Real_timer timer;
    timer.start();

    ReadLock l(polyhedron->mutex());
    SurfaceEventSPtr result = SurfaceEventSPtr();

    std::list<EdgeSPtr>::iterator it_e1 = polyhedron->edges().begin();
    while (it_e1 != polyhedron->edges().end()) {
        EdgeSPtr edge_1 = *it_e1++;
        FacetSPtr facet_1_src = getFacetSrc(edge_1);
        FacetSPtr facet_1_dst = getFacetDst(edge_1);
        std::list<EdgeSPtr> edges_2;
        edges_2.insert(edges_2.end(),
                facet_1_src->edges().begin(), facet_1_src->edges().end());
        edges_2.insert(edges_2.end(),
                facet_1_dst->edges().begin(), facet_1_dst->edges().end());
        std::list<EdgeSPtr>::iterator it_e2 = edges_2.begin();
        while (it_e2 != edges_2.end()) {
            EdgeSPtr edge_2 = *it_e2++;
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
            FacetSPtr facet_2_src = getFacetSrc(edge_2);
            FacetSPtr facet_2_dst = getFacetDst(edge_2);
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
            bool is_conv_split_event = false;
            std::list<EdgeSPtr> edges = edge_1->getFacetL()->findEdges(edge_1->getFacetR());
            std::list<EdgeSPtr>::iterator it_e = edges.begin();
            while (it_e != edges.end()) {
                EdgeSPtr edge = *it_e++;
                if (edge == edge_1) {
                    continue;
                }
                FacetSPtr facet_src = getFacetSrc(edge);
                FacetSPtr facet_dst = getFacetDst(edge);
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

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset_to_nearest_event);
            if (!point)
                continue;

            CGAL_assertion(offset_event > current_offset_to_nearest_event);

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = SurfaceEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setEdge1(edge_1);
            result->setEdge2(edge_2);

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

            current_offset_to_nearest_event = offset_event;
        }
    }

    timer.stop();
    std::cout << "Sought Surface Event in: " << timer.time() << std::endl;

    return result;
}

PolyhedronSplitEventSPtr SimpleStraightSkel::nextPolyhedronSplitEvent(PolyhedronSPtr polyhedron,
                                                                      const CGAL::FT current_offset,
                                                                      CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Polyhedron Split Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    PolyhedronSplitEventSPtr result = PolyhedronSplitEventSPtr();

    std::list<EdgeSPtr>::iterator it_e1 = polyhedron->edges().begin();
    while (it_e1 != polyhedron->edges().end()) {
        EdgeSPtr edge_1 = *it_e1++;
        if (!isReflex(edge_1)) {
            continue;
        }
        FacetSPtr facet_1_src = getFacetSrc(edge_1);
        FacetSPtr facet_1_dst = getFacetDst(edge_1);
        std::list<EdgeSPtr>::iterator it_e2 = facet_1_src->edges().begin();
        while (it_e2 != facet_1_src->edges().end()) {
            EdgeSPtr edge_2 = *it_e2++;
            if (edge_1->getVertexSrc()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
                    edge_1->getVertexSrc()->getPoint() == edge_2->getVertexDst()->getPoint() ||
                    edge_1->getVertexDst()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
                    edge_1->getVertexDst()->getPoint() == edge_2->getVertexDst()->getPoint()) {
                // share a vertex
                continue;
            }
            if (!((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() == facet_1_dst) ||
                    (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() == facet_1_src))) {
                // no polyhedron split event
                continue;
            }
            if (edge_1->getVertexSrc()->findEdge(edge_2->getVertexSrc()) ||
                    edge_1->getVertexSrc()->findEdge(edge_2->getVertexDst()) ||
                    edge_1->getVertexDst()->findEdge(edge_2->getVertexSrc()) ||
                    edge_1->getVertexDst()->findEdge(edge_2->getVertexDst())) {
                // does not work when there is only one edge in between
                continue;
            }

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion (offset_event > current_offset_to_nearest_event);

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = PolyhedronSplitEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setEdge1(edge_1);
            result->setEdge2(edge_2);

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

            current_offset_to_nearest_event = offset_event;
        }
    }
    return result;
}

SplitMergeEventSPtr SimpleStraightSkel::nextSplitMergeEvent(PolyhedronSPtr polyhedron,
                                                            const CGAL::FT current_offset,
                                                            CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Split Merge Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    SplitMergeEventSPtr result = SplitMergeEventSPtr();

    std::list<VertexSPtr>::iterator it_v1 = polyhedron->vertices().begin();
    while (it_v1 != polyhedron->vertices().end()) {
        VertexSPtr vertex_1 = *it_v1++;
        if (isConvex(vertex_1)) {
            continue;
        }

        std::list<VertexSPtr> vertices_2;
        std::list<FacetWPtr>::iterator it_f = vertex_1->facets().begin();
        while (it_f != vertex_1->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                vertices_2.insert(vertices_2.end(),
                        facet->vertices().begin(), facet->vertices().end());
            }
        }
        std::list<VertexSPtr>::iterator it_v2 = vertices_2.begin();
        while (it_v2 != vertices_2.end()) {
            VertexSPtr vertex_2 = *it_v2++;
            if (vertex_1 == vertex_2) {
                continue;
            }
            if (vertex_1->getPoint() == vertex_2->getPoint()) {
                continue;
            }
            if (isConvex(vertex_2)) {
                continue;
            }

            if (vertex_1->findEdge(vertex_2)) {
                // edge event
                continue;
            }

            FacetSPtr facet_1;
            FacetSPtr facet_2;
            int num_equal_facets = 0;
            std::list<FacetWPtr>::iterator it_f1 = vertex_1->facets().begin();
            while (it_f1 != vertex_1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = vertex_2->facets().begin();
                    while (it_f2 != vertex_2->facets().end()) {
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
            if (facet_1->next(vertex_1) != facet_2) {
                FacetSPtr facet_tmp = facet_1;
                facet_1 = facet_2;
                facet_2 = facet_tmp;
            }
            if (vertex_1->next(facet_1)->next(facet_1) == vertex_2 ||
                    vertex_1->next(facet_2)->next(facet_2) == vertex_2 ||
                    vertex_1->prev(facet_1)->prev(facet_1) == vertex_2 ||
                    vertex_1->prev(facet_2)->prev(facet_2) == vertex_2) {
                // edge merge event
                continue;
            }

            EdgeSPtr edge_11 = EdgeSPtr();
            EdgeSPtr edge_12 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
            while (it_e1 != vertex_1->edges().end()) {
                EdgeWPtr edge_1_wptr = *it_e1++;
                if (!edge_1_wptr.expired()) {
                    EdgeSPtr edge_1 = EdgeSPtr(edge_1_wptr);
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
            std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
            while (it_e2 != vertex_2->edges().end()) {
                EdgeWPtr edge_2_wptr = *it_e2++;
                if (!edge_2_wptr.expired()) {
                    EdgeSPtr edge_2 = EdgeSPtr(edge_2_wptr);
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
            bool conv_split_event = false;
            FacetSPtr facet_1b = facet_2->next(vertex_1);
            FacetSPtr facet_2b = facet_1->next(vertex_2);
            if (facet_2b == facet_2) {
                // flip vertex event
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
                continue;
            }

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion (offset_event > current_offset_to_nearest_event);

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = SplitMergeEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_2->getData());
            node->addArc(data_2->getArc());
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setVertex1(vertex_1);
            result->setVertex2(vertex_2);
            result->setFacet1(facet_1);
            result->setFacet2(facet_2);
            current_offset_to_nearest_event = offset_event;
        }
    }
    return result;
}

EdgeSplitEventSPtr SimpleStraightSkel::nextEdgeSplitEvent(PolyhedronSPtr polyhedron,
                                                          const CGAL::FT current_offset,
                                                          CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Edge Split Event" << std::endl;
    CGAL::Real_timer timer;
    timer.start();

    ReadLock l(polyhedron->mutex());
    EdgeSplitEventSPtr result = EdgeSplitEventSPtr();

    std::list<EdgeSPtr> edges_reflex;
    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (isReflex(edge)) {
            edges_reflex.push_back(edge);
        }
    }

    std::list<EdgeSPtr>::iterator it_e1 = edges_reflex.begin();
    while (it_e1 != edges_reflex.end()) {
        EdgeSPtr edge_1 = *it_e1++;
        FacetSPtr facet_1_src = getFacetSrc(edge_1);
        FacetSPtr facet_1_dst = getFacetDst(edge_1);

        std::list<EdgeSPtr>::iterator it_e2 = it_e1;
        while (it_e2 != edges_reflex.end()) {
            EdgeSPtr edge_2 = *it_e2++;
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
            FacetSPtr facet_2_src = getFacetSrc(edge_2);
            FacetSPtr facet_2_dst = getFacetDst(edge_2);
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

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion (offset_event > current_offset_to_nearest_event);

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = EdgeSplitEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setEdge1(edge_1);
            result->setEdge2(edge_2);

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

            current_offset_to_nearest_event = offset_event;
        }
    }

    timer.stop();
    std::cout << "Sought Edge Split Event in: " << timer.time() << std::endl;

    return result;
}

PierceEventSPtr SimpleStraightSkel::nextPierceEvent(PolyhedronSPtr polyhedron,
                                                    const CGAL::FT current_offset,
                                                    CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Pierce Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    PierceEventSPtr result = PierceEventSPtr();

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;

        if (isReflex(vertex)) {
            SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
            ArcSPtr arc = data->getArc();

            std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
            while (it_f != polyhedron->facets().end()) {
                FacetSPtr facet = *it_f++;
                bool contains_vertex = false;
                std::list<VertexSPtr>::iterator it_v2 = facet->vertices().begin();
                while (it_v2 != facet->vertices().end()) {
                    VertexSPtr vertex_2 = *it_v2++;
                    if (vertex_2->getPoint() == vertex->getPoint()) {
                        contains_vertex = true;
                        break;
                    }
                }
                if (contains_vertex) {
                    continue;
                }

                bool has_edge_to_facet = false;
                std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
                while (it_e != vertex->edges().end()) {
                    EdgeWPtr edge_wptr = *it_e++;
                    if (!edge_wptr.expired()) {
                        EdgeSPtr edge(edge_wptr);
                        // @todo checking both because we don't know if vertex is the src or the dst of the edge?
                        FacetSPtr facet_src = edge->getFacetL()->next(edge->getVertexSrc());
                        FacetSPtr facet_dst = edge->getFacetR()->next(edge->getVertexDst());
                        if (facet == facet_src || facet == facet_dst) {
                            has_edge_to_facet = true;
                            break;
                        }
                    }
                }
                if (has_edge_to_facet) {
                    continue;
                }

                // @fixme assumes positive weights
                if (KernelWrapper::side(facet->plane(), vertex->getPoint()) > 0) {
                    continue;
                }

                Point3SPtr point;
                CGAL::FT offset_event;

// #define CGAL_SS3_OLD_CODE_PIERCE_EVENT
#ifdef CGAL_SS3_OLD_CODE_PIERCE_EVENT
                // @fixme this filter seems overly restrictive?
                // Could not we have an intersection between the OFFSET face
                // and the arc's supporting line (but we don't know at what offset yet)?
                if (!IsLineInFacet(facet, arc->line())) {
                    continue;
                }

                FacetSPtr facet_vertex = FacetSPtr(vertex->facets().front());
                CGAL::FT facet_speed_vertex = 1.0;
                if (facet_vertex->hasData()) {
                    facet_speed_vertex = std::dynamic_pointer_cast<SkelFacetData>(
                            facet_vertex->getData())->getSpeed();
                }
                Plane3SPtr plane_vertex_offset = KernelWrapper::offsetPlane(facet_vertex->plane(), -facet_speed_vertex);
                Point3SPtr point_vertex_offset = KernelWrapper::intersection(plane_vertex_offset, arc->line());
                CGAL::FT speed_vertex = KernelWrapper::distance(vertex->getPoint(), point_vertex_offset);

                Point3SPtr point_facet = KernelWrapper::intersection(facet->plane(), arc->line());
                CGAL::FT facet_speed = 1.0;
                if (facet->hasData()) {
                    facet_speed = std::dynamic_pointer_cast<SkelFacetData>(
                            facet->getData())->getSpeed();
                }
                Plane3SPtr plane_facet_offset = KernelWrapper::offsetPlane(facet->plane(), -facet_speed);
                Point3SPtr point_facet_offset = KernelWrapper::intersection(plane_facet_offset, arc->line());
                CGAL::FT speed_facet = KernelWrapper::distance(point_facet, point_facet_offset);

                CGAL::FT distance = KernelWrapper::distance(vertex->getPoint(), point_facet);
                CGAL::FT dist_vertex = (distance * speed_vertex) / (speed_vertex + speed_facet);
                if (KernelWrapper::comparePoints(arc->getDirection(),
                        point_facet, point_facet_offset) < 0) {
                    // for weighted straight skeleton
                    // reflex vertex and facet move into same direction
                    if (speed_facet < speed_vertex) {
                        // facet to slow
                        continue;
                    }
                    dist_vertex = (distance * speed_vertex) / (speed_facet - speed_vertex);
                }

                offset_event = -dist_vertex / speed_vertex;
                point = KernelWrapper::offsetPoint(vertex->getPoint(), arc->getDirection(), dist_vertex);
                // std::cout << "old result: " << *point << " time: " << offset_event << std::endl;
#else
                CGAL_assertion(vertex->facets().size() == 3);

                FacetWPtr wf0 = *(std::next(vertex->facets().begin(), 0)); // @todo can be factorized
                FacetWPtr wf1 = *(std::next(vertex->facets().begin(), 1));
                FacetWPtr wf2 = *(std::next(vertex->facets().begin(), 2));

                CGAL_assertion(!wf0.expired());
                CGAL_assertion(!wf1.expired());
                CGAL_assertion(!wf2.expired());

                FacetSPtr f0(wf0);
                FacetSPtr f1(wf1);
                FacetSPtr f2(wf2);

                Plane3SPtr plane = facet->plane();
                Plane3SPtr plane_0 = f0->plane();
                Plane3SPtr plane_1 = f1->plane();
                Plane3SPtr plane_2 = f2->plane();

                CGAL::FT speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
                CGAL::FT speed_0 = std::dynamic_pointer_cast<SkelFacetData>(f0->getData())->getSpeed();
                CGAL::FT speed_1 = std::dynamic_pointer_cast<SkelFacetData>(f1->getData())->getSpeed();
                CGAL::FT speed_2 = std::dynamic_pointer_cast<SkelFacetData>(f2->getData())->getSpeed();

                std::tie(point, offset_event) = KernelWrapper::intersectionAndTimeOffsetPlanes(
                    plane, speed, plane_0, speed_0, plane_1, speed_1, plane_2, speed_2);

                // std::cout << " **" << std::endl;
                // std::cout << "facet = " << facet->getID() << std::endl;
                // std::cout << "plane = " << *plane << std::endl;
                // std::cout << "f0 = " << f0->getID() << std::endl;
                // std::cout << "f1 = " << f1->getID() << std::endl;
                // std::cout << "f2 = " << f2->getID() << std::endl;
                // std::cout << "ipoint: " << *point << " time: " << offset_event << std::endl;

                // Now, we might naively wish to filter using bisectors like in 2D SLS code,
                // but unlike a segment, a face in the 3D SS code has no reason to be convex,
                // which changes everything and can result in false positives.
                //
                // The bisector filter in 2D is equivalent to checking if the point is on the offset
                // face. We could check this here, but determining what is the offset face at this
                // point (i.e., while searching for events) is rough: plenty of other events
                // might modify the face before this particular pierce event appears, and so
                // we can't just do shift(facet) because the result might be a self-intersecting
                // polygon with holes.
                //
                // Instead, we do not filter here, but simply put it in the queue. When the event
                // will be popped, then we know it's the next event globally and nothing else
                // can mess up the face, and we can do the in-test check then.
                //
                // See HandlePierceEvent()

#endif // CGAL_SS3_OLD_CODE_PIERCE_EVENT
                std::cout << "Accepted event" << std::endl;

                if (offset_event == 0)
                {
                  // @fixme this should be <= like the other nextSomethingEvent
                  // but since we handle filtering at pop time, it will just loop indefinitely.
                  CGAL_assertion_msg(false, "this should have been accepted!!");
                }

                if (offset_event <= 0 && offset_event > current_offset_to_nearest_event) {
                    NodeSPtr node;
                    if (!result) {
                        node = Node::create(point);
                        result = PierceEvent::create(polyhedron);
                        result->setNode(node);
                    }
                    node = result->getNode();
                    node->clear();
                    node->addArc(arc);
                    node->setOffset(current_offset + offset_event);
                    node->setPoint(point);
                    result->setFacet(facet);
                    result->setVertex(vertex);
                    current_offset_to_nearest_event = offset_event;
                }
            }
        }
    }
    return result;
}
#endif // individual events of a certain type

void SimpleStraightSkel::collectEdgeEvents(PolyhedronSPtr polyhedron,
                                           const CGAL::FT current_offset,
                                           CGAL::FT& current_offset_to_nearest_event,
                                           PQ& queue)
{
    std::cout << ">>> Collect Next Edge Events" << std::endl;

    ReadLock l(polyhedron->mutex());

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;

        VertexSPtr vertex_src = edge->getVertexSrc();
        VertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->getPoint() == vertex_dst->getPoint()) {
            continue;
        }

        if (isLocked(edge)) { // one of the vertices cannot be rebuilt from its incident faces
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
        std::tie(point, offset_event) = vanishesAt(edge);
        if (!point) {
            continue;
        }

        FacetSPtr facet_src = getFacetSrc(edge);
        FacetSPtr facet_dst = getFacetDst(edge);

        // This does not work when there is more than one edge between both facets.
        // EdgeSPtr edge_2 = facet_src->findEdge(facet_dst);
        std::list<EdgeSPtr> edges_2 = facet_src->findEdges(facet_dst);

        bool split_event = false;
        std::list<EdgeSPtr>::iterator it_e2 = edges_2.begin();
        while (it_e2 != edges_2.end()) {
            EdgeSPtr edge_2 = *it_e2++;

            if (isLocked(edge_2)) {
                continue;
            }

            // @todo exit as soon as there is a single rejection (to be done once the old code is safe to remove)
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
            FacetSPtr facet_2_src = getFacetSrc(edge_2);
            FacetSPtr facet_2_dst = getFacetDst(edge_2);

            Plane3SPtr plane_l2 = facet_l2->plane();
            Plane3SPtr plane_r2 = facet_r2->plane();
            CGAL::FT speed_l2 = std::dynamic_pointer_cast<SkelFacetData>(facet_l2->getData())->getSpeed();
            CGAL::FT speed_r2 = std::dynamic_pointer_cast<SkelFacetData>(facet_r2->getData())->getSpeed();

            CGAL::FT l2a = plane_l2->a();
            CGAL::FT l2b = plane_l2->b();
            CGAL::FT l2c = plane_l2->c();
            CGAL::FT l2d = plane_l2->d();
            CGAL::FT r2a = plane_r2->a();
            CGAL::FT r2b = plane_r2->b();
            CGAL::FT r2c = plane_r2->c();
            CGAL::FT r2d = plane_r2->d();

            CGAL::FT lt2 = (l2a * point->x() + l2b * point->y() + l2c * point->z() + l2d) / speed_l2;
            CGAL::FT rt2 = (r2a * point->x() + r2b * point->y() + r2c * point->z() + r2d) / speed_r2;
            CGAL_assertion(lt2 == rt2);

            if ((lt2 > 0) || (rt2 > 0)) {
                split_event_current_1_b = false;
            }

            if (!check_bisector(edge_2, facet_r2, rt2, facet_2_src, point)) {
                split_event_current_2_b = false;
            }

            if (!check_bisector(edge_2, facet_l2, lt2, facet_2_dst, point)) {
                split_event_current_3_b = false;
            }

            const bool split_event_current_b = (split_event_current_1_b &&
                                                split_event_current_2_b &&
                                                split_event_current_3_b);
#endif

#ifdef CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS
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
        // @todo could move up this check
        if (offset_event <= 0 && offset_event >= current_offset_to_nearest_event)
        {
            NodeSPtr node = Node::create(point);
            EdgeEventSPtr event = EdgeEvent::create(polyhedron);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setEdge(edge);

#ifndef CGAL_SS3_NO_ARC_NO_SHEET
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

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
}

void SimpleStraightSkel::collectEdgeMergeEvents(PolyhedronSPtr polyhedron,
                                                const CGAL::FT current_offset,
                                                CGAL::FT& current_offset_to_nearest_event,
                                                PQ& queue)
{
    std::cout << ">>> Collect Edge Merge Events" << std::endl;

    ReadLock l(polyhedron->mutex());

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;

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

        // @todo do we need to test these other combinations if a previous one matched?
        EdgeSPtr edge_prev = edge->prev(facet_l);
        edge_next = edge->next(facet_l)->next(facet_l);
        if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
            facet = facet_l;
            edge_1 = edge_prev;
            edge_2 = edge_next;
        }

#ifndef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
        // @fixme not sure about that macro:
        // this is essentially the same as above but if 'edge' were 'edge->prev(facet_l)'
        // so when we pass by again wit this edge, then we will meet it in the check just above.
        // BUT: the dbl edge merge event filter above applies only to 'edge' so we might filter
        // for edge->prev(facet_l) which we wouldn't have filtered when seeing it from 'edge'?

        edge_prev = edge->prev(facet_l)->prev(facet_l);
        edge_next = edge->next(facet_l);
        if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
            facet = facet_l;
            edge_1 = edge_prev;
            edge_2 = edge_next;
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
        // @fixme not sure about that macro: see above
        edge_prev = edge->prev(facet_r)->prev(facet_r);
        edge_next = edge->next(facet_r);
        if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
            facet = facet_r;
            edge_1 = edge_prev;
            edge_2 = edge_next;
        }
#endif

        if (!(facet && edge_1 && edge_2)) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge);
        if (!point) {
            continue;
        }
        if (offset_event <= 0 && offset_event >= current_offset_to_nearest_event) {
            EdgeMergeEventSPtr event = EdgeMergeEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setFacet(facet);
            event->setEdge1(edge_1);
            event->setEdge2(edge_2);

#ifndef CGAL_SS3_NO_ARC_NO_SHEET
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

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
}

void SimpleStraightSkel::collectTriangleEvents(PolyhedronSPtr polyhedron,
                                               const CGAL::FT current_offset,
                                               CGAL::FT& current_offset_to_nearest_event,
                                               PQ& queue)
{
    std::cout << ">>> Collect Triangle Event" << std::endl;

    ReadLock l(polyhedron->mutex());

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;

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
        // for canonicity, only investigate this event if we are in the smallest edge of the face
        if (edge->getID() > edge->prev(facet)->getID() ||
            edge->getID() > edge->next(facet)->getID()) {
            continue;
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
        std::tie(point, offset_event) = vanishesAt(edge);
        if (!point) {
            continue;
        }

        if ((KernelWrapper::side(edge->getFacetL()->plane(), point) > 0) ||
                KernelWrapper::side(edge->getFacetR()->plane(), point) > 0) {
            // triangle may not be a hole
            // after pierce event
            continue;
        }

        if (offset_event <= 0 && offset_event >= current_offset_to_nearest_event) {
            TriangleEventSPtr event = TriangleEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setFacet(facet);
            event->setEdgeBegin(edge);

#ifndef CGAL_SS3_NO_ARC_NO_SHEET
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

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
}

void SimpleStraightSkel::collectDblEdgeMergeEvents(PolyhedronSPtr polyhedron,
                                                   const CGAL::FT current_offset,
                                                   CGAL::FT& current_offset_to_nearest_event,
                                                   PQ& queue)
{
    std::cout << ">>> Collect Dbl Edge Merge Events" << std::endl;

    ReadLock l(polyhedron->mutex());

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;

        if (isLocked(edge)) {
            continue;
        }

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
        std::tie(point, offset_event) = vanishesAt(edge);
        if (!point) {
            continue;
        }
        if (offset_event <= 0 && offset_event >= current_offset_to_nearest_event) {
            DblEdgeMergeEventSPtr event = DblEdgeMergeEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setFacet1(facet_1);
            event->setEdge11(edge_11);
            event->setEdge12(edge_12);
            event->setFacet2(facet_2);
            event->setEdge21(edge_21);
            event->setEdge22(edge_22);

#ifndef CGAL_SS3_NO_ARC_NO_SHEET
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

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
}

void SimpleStraightSkel::collectDblTriangleEvents(PolyhedronSPtr polyhedron,
                                                  const CGAL::FT current_offset,
                                                  CGAL::FT& current_offset_to_nearest_event,
                                                  PQ& queue)
{
    std::cout << ">>> Collect Dbl Triangle Events" << std::endl;

    ReadLock l(polyhedron->mutex());

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (isTetrahedron(edge)) {
            continue;
        }
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (!facet_l || !facet_r) {
            continue;
        }
        if (!(isTriangle(facet_l, edge) &&
                isTriangle(facet_r, edge))) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge);
        if (!point) {
            continue;
        }
        if (offset_event <= 0 && offset_event >= current_offset_to_nearest_event) {
            DblTriangleEventSPtr event = DblTriangleEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setEdge(edge);

#ifndef CGAL_SS3_NO_ARC_NO_SHEET
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

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
}

void SimpleStraightSkel::collectTetrahedronEvents(PolyhedronSPtr polyhedron,
                                                  const CGAL::FT current_offset,
                                                  CGAL::FT& current_offset_to_nearest_event,
                                                  PQ& queue)
{
    std::cout << ">>> Collect Tetrahedron Events" << std::endl;

    ReadLock l(polyhedron->mutex());

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (isTetrahedron(edge)) {

#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
            // for canonicity, only investigate this event if we are in the smallest edge
            // of the smallest face of the tetrahedron
            FacetSPtr facet = (edge->getFacetL()->getID() < edge->getFacetR()->getID()) ? edge->getFacetL()
                                                                                        : edge->getFacetR();

            // smallest incident face is the smallest of the tetrahedron
            if (facet->getID() > edge->prev(facet)->other(facet)->getID() ||
                facet->getID() > edge->other(facet)->getID() ||
                facet->getID() > edge->next(facet)->other(facet)->getID()) {
                continue;
            }

            // edge is the smallest within the face
            if (edge->getID() > edge->next(facet)->getID() || edge->getID() > edge->prev(facet)->getID()) {
                continue;
            }
#endif

            Point3SPtr point = Point3SPtr();
            CGAL::FT offset_event;
            std::tie(point, offset_event) = vanishesAt(edge);
            if (!point) {
                continue;
            }

            std::cout << "potential tet event @ " << *point << " (t=" << offset_event << ")" << std::endl;

            if (offset_event <= 0 && offset_event >= current_offset_to_nearest_event) {
                TetrahedronEventSPtr event = TetrahedronEvent::create(polyhedron);
                NodeSPtr node = Node::create(point);
                event->setNode(node);
                node->clear();
                node->setOffset(current_offset + offset_event);
                node->setPoint(point);
                event->setEdgeBegin(edge);

#ifndef CGAL_SS3_NO_ARC_NO_SHEET
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

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
                current_offset_to_nearest_event = offset_event;
#endif
            }
        }
    }
}

void SimpleStraightSkel::collectVertexEvents(PolyhedronSPtr polyhedron,
                                             const CGAL::FT current_offset,
                                             CGAL::FT& current_offset_to_nearest_event,
                                             PQ& queue)
{
    std::cout << ">>> Collect Vertex Events" << std::endl;

    ReadLock l(polyhedron->mutex());

    std::list<VertexSPtr>::iterator it_v1 = polyhedron->vertices().begin();
    while (it_v1 != polyhedron->vertices().end()) {
        VertexSPtr vertex_1 = *it_v1++;
        if (isConvex(vertex_1)) {
            continue;
        }

        std::set<VertexSPtr> vertices_2;
        std::list<FacetWPtr>::iterator it_f = vertex_1->facets().begin();
        while (it_f != vertex_1->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                vertices_2.insert(facet->vertices().begin(), facet->vertices().end());
            }
        }
        std::set<VertexSPtr>::iterator it_v2 = vertices_2.begin();
        while (it_v2 != vertices_2.end()) {
            VertexSPtr vertex_2 = *it_v2++;
            if (vertex_1 == vertex_2) {
                continue;
            }
#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
            if (vertex_1->getID() > vertex_2->getID()) {
                continue;
            }
#endif
            if (vertex_1->getPoint() == vertex_2->getPoint()) {
                continue;
            }
            if (isConvex(vertex_2)) {
                continue;
            }

            if (vertex_1->findEdge(vertex_2)) {
                // edge event
                continue;
            }

            FacetSPtr facet_1;
            FacetSPtr facet_2;
            int num_equal_facets = 0;
            std::list<FacetWPtr>::iterator it_f1 = vertex_1->facets().begin();
            while (it_f1 != vertex_1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = vertex_2->facets().begin();
                    while (it_f2 != vertex_2->facets().end()) {
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
            if (facet_1->next(vertex_1) != facet_2) {
                FacetSPtr facet_tmp = facet_1;
                facet_1 = facet_2;
                facet_2 = facet_tmp;
            }
            if (vertex_1->next(facet_1)->next(facet_1) == vertex_2 ||
                    vertex_1->next(facet_2)->next(facet_2) == vertex_2 ||
                    vertex_1->prev(facet_1)->prev(facet_1) == vertex_2 ||
                    vertex_1->prev(facet_2)->prev(facet_2) == vertex_2) {
                // edge merge event
                continue;
            }

            EdgeSPtr edge_11 = EdgeSPtr();
            EdgeSPtr edge_12 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
            while (it_e1 != vertex_1->edges().end()) {
                EdgeWPtr edge_1_wptr = *it_e1++;
                if (!edge_1_wptr.expired()) {
                    EdgeSPtr edge_1 = EdgeSPtr(edge_1_wptr);
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
            std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
            while (it_e2 != vertex_2->edges().end()) {
                EdgeWPtr edge_2_wptr = *it_e2++;
                if (!edge_2_wptr.expired()) {
                    EdgeSPtr edge_2 = EdgeSPtr(edge_2_wptr);
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
            if (!((edge_11->next(vertex_1) == edge_12 && edge_22->next(vertex_2) == edge_21) ||
                    (edge_12->next(vertex_1) == edge_11 && edge_21->next(vertex_2) == edge_22))) {
                // flip vertex event
                continue;
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
                continue;
            }

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion(offset_event >= current_offset_to_nearest_event);

            VertexEventSPtr event = VertexEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();

            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_2->getData());
            node->addArc(data_2->getArc());

            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setVertex1(vertex_1);
            event->setVertex2(vertex_2);
            event->setFacet1(facet_1);
            event->setFacet2(facet_2);

            queue.push(event);

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
}

void SimpleStraightSkel::collectFlipVertexEvents(PolyhedronSPtr polyhedron,
                                                 const CGAL::FT current_offset,
                                                 CGAL::FT& current_offset_to_nearest_event,
                                                 PQ& queue)
{
    std::cout << ">>> Collect Flip Vertex Events" << std::endl;

    ReadLock l(polyhedron->mutex());

    std::list<VertexSPtr>::iterator it_v1 = polyhedron->vertices().begin();
    while (it_v1 != polyhedron->vertices().end()) {
        VertexSPtr vertex_1 = *it_v1++;
        if (isConvex(vertex_1)) {
            continue;
        }

        std::set<VertexSPtr> vertices_2;
        std::list<FacetWPtr>::iterator it_f = vertex_1->facets().begin();
        while (it_f != vertex_1->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                vertices_2.insert(facet->vertices().begin(), facet->vertices().end());
            }
        }
        std::set<VertexSPtr>::iterator it_v2 = vertices_2.begin();
        while (it_v2 != vertices_2.end()) {
            VertexSPtr vertex_2 = *it_v2++;
            if (vertex_1 == vertex_2) {
                continue;
            }
#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
            if (vertex_1->getID() > vertex_2->getID()) {
                continue;
            }
#endif
            if (vertex_1->getPoint() == vertex_2->getPoint()) {
                continue;
            }
            if (isConvex(vertex_2)) {
                continue;
            }

            if (vertex_1->findEdge(vertex_2)) {
                // edge event
                continue;
            }

            FacetSPtr facet_1;
            FacetSPtr facet_2;
            int num_equal_facets = 0;
            std::list<FacetWPtr>::iterator it_f1 = vertex_1->facets().begin();
            while (it_f1 != vertex_1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = vertex_2->facets().begin();
                    while (it_f2 != vertex_2->facets().end()) {
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
            if (facet_1->getID() > facet_2->getID()) {
                continue;
            }
#endif
            if (facet_1->next(vertex_1) != facet_2) {
                FacetSPtr facet_tmp = facet_1;
                facet_1 = facet_2;
                facet_2 = facet_tmp;
            }
            if (vertex_1->next(facet_1)->next(facet_1) == vertex_2 ||
                    vertex_1->next(facet_2)->next(facet_2) == vertex_2 ||
                    vertex_1->prev(facet_1)->prev(facet_1) == vertex_2 ||
                    vertex_1->prev(facet_2)->prev(facet_2) == vertex_2) {
                // edge merge event
                continue;
            }

            EdgeSPtr edge_11 = EdgeSPtr();
            EdgeSPtr edge_12 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
            while (it_e1 != vertex_1->edges().end()) {
                EdgeWPtr edge_1_wptr = *it_e1++;
                if (!edge_1_wptr.expired()) {
                    EdgeSPtr edge_1 = EdgeSPtr(edge_1_wptr);
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
            std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
            while (it_e2 != vertex_2->edges().end()) {
                EdgeWPtr edge_2_wptr = *it_e2++;
                if (!edge_2_wptr.expired()) {
                    EdgeSPtr edge_2 = EdgeSPtr(edge_2_wptr);
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
            if (!(edge_12->next(vertex_1) == edge_11 && edge_22->next(vertex_2) == edge_21)) {
                // vertex event
                continue;
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
                continue;
            }

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion(offset_event >= current_offset_to_nearest_event);

            FlipVertexEventSPtr event = FlipVertexEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();

            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_2->getData());
            node->addArc(data_2->getArc());

            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setVertex1(vertex_1);
            event->setVertex2(vertex_2);
            event->setFacet1(facet_1);
            event->setFacet2(facet_2);

            queue.push(event);

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
}

void SimpleStraightSkel::collectSurfaceEvents(PolyhedronSPtr polyhedron,
                                              const CGAL::FT current_offset,
                                              CGAL::FT& current_offset_to_nearest_event,
                                              PQ& queue)
{
    std::cout << ">>> Collect Surface Events" << std::endl;
    CGAL::Real_timer timer;
    timer.start();

    ReadLock l(polyhedron->mutex());

    std::list<EdgeSPtr>::iterator it_e1 = polyhedron->edges().begin();
    while (it_e1 != polyhedron->edges().end()) {
        EdgeSPtr edge_1 = *it_e1++;
        FacetSPtr facet_1_src = getFacetSrc(edge_1);
        FacetSPtr facet_1_dst = getFacetDst(edge_1);
        std::list<EdgeSPtr> edges_2;
        edges_2.insert(edges_2.end(),
                facet_1_src->edges().begin(), facet_1_src->edges().end());
        edges_2.insert(edges_2.end(),
                facet_1_dst->edges().begin(), facet_1_dst->edges().end());
        std::list<EdgeSPtr>::iterator it_e2 = edges_2.begin();
        while (it_e2 != edges_2.end()) {
            EdgeSPtr edge_2 = *it_e2++;
            if (edge_1 == edge_2) {
                continue;
            }
// #ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
//             // @fixme? due to the convex check, it's not symmetric, but maybe it can be sometimes?
//             if (edge_1->getID() > edge_2->getID()) {
//                 continue;
//             }
// #endif
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
            FacetSPtr facet_2_src = getFacetSrc(edge_2);
            FacetSPtr facet_2_dst = getFacetDst(edge_2);
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
            bool is_conv_split_event = false;
            std::list<EdgeSPtr> edges = edge_1->getFacetL()->findEdges(edge_1->getFacetR());
            std::list<EdgeSPtr>::iterator it_e = edges.begin();
            while (it_e != edges.end()) {
                EdgeSPtr edge = *it_e++;
                if (edge == edge_1) {
                    continue;
                }
                FacetSPtr facet_src = getFacetSrc(edge);
                FacetSPtr facet_dst = getFacetDst(edge);
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

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion(offset_event >= current_offset_to_nearest_event);

            SurfaceEventSPtr event = SurfaceEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setEdge1(edge_1);
            event->setEdge2(edge_2);

#ifndef CGAL_SS3_NO_ARC_NO_SHEET
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

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }

    timer.stop();
    std::cout << "Sought Surface Events in: " << timer.time() << std::endl;
}

void SimpleStraightSkel::collectPolyhedronSplitEvents(PolyhedronSPtr polyhedron,
                                                      const CGAL::FT current_offset,
                                                      CGAL::FT& current_offset_to_nearest_event,
                                                      PQ& queue)
{
    std::cout << ">>> Collect Polyhedron Split Events" << std::endl;

    ReadLock l(polyhedron->mutex());

    std::list<EdgeSPtr>::iterator it_e1 = polyhedron->edges().begin();
    while (it_e1 != polyhedron->edges().end()) {
        EdgeSPtr edge_1 = *it_e1++;
        if (!isReflex(edge_1)) {
            continue;
        }

        FacetSPtr facet_1_src = getFacetSrc(edge_1);
        FacetSPtr facet_1_dst = getFacetDst(edge_1);
        std::list<EdgeSPtr>::iterator it_e2 = facet_1_src->edges().begin();
        while (it_e2 != facet_1_src->edges().end()) {
            EdgeSPtr edge_2 = *it_e2++;
            if (edge_1->getVertexSrc()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
                    edge_1->getVertexSrc()->getPoint() == edge_2->getVertexDst()->getPoint() ||
                    edge_1->getVertexDst()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
                    edge_1->getVertexDst()->getPoint() == edge_2->getVertexDst()->getPoint()) {
                // share a vertex
                continue;
            }
            if (!((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() == facet_1_dst) ||
                    (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() == facet_1_src))) {
                // no polyhedron split event
                continue;
            }
            if (edge_1->getVertexSrc()->findEdge(edge_2->getVertexSrc()) ||
                    edge_1->getVertexSrc()->findEdge(edge_2->getVertexDst()) ||
                    edge_1->getVertexDst()->findEdge(edge_2->getVertexSrc()) ||
                    edge_1->getVertexDst()->findEdge(edge_2->getVertexDst())) {
                // does not work when there is only one edge in between
                continue;
            }

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion (offset_event > current_offset_to_nearest_event);

            PolyhedronSplitEventSPtr event = PolyhedronSplitEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setEdge1(edge_1);
            event->setEdge2(edge_2);

#ifndef CGAL_SS3_NO_ARC_NO_SHEET
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

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
}

void SimpleStraightSkel::collectSplitMergeEvents(PolyhedronSPtr polyhedron,
                                                 const CGAL::FT current_offset,
                                                 CGAL::FT& current_offset_to_nearest_event,
                                                 PQ& queue)
{
    std::cout << ">>> Collect Split Merge Events" << std::endl;

    ReadLock l(polyhedron->mutex());

    std::list<VertexSPtr>::iterator it_v1 = polyhedron->vertices().begin();
    while (it_v1 != polyhedron->vertices().end()) {
        VertexSPtr vertex_1 = *it_v1++;
        if (isConvex(vertex_1)) {
            continue;
        }

        std::set<VertexSPtr> vertices_2;
        std::list<FacetWPtr>::iterator it_f = vertex_1->facets().begin();
        while (it_f != vertex_1->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                vertices_2.insert(facet->vertices().begin(), facet->vertices().end());
            }
        }
        std::set<VertexSPtr>::iterator it_v2 = vertices_2.begin();
        while (it_v2 != vertices_2.end()) {
            VertexSPtr vertex_2 = *it_v2++;
            if (vertex_1 == vertex_2) {
                continue;
            }
#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
            if (vertex_1->getID() > vertex_2->getID()) {
                continue;
            }
#endif
            if (vertex_1->getPoint() == vertex_2->getPoint()) {
                continue;
            }
            if (isConvex(vertex_2)) {
                continue;
            }

            if (vertex_1->findEdge(vertex_2)) {
                // edge event
                continue;
            }

            FacetSPtr facet_1;
            FacetSPtr facet_2;
            int num_equal_facets = 0;
            std::list<FacetWPtr>::iterator it_f1 = vertex_1->facets().begin();
            while (it_f1 != vertex_1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = vertex_2->facets().begin();
                    while (it_f2 != vertex_2->facets().end()) {
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
            if (facet_1->next(vertex_1) != facet_2) {
                FacetSPtr facet_tmp = facet_1;
                facet_1 = facet_2;
                facet_2 = facet_tmp;
            }
            if (vertex_1->next(facet_1)->next(facet_1) == vertex_2 ||
                    vertex_1->next(facet_2)->next(facet_2) == vertex_2 ||
                    vertex_1->prev(facet_1)->prev(facet_1) == vertex_2 ||
                    vertex_1->prev(facet_2)->prev(facet_2) == vertex_2) {
                // edge merge event
                continue;
            }

            EdgeSPtr edge_11 = EdgeSPtr();
            EdgeSPtr edge_12 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
            while (it_e1 != vertex_1->edges().end()) {
                EdgeWPtr edge_1_wptr = *it_e1++;
                if (!edge_1_wptr.expired()) {
                    EdgeSPtr edge_1 = EdgeSPtr(edge_1_wptr);
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
            std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
            while (it_e2 != vertex_2->edges().end()) {
                EdgeWPtr edge_2_wptr = *it_e2++;
                if (!edge_2_wptr.expired()) {
                    EdgeSPtr edge_2 = EdgeSPtr(edge_2_wptr);
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
            bool conv_split_event = false;
            FacetSPtr facet_1b = facet_2->next(vertex_1);
            FacetSPtr facet_2b = facet_1->next(vertex_2);
            if (facet_2b == facet_2) {
                // flip vertex event
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
                continue;
            }

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion (offset_event > current_offset_to_nearest_event);

            SplitMergeEventSPtr event = SplitMergeEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();

            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_2->getData());
            node->addArc(data_2->getArc());

            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setVertex1(vertex_1);
            event->setVertex2(vertex_2);
            event->setFacet1(facet_1);
            event->setFacet2(facet_2);

            queue.push(event);

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
}

void SimpleStraightSkel::collectEdgeSplitEvents(PolyhedronSPtr polyhedron,
                                                const CGAL::FT current_offset,
                                                CGAL::FT& current_offset_to_nearest_event,
                                                PQ& queue)
{
    std::cout << ">>> Collect Edge Split Events" << std::endl;

    CGAL::Real_timer timer;
    timer.start();

    ReadLock l(polyhedron->mutex());

    std::list<EdgeSPtr> edges_reflex;
    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (isReflex(edge)) {
            edges_reflex.push_back(edge);
        }
    }

    std::cout << edges_reflex.size() << " reflex edges" << std::endl;

    std::list<EdgeSPtr>::iterator it_e1 = edges_reflex.begin();
    while (it_e1 != edges_reflex.end()) {
        EdgeSPtr edge_1 = *it_e1++;
        FacetSPtr facet_1_src = getFacetSrc(edge_1);
        FacetSPtr facet_1_dst = getFacetDst(edge_1);

        std::list<EdgeSPtr>::iterator it_e2 = it_e1;
        while (it_e2 != edges_reflex.end()) {
            EdgeSPtr edge_2 = *it_e2++;
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
            FacetSPtr facet_2_src = getFacetSrc(edge_2);
            FacetSPtr facet_2_dst = getFacetDst(edge_2);
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

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion (offset_event > current_offset_to_nearest_event);

            EdgeSplitEventSPtr event = EdgeSplitEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setEdge1(edge_1);
            event->setEdge2(edge_2);

#ifndef CGAL_SS3_NO_ARC_NO_SHEET
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

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }

    timer.stop();
    std::cout << "Sought Edge Split Events in: " << timer.time() << std::endl;
}

void SimpleStraightSkel::collectPierceEvents(PolyhedronSPtr polyhedron,
                                            const CGAL::FT current_offset,
                                            CGAL::FT& current_offset_to_nearest_event,
                                            PQ& queue)
{
    std::cout << ">>> Collect Pierce Events" << std::endl;

    ReadLock l(polyhedron->mutex());

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (isLocked(vertex)) {
            continue;
        }

        if (isReflex(vertex)) {
            SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
            ArcSPtr arc = data->getArc();

            std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
            while (it_f != polyhedron->facets().end()) {
                FacetSPtr facet = *it_f++;
                bool contains_vertex = false;
                std::list<VertexSPtr>::iterator it_v2 = facet->vertices().begin();
                while (it_v2 != facet->vertices().end()) {
                    VertexSPtr vertex_2 = *it_v2++;
                    if (vertex_2->getPoint() == vertex->getPoint()) {
                        contains_vertex = true;
                        break;
                    }
                }
                if (contains_vertex) {
                    continue;
                }

                bool has_edge_to_facet = false;
                std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
                while (it_e != vertex->edges().end()) {
                    EdgeWPtr edge_wptr = *it_e++;
                    if (!edge_wptr.expired()) {
                        EdgeSPtr edge(edge_wptr);
                        // @todo checking both because we don't know if vertex is the src or the dst of the edge?
                        FacetSPtr facet_src = edge->getFacetL()->next(edge->getVertexSrc());
                        FacetSPtr facet_dst = edge->getFacetR()->next(edge->getVertexDst());
                        if (facet == facet_src || facet == facet_dst) {
                            has_edge_to_facet = true;
                            break;
                        }
                    }
                }
                if (has_edge_to_facet) {
                    continue;
                }

                // @fixme assumes positive weights
                if (KernelWrapper::side(facet->plane(), vertex->getPoint()) > 0) {
                    continue;
                }

                Point3SPtr point;
                CGAL::FT offset_event;

// #define CGAL_SS3_OLD_CODE_PIERCE_EVENT
#ifdef CGAL_SS3_OLD_CODE_PIERCE_EVENT
                // @fixme this filter seems overly restrictive?
                // Could not we have an intersection between the OFFSET face
                // and the arc's supporting line (but we don't know at what offset yet)?
                if (!IsLineInFacet(facet, arc->line())) {
                    continue;
                }

                FacetSPtr facet_vertex = FacetSPtr(vertex->facets().front());
                CGAL::FT facet_speed_vertex = 1.0;
                if (facet_vertex->hasData()) {
                    facet_speed_vertex = std::dynamic_pointer_cast<SkelFacetData>(
                            facet_vertex->getData())->getSpeed();
                }
                Plane3SPtr plane_vertex_offset = KernelWrapper::offsetPlane(facet_vertex->plane(), -facet_speed_vertex);
                Point3SPtr point_vertex_offset = KernelWrapper::intersection(plane_vertex_offset, arc->line());
                CGAL::FT speed_vertex = KernelWrapper::distance(vertex->getPoint(), point_vertex_offset);

                Point3SPtr point_facet = KernelWrapper::intersection(facet->plane(), arc->line());
                CGAL::FT facet_speed = 1.0;
                if (facet->hasData()) {
                    facet_speed = std::dynamic_pointer_cast<SkelFacetData>(
                            facet->getData())->getSpeed();
                }
                Plane3SPtr plane_facet_offset = KernelWrapper::offsetPlane(facet->plane(), -facet_speed);
                Point3SPtr point_facet_offset = KernelWrapper::intersection(plane_facet_offset, arc->line());
                CGAL::FT speed_facet = KernelWrapper::distance(point_facet, point_facet_offset);

                CGAL::FT distance = KernelWrapper::distance(vertex->getPoint(), point_facet);
                CGAL::FT dist_vertex = (distance * speed_vertex) / (speed_vertex + speed_facet);
                if (KernelWrapper::comparePoints(arc->getDirection(),
                        point_facet, point_facet_offset) < 0) {
                    // for weighted straight skeleton
                    // reflex vertex and facet move into same direction
                    if (speed_facet < speed_vertex) {
                        // facet to slow
                        continue;
                    }
                    dist_vertex = (distance * speed_vertex) / (speed_facet - speed_vertex);
                }

                offset_event = -dist_vertex / speed_vertex;
                point = KernelWrapper::offsetPoint(vertex->getPoint(), arc->getDirection(), dist_vertex);
                // std::cout << "old result: " << *point << " time: " << offset_event << std::endl;
#else
                CGAL_assertion(vertex->facets().size() == 3);

                FacetWPtr wf0 = *(std::next(vertex->facets().begin(), 0)); // @todo can be factorized
                FacetWPtr wf1 = *(std::next(vertex->facets().begin(), 1));
                FacetWPtr wf2 = *(std::next(vertex->facets().begin(), 2));

                CGAL_assertion(!wf0.expired());
                CGAL_assertion(!wf1.expired());
                CGAL_assertion(!wf2.expired());

                FacetSPtr f0(wf0);
                FacetSPtr f1(wf1);
                FacetSPtr f2(wf2);

                Plane3SPtr plane = facet->plane();
                Plane3SPtr plane_0 = f0->plane();
                Plane3SPtr plane_1 = f1->plane();
                Plane3SPtr plane_2 = f2->plane();

                CGAL::FT speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
                CGAL::FT speed_0 = std::dynamic_pointer_cast<SkelFacetData>(f0->getData())->getSpeed();
                CGAL::FT speed_1 = std::dynamic_pointer_cast<SkelFacetData>(f1->getData())->getSpeed();
                CGAL::FT speed_2 = std::dynamic_pointer_cast<SkelFacetData>(f2->getData())->getSpeed();

                std::tie(point, offset_event) = KernelWrapper::intersectionAndTimeOffsetPlanes(
                    plane, speed, plane_0, speed_0, plane_1, speed_1, plane_2, speed_2);

                if (!point) {
                    continue;
                }

                // std::cout << " **" << std::endl;
                // std::cout << "facet = " << facet->getID() << std::endl;
                // std::cout << "plane = " << *plane << std::endl;
                // std::cout << "f0 = " << f0->getID() << std::endl;
                // std::cout << "f1 = " << f1->getID() << std::endl;
                // std::cout << "f2 = " << f2->getID() << std::endl;
                // std::cout << "ipoint: " << *point << " time: " << offset_event << std::endl;

                // Now, we might naively wish to filter using bisectors like in 2D SLS code,
                // but unlike a segment, a face in the 3D SS code has no reason to be convex,
                // which changes everything and can result in false positives.
                //
                // The bisector filter in 2D is equivalent to checking if the point is on the offset
                // face. We could check this here, but determining what is the offset face at this
                // point (i.e., while searching for events) is rough: plenty of other events
                // might modify the face before this particular pierce event appears, and so
                // we can't just do shift(facet) because the result might be a self-intersecting
                // polygon with holes.
                //
                // Instead, we do not filter here, but simply put it in the queue. When the event
                // will be popped, then we know it's the next event globally and nothing else
                // can mess up the face, and we can do the in-test check then.
                //
                // See HandlePierceEvent()
                //
                // One exception: if offset_event is 0 (meaning, pierce event at the current time),
                // then we can and should do the filtering at this particular point for two reasons:
                // - if we put it in the queue and filter at pop, we might get an endless loop
                // - it is safe: there can't be topological changes like described above

                if (offset_event == 0)
                {
                    if (!SelfIntersection::isInsideWithRayShooting(point, facet)) {
                        std::cout << "Filtere T=0 pierce event" << std::endl;
                        continue;
                    }
                }
#endif // CGAL_SS3_OLD_CODE_PIERCE_EVENT
                if (offset_event <= 0 && offset_event >= current_offset_to_nearest_event) {
                    PierceEventSPtr event = PierceEvent::create(polyhedron);
                    NodeSPtr node = Node::create(point);
                    event->setNode(node);
                    node->clear();
                    node->addArc(arc);
                    node->setOffset(current_offset + offset_event);
                    node->setPoint(point);
                    event->setFacet(facet);
                    event->setVertex(vertex);

                    queue.push(event);

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
                    current_offset_to_nearest_event = offset_event;
#endif
                }
            }
        }
    }
}

void SimpleStraightSkel::collectEvents(PolyhedronSPtr polyhedron,
                                       const CGAL::FT current_offset,
                                       PQ& queue)
{
    AbstractEventSPtr result = AbstractEventSPtr();
    if (!polyhedron || polyhedron->facets().size() == 0) {
        return;
    }

    if (!save_offsets_.empty()) {
        queue.push(SaveOffsetEvent::create(save_offsets_.front()));
    }

    CGAL::FT const_offset = util::Configuration::getInstance()->getDouble(
            "algo_3d_SimpleStraightSkel", "const_offset");
    if (const_offset != 0.0) {
        CGAL::FT next_offset = floor(CGAL::to_double(current_offset/const_offset) + 1.0) * const_offset; // @fixme to interval?
        if (next_offset >= current_offset) {
            next_offset += const_offset;
        }

        queue.push(ConstOffsetEvent::create(next_offset));
    }
    // two types of useless events:
    // - events that are in the past (i.e., offset > current_offset) (values are negative!)
    // - events that are (stricly) later than the current next tentative offset (i.e., offset < curr_next_offset)
    CGAL::FT current_offset_to_nearest_event = - (std::numeric_limits<double>::max());

    // --- Vanish Events

    collectEdgeEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectEdgeMergeEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectTriangleEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectDblEdgeMergeEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectDblTriangleEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectTetrahedronEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);

    // --- Contact Event

    collectVertexEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectFlipVertexEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectPolyhedronSplitEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectSplitMergeEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectPierceEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);

    // the next two are particularly slow, so reduce the bound by going through other events first
    collectSurfaceEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectEdgeSplitEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
}

// Note that this takes care of the pop
// @fixme can handling an event reveal another event at the same times?
std::pair<AbstractEventSPtr, bool> SimpleStraightSkel::nextEvent(PQ& queue) {
    if (queue.empty()) {
        return { };
    }

    AbstractEventSPtr event = queue.top();
    queue.pop();

    if (queue.empty()) {
        return { event, false /*not simultaneous*/};
    }

    // we could have simultaneous events at the save event but treat them *after*
    //
    // @fixme we expect the save output to be what is obtained *AFTER* the events
    // (i.e., surfaces merging) not *BEFORE* the events...
    if (event->getType() == AbstractEvent::SAVE_OFFSET_EVENT) {
        return { event, false  /*not simultaneous*/ };
    }

    AbstractEventSPtr next_event = queue.top();
    bool simultaneousEvents = (event->getOffset() == next_event->getOffset());
    if (event->getType() == AbstractEvent::CONST_OFFSET_EVENT && simultaneousEvents) {
        return nextEvent(queue); // the const event has been popped and is ignored
    }

    return { event, simultaneousEvents };
}

// normalize coefficients and make sure that coplanar facets have the same coefficients (@todo for the second part)
void SimpleStraightSkel::normalizeFacetPlanes(PolyhedronSPtr polyhedron)
{
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;

        Plane3SPtr plane = facet->plane(); // calls initPlane() if needed

#ifdef USE_CGAL
        const CGAL::FT a = plane->a();
        const CGAL::FT b = plane->b();
        const CGAL::FT c = plane->c();
        const CGAL::FT d = plane->d();
        // this should be the only place with unavoidable SQRTs
        const CGAL::FT n = CGAL::approximate_sqrt(CGAL::square(a) + CGAL::square(b) + CGAL::square(c));
#else
        const double a = plane->getA();
        const double b = plane->getB();
        const double c = plane->getC();
        const double d = plane->getD();
        const double n = sqrt(a*a + b*b + c*c);
#endif

        if(!is_zero(n))
        {
          Plane3SPtr normalized_plane = KernelFactory::createPlane3(a/n, b/n, c/n, d/n);
          facet->setPlane(normalized_plane);
        }
    }
}

// normalize coefficients and make sure that coplanar facets have the same coefficients (@todo for the second part)
void SimpleStraightSkel::harmonizeFacetPlanes(PolyhedronSPtr polyhedron)
{
    std::cout << "\n> harmonizeFacetPlanes()" << std::endl;

    // @todo this needs to be rewritten to use kernel predicates to sort the planes
    // Order all the supporting planes in a global order
    // The order is completely arbitrary, the only thing that we care about
    // is to quickly find which planes are identical

    // @todo also harmonize parallel but non equal planes (including opposite directions)

    // for now, abuse unordered containers because it's less code to write than above
    struct FEqual
    {
        bool operator()(FacetSPtr facet_1, FacetSPtr facet_2) const {
            Plane3SPtr plane_1 = facet_1->plane(); // calls initPlane() if needed
            Plane3SPtr plane_2 = facet_2->plane();
            std::cout << "facet #" << facet_1->getID() << " VS facet #" << facet_2->getID()
                      << " --> " << (*plane_1 == *plane_2) << std::endl;
            return CGAL::parallel(*plane_1, *plane_2);
        }
    };

    struct FHash
    {
        std::size_t operator()(FacetSPtr facet) const {
            return 1; // equality is only checked if the hash is the same
        }
    };

    std::unordered_map<FacetSPtr, Plane3SPtr, FHash, FEqual> facet_coefficients;

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;

        auto res = facet_coefficients.emplace(facet, Plane3SPtr{});
        if(res.second) { // never seen this direction of planes before
            Plane3SPtr plane = facet->plane(); // calls initPlane() if needed

#ifdef USE_CGAL
            const CGAL::FT a = plane->a();
            const CGAL::FT b = plane->b();
            const CGAL::FT c = plane->c();
            const CGAL::FT d = plane->d();
            // this should be the only place with unavoidable SQRTs
            const CGAL::FT n = CGAL::approximate_sqrt(CGAL::square(a) + CGAL::square(b) + CGAL::square(c));
#else
            const double a = plane->getA();
            const double b = plane->getB();
            const double c = plane->getC();
            const double d = plane->getD();
            const double n = sqrt(a*a + b*b + c*c);
#endif

            CGAL_assertion(!is_zero(n)); // degenerate faces are not allowed

            Plane3SPtr normalized_plane = KernelFactory::createPlane3(a/n, b/n, c/n, d/n);
            facet->setPlane(normalized_plane);
            res.first->second = normalized_plane;
        } else {
            std::cout << "Facet #" << facet->getID() << " is reusing coefficients from Facet #" << res.first->first->getID() << std::endl;

            // a parallel plane already exists, so we re-use its a, b, c coordinates
            // but need the 'd' to be shifted so the plane goes through the facet's vertices
            Plane3SPtr plane = facet->plane(); // calls initPlane() if needed
            Plane3SPtr parallel_plane = res.first->second;

            Vector3SPtr n = KernelFactory::createVector3(plane);
            Vector3SPtr n_p = KernelFactory::createVector3(parallel_plane);
            CGAL::FT sign = (CGAL::scalar_product(*n, *n_p) > 0) ? 1 : -1;

            VertexSPtr v = facet->vertices().front();
            Point3SPtr p = v->getPoint();
            const CGAL::FT a = sign * parallel_plane->a();
            const CGAL::FT b = sign * parallel_plane->b();
            const CGAL::FT c = sign * parallel_plane->c();
            const CGAL::FT d = - (a * p->x() + b * p->y() + c * p->z());
            Plane3SPtr normalized_plane = KernelFactory::createPlane3(a, b, c, d);
            CGAL_assertion(normalized_plane->has_on(*p));
            facet->setPlane(normalized_plane);
        }

        const CGAL::FT a = facet->plane()->a();
        const CGAL::FT b = facet->plane()->b();
        const CGAL::FT c = facet->plane()->c();
        const CGAL::FT d = facet->plane()->d();
        CGAL_assertion_code(CGAL::FT sq_n = a*a + b*b + c*c;)
        // std::cout << "normalization check (post harmonize): " << sq_n
        //           << " (" << CGAL::to_interval(sq_n).first << "; "
        //           << CGAL::to_interval(sq_n).second << ")" << std::endl;
        // std::cout << "a: " << a << std::endl;
        // std::cout << "b: " << b << std::endl;
        // std::cout << "c: " << c << std::endl;
        // std::cout << "d: " << d << std::endl;

        // inaccuracies during normalization since the sqrt is (usually) not exact
        CGAL_assertion((sq_n - 1) < 1e-5);
    }
}

void SimpleStraightSkel::appendEventNode(NodeSPtr node) {
    std::list<ArcWPtr>::iterator it_a = node->arcs().begin();
    while (it_a != node->arcs().end()) {
        ArcWPtr arc_wptr = *it_a++;
        if (!arc_wptr.expired()) {
            ArcSPtr arc = ArcSPtr(arc_wptr);
            arc->setNodeDst(node);
            arc->setNodeDstListIt(
                    std::find(node->arcs().begin(), node->arcs().end(), arc_wptr));
        }
    }
    std::list<SheetWPtr>::iterator it_s = node->sheets().begin();
    while (it_s != node->sheets().end()) {
        SheetWPtr sheet_wptr = *it_s++;
        if (!sheet_wptr.expired()) {
            SheetSPtr sheet = SheetSPtr(sheet_wptr);
            sheet->addNode(node);
        }
    }
    skel_result_->addNode(node);
}

void SimpleStraightSkel::handleEdgeEvent(EdgeEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());

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

    std::cout << "########################################" << std::endl;
    std::cout << "#########  Handle Edge Event  ##########" << std::endl;
    std::cout << "########################################" << std::endl;
    std::cout << node->toString() << std::endl;
    std::cout << edge->toString() << std::endl;
    std::cout << "Edge Offset: " << edge_offset->toString() << std::endl;

    // grab vertices and facets counter clockwise around node
    VertexSPtr vertices[4];
    for (unsigned int i = 0; i < 4; i++) {
        vertices[i] = VertexSPtr();
    }

    // @fixme? should seek the first non-collinear vertex? Or is a different type of event?
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

    std::cout << "vertices:\n"
              << "     " << *(vertices[0]->getPoint()) << "\n"
              << "     " << *(vertices[1]->getPoint()) << "\n"
              << "     " << *(vertices[2]->getPoint()) << "\n"
              << "     " << *(vertices[3]->getPoint()) << std::endl;

    for(std::size_t i=0; i<4; ++i)
    {
        std::cout << "facet #" << i << "\n";
        for(auto p : facets[i]->vertices_)
            std::cout << "  " << *(p->getPoint()) << "\n";
    }

    // check if edge should be flipped
    // @todo if one fails, just take the other immediately without SI checks
    bool flip_edge = true;

    bool not_flipped_valid = false;
    {
        std::cout << "== Trying WITHOUT a flip ==" << std::endl;

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
        FacetSPtr facets_clone[4];
        for (unsigned int i = 0; i < 4; i++) {
            facets_clone[i] = Facet::create();
            facets_clone[i]->setPlane(facets[i]->getPlane());
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

        PolyhedronSPtr polyhedron_no_flip = Polyhedron::create(4, facets_clone);
        PolyhedronSPtr polyhedron_no_flip_offset = PolyhedronTransformation::shiftFacets(polyhedron_no_flip, -1.0);

        if (!SelfIntersection::hasSelfIntersectingSurface(polyhedron_no_flip_offset)) {
            not_flipped_valid = true;
        }
    }

    std::cout << "not_flipped_valid = " << not_flipped_valid << std::endl;

    bool flipped_valid = false;
    {
        std::cout << "== Trying WITH a flip ==" << std::endl;

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
        FacetSPtr facets_clone[4];
        for (unsigned int i = 0; i < 4; i++) {
            facets_clone[i] = Facet::create();
            facets_clone[i]->setPlane(facets[i]->getPlane());
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

        PolyhedronSPtr polyhedron_flipped = Polyhedron::create(4, facets_clone);
        PolyhedronSPtr polyhedron_flipped_offset = PolyhedronTransformation::shiftFacets(polyhedron_flipped, -1.0);

        if (!SelfIntersection::hasSelfIntersectingSurface(polyhedron_flipped_offset)) {
            flipped_valid = true;
        }
    }

    std::cout << "flipped_valid = " << flipped_valid << std::endl;

    if (flipped_valid && !not_flipped_valid) {
        flip_edge = true;
    } else if (not_flipped_valid && !flipped_valid) {
        flip_edge = false;
    } else if (flipped_valid && not_flipped_valid) {

// #define CGAL_SS3_USE_APPROXIMATE_ANGLES // old code
#ifdef CGAL_SS3_USE_APPROXIMATE_ANGLES
        double angle_flipped = KernelWrapper::angle(facets[0]->plane(), facets[2]->plane());
        double angle_no_flip = KernelWrapper::angle(facets[1]->plane(), facets[3]->plane());
        DEBUG_VAR(angle_no_flip);
        DEBUG_VAR(angle_flipped);
        if (edge_event_ == 0) {
            // convex
            flip_edge = (angle_flipped <= angle_no_flip);
        } else if (edge_event_ == 1) {
            // reflex
            // choose edge that moves faster
            flip_edge = (angle_flipped >= angle_no_flip);
        }
#else
# ifndef USE_CGAL
#  error
# endif
        Vector3SPtr n0 = KernelFactory::createVector3(facets[0]->plane());
        Vector3SPtr n2 = KernelFactory::createVector3(facets[2]->plane());
        Vector3SPtr n1 = KernelFactory::createVector3(facets[1]->plane());
        Vector3SPtr n3 = KernelFactory::createVector3(facets[3]->plane());
        CGAL::Comparison_result dac = CGAL::compare_angle(*n0, *n2, *n1, *n3);

        if (edge_event_ == 0) {
            // convex
            flip_edge = (dac != CGAL::LARGER); // (angle_flipped <= angle_no_flip);
        } else if (edge_event_ == 1) {
            // reflex
            // choose edge that moves faster
            flip_edge = (dac != CGAL::SMALLER); // (angle_flipped >= angle_no_flip);
        }
#endif
        else if (edge_event_ == 2) {
            // flip_when_possible
            flip_edge = true;
        } else if (edge_event_ == 3) {
            // sphere
            VertexSPtr vertex_c = Vertex::create(node->getPoint());
            EdgeSPtr edges[4];
            for (unsigned int i = 0; i < 4; i++) {
                edges[i] = Edge::create(vertex_c, vertices[i]->clone());
            }
            FacetSPtr facets_c[4];
            for (unsigned int i = 0; i < 4; i++) {
                facets_c[i] = Facet::create();
                facets_c[i]->setPlane(facets[i]->getPlane());
            }
            for (unsigned int i = 0; i < 4; i++) {
                edges[i]->setFacetR(facets_c[(i+3)%4]);
                edges[i]->setFacetL(facets_c[i]);
                facets_c[(i+3)%4]->addEdge(edges[i]);
                facets_c[i]->addEdge(edges[i]);
            }
            PolyhedronSPtr polyhedron_sphere = Polyhedron::create(4, facets_c);
            SphereVertexSplitterSPtr splitter = SphereVertexSplitter::create();
            splitter->splitVertex(vertex_c);
            if (facets_c[0]->findEdge(facets_c[2])) {
                flip_edge = true;
            } else if (facets_c[1]->findEdge(facets_c[3])) {
                flip_edge = false;
            } else {
                throw std::runtime_error("Not able to handle EdgeEvent (1).");
            }
        }
    } else {
        throw std::runtime_error("Not able to handle EdgeEvent (2).");
    }
    DEBUG_VAR(flip_edge);

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
    }

    // update arcs and sheets
    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
            edge_offset->getData());
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
    skel_result_->addEvent(event);
}

void SimpleStraightSkel::handleEdgeMergeEvent(EdgeMergeEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());
    appendEventNode(event->getNode());

    std::cout << "########################################" << std::endl;
    std::cout << "######  Handle Edge Merge Event  #######" << std::endl;
    std::cout << "########################################" << std::endl;

    SkelFacetDataSPtr facet_data = std::dynamic_pointer_cast<SkelFacetData>(
            event->getFacet()->getData());
    FacetSPtr facet = facet_data->getOffsetFacet();
    SkelEdgeDataSPtr edge_data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
                event->getEdge1()->getData());
    SkelEdgeDataSPtr edge_data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
                event->getEdge2()->getData());
    EdgeSPtr edge_1 = edge_data_1->getOffsetEdge();
    EdgeSPtr edge_2 = edge_data_2->getOffsetEdge();

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
    std::list<FacetWPtr>::iterator it_f = vertex_1->facets().begin();
    while (it_f != vertex_1->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (!facet_wptr.expired()) {
            FacetSPtr facet = FacetSPtr(facet_wptr);
            facet->removeVertex(vertex_1);
        }
    }
    polyhedron->removeVertex(vertex_1);
    it_f = vertex_2->facets().begin();
    while (it_f != vertex_2->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (!facet_wptr.expired()) {
            FacetSPtr facet = FacetSPtr(facet_wptr);
            facet->removeVertex(vertex_2);
        }
    }
    polyhedron->removeVertex(vertex_2);

    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
            vertex->getData());
    vertex_data->setNode(event->getNode());
    ArcSPtr arc = createArc(vertex);
    skel_result_->addArc(arc);
    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);

    if (crashAt(edge_1, edge_b).first || crashAt(edge_b, edge_1).first) {
        // crashAt(...) should be commutative,
        // but because of rounding errors it's not in certain cases.
        // TODO: this may need an improvement
        // fix coordinates because point is exactly on merged edge
        DEBUG_VAL("Fixing coordinates of " << vertex->toString());
        Point3SPtr p = vertex->getPoint();
        Vector3SPtr dir = KernelWrapper::normalize(arc->getDirection());
        Point3SPtr p_fixed = KernelFactory::createPoint3((*p) + (*dir)*0.0001);
        event->getNode()->setPoint(p_fixed);
        vertex->setPoint(p_fixed);
    }
}

void SimpleStraightSkel::handleTriangleEvent(TriangleEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());
    appendEventNode(event->getNode());

    std::cout << "########################################" << std::endl;
    std::cout << "#######  Handle Triangle Event  ########" << std::endl;
    std::cout << "########################################" << std::endl;

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

    std::cout << "VS:\n"
              << vertices[0]->toString() << "\n"
              << vertices[1]->toString() << "\n"
              << vertices[2]->toString() << std::endl;
    std::cout << "VSO:\n"
              << vertices_offset[0]->toString() << "\n"
              << vertices_offset[1]->toString() << "\n"
              << vertices_offset[2]->toString() << std::endl;

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
        EdgeSPtr edge_offset = EdgeSPtr(vertices_offset[i]->edges().front());
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

    SkelVertexDataSPtr data_offset = SkelVertexData::create(vertex_offset);
    data_offset->setNode(event->getNode());
    ArcSPtr arc = createArc(vertex_offset);
    skel_result_->addArc(arc);
    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
}

void SimpleStraightSkel::handleDblEdgeMergeEvent(DblEdgeMergeEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());
    appendEventNode(event->getNode());

    std::cout << "########################################" << std::endl;
    std::cout << "#######  Handle Dbl Edge Event  ########" << std::endl;
    std::cout << "########################################" << std::endl;

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
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet = FacetSPtr(facet_wptr);
                facet->removeVertex(vertex);
            }
        }
        polyhedron->removeVertex(vertex);
    }

    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
}

void SimpleStraightSkel::handleDblTriangleEvent(DblTriangleEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());
    appendEventNode(event->getNode());

    std::cout << "########################################" << std::endl;
    std::cout << "#####  Handle Dbl Triangle Event  ######" << std::endl;
    std::cout << "########################################" << std::endl;

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
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet = FacetSPtr(facet_wptr);
                facet->removeVertex(vertex);
            }
        }
        polyhedron->removeVertex(vertex);
    }

    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
}

void SimpleStraightSkel::handleTetrahedronEvent(TetrahedronEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());
    appendEventNode(event->getNode());

    std::cout << "########################################" << std::endl;
    std::cout << "######  Handle Tetrahedron Event  ######" << std::endl;
    std::cout << "########################################" << std::endl;

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
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet = FacetSPtr(facet_wptr);
                facet->removeVertex(vertex);
            }
        }
        polyhedron->removeVertex(vertex);
    }

    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
}

void SimpleStraightSkel::handleVertexEvent(VertexEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());
    appendEventNode(event->getNode());

    std::cout << "########################################" << std::endl;
    std::cout << "#######  Handle Vertex Event  #########" << std::endl;
    std::cout << "########################################" << std::endl;

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
    std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
    while (it_e1 != vertex_1->edges().end()) {
        EdgeWPtr edge_wptr = *it_e1++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
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
    std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
    while (it_e2 != vertex_2->edges().end()) {
        EdgeWPtr edge_wptr = *it_e2++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
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

    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
            edge_tomerge_2->getData());
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
    skel_result_->addEvent(event);
}

void SimpleStraightSkel::handleFlipVertexEvent(FlipVertexEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());
    appendEventNode(event->getNode());

    std::cout << "########################################" << std::endl;
    std::cout << "######  Handle Flip Vertex Event  ######" << std::endl;
    std::cout << "########################################" << std::endl;

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
    std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
    while (it_e1 != vertex_1->edges().end()) {
        EdgeWPtr edge_wptr = *it_e1++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
            if ((edge->getFacetL() == facet_1 && edge->getFacetR() == facet_2) ||
                    (edge->getFacetL() == facet_2 && edge->getFacetR() == facet_1)) {
                edge_1 = edge;
                break;
            }
        }
    }
    EdgeSPtr edge_2;
    std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
    while (it_e2 != vertex_2->edges().end()) {
        EdgeWPtr edge_wptr = *it_e2++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
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

    vertex_data_1 = std::dynamic_pointer_cast<SkelVertexData>(vertex_1->getData());
    vertex_data_1->setNode(event->getNode());
    vertex_data_2 = std::dynamic_pointer_cast<SkelVertexData>(vertex_2->getData());
    vertex_data_2->setNode(event->getNode());
    ArcSPtr arc_1 = createArc(vertex_1);
    skel_result_->addArc(arc_1);
    ArcSPtr arc_2 = createArc(vertex_2);
    skel_result_->addArc(arc_2);

    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
}

void SimpleStraightSkel::handleSurfaceEvent(SurfaceEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());

    std::cout << "########################################" << std::endl;
    std::cout << "#######  Handle Surface Event  #########" << std::endl;
    std::cout << "########################################" << std::endl;

    NodeSPtr node = event->getNode();
    std::cout << "N: " << node->toString() << std::endl;

    appendEventNode(node);

    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge1()->getData());
    EdgeSPtr edge_1 = data_1->getOffsetEdge();
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge2()->getData());
    EdgeSPtr edge_2 = data_2->getOffsetEdge();

    FacetSPtr facet_1_src = getFacetSrc(edge_1);
    FacetSPtr facet_1_dst = getFacetDst(edge_1);

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

    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
            edge_2->getData());
    SheetSPtr sheet = edge_data->getSheet();
    edge_data = SkelEdgeData::create(edge_tmp);
    edge_data->setSheet(sheet);

    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
            vertex->getData());
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
    skel_result_->addEvent(event);

    std::cout << "########################################" << std::endl;
    std::cout << "#####  Handle Surface Event (END)  #####" << std::endl;
    std::cout << "########################################" << std::endl;
}

void SimpleStraightSkel::handlePolyhedronSplitEvent(PolyhedronSplitEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());

    std::cout << "########################################" << std::endl;
    std::cout << "####  Handle Polyhedron Split Event  ###" << std::endl;
    std::cout << "########################################" << std::endl;

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge1()->getData());
    EdgeSPtr edge_1 = data_1->getOffsetEdge();
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge2()->getData());
    EdgeSPtr edge_2 = data_2->getOffsetEdge();

    FacetSPtr facet_1_src = getFacetSrc(edge_1);
    FacetSPtr facet_1_dst = getFacetDst(edge_1);

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

    SkelVertexDataSPtr data_l = std::dynamic_pointer_cast<SkelVertexData>(
            vertex_l->getData());
    data_l->setNode(node);
    SkelVertexDataSPtr data_r = std::dynamic_pointer_cast<SkelVertexData>(
            vertex_r->getData());
    data_r->setNode(node);
    SkelEdgeDataSPtr data_22 = std::dynamic_pointer_cast<SkelEdgeData>(
            edge_22->getData());
    data_22->setSheet(data_2->getSheet());
    ArcSPtr arc_l = createArc(vertex_l);
    ArcSPtr arc_r = createArc(vertex_r);
    skel_result_->addArc(arc_l);
    skel_result_->addArc(arc_r);

    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
}

void SimpleStraightSkel::handleSplitMergeEvent(SplitMergeEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());
    appendEventNode(event->getNode());

    std::cout << "########################################" << std::endl;
    std::cout << "#####  Handle Split Merge Event  #######" << std::endl;
    std::cout << "########################################" << std::endl;

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
    std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
    while (it_e1 != vertex_1->edges().end()) {
        EdgeWPtr edge_wptr = *it_e1++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
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
    std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
    while (it_e2 != vertex_2->edges().end()) {
        EdgeWPtr edge_wptr = *it_e2++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
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
    DEBUG_SPTR(edge_tosplit);

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

    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
            edge_tosplit->getData());
    SheetSPtr sheet = edge_data->getSheet();
    edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
            edge_tomerge_2->getData());
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
    skel_result_->addEvent(event);

    std::cout << "########################################" << std::endl;
    std::cout << "###  Handle Split Merge Event (END)  ###" << std::endl;
    std::cout << "########################################" << std::endl;
}

void SimpleStraightSkel::handleEdgeSplitEvent(EdgeSplitEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());

    std::cout << "########################################" << std::endl;
    std::cout << "######  Handle Edge Split Event  #######" << std::endl;
    std::cout << "########################################" << std::endl;

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge1()->getData());
    EdgeSPtr edge_1 = data_1->getOffsetEdge();
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge2()->getData());
    EdgeSPtr edge_2 = data_2->getOffsetEdge();

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
    if (KernelWrapper::orientation(line(event->getEdge1()),
            line(event->getEdge2())) > 0) {
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
    skel_result_->addEvent(event);
}

void SimpleStraightSkel::handlePierceEvent(PierceEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());

    std::cout << "########################################" << std::endl;
    std::cout << "########  Handle Pierce Event  #########" << std::endl;
    std::cout << "########################################" << std::endl;

    NodeSPtr node = event->getNode();

    std::cout << "Node: " << node->toString() << std::endl;
    std::cout << "V: " << event->getVertex()->toString() << std::endl;
    std::cout << "F: " << event->getFacet()->toString() << std::endl;

    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
            event->getVertex()->getData());
    VertexSPtr vertex_offset = vertex_data->getOffsetVertex();
    SkelFacetDataSPtr facet_data = std::dynamic_pointer_cast<SkelFacetData>(
            event->getFacet()->getData());
    FacetSPtr facet_offset = facet_data->getOffsetFacet();

#ifndef CGAL_SS3_OLD_CODE_PIERCE_EVENT
    // Check if it is a valid pierce event: the piercing point must be in the face
    bool reject_event = false;

# if 1 // not sure which one is better
    // with ray shooting
    reject_event = !(SelfIntersection::isInsideWithRayShooting(node->getPoint(), facet_offset));
# else
    using Itag = CGAL::No_constraint_intersection_tag;
    using PK = CGAL::Projection_traits_3<CGAL::K>;
    using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>; // @todo not needed
    using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
    using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
    using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
    using PCDT_VH = PCDT::Vertex_handle;
    using PCDT_FH = PCDT::Face_handle;

    Vector3SPtr n = KernelFactory::createVector3(facet_offset->plane());
    CGAL_assertion(*n != CGAL::NULL_VECTOR);
    PK traits(*n);
    PCDT pcdt(traits);

    std::map<VertexSPtr, PCDT_VH> face_vhs;
    std::list<VertexSPtr>::iterator it_v = facet_offset->vertices().begin();
    while (it_v != facet_offset->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        auto res = face_vhs.emplace(vertex, PCDT_VH());
        if(res.second) // first time seeing this point
        {
            PCDT_VH vh = pcdt.insert(*(vertex->getPoint()));
            res.first->second = vh;
        }
    }

    std::list<EdgeSPtr>::iterator it_e = facet_offset->edges().begin();
    while (it_e != facet_offset->edges().end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr v0 = edge->src(facet_offset);
        VertexSPtr v1 = edge->dst(facet_offset);

        if(*(v0->getPoint()) == *(v1->getPoint()))
        {
            std::cerr << "W: encountered degenerate edge @ " << *(v0->getPoint()) << std::endl;
        }
        else
        {
            PCDT_VH vh0 = face_vhs.at(v0);
            PCDT_VH vh1 = face_vhs.at(v1);

            try
            {
                pcdt.insert_constraint(vh0, vh1);
            }
            catch(const typename PCDT::Intersection_of_constraints_exception&)
            {
                std::cerr << "Error: Intersection of constraints" << std::endl;
                DEBUG_VAR(facet_offset->toString());
                CGAL_assertion_msg(false, "Intersections in CDT2 not allowed");
                break;
            }
        }
    }

    std::unordered_map<PCDT_FH, bool> in_domain_map;
    boost::associative_property_map< std::unordered_map<PCDT_FH, bool> > in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(pcdt, in_domain);

    int li;
    PCDT::Locate_type lt;
    PCDT_FH fh = pcdt.locate(*(node->getPoint()), lt, li);

    std::cout << "lt = " << lt << std::endl;

    if (lt == PCDT::VERTEX || lt == PCDT::EDGE) {
        CGAL_assertion(false); // @todo handle this by looping over incident faces
    } else if (lt == PCDT::FACE) {
        reject_event = !(get(in_domain, fh));
    } else { // outside of convex hull
        reject_event = true;
    }
# endif // ray shooting or CDT
    if (reject_event) {
        std::cout << "Pierce rejected" << std::endl;
        event->setPolyhedronResult(polyhedron);
        // skel_result_->addEvent(event);
        return;
    }
#endif

    std::cout << "Pierce accepted" << std::endl;

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
    skel_result_->addEvent(event);
}


StraightSkeletonSPtr SimpleStraightSkel::getResult() const {
    return this->skel_result_;
}

} }
