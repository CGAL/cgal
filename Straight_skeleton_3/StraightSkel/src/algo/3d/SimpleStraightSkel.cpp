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

// @todo figure out a proper stack so shared pointers do not contaminate everything
// currently, we have:
// - event -> polyhedron
// - skeleton -> event
// - skeleton -> polyhedron
// etc.
#define CGAL_SS3_NO_SKELETON_DS

// #define CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
#define CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
// #define CGAL_SS3_USE_AUTOREF_FOR_ALL_EVENTS
#define CGAL_SS3_USE_AUTOREF_FOR_SIMULTANEOUS_EVENTS
#define CGAL_SS3_FILTER_VOLUMES_WITH_ONLY_REACHABLE_FACES
// #define CGAL_SS3_DO_NOT_ADD_UNREACHED_TRIANGLES_TO_CONTRIBUTIONS

#if defined(CGAL_SS3_USE_AUTOREF_FOR_ALL_EVENTS) || defined(CGAL_SS3_USE_AUTOREF_FOR_SIMULTANEOUS_EVENTS)
# ifndef CGAL_SS3_NO_SKELETON_DS
#  error "You must disable Skeleton DS" // can't maintain sheets and arcs in that mode for now
# endif
#endif

#define CGAL_SS3_PROFILE_FILTERING_MECHANISMS

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
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <limits>
#include <list>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>

namespace algo { namespace _3d {

template <typename ValueType,
          typename VisitorBase = CGAL::Polygon_mesh_processing::Autorefinement::Default_visitor>
struct Range_updating_autoref_visitor : public VisitorBase {
    Range_updating_autoref_visitor(const std::vector<ValueType>& old_range,
                                   std::vector<ValueType>& new_range,
                                   const VisitorBase& base = VisitorBase{})
        : VisitorBase(base), old_range_(old_range), new_range_(new_range) {
        new_range.reserve(old_range.size());
    }

    void verbatim_triangle_copy(std::size_t tgt_id, std::size_t src_id) {
        // std::cout << "verbatim_triangle_copy " << tgt_id << " from " << src_id << std::endl;
        VisitorBase::verbatim_triangle_copy(tgt_id, src_id);
        new_range_.resize(tgt_id + 1);
        new_range_[tgt_id] = old_range_[src_id];
    }

    void new_subtriangle (std::size_t tgt_id, std::size_t src_id) {
        // std::cout << "new_subtriangle " << tgt_id << " from " << src_id << std::endl;
        VisitorBase::new_subtriangle(tgt_id, src_id);
        new_range_.resize(tgt_id + 1);
        new_range_[tgt_id] = old_range_[src_id];
    }

private:
    const std::vector<ValueType>& old_range_;
    std::vector<ValueType>& new_range_;
};

template <typename ValueType,
          typename BaseVisitor = CGAL::Polygon_mesh_processing::internal::Default_repair_PS_visitor>
struct Range_updating_repair_PS_visitor : public BaseVisitor {
    Range_updating_repair_PS_visitor(std::vector<ValueType>& range,
                                     const BaseVisitor& base_visitor = BaseVisitor{})
        : BaseVisitor(base_visitor), range_(range) { }

    void swap(std::size_t pos_1, std::size_t pos_2) {
        BaseVisitor::swap(pos_1, pos_2);
        std::swap(range_[pos_1], range_[pos_2]);
    }
    void resize(std::size_t new_size) {
        BaseVisitor::resize(new_size);
        range_.resize(new_size);
    }
private:
    std::vector<ValueType>& range_;
};

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
                                       const std::list<CGAL::FT>& save_offsets,
                                       const std::filesystem::path& save_path) {
    polyhedron_ = polyhedron;
    controller_ = controller;
    save_offsets_ = save_offsets;
    save_path_ = save_path;
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
                                                  const std::list<CGAL::FT>& save_offsets,
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

bool SimpleStraightSkel::isLocked(VertexSPtr vertex) {
    return false;

    Point3SPtr recomputed_point = PolyhedronTransformation::shiftPoint(vertex, 0);

    if (!recomputed_point) {
        std::cout << vertex->getID() << " is locked" << std::endl;
    }

    return (recomputed_point == Point3SPtr());
}

bool SimpleStraightSkel::isLocked(EdgeSPtr edge) {
    // only lock degenerate edges, as we should not have to recompute the position otherwise
    return (edge->getVertexSrc()->getPoint() == edge->getVertexDst()->getPoint()) &&
            (isLocked(edge->getVertexSrc()) || isLocked(edge->getVertexDst()));
}

// @fixme return three states to handle flat cases?
bool SimpleStraightSkel::isReflex(EdgeSPtr edge) {
    bool result = false;

    VertexSPtr vertex_src = edge->getVertexSrc();
    VertexSPtr vertex_dst = edge->getVertexDst();

    if (*(vertex_src->getPoint()) == *(vertex_dst->getPoint())) {
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        FacetSPtr facet_src = getFacetSrc(edge);
        FacetSPtr facet_dst = getFacetDst(edge);

        // std::cout << "Facet L = " << facet_l->getID() << std::endl;
        // std::cout << "Facet R = " << facet_r->getID() << std::endl;
        // std::cout << "Facet SRC = " << facet_src->getID() << std::endl;
        // std::cout << "Facet DST = " << facet_dst->getID() << std::endl;

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

        Plane3SPtr offset_plane_l = KernelWrapper::offsetPlane(facet_l->plane(), - speed_l);
        Plane3SPtr offset_plane_r = KernelWrapper::offsetPlane(facet_r->plane(), - speed_r);
        Plane3SPtr offset_plane_src = KernelWrapper::offsetPlane(facet_src->plane(), - speed_src);
        Plane3SPtr offset_plane_dst = KernelWrapper::offsetPlane(facet_dst->plane(), - speed_dst);

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
    std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
    while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge = EdgeSPtr(edge_wptr);
            if (isLocked(edge)) {
                // @fixme I have no idea what effects this 'continue' has...
                //
                // Even if the edge is degen and locked, it still has geometry
                // and could be convex or reflex?
                continue;
            }
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
            if (isLocked(edge)) {
                continue; // I have no idea what effects this has; see isReflex()
            }
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

std::pair<PolyhedronSPtr, CGAL::FT> SimpleStraightSkel::enablePerturbedMode(PolyhedronSPtr polyhedron,
                                                                            CGAL::FT currentOffset,
                                                                            CGAL::FT simultaneousOffset) {
    std::cout << "Enabling perturbed mode..." << std::endl;

    // @fixme the save offset should be taken into account in the limitation of the perturbation:
    // we need perturbed mode to end before the next save offset

    // shift a little bit the polyhedron inwards to avoid self intersections
    CGAL::FT shift = (currentOffset + simultaneousOffset) / 2 - currentOffset; // @fixme complete hack for now
    std::cout << "pre-enable shift = " << shift << std::endl;
    CGAL_assertion(!is_zero(shift));
    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, shift);

    db::_3d::OBJFile::save("results/pertubation_pre_enable.obj", polyhedron, false /*do not triangulate*/);

    perturbationOffset_ = currentOffset + shift;

    // perturb
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

    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, 0.0);
    db::_3d::OBJFile::save("results/pertubation_post_enable.obj", polyhedron, false /*do not triangulate*/);
    db::_3d::OBJFile::save("results/pertubation_post_enable_triangulated.obj", polyhedron, true /*do not triangulate*/);
    CGAL_assertion(bool(polyhedron));
    CGAL_assertion(polyhedron->isConsistent());

    usingTemporaryPerturbedMode_ = true;
    simultaneousOffset_ = simultaneousOffset;

    return { polyhedron, perturbationOffset_ };
}

std::pair<PolyhedronSPtr, CGAL::FT> SimpleStraightSkel::disablePerturbedMode(PolyhedronSPtr polyhedron,
                                                                             CGAL::FT currentOffset,
                                                                             CGAL::FT nextEventOffset) {
    std::cout << "Disabling perturbed mode..." << std::endl;
    std::cout << "  Current offset: " << currentOffset << std::endl;
    std::cout << "  Next offset: " << nextEventOffset << std::endl;
    CGAL_precondition(usingTemporaryPerturbedMode_);

    // shift a little more to get some waylay
    CGAL::FT shift = (currentOffset + nextEventOffset) / 2 - currentOffset; // @fixme complete hack for now
    std::cout << "pre-disable shift = " << shift << std::endl;
    CGAL_assertion(!is_zero(shift));
    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, shift);

    db::_3d::OBJFile::save("results/pertubation_pre_disable.obj", polyhedron, false /*do not triangulate*/);

    CGAL::FT perturbationEndOffset = (currentOffset + shift);
    std::cout << "  Perturbation Start: " << perturbationOffset_ << std::endl;
    std::cout << "  Perturbation End: " << perturbationEndOffset << std::endl;

    // un-perturb
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
#ifdef CGAL_SS3_PERTURB_PLANE_COEFFICIENTS
        facet->restorePlaneCoefficients(perturbationOffset_, perturbationEndOffset);
#else
        CGAL_assertion(facet->hasData());
        SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
        data->setSpeed(facet->cachedSpeed_);
        std::cout << "speed is back to: " << data->getSpeed() << std::endl;
#endif
    }

    // @fixme what if faces incident to a vertex become a degenerate configuration and we can't
    // recompute the vertex? ==> just continue another iteration in perturbed mode?
    //
    // In this case, need to handle the case of the next offset event being the desired save offset

    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, 0.0);
    db::_3d::OBJFile::save("results/pertubation_post_disable.obj", polyhedron, false /*do not triangulate*/);
    CGAL_assertion(bool(polyhedron));

    usingTemporaryPerturbedMode_ = false;
    simultaneousOffset_ = 0;

    return { polyhedron, perturbationEndOffset };
}

bool SimpleStraightSkel::handleSaveEventAtSimultaneity(PolyhedronSPtr polyhedron,
                                                       CGAL::FT current_offset,
                                                       CGAL::FT simultaneity_offset) {
    namespace PMP = CGAL::Polygon_mesh_processing;

    std::cout << "handleSaveEventAtSimultaneity()" << std::endl;

    CGAL_precondition(current_offset != simultaneity_offset);

    const CGAL::FT shift = simultaneity_offset - current_offset;
    polyhedron = PolyhedronTransformation::shiftFacets(polyhedron, shift);

    std::stringstream ss_filename;
    ss_filename << "results/offset_" << simultaneity_offset << ".obj";
    db::_3d::OBJFile::save(ss_filename.str(), polyhedron);

    polyhedron->initializeAllIDs();

    // convert to triangulated soup
    std::vector<Point3> points;
    std::vector<std::vector<std::size_t> > triangles;

    // @todo factorize CDT usages but mind the tags
    using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
    using PK = CGAL::Projection_traits_3<CGAL::K>;
    using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
    using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
    using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
    using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
    using PCDT_VH = PCDT::Vertex_handle;
    using PCDT_FH = PCDT::Face_handle;

    // not really required because the IDs are initialized in the order of vertex iteration
    std::map<VertexSPtr, std::size_t> v_ids;

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        unsigned int id = vertex->getID();
        points.emplace_back(vertex->getX(), vertex->getY(), vertex->getZ());
        v_ids[vertex] = id;
    }

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        facet->makeFirstConvex();

        if (facet->edges().size() < 3) {
            std::cerr << "Warning: face with < 3 edges" << std::endl;
            continue;
        } else {
            Vector3SPtr n = KernelFactory::createVector3(facet->plane());
            CGAL_assertion(*n != CGAL::NULL_VECTOR);

            PK traits(*n);
            PCDT pcdt(traits);

            std::map<VertexSPtr, PCDT_VH> face_vhs; // might have multiple vertices at the same position

            std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
            while (it_v != facet->vertices().end()) {
                VertexSPtr vertex = *it_v++;
                auto res = face_vhs.emplace(vertex, PCDT_VH());
                if(res.second) // first time seeing this point
                {
                    PCDT_VH vh = pcdt.insert(*(vertex->getPoint()));
                    vh->info() = vertex;
                    res.first->second = vh;
                }
            }

            auto ne = 0;
            std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
            while (it_e != facet->edges().end()) {
                EdgeSPtr edge = *it_e++;
                VertexSPtr v0 = edge->src(facet);
                VertexSPtr v1 = edge->dst(facet);

                if(*(v0->getPoint()) == *(v1->getPoint()))
                {
                    std::cerr << "W: encountered degenerate edge @ " << *(v0->getPoint()) << std::endl;

                    CGAL_assertion(v0->degree() != 1); // @todo handle that...
                    VertexSPtr vm1 = edge->prev(facet)->src(facet);

                    // create a degenerate face, commented here because we will purge it anyway
                    // triangles.emplace_back({v_ids.at(vm1), v_ids.at(v0), v_ids.at(v1)});
                } else {
                    PCDT_VH vh0 = face_vhs.at(v0);
                    PCDT_VH vh1 = face_vhs.at(v1);

                    try {
                        pcdt.insert_constraint(vh0, vh1);
                    } catch(const typename PCDT::Intersection_of_constraints_exception&) {
                        std::cerr << "Error: Intersection of constraint w/ " << vh0->point() << " " << vh1->point() << std::endl;
                        DEBUG_VAR(facet->toString());
                        CGAL_warning_msg(false, "Intersections in CDT2 not allowed");
                        return false;
                    }
                    ++ne;
                }
            }

            if(ne < 3) { // degenerate face
                std::cerr << "Warning: skipping degenerate face" << std::endl;
                continue;
            }

            std::unordered_map<PCDT_FH, bool> in_domain_map;
            boost::associative_property_map<std::unordered_map<PCDT_FH, bool> > in_domain(in_domain_map);
            CGAL::mark_domain_in_triangulation(pcdt, in_domain);

            for(auto fh : pcdt.finite_face_handles()) {
                if(!get(in_domain, fh)) {
                    continue;
                }

                triangles.push_back({v_ids.at(fh->vertex(0)->info()),
                                     v_ids.at(fh->vertex(1)->info()),
                                     v_ids.at(fh->vertex(2)->info())});

            }
        }
    }

    CGAL::IO::write_OFF("results/simultaneous_save.off", points, triangles);

    // remove degenerate faces
    auto removal_predicate = [&](const std::vector<std::size_t>& polygon) {
        CGAL_precondition(polygon.size() == 3);
        return CGAL::collinear(points[polygon[0]], points[polygon[1]], points[polygon[2]]);
    };
    std::size_t size_pre = triangles.size();
    triangles.erase(std::remove_if(std::begin(triangles), std::end(triangles), removal_predicate),
                   std::end(triangles));
    std::cout << "Removed " << size_pre - triangles.size() << " degenerate faces" << std::endl;

    PMP::autorefine_triangle_soup(points, triangles,
                                  CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag()));
    CGAL::IO::write_OFF("results/simultaneous_save_autoref.off", points, triangles);

    PMP::repair_polygon_soup(points, triangles, CGAL::parameters::erase_all_duplicates(true)
                                                                .require_same_orientation(false)
                                                                .verbose(true));
    CGAL::IO::write_OFF("results/simultaneous_save_autoref_repaired.off", points, triangles);

    PMP::orient_polygon_soup(points, triangles);
    CGAL::IO::write_OFF("results/simultaneous_save_autoref_repaired_oriented.off", points, triangles);

    CGAL_assertion(PMP::does_triangle_soup_self_intersect(points, triangles));

    if (!PMP::is_polygon_soup_a_polygon_mesh(triangles)) {
        std::cerr << "Error: PS not a PM" << std::endl;
        return false;
    }

    return true;
}

bool SimpleStraightSkel::savePolyhedron(PolyhedronSPtr polyhedron,
                                        CGAL::FT current_offset,
                                        const bool do_triangulate,
                                        const bool dump_exact,
                                        const bool attempt_untilting)
{
    bool result;

    // attempt naive un-tilting
    if (attempt_untilting) {
        std::cout << "Attempting to un-tilt polyhedron..." << std::endl;

        db::_3d::OBJFile::save("results/pre-untilt_attempt.obj", polyhedron,
                               false /*do_triangulate*/,
                               true /*convert_to_double*/);

        PolyhedronSPtr polyhedron_cpy = polyhedron->clone();

        std::list<FacetSPtr>::iterator it_f = polyhedron_cpy->facets().begin();
        while (it_f != polyhedron_cpy->facets().end()) {
            FacetSPtr facet = *it_f++;

            // this assumes we have perturbed at the start ('0')
            facet->restorePlaneCoefficients(0, current_offset);
        }

        // As to avoid having a vertex not be defined by e.g. 3 planes with 2 being equal post un-tilt
        // That vertex is then useless, so just remove it
        // Do it here, before we recompute point positions
        db::_3d::AbstractFile::mergeCoplanarFacets(polyhedron_cpy, 0.0);

        db::_3d::OBJFile::save("results/restored_merged.obj", polyhedron_cpy,
                               false /*do_triangulate*/,
                               true /*convert_to_double*/);

        CGAL_assertion(polyhedron_cpy->isConsistent());

        db::_3d::AbstractFile::removeVerticesDegLt3(polyhedron_cpy);

        db::_3d::OBJFile::save("results/restored_final.obj", polyhedron_cpy,
                               false /*do_triangulate*/,
                               true /*convert_to_double*/);

        polyhedron_cpy = PolyhedronTransformation::shiftFacets(polyhedron_cpy, 0.0);
        if(!polyhedron_cpy || !polyhedron_cpy->isConsistent()) {
            std::cerr << "Warning: failed to un-tilt polyhedron" << std::endl;
            result = false;
        } else {
            polyhedron = polyhedron_cpy;
            if (SelfIntersection::hasSelfIntersectingSurface(polyhedron)) {
                std::cerr << "Warning: self-intersections after un-tilting" << std::endl;
            } else {
                std::cerr << "Successfully un-tilted polyhedron" << std::endl;
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

    DEBUG_PRINT("== Straight Skeleton 3D started ==");

    double t_start = util::Timer::now();

    PolyhedronSPtr polyhedron = polyhedron_->clone();

    basePlanes_.reserve(polyhedron->facets().size());
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        CGAL_assertion(bool(facet->getPlane()));
        facet->setBasePlaneID(basePlanes_.size());
        basePlanes_.push_back(facet->getPlane());
    }

    db::_3d::OBJFile::save("results/init_pre.obj", polyhedron, false /*do not triangulate*/);

// #define CGAL_SSE_ACUTE_WEIGHTS
// #define CGAL_SSE_MERGING_WEIGHTS
// #define CGAL_SSE_PERFORMANCE_WEIGHTS

#if defined(CGAL_SSE_ACUTE_WEIGHTS) || defined(CGAL_SSE_MERGING_WEIGHTS) || defined(CGAL_SSE_PERFORMANCE_WEIGHTS)
# ifdef CGAL_SSE_ACUTE_WEIGHTS
    const CGAL::FT x_speed = 20;
    const CGAL::FT y_speed = 20;
    const CGAL::FT z_speed = 20;
    const CGAL::FT other_speed = 18.7939;
# elif defined(CGAL_SSE_MERGING_WEIGHTS)
    const CGAL::FT x_speed = 20;
    const CGAL::FT y_speed = 20;
    const CGAL::FT z_speed = 20;
    const CGAL::FT other_speed = 19.8777;
# elif defined(CGAL_SSE_PERFORMANCE_WEIGHTS)
    const CGAL::FT x_speed = 5;
    const CGAL::FT y_speed = 5;
    const CGAL::FT z_speed = 2;
    const CGAL::FT other_speed = 5;
# else
#  error
# endif

    std::size_t fi = 0;
    it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        CGAL::FT speed = other_speed;
        const auto pl = facet->plane();
        const auto normal = KernelFactory::createVector3(pl);
        std::cout << "SP X " << CGAL::scalar_product(*normal, Vector3(1,0,0)) << std::endl;
        std::cout << "SP Y " << CGAL::scalar_product(*normal, Vector3(0,1,0)) << std::endl;
        std::cout << "SP Z " << CGAL::scalar_product(*normal, Vector3(0,0,1)) << std::endl;
        if(CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector3(1,0,0))) - 1) < 1e-3)
          speed = x_speed;
        if(CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector3(0,1,0))) - 1) < 1e-3)
          speed = y_speed;
        if(CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector3(0,0,1))) - 1) < 1e-3)
          speed = z_speed;

        SkelFacetDataSPtr data = SkelFacetData::create(facet);
        data->setSpeed(speed);
        std::cout << "speed to " << speed << std::endl;

        // for visualization, use color_ply_inputs.cpp to compile it into a single colored ply file
        std::vector<CGAL::EPICK::Point_3> points;
        std::vector<std::vector<std::size_t> > faces;

        CGAL::Cartesian_converter<CGAL::K, CGAL::EPICK> to_epick;

        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
          points.push_back(to_epick(*((*it_v++)->getPoint())));
        }

        std::vector<std::size_t> f(points.size());
        std::iota(f.begin(), f.end(), 0);
        faces.push_back(f);

        std::cout << points.size() << std::endl;
        for(auto p : points)
          std::cout << "  " << p << std::endl;

        std::cout << faces.size() << std::endl;
        std::cout << faces[0].size() << std::endl;
        for(std::size_t i : faces[0])
          std::cout << "  " << i << std::endl;

        ++fi;
    }
#endif

    DEBUG_VAL("Using " << vertex_splitter_->toString() << " to initialize polyhedron.");
    if (init(polyhedron)) {

        if (controller_) {
            controller_->wait();
        }

        CGAL::FT offset = 0.0;
        CGAL::FT offset_prev = 0.0;
        CGAL::FT offset_next = 0.0;

        db::_3d::OBJFile::save("results/init_post.obj", polyhedron, false /*do not triangulate*/);

        for(;;) {
            static int event_id = 0;

            std::cout << " =========== ITERATION #" << event_id << " AT OFFSET " << offset << std::endl;

            // db::_3d::OBJFile::save("results/input_" + std::to_string(event_id) + ".obj", polyhedron, false /*do triangulate*/);
            // db::_3d::OBJFile::save("results/input_" + std::to_string(event_id) + "_triangulated.obj", polyhedron);

            PQ queue;
            collectEvents(polyhedron, offset, queue);

            // @debug +
            {
                std::cout << "------------------------------" << std::endl;
                std::cout << "--- Event queue (" << event_id << ") ---" << std::endl;
                std::cout << "------------------------------" << std::endl;
                PQ duplicate_queue = queue;
                while (!duplicate_queue.empty()) {
                    AbstractEventSPtr event = duplicate_queue.top();
                    std::cout << event->toString() << std::endl;
                    duplicate_queue.pop();
                }
                std::cout << "Saves:";
                for (CGAL::FT save_offset : save_offsets_) {
                    std::cout << " " << save_offset;
                }
                std::cout << std::endl;
                std::cout << "-------------------" << std::endl;
                std::cout << "-------------------" << std::endl;
            }
            // @ debug -

            if (queue.empty()) {
                break; // we are done
            }

            // treat the next event
            AbstractEventSPtr event;
            bool doSave;
            bool simultaneousEvents;

            if (!save_offsets_.empty()) {
                CGAL::FT next_save_offset = save_offsets_.front();
                if (next_save_offset > queue.top()->getOffset()) { // save is strictly earlier
                    simultaneousEvents = false;
                    doSave = true;
                    event = SaveOffsetEvent::create(save_offsets_.front());
                } else {
                    std::tie(event, simultaneousEvents) = nextEvent(queue);
                    doSave = (next_save_offset == event->getOffset());
                }
            } else {
                std::tie(event, simultaneousEvents) = nextEvent(queue);
                doSave = false;
            }

            // @tmp +
            if (simultaneousEvents) {
              std::cerr << "Error: you should not be encountering simultaneous events these days" << std::endl;
              std::exit(1);
            }
            // @tmp -

            offset_next = event->getOffset();

            std::cout << " current offset: " << offset << "\n"
                      << " next offset: " << offset_next << " (type " << event->getType() << ")\n"
                      << " simultaneous? " << simultaneousEvents << "\n"
                      << " save? " << doSave << "\n"
                      << " in perturbed mode? " << usingTemporaryPerturbedMode_
                      << " @ " << simultaneousOffset_ << std::endl;
            std::cout << "current elapsed time: " << util::Timer::now() - t_start << std::endl;

#ifdef CGAL_SS3_USE_AUTOREF_FOR_ALL_EVENTS
            std::cout << "============== USING AUTOREF APPROACH =============" << std::endl;
            simultaneousEvents = false;
#endif

            // switch ON perturbation if needed
            if (simultaneousEvents) {
                if (doSave) {
                    std::cout << "!!!! SAVE EVENT AT SIMULTANEITY !!!!" << std::endl;
                    return handleSaveEventAtSimultaneity(polyhedron, offset, offset_next);
                }

#ifndef CGAL_SS3_USE_AUTOREF_FOR_SIMULTANEOUS_EVENTS
                if (offset_next == 0) {
                    std::cerr << "NYI: simultaneous event at offset=0" << std::endl;
                    std::exit(1);
                }

                if (usingTemporaryPerturbedMode_) {
                    std::cerr << "Error: met a simultaneous event in perturbed mode!" << std::endl;
                    std::exit(1); // nuke it
                } else {
                    std::cout << "simultaneous events ahead, at offset: " << offset_next << std::endl;
                    std::tie(polyhedron, offset) = enablePerturbedMode(polyhedron, offset, offset_next);
                    CGAL_assertion(offset > offset_next);
                    continue; // recompute events
                }
#endif
            }

            // switch OFF perturbation if needed
            if (usingTemporaryPerturbedMode_) {
                CGAL::FT delta = 1e-5; // @fixme don't hardcode this value
                // if the next event is far and we are in perturbed mode, disable perturbation
                if (simultaneousOffset_ - delta > offset_next) { // offsets are negative
                    std::tie(polyhedron, offset) = disablePerturbedMode(polyhedron, offset, offset_next);
                    CGAL_assertion(offset > offset_next);
                    continue; // recompute the polyhedron and its events
                }
            }

            std::cout << "\n-----------------------------------------------------" << std::endl;
            std::cout << "-- Event #" << event_id << " " << event->toString() << " --" << std::endl;
            std::cout << "-----------------------------------------------------\n" << std::endl;

            if (controller_) {
                controller_->wait();
            }

            Point3SPtr p_box_min = PolyhedronTransformation::boundingBoxMin(polyhedron);
            Point3SPtr p_box_max = PolyhedronTransformation::boundingBoxMax(polyhedron);

#ifdef CGAL_SS3_USE_AUTOREF_FOR_ALL_EVENTS
            // @fixme delta should depend on the next save and next constant events
            // @fixme delta should be small enough so that there are no other events
            CGAL::FT nudged_offset;
            std::tie(polyhedron, nudged_offset) = handleEventWithAutoref(event, offset, polyhedron);

            offset_prev = offset;
            offset = nudged_offset;
#else
# ifdef CGAL_SS3_USE_AUTOREF_FOR_SIMULTANEOUS_EVENTS
            if (simultaneousEvents) {
                if (doSave) {
                    return handleSaveEventAtSimultaneity(polyhedron, offset, offset_next);
                } else {
                    CGAL::FT nudged_offset;
                    std::tie(polyhedron, nudged_offset) = handleEventWithAutoref(event, offset, polyhedron);

                    offset_prev = offset;
                    offset = nudged_offset;
                }
            } else
# endif
            {
                offset_prev = offset;
                offset = event->getOffset();
                const CGAL::FT shift = offset - offset_prev;
                CGAL_warning(!usingTemporaryPerturbedMode_ || !is_zero(shift));
                const bool recompute_positions = (shift != 0);
                polyhedron = PolyhedronTransformation::shiftFacets(polyhedron, shift, recompute_positions);

                // will have degeneracies since we haven't treated the event yet
                db::_3d::OBJFile::save("results/shift_" + std::to_string(event_id) + ".obj", polyhedron, false /*do not triangulate*/);

                if (event->getType() == AbstractEvent::CONST_OFFSET_EVENT) {
#ifndef CGAL_SS3_NO_SKELETON_DS
                    event->setPolyhedronResult(polyhedron);
                    skel_result_->addEvent(event);
#endif
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
            }
#endif // CGAL_SS3_USE_AUTOREF_FOR_ALL_EVENTS

            std::cout << "-- Finished handling Event --" << std::endl;

            db::_3d::OBJFile::save("results/iter_" + std::to_string(event_id) + ".obj", polyhedron, false /*do triangulate*/);
            db::_3d::OBJFile::save("results/iter_" + std::to_string(event_id) + "_triangulated.obj", polyhedron);

            // @debug +
            DEBUG_PRINT("-- Degen count --");
            std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
            while (it_f != polyhedron->facets().end()) {
                FacetSPtr facet = *it_f++;
                auto it_v = facet->vertices().begin();
                Point3SPtr p0 = (*(it_v++))->getPoint();
                Point3SPtr p1 = (*(it_v++))->getPoint();
                Point3SPtr p2 = (*it_v)->getPoint();
                if (CGAL::collinear(*p0, *p1, *p2)) {
                    std::cout << *p0 << " " << *p1 << " " << *p2 << " is degen" << std::endl;
                }
            }
            // @debug -

            CGAL_assertion(polyhedron->isConsistent());
#ifndef CGAL_SS3_NO_SKELETON_DS
            CGAL_assertion(skel_result_->isConsistent());
#endif
            CGAL_assertion(p_box_min && p_box_max);
            CGAL_assertion(PolyhedronTransformation::isInsideBox(polyhedron, p_box_min, p_box_max));

            // this is tempting, but the mesh is usually not in a nice state here
            // CGAL_assertion(!SelfIntersection::hasSelfIntersectingSurface(polyhedron));

            if (controller_) {
                controller_->wait();
            }

            if (doSave) {
                savePolyhedron(polyhedron, offset,
                               true /*triangulate*/,
                               true /*convert to double*/,
                               false /*attempt untilting*/);

                save_offsets_.pop_front();
                if (save_offsets_.empty()) {
                    break; // @todo ought to be a config flag
                }
            }

            // Can't perturb to treat simultaneous events AND get a meaningful "SAVE" result
            if (simultaneousEvents && event->getType() == AbstractEvent::SAVE_OFFSET_EVENT) {
                std::cout << "WARNING: simultaneous event @ save time, not proceeding farther" << std::endl;
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
        std::cout << skel_result_->toString() << std::endl;
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

            // @fixme is division by "speed" correct?
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

bool SimpleStraightSkel::init(PolyhedronSPtr polyhedron) {
    std::cout << "Simple Skeleton -- Initialization" << std::endl;

    WriteLock l(polyhedron->mutex());
    bool result = true;

    std::list<VertexSPtr>::iterator it_v;

    it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (vertex->degree() < 3) {
            continue;
        }
        if (!vertex->hasData()) {
            SkelVertexData::create(vertex);
        }
        // leaving this regardless of the SKELETON_DS macro because of splitters
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
        std::cout << "Split #" << v2s_i++ << ": " << vertex->toString() << std::endl;

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

#ifndef CGAL_SS3_NO_SKELETON_DS
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
#endif

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        if (!facet->hasData()) {
            SkelFacetData::create(facet);
        }
    }

    CGAL_postcondition(skel_result_->isConsistent());

    polyhedron->initializeAllIDs();

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

std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::intersectionAndTimeOffsetPlanesWithCache(FacetSPtr facet_0,
                                                                                             FacetSPtr facet_1,
                                                                                             FacetSPtr facet_2,
                                                                                             FacetSPtr facet_3)
{
// #define CGAL_SS3_NO_CACHING
#ifdef CGAL_SS3_NO_CACHING
      Plane3SPtr plane_0 = basePlanes_.at(facet_0->getBasePlaneID());
      Plane3SPtr plane_1 = basePlanes_.at(facet_1->getBasePlaneID());
      Plane3SPtr plane_2 = basePlanes_.at(facet_2->getBasePlaneID());
      Plane3SPtr plane_3 = basePlanes_.at(facet_3->getBasePlaneID());
      CGAL::FT speed_0 = std::dynamic_pointer_cast<SkelFacetData>(facet_0->getData())->getSpeed();
      CGAL::FT speed_1 = std::dynamic_pointer_cast<SkelFacetData>(facet_1->getData())->getSpeed();
      CGAL::FT speed_2 = std::dynamic_pointer_cast<SkelFacetData>(facet_2->getData())->getSpeed();
      CGAL::FT speed_3 = std::dynamic_pointer_cast<SkelFacetData>(facet_3->getData())->getSpeed();
      return KernelWrapper::intersectionAndTimeOffsetPlanes(plane_0, speed_0, plane_1, speed_1, plane_2, speed_2, plane_3, speed_3);
#endif

    std::size_t ids[] = {facet_0->getBasePlaneID(), facet_1->getBasePlaneID(), facet_2->getBasePlaneID(), facet_3->getBasePlaneID()};
    std::sort(std::begin(ids), std::end(ids));
    auto canonical_ids = CGAL::make_array(ids[0], ids[1], ids[2], ids[3]);
    std::pair<Point3SPtr, CGAL::FT> dummy_intersection;
    auto res = intersectionCache_.emplace(canonical_ids, dummy_intersection);
    if (res.second) {
      // successful insertion, first time seeing it so actual computation is required
      std::pair<Point3SPtr, CGAL::FT>& point_and_time = res.first->second;

      Plane3SPtr plane_0 = basePlanes_.at(facet_0->getBasePlaneID()); // canonical order doesn't matter
      Plane3SPtr plane_1 = basePlanes_.at(facet_1->getBasePlaneID());
      Plane3SPtr plane_2 = basePlanes_.at(facet_2->getBasePlaneID());
      Plane3SPtr plane_3 = basePlanes_.at(facet_3->getBasePlaneID());

      // @fixme in theory, one would need to check if the data exists etc. etc.
      // for now, it's a good to check if speeds are properly set
      CGAL::FT speed_0 = std::dynamic_pointer_cast<SkelFacetData>(facet_0->getData())->getSpeed();
      CGAL::FT speed_1 = std::dynamic_pointer_cast<SkelFacetData>(facet_1->getData())->getSpeed();
      CGAL::FT speed_2 = std::dynamic_pointer_cast<SkelFacetData>(facet_2->getData())->getSpeed();
      CGAL::FT speed_3 = std::dynamic_pointer_cast<SkelFacetData>(facet_3->getData())->getSpeed();

      std::tie(point_and_time.first, point_and_time.second) = KernelWrapper::intersectionAndTimeOffsetPlanes(
        plane_0, speed_0, plane_1, speed_1, plane_2, speed_2, plane_3, speed_3);

      CGAL_assertion(*(point_and_time.first) == *(KernelWrapper::intersectionOffsetPlanes(
          facet_0->getPlane(), speed_0, facet_1->getPlane(), speed_1,
          facet_2->getPlane(), speed_2, facet_3->getPlane(), speed_3)));

      // std::cout << "recompute needed: " << canonical_ids[0] << " " << canonical_ids[1] << " " << canonical_ids[2] << " " << canonical_ids[3] << std::endl;
    } else {
      // std::cout << "used cache value: " << canonical_ids[0] << " " << canonical_ids[1] << " " << canonical_ids[2] << " " << canonical_ids[3] << std::endl;
    }

    return res.first->second;
}

std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::vanishesAtGeneric(FacetSPtr facet_0,
                                                                      FacetSPtr facet_1,
                                                                      FacetSPtr facet_2,
                                                                      FacetSPtr facet_3,
                                                                      CGAL::FT current_offset)
{
    Point3SPtr point = Point3SPtr();
    CGAL::FT event_offset;
    std::tie(point, event_offset) = intersectionAndTimeOffsetPlanesWithCache(facet_0, facet_1, facet_2, facet_3);

    if(point && event_offset <= current_offset) {
        // filtering with 'current_event' is done within the functions calling vanishesAt
        return { point, event_offset - current_offset };
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

std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::vanishesAt(EdgeSPtr edge,
                                                               CGAL::FT current_offset)
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

        util::ConfigurationSPtr config = util::Configuration::getInstance();
        bool usePerturbations = false;
        if (config->isLoaded()) {
            if ((config->contains("main", "rand_move_points") &&
                config->getBool("main", "rand_move_points")) ||
                (config->contains("main", "rand_move_points_when_degenerated") &&
                config->getBool("main", "rand_move_points_when_degenerated"))) {
                usePerturbations = true;
            }
        }

        Quadplane_parallelism parallelism;
        if (usePerturbations) {
            parallelism = PLANES_PARALLELISM_NONE;
        } else {
            parallelism = quadplane_parallelism(facetL, facetP, facetR, facetN);
        }

        if(parallelism != PLANES_PARALLELISM_NONE)
          std::cout << "parallelism = " << parallelism << std::endl;

        if (parallelism == PLANES_PARALLELISM_NONE) {
            // generic case
            std::tie(point, offset_event) = vanishesAtGeneric(facetL, facetP, facetR, facetN, current_offset);
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
// @speed, should be able to not solve the system but just exit early if the 4 planes are clearly not intersecting (diametral spheres around the edges of size something?)
std::pair<Point3SPtr, CGAL::FT>
SimpleStraightSkel::crashAt(EdgeSPtr edge_1, EdgeSPtr edge_2,
                            const CGAL::FT current_offset,
                            const CGAL::FT offset_to_farthest_event)
{
    FacetSPtr facet_l1 = edge_1->getFacetL();
    FacetSPtr facet_r1 = edge_1->getFacetR();
    FacetSPtr facet_l2 = edge_2->getFacetL();
    FacetSPtr facet_r2 = edge_2->getFacetR();
    Plane3SPtr plane_l1 = facet_l1->getPlane();
    Plane3SPtr plane_r1 = facet_r1->getPlane();
    Plane3SPtr plane_l2 = facet_l2->getPlane();
    Plane3SPtr plane_r2 = facet_r2->getPlane();
    CGAL::FT speed_l1 = std::dynamic_pointer_cast<SkelFacetData>(facet_l1->getData())->getSpeed();
    CGAL::FT speed_r1 = std::dynamic_pointer_cast<SkelFacetData>(facet_r1->getData())->getSpeed();
    CGAL::FT speed_l2 = std::dynamic_pointer_cast<SkelFacetData>(facet_l2->getData())->getSpeed();
    CGAL::FT speed_r2 = std::dynamic_pointer_cast<SkelFacetData>(facet_r2->getData())->getSpeed();

    util::ConfigurationSPtr config = util::Configuration::getInstance();
    bool usePerturbations = false;
    if (config->isLoaded()) {
        if ((config->contains("main", "rand_move_points") &&
             config->getBool("main", "rand_move_points")) ||
            (config->contains("main", "rand_move_points_when_degenerated") &&
             config->getBool("main", "rand_move_points_when_degenerated"))) {
            usePerturbations = true;
        }
    }

    Quadplane_parallelism parallelism;
    if (usePerturbations) {
        parallelism = PLANES_PARALLELISM_NONE;
    } else {
        parallelism = quadplane_parallelism(facet_l1, facet_r1, facet_l2, facet_r2);
    }

// #define CGAL_SS3_DEBUG_CRASH_AT
#ifdef CGAL_SS3_DEBUG_CRASH_AT
    std::cout << "-- Crash At\n    " << edge_1->toString() << "\n    " << edge_2->toString() << std::endl;

    std::cout << "Facet L1 = " << facet_l1->getID() << std::endl;
    std::cout << "Facet R1 = " << facet_r1->getID() << std::endl;
    std::cout << "Facet L2 = " << facet_l2->getID() << std::endl;
    std::cout << "Facet R2 = " << facet_r2->getID() << std::endl;

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
    std::tie(point, offset_event) = intersectionAndTimeOffsetPlanesWithCache(facet_l1, facet_r1, facet_l2, facet_r2);

    // @todo modify the rest of the code so that it uses 'offset_event' directly without
    // having to intermediary incremental offsets
    // Be careful about shifts, for example (shift_point(p, offset_event))
    offset_event = offset_event - current_offset;

    if (!point) {
        std::cerr << "Warning: no crashing?" << std::endl;
        return { };
    }

#ifdef CGAL_SS3_DEBUG_CRASH_AT
    std::cout << "Intersection: " << *point << " @ " << offset_event << std::endl;
#endif

    CGAL_assertion(*point == *(KernelWrapper::intersectionOffsetPlanes(
        plane_l1, speed_l1, plane_r1, speed_r1, plane_l2, speed_l2, plane_r2, speed_r2)));

    // Ignore the intersection if it is in the past
    // This is equivalent to the check that the point lives
    // on the inward half plane (split by the edge's supporting line) of the bisector plane
    //
    // @fixme fix the == behavior... assuming no simultaneous events means we don't care for ==
    if (!(offset_event < 0)) { // written that way so that it matches when using CTRL+F
#ifdef CGAL_SS3_DEBUG_CRASH_AT
        std::cout << "event is strictly in the past" << std::endl;
#endif
      return { };
    }

    // Ignore the event if it happens further the future compared to the soonest top event.
    if(offset_event <= offset_to_farthest_event) {
#ifdef CGAL_SS3_DEBUG_CRASH_AT
        std::cout << "event is too far in the future" << std::endl;
#endif
        return { };
    }

    // Check that the point is inside bounds
    FacetSPtr facet_1_src = getFacetSrc(edge_1);
    FacetSPtr facet_1_dst = getFacetDst(edge_1);
    FacetSPtr facet_2_src = getFacetSrc(edge_2);
    FacetSPtr facet_2_dst = getFacetDst(edge_2);

#ifdef CGAL_SS3_DEBUG_CRASH_AT
    std::cout << "Facet 1 SRC = " << facet_1_src->getID() << std::endl;
    std::cout << "Facet 1 DST = " << facet_1_dst->getID() << std::endl;
    std::cout << "Facet 2 SRC = " << facet_2_src->getID() << std::endl;
    std::cout << "Facet 2 DST = " << facet_2_dst->getID() << std::endl;

    std::cout << "Point: " << *point << std::endl;
#endif

// #define CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS
#if defined(CGAL_SS3_OLD_CODE_BOUND_CHECKS) || defined(CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS)
    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(edge_1->getData());
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(edge_2->getData());

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

    // this assumes positive weights
    CGAL_assertion(!(KernelWrapper::side(plane_l1, point) > 0 ||
                     KernelWrapper::side(plane_r1, point) > 0));
    CGAL_assertion(!(KernelWrapper::side(plane_l2, point) > 0 ||
                     KernelWrapper::side(plane_r2, point) > 0));

    // @todo this should be a predicate (oriented_side_of_event_point_wrt_bisectorC2)
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

    if (!(facet_1_src == facet_l2 || // @fixme are these facet checks required?
            facet_1_src == facet_r2)) {
        // std::cout << "-- check 1_src" << std::endl;
        // src is the target of the edge when in the right face
        if (!check_bisector(edge_1, facet_r1, rt1, facet_1_src, point)) {
            // std::cout << "reject #2b" << std::endl;
// #define CGAL_SS3_EXIT_ASAP
#ifdef CGAL_SS3_EXIT_ASAP
            return { };
#else
            reject_2b = true;
#endif
        }
    }

    if (!(facet_1_dst == facet_l2 ||
            facet_1_dst == facet_r2)) {
        // std::cout << "-- check 1_dst" << std::endl;
        if (!check_bisector(edge_1, facet_l1, lt1, facet_1_dst, point)) {
            // std::cout << "reject #3b" << std::endl;
#ifdef CGAL_SS3_EXIT_ASAP
            return { };
#else
            reject_3b = true;
#endif
        }
    }

    if (!(facet_2_src == facet_l1 ||
            facet_2_src == facet_r1)) {
        // std::cout << "-- check 2_src" << std::endl;
        if (!check_bisector(edge_2, facet_r2, rt2, facet_2_src, point)) {
            // std::cout << "reject #5b" << std::endl;
#ifdef CGAL_SS3_EXIT_ASAP
            return { };
#else
            reject_5b = true;
#endif
        }
    }

    if (!(facet_2_dst == facet_l1 ||
            facet_2_dst == facet_r1)) {
          // std::cout << "-- check 2_dst" << std::endl;
          if (!check_bisector(edge_2, facet_l2, lt2, facet_2_dst, point)) {
              // std::cout << "reject #6b" << std::endl;
#ifdef CGAL_SS3_EXIT_ASAP
              return { };
#else
              reject_6b = true;
#endif
        }
    }


    const bool reject_new = (reject_2b || reject_3b || reject_5b || reject_6b);
#endif

#ifdef CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS
# ifdef CGAL_SS3_EXIT_ASAP
#  error // can't use ASAP-exits if we want to compare
# endif
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

void SimpleStraightSkel::collectEdgeEvents(PolyhedronSPtr polyhedron,
                                           const CGAL::FT current_offset,
                                           CGAL::FT& current_offset_to_nearest_event,
                                           PQ& queue)
{
    std::cout << ">>> Collect Edge Events" << std::endl;

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
        std::tie(point, offset_event) = vanishesAt(edge, current_offset);
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
        if (offset_event < 0 && offset_event >= current_offset_to_nearest_event)
        {
            NodeSPtr node = Node::create(point);
            EdgeEventSPtr event = EdgeEvent::create(polyhedron);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
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
        std::tie(point, offset_event) = vanishesAt(edge, current_offset);
        if (!point) {
            continue;
        }
        if (offset_event < 0 && offset_event >= current_offset_to_nearest_event) {
            EdgeMergeEventSPtr event = EdgeMergeEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
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
        std::tie(point, offset_event) = vanishesAt(edge, current_offset);
        if (!point) {
            continue;
        }

        if ((KernelWrapper::side(edge->getFacetL()->plane(), point) > 0) ||
                KernelWrapper::side(edge->getFacetR()->plane(), point) > 0) {
            // triangle may not be a hole
            // after pierce event
            continue;
        }

        if (offset_event < 0 && offset_event >= current_offset_to_nearest_event) {
            TriangleEventSPtr event = TriangleEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
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
        std::tie(point, offset_event) = vanishesAt(edge, current_offset);
        if (!point) {
            continue;
        }
        if (offset_event < 0 && offset_event >= current_offset_to_nearest_event) {
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
        std::tie(point, offset_event) = vanishesAt(edge, current_offset);
        if (!point) {
            continue;
        }
        if (offset_event < 0 && offset_event >= current_offset_to_nearest_event) {
            DblTriangleEventSPtr event = DblTriangleEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
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
            std::tie(point, offset_event) = vanishesAt(edge, current_offset);
            if (!point) {
                continue;
            }

            std::cout << "potential tet event @ " << *point << " (t=" << offset_event << ")" << std::endl;

            if (offset_event < 0 && offset_event >= current_offset_to_nearest_event) {
                TetrahedronEventSPtr event = TetrahedronEvent::create(polyhedron);
                NodeSPtr node = Node::create(point);
                event->setNode(node);
                node->clear();
                node->setOffset(current_offset + offset_event);
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
            CGAL_assertion(vertex_1->getID() != std::size_t(-1) && vertex_2->getID() != std::size_t(-1));
            if (vertex_1->getID() > vertex_2->getID()) {
                continue;
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
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset, current_offset_to_nearest_event);
            if (!point)
                continue;

            VertexEventSPtr event = VertexEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();

#ifndef CGAL_SS3_NO_SKELETON_DS
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(->getData());
            node->addArc(data_2->getArc());
#endif
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            event->setVertex1(vertex_1);
            event->setVertex2(vertex_2);
            event->setFacet1(facet_1);
            event->setFacet2(facet_2);

            queue.push(event);

            std::cout << "Accepted vertex event: " << event->toString() << std::endl;
            std::cout << "point at zero x ? " << is_zero(point->x()) << std::endl;

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
            if (vertex_1->findEdge(vertex_2)) {
                // edge event
                continue;
            }
            if (isConvex(vertex_2)) {
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
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset, current_offset_to_nearest_event);
            if (!point)
                continue;

            FlipVertexEventSPtr event = FlipVertexEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();

#ifndef CGAL_SS3_NO_SKELETON_DS
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(vertex_2->getData());
            node->addArc(data_2->getArc());
#endif

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
//             // @fixme? it's not symmetric due to the convex check, but maybe it can be sometimes?
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
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset, current_offset_to_nearest_event);
            if (!point)
                continue;

            SurfaceEventSPtr event = SurfaceEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
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

            std::cout << "accepted surface event " << event->toString() << std::endl;

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
        CGAL_assertion(facet_1_src && facet_1_dst);
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
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset, current_offset_to_nearest_event);
            if (!point)
                continue;

            Point3SPtr e1so = PolyhedronTransformation::shiftPoint(edge_1->getVertexSrc(), offset_event);
            Point3SPtr e1to = PolyhedronTransformation::shiftPoint(edge_1->getVertexDst(), offset_event);
            if(*(e1so) != *(e1to)) {
                // not a polyhedron split (edge split?)
                continue;
            }

            PolyhedronSplitEventSPtr event = PolyhedronSplitEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
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

            std::cout << "accepted polyhedron split: " << event->toString() << std::endl;

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
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset, current_offset_to_nearest_event);
            if (!point)
                continue;

            SplitMergeEventSPtr event = SplitMergeEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();

#ifndef CGAL_SS3_NO_SKELETON_DS
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(vertex_2->getData());
            node->addArc(data_2->getArc());
#endif

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

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    unsigned int filtered_candidates = 0;
    unsigned int tested_candidates = 0;
#endif

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

        std::list<EdgeSPtr>::iterator it_e2 = it_e1; // starts at 'it_e1' => unique rep.
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

            // Filter if edges are too far
            //
            // @speed the best filtering would be to first have a loop over all combinations
            // to tighten the bound, and only then compute the crashAt, but below is
            // to get an easy speed-up
            //
            //
#if 0
            this does not filter enough because edges can still see very far

            // check shifted planes and edges orientations: if a shifted edge is still in the negative
            // side of the planes, we can 'continue'
            auto is_shifted_edge_definitely_on_negative_side_of_planes = [&current_offset_to_nearest_event](EdgeSPtr lhs, EdgeSPtr rhs) {
                // @speed note that here and in other places, we could speed up shifting computations
                // because there is a lot of redundant computations: for example here, the two
                // points are the intersection of 2 planes with the same shifted edge line
                // but we compute it from scratch for both
                Point3SPtr offset_e1so = PolyhedronTransformation::shiftPoint(lhs->getVertexSrc(), current_offset_to_nearest_event);
                Point3SPtr offset_e1to = PolyhedronTransformation::shiftPoint(lhs->getVertexDst(), current_offset_to_nearest_event);

                Plane3SPtr offset_pl0 = PolyhedronTransformation::shiftPlane(rhs->getFacetL(), current_offset_to_nearest_event);
                Plane3SPtr offset_pl1 = PolyhedronTransformation::shiftPlane(rhs->getFacetR(), current_offset_to_nearest_event);

                auto is_shifted_edge_definitely_on_negative_side_of_plane = [](Point3SPtr os, Point3SPtr od, Plane3SPtr opl, Plane3SPtr other_opl) {
                    auto orient_s = KernelWrapper::side(opl, os);
                    auto orient_d = KernelWrapper::side(opl, od);
                    if (orient_s * orient_d < 0) {
                        // the offset edge crosses opl; is it entirely on the negative side of the other plane?
                        //
                        // @todo we can still have the segment be on the union of the negative sides
                        // even if it crosses both planes
                        return (KernelWrapper::side(other_opl, os) < 0 && KernelWrapper::side(other_opl, od) < 0);
                    } else {
                        return (orient_s < 0);
                    }
                };

                return (is_shifted_edge_definitely_on_negative_side_of_plane(offset_e1so, offset_e1to, offset_pl0, offset_pl1) ||
                        is_shifted_edge_definitely_on_negative_side_of_plane(offset_e1so, offset_e1to, offset_pl1, offset_pl0));
            };

            if (is_shifted_edge_definitely_on_negative_side_of_planes(edge_1, edge_2) ||
                is_shifted_edge_definitely_on_negative_side_of_planes(edge_2, edge_1)) {
                ++filtered_candidates;
                continue;
            } else {
                // std::cout << "Checking possible edge split event\n\t"
                //           << edge_1->toString() << "\n\t"
                //           << edge_2->toString() << std::endl;
            }
#else
            // build 2 pairs of quads from each edge + shifted edge and check intersections
            //
            // a very important point is that the shifted edge could definitely be different due
            // to other events... but then another event will come first to modify the shifted edge!

            // let's just check bbox overlaps first
            Point3SPtr offset_e1so = PolyhedronTransformation::shiftPoint(edge_1->getVertexSrc(), current_offset_to_nearest_event);
            Point3SPtr offset_e1to = PolyhedronTransformation::shiftPoint(edge_1->getVertexDst(), current_offset_to_nearest_event);
            CGAL::Bbox_3 b1 = edge_1->getVertexSrc()->getPoint()->bbox();
            b1 += edge_1->getVertexDst()->getPoint()->bbox();
            b1 += offset_e1so->bbox();
            b1 += offset_e1to->bbox();

            Point3SPtr offset_e2so = PolyhedronTransformation::shiftPoint(edge_2->getVertexSrc(), current_offset_to_nearest_event);
            Point3SPtr offset_e2to = PolyhedronTransformation::shiftPoint(edge_2->getVertexDst(), current_offset_to_nearest_event);
            CGAL::Bbox_3 b2 = edge_2->getVertexSrc()->getPoint()->bbox();
            b2 += edge_2->getVertexDst()->getPoint()->bbox();
            b2 += offset_e2so->bbox();
            b2 += offset_e2to->bbox();
            if (!CGAL::do_overlap(b1, b2)) {
# ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
                ++filtered_candidates;
# endif
                // std::cout << "Filtered edge split candidates\n\t"
                //           << edge_1->toString() << "\n\t"
                //           << edge_2->toString() << std::endl;
                continue;
            } else {
# ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
                ++tested_candidates;
# endif
                // std::cout << "Checking possible edge split event\n\t"
                //           << edge_1->toString() << "\n\t"
                //           << edge_2->toString() << std::endl;
            }
#endif

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset, current_offset_to_nearest_event);
            if (!point) {
                continue;
            }

            // @fixme below needs to be double checked
#ifdef CGAL_SS3_ENFORCE_UNIQUE_EVENT_REPRESENTATIONS
            Point3SPtr e1so = PolyhedronTransformation::shiftPoint(edge_1->getVertexSrc(), offset_event);
            Point3SPtr e1to = PolyhedronTransformation::shiftPoint(edge_1->getVertexDst(), offset_event);
            if(*(e1so) == *(e1to)) {
                // polyhedron split
                continue;
            }
            Point3SPtr e2so = PolyhedronTransformation::shiftPoint(edge_2->getVertexSrc(), offset_event);
            Point3SPtr e2to = PolyhedronTransformation::shiftPoint(edge_2->getVertexDst(), offset_event);
            if(*(e2so) == *(e2to)) {
                // polyhedron split
                continue;
            }
#endif

            EdgeSplitEventSPtr event = EdgeSplitEvent::create(polyhedron);
            NodeSPtr node = Node::create(point);
            event->setNode(node);
            node->clear();
            node->setOffset(current_offset + offset_event);
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

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }

    timer.stop();
    std::cout << "Sought Edge Split Events in: " << timer.time() << std::endl;
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    std::cout << "  " << filtered_candidates << " filtered, " << tested_candidates << " tests" << std::endl;
#endif
}

void SimpleStraightSkel::collectPierceEvents(PolyhedronSPtr polyhedron,
                                             const CGAL::FT current_offset,
                                             CGAL::FT& current_offset_to_nearest_event,
                                             PQ& queue)
{
    std::cout << ">>> Collect Pierce Events" << std::endl;

    CGAL::Real_timer timer;
    timer.start();

#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    unsigned int pierce_vertex_counter = 0;
    unsigned int filtered_candidates = 0;
    unsigned int tested_candidates = 0;
#endif

    ReadLock l(polyhedron->mutex());

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (isLocked(vertex)) {
            continue;
        }

        // actual check
        if (isReflex(vertex)) {
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
            ++pierce_vertex_counter;
#endif

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

                // shrinking, so the vertex must be on the backside of the plane
                if (KernelWrapper::side(facet->plane(), vertex->getPoint()) > 0) {
                    continue;
                }

                // if the face is so far that even when shifting point and plane by the current
                // best lower bound on offset delta, the vertex has not crossed it yet, then we are done
                Point3SPtr shifted_pt = PolyhedronTransformation::shiftPoint(vertex, current_offset_to_nearest_event);
                Plane3SPtr shifted_plane = PolyhedronTransformation::shiftPlane(facet, current_offset_to_nearest_event);
                if (KernelWrapper::side(shifted_plane, shifted_pt) < 0) {
                    // std::cerr << "Filtering " << facet->getID() << " & " << *(vertex->getPoint()) << std::endl;
                    // std::cerr << "current_offset_to_nearest_event: " << current_offset_to_nearest_event << std::endl;
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
                    ++filtered_candidates;
#endif
                    continue;
                } else {
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
                    ++tested_candidates;
#endif
                }

                Point3SPtr point;
                CGAL::FT offset_event;

// #define CGAL_SS3_OLD_CODE_PIERCE_EVENT
#ifdef CGAL_SS3_OLD_CODE_PIERCE_EVENT
                SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
                ArcSPtr arc = data->getArc();

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
                        // facet too slow
                        continue;
                    }
                    dist_vertex = (distance * speed_vertex) / (speed_facet - speed_vertex);
                }

                offset_event = -dist_vertex / speed_vertex;
                point = KernelWrapper::offsetPoint(vertex->getPoint(), arc->getDirection(), dist_vertex);
                // std::cout << "old result: " << *point << " time: " << offset_event << std::endl;
#else
                CGAL_assertion(vertex->facets().size() >= 3);

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

                // @todo ugly squatting of unrelated function 'vanishesAtGeneric'
                std::tie(point, offset_event) = vanishesAtGeneric(facet, f0, f1, f2, current_offset);
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
#endif // CGAL_SS3_OLD_CODE_PIERCE_EVENT

                if (offset_event < 0 && offset_event >= current_offset_to_nearest_event) {
                    // Filter if the event point is on an edge (and a fortiori on a vertex)
                    // as it will be a different kind of event
                    FacetSPtr facet_offset = facet->clone();

                    Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->plane(), offset_event*speed);
                    facet_offset->setPlane(offset_plane);

                    // abusing the fact that vertices will have the same order in both facets
                    std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
                    std::list<VertexSPtr>::iterator it_v_offset = facet_offset->vertices().begin();
                    while (it_v != facet->vertices().end()) {
                        VertexSPtr vertex = *it_v++;
                        VertexSPtr offset_vertex = *it_v_offset++;
                        Point3SPtr point_offset = PolyhedronTransformation::shiftPoint(vertex, offset_event);
                        offset_vertex->setPoint(point_offset);
                    }

// #define CGAL_SS3_FILTER_PIERCE_EVENTS_AT_POP_TIME
# ifdef CGAL_SS3_FILTER_PIERCE_EVENTS_AT_POP_TIME
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
                    //
                    // Note that it'll avoid the infinite loop as such:
                    // - pierce event at t = t_0 is put in the queue
                    // - pierce event is top of the queue
                    // - polyhedron shifts
                    // - the pierce is rejected
                    // - pierce event is detected (again) and rejected (again)
                    // - pierce event is not put in the queue
                    // @todo needless computations on the penultimate step
                    if (offset_event == 0)
                    {
                        if (!SelfIntersection::isInsideWithRayShooting(point, facet)) {
                            std::cout << "Filtered T=0 pierce event" << std::endl;
                            continue;
                        }
                    }
# else

#  if 0
                    // @fixme actually bugged? (construction of the offset of the facet was broken when I put the other function...)
                    if (!IsLineInFacet(facet_offset, arc->line())) {
                        // std::cout << "IsLineInFacet rejects" << std::endl;
                        continue;
                    }
#  else
                    // Note that this result could be meaningless if the offset face
                    // is not a simple polygon. However, if it's not simple, then some event
                    // has happened before the pierce, and the pierce event - if whitelisted -
                    // would be checked again later, thus it's safe to call.
                    if (!SelfIntersection::isInsideWithRayShooting(point, facet_offset)) {
                        // std::cout << "isInsideWithRayShooting rejects" << std::endl;
                        continue;
                    }
#  endif
# endif
                    // @todo would be good if it could be merged with the function above...
                    bool boundary_rejection = false;
                    std::list<EdgeSPtr>::iterator it_fe = facet_offset->edges().begin();
                    while (it_fe != facet_offset->edges().end()) {
                        EdgeSPtr edge = *it_fe++;
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
                        std::cout << "Pierce event on boundary --> rejected" << std::endl;
                        continue;
                    }

                    PierceEventSPtr event = PierceEvent::create(polyhedron);
                    NodeSPtr node = Node::create(point);
                    event->setNode(node);
                    node->clear();
#ifndef CGAL_SS3_NO_SKELETON_DS
                    node->addArc(arc);
#endif
                    node->setOffset(current_offset + offset_event);
                    node->setPoint(point);
                    event->setFacet(facet);
                    event->setVertex(vertex);

                    // std::cout << "Accepted pierce event " << event->toString() << std::endl;

                    queue.push(event);

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
                    current_offset_to_nearest_event = offset_event;
#endif
                }
            }
        }
    }

    timer.stop();
    std::cout << "Sought Pierce Events in: " << timer.time() << std::endl;
#ifdef CGAL_SS3_PROFILE_FILTERING_MECHANISMS
    std::cout << "  " << pierce_vertex_counter << " reflex vertices" << std::endl;
    std::cout << "  " << filtered_candidates << " filtered, " << tested_candidates << " tests" << std::endl;
#endif
}

void SimpleStraightSkel::collectEvents(PolyhedronSPtr polyhedron,
                                       const CGAL::FT current_offset,
                                       PQ& queue)
{
    AbstractEventSPtr result = AbstractEventSPtr();
    if (!polyhedron || polyhedron->facets().size() == 0) {
        return;
    }

    CGAL::Real_timer timer;
    timer.start();

    // To handle save events at simultaneous events, save events are dealt with directly in the loop
    // if (!save_offsets_.empty()) {
    //     queue.push(SaveOffsetEvent::create(save_offsets_.front()));
    // }

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
    // - events that are in the past (i.e., offset > current_offset) (values are negative and decreasing!)
    // - events that are (stricly) later than the current next tentative offset (i.e., offset < curr_earliest_next_offset)
    CGAL::FT current_offset_to_nearest_event = queue.empty() ? - (std::numeric_limits<double>::max())
                                                             : (queue.top()->getOffset() - current_offset);

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

    // the next event types are particularly slow, so reduce the bound by doing them last
    // so other events lower the bound
    collectPierceEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectSurfaceEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);
    collectEdgeSplitEvents(polyhedron, current_offset, current_offset_to_nearest_event, queue);

    timer.stop();
    std::cout << "Sought All Events in: " << timer.time() << std::endl;
}

// Note that this takes care of the pop
// @fixme can handling an event reveal another event at the same times?
std::pair<AbstractEventSPtr, bool> SimpleStraightSkel::nextEvent(PQ& queue) {
    if (queue.empty()) {
        return { };
    }

    AbstractEventSPtr event = queue.top();
    CGAL_assertion(event->getType() != AbstractEvent::SAVE_OFFSET_EVENT);

    queue.pop();
    if (queue.empty()) {
        return { event, false /*not simultaneous*/};
    }

    AbstractEventSPtr next_event = queue.top();
    bool simultaneousEvents = (event->getOffset() == next_event->getOffset());
    if (event->getType() == AbstractEvent::CONST_OFFSET_EVENT && simultaneousEvents) {
        return nextEvent(queue); // the const event is popped and ignored
    }

    return { event, simultaneousEvents };
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

// belongs somewhere else, probably
PolyhedronSPtr SimpleStraightSkel::soup_to_polyhedron(const std::vector<Point3>& points,
                                                      const std::vector<std::vector<std::size_t> >& triangles,
                                                      const std::vector<Plane3SPtr>& planes,
                                                      const std::vector<CGAL::FT>& speeds) {
    PolyhedronSPtr result = Polyhedron::create();

    unsigned int vertex_id_new = 0;
    for (const Point3& p : points) {
        Point3SPtr point = KernelFactory::createPoint3(p);
        VertexSPtr vertex = Vertex::create(point);
        vertex->setID(vertex_id_new++);
        result->addVertex(vertex);
    }

    unsigned int facet_id_new = 0;
    for (std::size_t tid=0; tid<triangles.size(); ++tid) {
        const std::vector<std::size_t>& t = triangles[tid];
        CGAL_assertion(t.size() == 3); // @tmp could handle more, but no need right now
        VertexSPtr poly_vertices[3];
        for (unsigned int i = 0; i < 3; i++) {
            poly_vertices[i] = *(std::next(result->vertices().begin(), t[i]));
        }

        FacetSPtr facet = Facet::create(3, poly_vertices);
        facet->setID(facet_id_new++);
        facet->setPlane(planes[tid]);
        Triangle::create(facet, poly_vertices);

        data::_3d::skel::SkelFacetDataSPtr data;
        if (facet->hasData()) {
            data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
        } else {
            data = data::_3d::skel::SkelFacetData::create(facet);
        }
        data->setSpeed(speeds[tid]);

        result->addFacet(facet);
    }

    CGAL_postcondition(result->vertices().size() == points.size());
    CGAL_postcondition(result->facets().size() == triangles.size());

    return result;
}

std::pair<PolyhedronSPtr, CGAL::FT> SimpleStraightSkel::handleEventWithAutoref(AbstractEventSPtr event,
                                                                               CGAL::FT currentOffset,
                                                                               PolyhedronSPtr polyhedron) {
    namespace PMP = CGAL::Polygon_mesh_processing;

    using PID = std::size_t;
    using TID = std::size_t;
    using VID = std::size_t;

    std::cout << "handleEventWithAutoref()" << std::endl;

    CGAL_precondition(!usingTemporaryPerturbedMode_);

    // appendEventNode(event->getNode()); // @todo, if there is a point

    const CGAL::FT shift = event->getOffset() - currentOffset;
    const CGAL::FT delta = 0.2 * shift; // @fixme don't hardcode this value
    const CGAL::FT nudged_shift = shift + delta;
    std::cout << "currentOffset: " << currentOffset << std::endl;
    std::cout << "event->getOffset(): " << event->getOffset() << std::endl;
    std::cout << "shift: " << shift << std::endl;
    std::cout << "delta: " << delta << std::endl;
    std::cout << "total shift: " << nudged_shift << std::endl;

    static int autoref_event_id = -1;
    ++autoref_event_id;

    // {
    //     PolyhedronSPtr polyhedron_tmp = PolyhedronTransformation::shiftFacets(polyhedron, nudged_shift);
    //     db::_3d::OBJFile::save("results/autoref_event-shifted.obj", polyhedron_tmp, false /*do not triangulate*/);
    // }

    auto triangulate_quad_with_CDT2 = [](VertexSPtr vertex, VertexSPtr vertex_offset,
                                         VertexSPtr next_vertex, VertexSPtr next_vertex_offset,
                                         std::map<Point3, std::size_t>& pids,
                                         std::vector<Point3>& points,
                                         std::vector<std::vector<PID> >& triangles) {

        // @todo factorize CDT usages, but mind the tags
        using Itag = CGAL::Exact_intersections_tag;
        using PK = CGAL::Projection_traits_3<CGAL::K>;
        using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
        using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
        using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
        using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
        using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
        using PCDT_VH = PCDT::Vertex_handle;
        using PCDT_FH = PCDT::Face_handle;

        std::cout << "CDT2 QUAD WITH\n";
        std::cout << "vertex = " << *(vertex->getPoint()) << std::endl;
        std::cout << "vertex offset = " << *(vertex_offset->getPoint()) << std::endl;
        std::cout << "next vertex = " << *(next_vertex->getPoint()) << std::endl;
        std::cout << "next vertex offset = " << *(next_vertex_offset->getPoint()) << std::endl;

        CGAL_precondition(*(vertex->getPoint()) != *(vertex_offset->getPoint()));
        CGAL_precondition(*(vertex_offset->getPoint()) != *(next_vertex_offset->getPoint()));
        CGAL_precondition(*(vertex->getPoint()) != *(next_vertex_offset->getPoint()));

        // currently the orientation doesn't matter because we re-order for compatibility
        // with base faces later
        Vector3 n = CGAL::cross_product(Vector3(*(vertex->getPoint()), *(vertex_offset->getPoint())),
                                        Vector3(*(vertex->getPoint()), *(next_vertex_offset->getPoint())));

        CGAL_assertion(n != CGAL::NULL_VECTOR);

        PK projection_traits(n);
        PCDT pcdt(projection_traits);

        Point3SPtr vertex_pt = vertex->getPoint();
        PCDT_VH vertex_vh = pcdt.insert(*vertex_pt);
        vertex_vh->info() = vertex;

        Point3SPtr next_vertex_pt = next_vertex->getPoint();
        PCDT_VH next_vertex_vh = pcdt.insert(*next_vertex_pt);
        next_vertex_vh->info() = next_vertex;

        Point3SPtr vertex_offset_pt = vertex_offset->getPoint();
        PCDT_VH vertex_offset_vh = pcdt.insert(*vertex_offset_pt);
        vertex_offset_vh->info() = vertex_offset;

        Point3SPtr next_vertex_offset_pt = next_vertex_offset->getPoint();
        PCDT_VH next_vertex_offset_vh = pcdt.insert(*next_vertex_offset_pt);
        next_vertex_offset_vh->info() = next_vertex_offset;

        CGAL_precondition(CGAL::orientation(*vertex_pt, *vertex_offset_pt,
                                            *next_vertex_pt, *next_vertex_offset_pt) == CGAL::COPLANAR);

        pcdt.insert_constraint(vertex_vh, next_vertex_vh);
        pcdt.insert_constraint(next_vertex_vh, next_vertex_offset_vh);
        pcdt.insert_constraint(next_vertex_offset_vh, vertex_offset_vh);
        pcdt.insert_constraint(vertex_offset_vh, vertex_vh);

        std::unordered_map<PCDT_FH, bool> in_domain_map;
        boost::associative_property_map<std::unordered_map<PCDT_FH, bool> > in_domain(in_domain_map);
        CGAL::mark_domain_in_triangulation(pcdt, in_domain);

        auto get_pid = [&pids, &points] (PCDT_VH vh) -> PID {
            auto res = pids.emplace(vh->point(), points.size());
            if(res.second) { // first time seeing the vertex handle
                points.push_back(vh->point());
            }
            return res.first->second;
        };

        std::cout << "new triangles from CDT2" << std::endl;

        for(auto fh : pcdt.finite_face_handles()) {
            if(!get(in_domain, fh)) {
                continue;
            }

            triangles.push_back({get_pid(fh->vertex(0)),
                                 get_pid(fh->vertex(1)),
                                 get_pid(fh->vertex(2))});

            std::cout << points[triangles.back()[0]] << " " << points[triangles.back()[1]] << " " << points[triangles.back()[2]] << std::endl;
        }
    }; // lambda 'triangulate_quad_with_CDT2'

    enum class CDT2_Filtering {
        ODD_EVEN = 0, // in domain := odd nesting levels, see mark_domain_in_triangulation()
        NOT_OUT // in domain := everything that is not in the connected component incident to the infinite faces
    };

    auto triangulate_facet_with_CDT2 = [](FacetSPtr facet,
                                          CDT2_Filtering filtering_policy,
                                          std::map<Point3, PID>& pids,
                                          std::vector<Point3>& points,
                                          std::vector<std::vector<PID> >& triangles) {
        std::cout << "triangulate_facet_with_CDT2()" << std::endl;

        if (facet->edges().size() < 3) {
            std::cerr << "Warning: face with < 3 edges" << std::endl;
            return;
        } else {
            // @todo factorize CDT usages, but mind the tags
            using Itag = CGAL::Exact_intersections_tag; // since we do shift + epsilon, we could have intersections within a face
            using PK = CGAL::Projection_traits_3<CGAL::K>;
            using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
            using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
            using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
            using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
            using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
            using PCDT_VH = PCDT::Vertex_handle;
            using PCDT_FH = PCDT::Face_handle;

            Vector3SPtr n = KernelFactory::createVector3(facet->plane());
            CGAL_assertion(*n != CGAL::NULL_VECTOR);

            PK projection_traits(*n);
            PCDT pcdt(projection_traits);

            std::map<VertexSPtr, PCDT_VH> face_vhs; // might have multiple vertices at the same position

            std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
            while (it_v != facet->vertices().end()) {
                VertexSPtr vertex = *it_v++;
                auto res = face_vhs.emplace(vertex, PCDT_VH());
                if(res.second) // first time seeing this point
                {
                    PCDT_VH vh = pcdt.insert(*(vertex->getPoint()));
                    res.first->second = vh;
                }
            }

            auto ne = 0;
            std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
            while (it_e != facet->edges().end()) {
                EdgeSPtr edge = *it_e++;
                VertexSPtr v0 = edge->src(facet);
                VertexSPtr v1 = edge->dst(facet);
                std::cout << "CDT2 constraint: " << *(v0->getPoint()) << " || " << *(v1->getPoint()) << std::endl;

                // degenerate edges can happen e.g. at t=0 for vertices with deg > 3
                if(*(v0->getPoint()) == *(v1->getPoint())) {
                    std::cerr << "Warning: encountered degenerate edge" << std::endl;
                } else {
                    PCDT_VH vh0 = face_vhs.at(v0);
                    PCDT_VH vh1 = face_vhs.at(v1);

                    try {
                        pcdt.insert_constraint(vh0, vh1);
                    } catch(const typename PCDT::Intersection_of_constraints_exception&) {
                        std::cerr << "Error: Intersection of constraint w/ " << vh0->point() << " " << vh1->point() << std::endl;
                        DEBUG_VAR(facet->toString());
                        CGAL_assertion_msg(false, "Intersections are allowed, throw shouldn't happen");
                        std::exit(1);
                    }
                    ++ne;
                }
            }

            if(ne < 3) { // degenerate face
                std::cerr << "Warning: skipping degenerate face (ne < 3)" << std::endl;
                return;
            }

            // @fixme there can be self-intersections, so doing the domain thing is wrong?
#ifndef CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER
            std::unordered_map<PCDT_FH, bool> in_domain_map;
            boost::associative_property_map<std::unordered_map<PCDT_FH, bool> > in_domain(in_domain_map);

            if (filtering_policy == CDT2_Filtering::ODD_EVEN) {
                CGAL::mark_domain_in_triangulation(pcdt, in_domain);
            } else {
                for (PCDT_FH fh : pcdt.all_face_handles()) {
                    put(in_domain, fh, true);
                }

                PCDT_FH seed_fh = pcdt.infinite_vertex()->face();
                std::stack<PCDT_FH> to_explore;
                to_explore.push(seed_fh);
                while (!to_explore.empty()) {
                    PCDT_FH fh = to_explore.top();
                    to_explore.pop();
                    if (!get(in_domain, fh)) {
                        continue;
                    } else {
                        put(in_domain, fh, false);
                    }
                    for (std::size_t j=0; j<3; ++j) {
                        if (!fh->is_constrained(j)) {
                            to_explore.push(fh->neighbor(j));
                        }
                    }
                }
            }
#endif // CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER

            for(auto fh : pcdt.finite_face_handles()) {
            std::cout << "Face handle " << fh->vertex(0)->point() << " " << fh->vertex(1)->point() << " " << fh->vertex(2)->point() << " in domain? " << get(in_domain, fh) << std::endl;
#ifndef CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER
                if(!get(in_domain, fh)) {
                    continue;
                }
#endif // CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER

                // purge degenerate faces as we will things conformal with autoref anyway
                if (CGAL::collinear(fh->vertex(0)->point(),
                                    fh->vertex(1)->point(),
                                    fh->vertex(2)->point())) {
                    std::cerr << "Warning: degenerate face in autoref event handling" << std::endl;
                    continue;
                }

                auto get_pid = [&pids, &points] (PCDT_VH vh) -> PID {
                    auto res = pids.emplace(vh->point(), points.size());
                    if(res.second) { // first time seeing the vertex handle
                        points.push_back(vh->point());
                    }
                    return res.first->second;
                };

                // trying something to avoid inverted subparts of a face from surviving an event:
                // compute the winding number with angles (obviously terrible both in complexity AND robustness!)
#ifdef CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER

                // this doesn't work because we can't distinguish between doubly inverted faces pointing
                // in the wrong direction, and non-inverted faces
#               error

                Point3 centroid = CGAL::centroid(fh->vertex(0)->point(),
                                                 fh->vertex(1)->point(),
                                                 fh->vertex(2)->point());

                CGAL::FT cumulative_angle = 0;
                CGAL::internal::Evaluate<CGAL::FT> evaluate;

                std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
                while (it_e != facet->edges().end()) {
                    EdgeSPtr edge = *it_e++;
                    VertexSPtr v_src = edge->src(facet);
                    VertexSPtr v_dst = edge->dst(facet);
                    CGAL_assertion(*(v_src->getPoint()) != *(v_dst->getPoint())); // @todo handle this, if it can happen

                    Vector3 ln = CGAL::cross_product(Vector3(centroid, *(v_src->getPoint())),
                                                     Vector3(centroid, *(v_dst->getPoint())));
                    CGAL_assertion(ln != CGAL::NULL_VECTOR);
                    const CGAL::FT s = (CGAL::scalar_product(*n, ln) > 0) ? 1 : -1;

                    cumulative_angle += s * CGAL::approximate_angle(*(v_src->getPoint()), centroid, *(v_dst->getPoint()));
                    evaluate(cumulative_angle);

                    std::cout << centroid << " || " << *(v_src->getPoint()) << " || " << *(v_dst->getPoint()) << std::endl;
                    std::cout << "sign: " << s << std::endl;
                    std::cout << "angle: " << CGAL::approximate_angle(*(v_src->getPoint()), centroid, *(v_dst->getPoint())) << std::endl;
                }

                std::cout << "cumulative angle = " << cumulative_angle
                          << " (factor = " << cumulative_angle / 360 << ")" << std::endl;
                if (cumulative_angle > 180) { // && < 540
                    std::cout << "1" << std::endl;
                    triangles.push_back({get_pid(fh->vertex(0)),
                                         get_pid(fh->vertex(1)),
                                         get_pid(fh->vertex(2))});
                } else if (cumulative_angle < -180) {
                    std::cout << "2" << std::endl;
                    triangles.push_back({get_pid(fh->vertex(0)),
                                         get_pid(fh->vertex(2)), // invert
                                         get_pid(fh->vertex(1))});
                }
#else
                triangles.push_back({get_pid(fh->vertex(0)),
                                     get_pid(fh->vertex(1)),
                                     get_pid(fh->vertex(2))});
#endif // CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER
            }
        }
    }; // lambda 'triangulate_facet_with_CDT2'

    auto fill_edge_map = [](const std::vector<Point3>& points,
                            const std::vector<std::vector<PID> >& triangles,
                            std::vector<std::unordered_map<PID, std::vector<TID> > >& edge_map) {
        CGAL_precondition(edge_map.size() == points.size());

        // collect duplicated edges
        for (TID ti=0; ti<triangles.size(); ++ti) {
            for (std::size_t j=0; j<3; ++j) {
                std::pair<PID, PID> e_pids = CGAL::make_sorted_pair(triangles[ti][j],
                                                                    triangles[ti][(j+1)%3]);
                edge_map[e_pids.first][e_pids.second].push_back(ti);
            }
        }

        namespace pred = CGAL::Polygon_mesh_processing::Corefinement;

        for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
            for (auto& pid1_and_edges : edge_map[pid0]) {
                std::vector<TID>& inc_triangles = pid1_and_edges.second;

                if (inc_triangles.size() == 2) { // edge is only incident to a single SS3 face
                    continue;
                }

                const PID pid1 = pid1_and_edges.first;
                CGAL_assertion(pid0 != pid1);

                std::cout << "processing a non-manifold edge: " << points[pid0] << " -- " << points[pid1] << std::endl;

                auto get_third_point_id = [&triangles, pid0, pid1](TID tid) -> PID
                {
                    std::size_t third;

                    // need to be careful that the orientation of the edge might not match the orientation of the triangle
                    if (triangles[tid][0] == pid0 || triangles[tid][0] == pid1) {
                        if (triangles[tid][1] == pid0 || triangles[tid][1] == pid1) {
                            third = triangles[tid][2];
                        } else {
                            third = triangles[tid][1];
                        }
                    } else {
                        third = triangles[tid][0];
                    }

                    CGAL_postcondition(third != pid0 && third != pid1);
                    return third;
                };

                const Point3& ref_pt = points.at(get_third_point_id(inc_triangles[0]));
                auto less = [&ref_pt, &points, pid0, pid1, get_third_point_id](TID tid1, TID tid2)
                {
                    return pred::sorted_around_edge<CGAL::K>(points.at(pid0), points.at(pid1),
                                                             ref_pt,
                                                             points.at(get_third_point_id(tid1)),
                                                             points.at(get_third_point_id(tid2)));
                };

                std::sort(inc_triangles.begin()+1, inc_triangles.end(), less);

                // std::cout << "Around edge [" << pid0 << " " << pid1 << "], faces are sorted: ";
                // for(TID tid : inc_triangles)
                //   std::cout << " " << tid;
                // std::cout << std::endl;
            }
        }
    }; // lambda 'fill_edge_map'

    enum class Volume_orientation {
        UNKNOWN = 0,
        UNREACHABLE, // 1
        INWARD, // 2
        OUTWARD, // 3
        INCONSISTENT // 4
    };

    auto build_volume_CC = [](const TID seed_tid,
                              const VID CC_ID,
                              const bool start_from_inverted_face,
                              const bool ignore_facet_with_incompatible_orientations,
                              const bool ignore_dangling_outside_facets,
                              const std::vector<Point3>& points,
                              const std::vector<std::vector<PID> >& triangles,
                              const auto& edge_map,
                              auto& volume_CCs,
                              auto& volume_orientations,
                              auto& face_volume_IDs) {

        std::cout << "Building volume #" << CC_ID << " from seed face " << seed_tid << std::endl;

        volume_CCs.emplace_back();
        volume_orientations.emplace_back();

        std::stack<std::pair<TID, bool> > to_visit;
        to_visit.emplace(seed_tid, start_from_inverted_face);

        while (!to_visit.empty())
        {
            TID current_tid;
            bool invert_face;
            std::tie(current_tid, invert_face) = to_visit.top();
            to_visit.pop();

            std::cout << "At face " << current_tid << " [" << triangles[current_tid][0]
                                                   << ", " << triangles[current_tid][1]
                                                   << ", " << triangles[current_tid][2] << "], ";
            std::cout << "invert: " << invert_face << ", ";
            std::cout << "VIDS: " << face_volume_IDs[current_tid][0] << " " << face_volume_IDs[current_tid][1] << std::endl;

            std::size_t pos = invert_face ? 0 : 1;
            if (face_volume_IDs[current_tid][pos] == CC_ID) {
                // already visited this facet during the flooding of this volume's boundary
                continue;
            }

            CGAL_assertion(face_volume_IDs[current_tid][pos] == VID(-1)); // triangle should only be encountered once
            CGAL_warning(face_volume_IDs[current_tid][(pos+1)%2] != CC_ID); // Moebius shenanigans should be an instance of a bug

            volume_CCs.back().push_back(current_tid);

            // mark face as visited
            face_volume_IDs[current_tid][pos] = CC_ID;

            // flood through the edges
            for (int j=0; j<3; ++j) {
                std::pair<PID, PID> e_pids = CGAL::make_sorted_pair(triangles[current_tid][j],
                                                                    triangles[current_tid][(j+1)%3]);
                const std::vector<TID>& inc_triangles = edge_map.at(e_pids.first).at(e_pids.second);
                CGAL_assertion(!inc_triangles.empty());

                std::cout << "  ~~ Crossing edge [" << e_pids.first << ", " << e_pids.second << "]" << std::endl;
                std::cout << "    pos: " << points[e_pids.first] << " " << points[e_pids.second] << std::endl;

                // The faces are ordered CCW while looking from pid0.
                // So the walking while looking from [j] depends on whether [j] is pid0 or not
                int iter_direction = (e_pids.first == triangles[current_tid][j]) ? 1 : -1;
                std::cout << "    iter_direction = " << iter_direction << std::endl;

                // and it also depends on whether we are walking above or below the face
                iter_direction *= invert_face ? 1 : -1;
                std::cout << "    invert_face = " << invert_face << std::endl;

                TID next_tid = current_tid;
                for (;;) {
                    if (inc_triangles.size() == 1) {
                        std::cerr << "Warning: dangling triangle..." << std::endl;
                        std::cout << "    over the edge, the triangle is ITSELF " << current_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                        to_visit.emplace(current_tid, !invert_face);
                        break;
                    } else if (inc_triangles.size() == 2) {
                        // we should only be there once, meaning if we do not ignore orientations,
                        // then the faces MUST be compatible
                        CGAL_assertion(next_tid == current_tid);

                        next_tid = (inc_triangles[0] == current_tid) ? inc_triangles[1] : inc_triangles[0];
                        std::cout << "    over the edge, the triangle is TRIVIALLY " << next_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                        std::cout << "VIDS " << face_volume_IDs[next_tid][0] << " " << face_volume_IDs[next_tid][1] << std::endl;
                        CGAL_assertion(next_tid != current_tid);
                    } else {
                        // tricky part, now
                        auto tid_it = std::find(std::begin(inc_triangles), std::end(inc_triangles), next_tid /*updates on every iteration*/);
                        CGAL_assertion(tid_it != inc_triangles.end());

                        if (iter_direction == 1) { // CCW
                            std::cout << "    CCW walk" << std::endl;
                            auto next_it = std::next(tid_it);
                            next_tid = (next_it == inc_triangles.end()) ? inc_triangles[0] : *next_it;
                        } else { // CW
                            std::cout << "    CW walk" << std::endl;
                            next_tid = (tid_it == inc_triangles.begin()) ? inc_triangles.back() : *(std::prev(tid_it));
                        }

                        std::cout << "    over the edge, the triangle is " << next_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                        std::cout << "VIDS " << face_volume_IDs[next_tid][0] << " " << face_volume_IDs[next_tid][1] << std::endl;
                        CGAL_assertion(next_tid != current_tid);
                    }

                    // If the next face is incident to the outside (CC_ID == 0) on both sides, ignore it
                    if (ignore_dangling_outside_facets) {
                        // CC id #0 is the outside marker for volume building of facet prisms
                        bool is_dangling = (face_volume_IDs[next_tid][0] == 0 && face_volume_IDs[next_tid][1] == 0);
                        if (is_dangling) {
                            CGAL_assertion(inc_triangles.size() != 2);
                            std::cout << "  ignoring next TID because it is dangling" << std::endl;
                            continue;
                        }
                    }

                    // If the edge has the same direction in both faces (aka, the orientation changes),
                    // then we have to flip the direction of turning around the edge)
                    const auto j_it = std::find(std::begin(triangles[next_tid]),
                                              std::end(triangles[next_tid]),
                                              triangles[current_tid][j]);
                    CGAL_assertion(j_it != std::end(triangles[next_tid]));
                    const std::size_t pos = std::distance(std::begin(triangles[next_tid]), j_it);
                    CGAL_assertion(triangles[next_tid][pos] == triangles[current_tid][j]);

                    const bool flip_side = (triangles[next_tid][(pos+1)%3] == triangles[current_tid][(j+1)%3]);
                    std::cout << "    flipping? " << flip_side
                              << " (N: " << triangles[next_tid][(pos+1)%3] << " C: " << triangles[current_tid][(j+1)%3] << ")" << std::endl;
                    if (ignore_facet_with_incompatible_orientations && flip_side) {
                        if (inc_triangles.size() != 2) { // @tmp
                            std::cout << "  ignoring next TID because of incompatible orientation" << std::endl;
                            continue; // keep turning
                        } else {
                            std::cerr << "Warning: we should not have ignored incompatible face orientations, but there are only 2 incident triangles" << std::endl;
                        }
                    }

                    std::cout << "Final TID = " << next_tid << std::endl;
                    if (flip_side) {
                        CGAL_warning(!ignore_facet_with_incompatible_orientations);
                        volume_orientations[CC_ID] = Volume_orientation::INCONSISTENT;
                        to_visit.emplace(next_tid, !invert_face);
                    } else {
                        to_visit.emplace(next_tid, invert_face);
                    }

                    break;
                }
            }
        }
    }; // lambda 'build_volume_CC'

    // -- ACTUAL START

    // polyhedron offset, in a triangle soup form
    std::vector<Point3> points;
    std::vector<std::vector<PID> > triangles;

#ifdef CGAL_SS3_NO_INVERTED_FACE_FILTERING
    obsolete code, purge it
    polyhedron = PolyhedronTransformation::shiftFacets(polyhedron, nudged_shift);

    db::_3d::OBJFile::save("results/autoref_event-shifted.obj", polyhedron, false /*do not triangulate*/);

    polyhedron->initializeAllIDs();

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        triangulate_facet_with_CDT2(facet, triangles);
    }

    std::cout << points.size() << " points (pre auto-refine)" << std::endl;
    std::cout << triangles.size() << " triangles (pref auto-refine)" << std::endl;

    CGAL::IO::write_OFF("results/autoref_event-pre_autoref.off", points, triangles);

    PMP::autorefine_triangle_soup(points, triangles,
                                  CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag()));

    CGAL::IO::write_OFF("results/autoref_event-post_autoref.off", points, triangles);
#else // CGAL_SS3_NO_INVERTED_FACE_FILTERING

    // We need to offset each face, but we only output as correctly oriented the parts that
    // have never been inverted.
    //
    // To detect which parts are inverted, we first build the prism 3D polyhedron going from the base
    // face to the offset face (quadrangular faces to stitch those).
    // Then, we autorefine this triangle soup and identify the volume(s) incident to the base face.
    // Each top face that gets flagged is a face that hasn't been inverted.

    CGAL_assertion(triangles.empty());

    std::vector<int> is_unreachable_triangle; // because we want to swap, and vector<bool> is incompatible
    std::vector<Plane3SPtr> supporting_planes;
    std::vector<CGAL::FT> speeds;

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        std::cout << "\n-- Prism of Facet " << facet->getID() << std::endl;

        facet->makeFirstConvex();

        // Build the offset face
        FacetSPtr facet_offset = facet->clone();
        CGAL::FT speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
        Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->plane(), nudged_shift*speed);
        facet_offset->setPlane(offset_plane);

        CGAL_assertion(facet->vertices().size() == facet_offset->vertices().size());

        // abusing the fact that vertices will have the same order in both facets
        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        std::list<VertexSPtr>::iterator it_v_offset = facet_offset->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            VertexSPtr offset_vertex = *it_v_offset++;
            CGAL_assertion(vertex->getPoint() == offset_vertex->getPoint()); // haven't offset yet
            Point3SPtr point_offset = PolyhedronTransformation::shiftPoint(vertex, nudged_shift);
            offset_vertex->setPoint(point_offset);
            std::cout << *(vertex->getPoint()) << " offsets to " << *point_offset << std::endl;
        }

        // Triangulate the base & offset faces
        std::vector<Point3> prism_points;
        std::map<Point3, PID> pids;

        std::vector<std::vector<PID> > prism_base_triangles, prism_offset_triangles;
        // odd even is good enough since we know the face is sane
        triangulate_facet_with_CDT2(facet, CDT2_Filtering::ODD_EVEN, pids, prism_points, prism_base_triangles);
        // here, we want all the faces initially, and volume filtering knows which are His own
        triangulate_facet_with_CDT2(facet_offset, CDT2_Filtering::NOT_OUT, pids, prism_points, prism_offset_triangles);

        // invert the offset triangles such that the orientation is opposite that of the base,
        // which will be used later on because we want to enable traversing faces with opposite
        // direction as to avoid erroneous clipping
        for (std::vector<PID>& t : prism_offset_triangles) {
            std::swap(t[0], t[1]);
        }

        std::cout << prism_base_triangles.size() << " triangles (base)" << std::endl;
        std::cout << prism_offset_triangles.size() << " triangles (offset)" << std::endl;

        // lateral faces link the base and the offset triangles
        std::vector<std::vector<PID> > prism_lateral_triangles;

        std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
        std::list<EdgeSPtr>::iterator it_e_offset = facet_offset->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            EdgeSPtr edge_offset = *it_e_offset++;
            VertexSPtr vertex = edge->src(facet);
            VertexSPtr vertex_offset = edge_offset->src(facet_offset);
            VertexSPtr next_vertex = edge->dst(facet);
            VertexSPtr next_vertex_offset = edge_offset->dst(facet_offset);

            // whatever the diagonal, we might have a bow tie, and we don't want a fold
            // because the sort_around_edge() predicate cannot handle folds
            std::vector<std::vector<PID> > local_triangles;
            triangulate_quad_with_CDT2(vertex, vertex_offset, next_vertex, next_vertex_offset,
                                       pids, prism_points, local_triangles);

            // Re-orient triangles so that they match the boundaries of base & offset faces.
            // Orientations can be flipped because the edges of inner holes can be badly oriented.
            bool flip_orientations = false;

            // 'prism_base_triangles' are correctly oriented, and we need the lateral triangle
            // to have compatible orientations.
            // @todo do something less brute force
            for (const std::vector<PID>& lt : local_triangles) {
                for (std::size_t lj=0; lj<3; ++lj) {
                    for (const std::vector<PID>& bt : prism_base_triangles) {
                        for (std::size_t bj=0; bj<3; ++bj) {
                            if (prism_points[lt[lj]] == prism_points[bt[bj]] &&
                                prism_points[lt[(lj+1)%3]] == prism_points[bt[(bj+1)%3]]) {
                                std::cout << " incompatible orientations on " << prism_points[lt[lj]] << " " << prism_points[lt[(lj+1)%3]] << std::endl;
                                flip_orientations = true;
                                break;
                            }
                        }
                        if (flip_orientations) { break; }
                    }
                    if (flip_orientations) { break; }
                }
                if (flip_orientations) { break; }
            }

            if (flip_orientations) {
                for (std::vector<PID>& lt : local_triangles) {
                    std::swap(lt[0], lt[1]);
                }
            }

            prism_lateral_triangles.insert(std::end(prism_lateral_triangles),
                                           std::cbegin(local_triangles), std::cend(local_triangles));
        }

        std::cout << prism_lateral_triangles.size() << " triangles (lateral)" << std::endl;

        enum class Triangle_prism_location {
            UNKNOWN = 0,
            LATERAL, // = 1
            BASE, // = 2
            OFFSET, // = 3
        };

        // these are the ranges for the merge of base + offset + lateral
        std::vector<std::vector<PID> > prism_triangles;
        std::vector<Triangle_prism_location> triangle_locations;

        prism_triangles.insert(std::end(prism_triangles),
                               std::begin(prism_offset_triangles), std::end(prism_offset_triangles));
        triangle_locations.resize(prism_triangles.size(), Triangle_prism_location::OFFSET);

        CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) + "-top.off",
                            prism_points, prism_triangles);

        prism_triangles.insert(std::end(prism_triangles),
                               std::begin(prism_base_triangles), std::end(prism_base_triangles));
        triangle_locations.resize(prism_triangles.size(), Triangle_prism_location::BASE);

        // ESSENTIAL: lateral faces are guaranted to be planar quads because the vertices
        // move in the bisecting plane of the two planes incident to the edge incident to the two vertices.
        prism_triangles.insert(std::end(prism_triangles),
                               std::begin(prism_lateral_triangles), std::end(prism_lateral_triangles));
        triangle_locations.resize(prism_triangles.size(), Triangle_prism_location::LATERAL);

        // @tmp debug
        {
            std::cout << prism_points.size() << " points in shifting prism (pre autorefine)" << std::endl;
            std::cout << prism_triangles.size() << " triangles in shifting prism (pre autorefine)" << std::endl;

            std::cout << prism_points.size() << " prism points DEBUG" << std::endl;
            std::cout << prism_triangles.size() << " prism triangles DEBUG" << std::endl;

            for (const auto& p : prism_points) {
                std::cout << p << std::endl;
            }
            for (const auto& t : prism_triangles) {
                std::cout << t[0] << " " << t[1] << " " << t[2] << std::endl;
            }

            for (std::size_t i=0; i<prism_triangles.size(); ++i) {
                CGAL_assertion(triangle_locations[i] != Triangle_prism_location::UNKNOWN);
                std::cout << "triangle location[" << i << "] = " << int(triangle_locations[i]) << std::endl;
            }

            CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) + "-pre_autoref.off",
                                prism_points, prism_triangles);
        }

        // normally we shouldn't need to autoref the base & offset, but if autoref merges points (does it?),
        // it would mess up the other ranges
        //
        // And because we autorefine everything, we need a visitor because autoref might re-order
        // bottom and top faces even though they will not be refined (since we have already
        // computed a CDT2 and inserted intersections.)
        std::vector<Triangle_prism_location> updated_triangle_locations;
        Range_updating_autoref_visitor<Triangle_prism_location> autoref_visitor(triangle_locations, updated_triangle_locations);

        PMP::autorefine_triangle_soup(prism_points, prism_triangles,
                                      CGAL::parameters::visitor(autoref_visitor)
                                                       .concurrency_tag(CGAL::Parallel_if_available_tag()));

        triangle_locations = std::move(updated_triangle_locations);

        CGAL_assertion(triangle_locations.size() == prism_triangles.size());
        for (std::size_t i=0; i<prism_triangles.size(); ++i) {
            CGAL_assertion(triangle_locations[i] != Triangle_prism_location::UNKNOWN);
        }

        // because we could have overlaps and that's troublesome to determine volumes
        // @fixme does it need to be erase_all_duplicates(*TRUE*) sometimes?
        // @fixme the visitor needs to be smart about which face should be purged or not to have
        //        the correct value of reachable
        Range_updating_repair_PS_visitor<Triangle_prism_location> repair_ps_visitor(triangle_locations);
        PMP::merge_duplicate_polygons_in_polygon_soup(prism_points, prism_triangles,
                                                      CGAL::parameters::visitor(repair_ps_visitor)
                                                                       .erase_all_duplicates(false) /*keep one*/
                                                                       .require_same_orientation(false));

        std::cout << prism_points.size() << " points in shifting prism (post autorefine)" << std::endl;
        std::cout << prism_triangles.size() << " triangles in shifting prism (post autorefine)" << std::endl;

        CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) + "-post_autoref.off",
                            prism_points, prism_triangles);

        // ----

        std::cout << "Walking into the prism offset of facet " << facet->getID() << std::endl;

        // extremity --> other extremity --> range of incident faces
        std::vector<std::unordered_map<PID, std::vector<TID> > > edge_map(prism_points.size());
        fill_edge_map(prism_points, prism_triangles, edge_map);

        // identify volumes in the shifting faces soup, and tag faces of the volumes
        // that are incident to the base face(s)
        std::vector<std::vector<TID> > unused_volume_CCs; // range of range (volume) of triangle IDs
        std::vector<Volume_orientation> unused_volume_orientations;
        std::vector<std::array<VID, 2> > face_volume_IDs(prism_triangles.size(),
                                                         // [0] is down, [1] is up
                                                         std::array<VID, 2>{VID(-1), VID(-1)});

        // need to flood from all possible base faces because the faces is not necessarily a single CC
        //
        // Two independent loops for clarity and safety:
        // - Flood outside
        // - Flood inner volumes

        // first loop, outside
        for(std::size_t i=0; i<prism_triangles.size(); ++i) {
            if (triangle_locations[i] != Triangle_prism_location::BASE ||
                face_volume_IDs[i][1] != VID(-1)) {
                continue;
            }

            // outer flooding starts 'up' because the base face points out
            build_volume_CC(i /*tid*/, 0 /*CC ID*/, false /*do not use inverted orientation*/,
                            false /*do not traverse boundaries with wrong orientations*/,
                            false /*ignore dangling faces while walking*/,
                            prism_points, prism_triangles, edge_map,
                            unused_volume_CCs, unused_volume_orientations, face_volume_IDs);
        }

        // second loop, inner volumes
        for(std::size_t i=0; i<prism_triangles.size(); ++i) {
            if (triangle_locations[i] != Triangle_prism_location::BASE ||
                face_volume_IDs[i][0] != VID(-1)) {
                continue;
            }

            // inner flooding starts 'down' because the base face points out
            build_volume_CC(i /*tid*/, 1 /*CC ID*/, true /*use inverted orientation*/,
                            true /*traverse boundaries with wrong orientations*/,
                            true /*ignore dangling faces while walking*/,
                            prism_points, prism_triangles, edge_map,
                            unused_volume_CCs, unused_volume_orientations, face_volume_IDs);

            CGAL_assertion(face_volume_IDs[i][0] == VID(1)); // inward
            CGAL_assertion(face_volume_IDs[i][1] == VID(0)); // outward
        }


        // just for debugging
        // {
        //     for(VID cc_id=0; cc_id<unused_volume_CCs.size(); ++cc_id) {
        //         std::cout << "Volume #" << cc_id
        //                   << ", size: " << unused_volume_CCs[cc_id].size() << std::endl;

        //         std::vector<std::vector<PID> > CC_triangles;
        //         for (TID tid : unused_volume_CCs[cc_id]) {
        //             std::cout << tid << " ";
        //             CC_triangles.push_back(prism_triangles[tid]);
        //         }
        //         std::cout << std::endl;
        //         CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) +  "-volume_CC-" + std::to_string(cc_id) + ".off", prism_points, CC_triangles);
        //     }
        // }

        std::vector<std::vector<PID> > prism_triangle_contributions; // @tmp debug only
        std::vector<std::vector<PID> > prism_triangle_reachable_contributions; // @tmp debug only

        // Now, flag the offset faces that are not reachable
        for (std::size_t i=0; i<prism_triangles.size(); ++i) {
            CGAL_assertion(triangle_locations[i] != Triangle_prism_location::UNKNOWN);
            // std::cout << "triangle [" << i << "] position = " << prism_points[prism_triangles[i][0]] << " " << prism_points[prism_triangles[i][1]] << " " << prism_points[prism_triangles[i][2]] << std::endl;
            // std::cout << "triangle location[" << i << "] = " << int(triangle_locations[i]) << std::endl;

            if (triangle_locations[i] != Triangle_prism_location::OFFSET) {
                continue;
            }

            const bool is_dangling = (face_volume_IDs[i][0] == 0 /*exterior*/ &&
                                      face_volume_IDs[i][1] == 0 /*exterior*/);

            // we have walked inner volumes only, so check if the _lower_ side of offset faces is reached
            const bool is_unreached = (face_volume_IDs[i][0] == VID(-1));
            std::cout << "Triangle [" << i << "] reached? " << !is_unreached << std::endl;

#ifdef CGAL_SS3_DO_NOT_ADD_UNREACHED_TRIANGLES_TO_CONTRIBUTIONS
            this is too aggressive because we do not get volumes in the union of contributions
            and it is not clear if dangling faces must be extended or purged
            if (!is_unreached) {
#endif
                if(!is_dangling) {
                    // top's orientation was flipped to create a real prism volume, restore the real orientation
                    triangles.push_back({points.size() + prism_triangles[i][0],
                                         points.size() + prism_triangles[i][2],
                                         points.size() + prism_triangles[i][1]});

                    speeds.push_back(speed);
                    supporting_planes.push_back(offset_plane);
                    is_unreachable_triangle.push_back(is_unreached);
                    prism_triangle_contributions.push_back(prism_triangles[i]);
                }

                if (!is_unreached) {
                    prism_triangle_reachable_contributions.push_back(prism_triangles[i]);
                }

#ifdef CGAL_SS3_DO_NOT_ADD_UNREACHED_TRIANGLES_TO_CONTRIBUTIONS
            }
#endif
        }

        CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) + "-contributions.off",
                            prism_points, prism_triangle_contributions);

        CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) + "-reachable_contributions.off",
                            prism_points, prism_triangle_reachable_contributions);

        points.insert(std::end(points), std::cbegin(prism_points), std::cend(prism_points));

        std::cout << "Now, " << triangles.size() << " triangles in the polyhedron offset" << std::endl;
    }

    std::cout << " == MAIN POLYHEDRON OFFSET AUTOREF / VOLUME CHECKS ==" << std::endl;

    CGAL_postcondition(triangles.size() == speeds.size());
    CGAL_postcondition(triangles.size() == supporting_planes.size());
    CGAL_postcondition(triangles.size() == is_unreachable_triangle.size());
    for (std::size_t i=0; i<triangles.size(); ++i) {
        std::cout << "triangle [" << i << "] positions = " << points[triangles[i][0]] << " " << points[triangles[i][1]] << " " << points[triangles[i][2]] << std::endl;
        std::cout << "is_unreached[" << i << "] = " << is_unreachable_triangle[i] << " (before autoref)" << std::endl;
    }

    std::cout << points.size() << " points (pre auto-refine)" << std::endl;
    std::cout << triangles.size() << " triangles (pref auto-refine)" << std::endl;

    CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-pre_autoref.off", points, triangles);

    std::vector<int> updated_is_unreachable_triangle;
    std::vector<Plane3SPtr> updated_supporting_planes;
    std::vector<CGAL::FT> updated_speeds;
    Range_updating_autoref_visitor<int> bbvis(is_unreachable_triangle, updated_is_unreachable_triangle);
    Range_updating_autoref_visitor<Plane3SPtr, decltype(bbvis)> bvis(supporting_planes, updated_supporting_planes, bbvis);
    Range_updating_autoref_visitor<CGAL::FT, decltype(bvis)> vis(speeds, updated_speeds, bvis);

    PMP::autorefine_triangle_soup(points, triangles,
                                  CGAL::parameters::visitor(vis)
                                                   .concurrency_tag(CGAL::Parallel_if_available_tag()));

    is_unreachable_triangle = std::move(updated_is_unreachable_triangle);
    supporting_planes = std::move(updated_supporting_planes);
    speeds = std::move(updated_speeds);

    Range_updating_repair_PS_visitor<int> bbvis2(is_unreachable_triangle);
    Range_updating_repair_PS_visitor<Plane3SPtr, decltype(bbvis2)> bvis2(supporting_planes, bbvis2);
    Range_updating_repair_PS_visitor<CGAL::FT, decltype(bvis2)> vis2(speeds, bvis2);

    PMP::merge_duplicate_polygons_in_polygon_soup(points, triangles,
                                                  CGAL::parameters::visitor(vis2)
                                                                   .erase_all_duplicates(false) /*keep one*/
                                                                   .require_same_orientation(false));

    CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-post_autoref.off", points, triangles);
#endif // CGAL_SS3_NO_INVERTED_FACE_FILTERING

    std::cout << points.size() << " points (post auto-refine)" << std::endl;
    std::cout << triangles.size() << " triangles (post auto-refine)" << std::endl;

    CGAL_postcondition(triangles.size() == speeds.size());
    CGAL_postcondition(triangles.size() == supporting_planes.size());
    CGAL_postcondition(triangles.size() == is_unreachable_triangle.size());
    for (std::size_t i=0; i<triangles.size(); ++i) {
        std::cout << "triangle [" << i << "] position = " << points[triangles[i][0]] << " " << points[triangles[i][1]] << " " << points[triangles[i][2]] << std::endl;
        std::cout << "is_unreached[" << i << "] = " << is_unreachable_triangle[i] << " (after autoref)" << std::endl;
    }

    // ----------
    // Now, we have built a triangulation of all the offset faces of the polyhedron,
    // with proper orientation and flags
    //
    // What's next is identifying the volumes within this offset polyhedron, and getting rid
    // of the volumes that are not interesting

    // extremity --> other extremity --> range of incident faces
    std::vector<std::unordered_map<PID, std::vector<TID> > > edge_map(points.size());
    fill_edge_map(points, triangles, edge_map);

    // split the volumetric partition into volume connected components
    std::vector<std::vector<TID> > volume_CCs; // range of range (volume) of triangle IDs
    std::vector<Volume_orientation> volume_orientations;
    std::vector<std::array<VID, 2> > face_volume_IDs(triangles.size(),
                                                     // [0] is down, [1] is up
                                                     std::array<VID, 2>{VID(-1), VID(-1)});

#ifdef CGAL_SS3_FILTER_VOLUMES_WITH_ONLY_REACHABLE_FACES
    VID CC_ID = 0;
    for (std::size_t i=0; i<triangles.size(); ++i) {
        if (face_volume_IDs[i][0] == VID(-1) && // 'down' as we have inverted the offset triangles
            !is_unreachable_triangle[i]) {
            build_volume_CC(i /*seed ID*/, CC_ID++, true /*inverted orientation (down)*/,
                            false /*traverse boundaries with wrong orientations*/,
                            false /*ignore dangling faces while walking*/,
                            points, triangles, edge_map,
                            volume_CCs, volume_orientations, face_volume_IDs);
        }
    }

    for(VID cc_id=0; cc_id<volume_CCs.size(); ++cc_id) {
        for (TID tid : volume_CCs[cc_id]) {
            if (is_unreachable_triangle[tid]) {
                volume_orientations[cc_id] = Volume_orientation::UNREACHABLE;
                break;
            }
        }
    }

    // debug-only loop
    for(VID cc_id=0; cc_id<volume_CCs.size(); ++cc_id) {
        std::cout << "Volume #" << cc_id
                  << ", size: " << volume_CCs[cc_id].size()
                  << ", flag: " << int(volume_orientations[cc_id]) << std::endl;

        std::vector<std::vector<PID> > CC_triangles;
        for (TID tid : volume_CCs[cc_id]) {
            CC_triangles.push_back(triangles[tid]);
        }

        CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-volume_CC-" + std::to_string(cc_id) + ".off", points, CC_triangles);
    }

    if (true || event->getType() == AbstractEvent::SURFACE_EVENT) {
        for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
            for (auto& pid1_and_edges : edge_map[pid0]) {
                const PID pid1 = pid1_and_edges.first;

                std::vector<TID>& inc_triangles = pid1_and_edges.second;
                if (inc_triangles.size() != 2) {
                    continue;
                }

                // convex | concave angle
                TID tid0 = inc_triangles[0], tid1 = inc_triangles[1];

                auto third_pos = [&triangles](const TID tid, const TID other) -> std::size_t {
                    for (std::size_t i=0; i<3; ++i)
                        if (std::find(std::cbegin(triangles[other]), std::cend(triangles[other]),
                                      triangles[tid][i]) == std::cend(triangles[other]))
                            return i;
                    return -1; // In case no such element is found
                };

                std::size_t third_pos0 = third_pos(tid0, tid1), third_pos1 = third_pos(tid1, tid0);
                PID pid = triangles[tid0][(third_pos0 + 1)%3];
                PID qid = triangles[tid0][(third_pos0 + 2)%3];
                PID rid = triangles[tid0][third_pos0];
                PID sid = triangles[tid1][third_pos1];
                CGAL_assertion((pid != qid) && (pid != rid) && (pid != sid) && (qid != rid) && (qid != sid) && (rid != sid));
                CGAL::FT dh = CGAL::approximate_dihedral_angle(points[pid], points[qid], points[rid], points[sid]);

                std::cout << "checking checking" << std::endl;
                if (is_unreachable_triangle[tid0] != is_unreachable_triangle[tid1]) {
                    std::cerr << "deg 2 edge with reachable/unreachable triangles @ edge ";
                    std::cerr << points[pid0] << " -- " << points[pid1] << std::endl;
                }
            }
        }
    }

    // end debug

    // Purge faces incident to volumes either unreachable or inconsistent
    auto removal_predicate = [&](const TID tid) {
        Volume_orientation down_flag = (face_volume_IDs[tid][0] != VID(-1)) ? volume_orientations.at(face_volume_IDs[tid][0]) : Volume_orientation::UNREACHABLE;
        Volume_orientation up_flag = (face_volume_IDs[tid][1] != VID(-1)) ? volume_orientations.at(face_volume_IDs[tid][1]) : Volume_orientation::UNREACHABLE;

        const bool down_rejected = (down_flag == Volume_orientation::UNREACHABLE || down_flag == Volume_orientation::INCONSISTENT);
        const bool up_rejected = (up_flag == Volume_orientation::UNREACHABLE || up_flag == Volume_orientation::INCONSISTENT);

        return (down_rejected && up_rejected);
    };
#else
    this is bad because if there are holes, the volume will be negative
    but we do not want to discard it

    VID CC_ID = 0;
    for (std::size_t i=0; i<triangles.size(); ++i) {
        std::cout << "face " << i << " VIDS: " << face_volume_IDs[i][0] << " " << face_volume_IDs[i][1] << std::endl;
        if (face_volume_IDs[i][0] == VID(-1)) {
            build_volume_CC(i /*seed ID*/, CC_ID++, true /*inverted orientation (down)*/,
                            false /*traverse boundaries with wrong orientations*/,
                            false /*ignore dangling faces while walking*/,
                            points, triangles, edge_map,
                            volume_CCs, volume_orientations, face_volume_IDs);
        }

        if (face_volume_IDs[i][1] == VID(-1)) {
            build_volume_CC(i /*seed ID*/, CC_ID++, false /*normal orientation (up)*/,
                            false /*traverse boundaries with wrong orientations*/,
                            false /*ignore dangling faces while walking*/,
                            points, triangles, edge_map,
                            volume_CCs, volume_orientations, face_volume_IDs);
        }


        CGAL_postcondition(face_volume_IDs[i][0] != VID(-1) && face_volume_IDs[i][1] != VID(-1));
    }

    // @fixme some volumes are computed twice (take, e.g., a sole tetrahedron)
    std::cout << CC_ID << " volumes in the offset polyhedron" << std::endl;

    // ----------

    // determine the orientation of consistent CCs
    for(VID cc_id=0; cc_id<volume_CCs.size(); ++cc_id) {
        // just for debugging
        {
            std::cout << "Volume #" << cc_id
                      << ", size: " << volume_CCs[cc_id].size()
                      << ", orientation: " << int(volume_orientations[cc_id]) << std::endl;

            std::vector<std::vector<PID> > CC_triangles;
            for (TID tid : volume_CCs[cc_id]) {
                CC_triangles.push_back(triangles[tid]);
            }
            CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-volume_CC-" + std::to_string(cc_id) + ".off", points, CC_triangles);
        }

        if (volume_orientations[cc_id] != Volume_orientation::UNKNOWN) {
            continue;
        }

#ifndef CGAL_SS3_NO_INVERTED_FACE_FILTERING
        bool reject_due_to_incident_unreachable_triangle = false;
        for (TID tid : volume_CCs[cc_id]) {
#ifdef CGAL_SS3_DO_NOT_ADD_UNREACHED_TRIANGLES_TO_CONTRIBUTIONS
            CGAL_assertion(!is_unreachable_triangle[tid]);
# else
            if (is_unreachable_triangle[tid]) {
                std::cout << "Rejecting because " << tid << " is NOT reachable" << std::endl;
                volume_orientations[cc_id] = Volume_orientation::UNREACHABLE;
                reject_due_to_incident_unreachable_triangle = true;
            }

            if (reject_due_to_incident_unreachable_triangle) {
                break;
            }
# endif
        }

        if (reject_due_to_incident_unreachable_triangle) {
            continue;
        }
#endif

        // the volume CC has consistent orientations, but we need to know if it's inward or outward oriented

        // @fixme for now, evaluating volumes because finding the equivalent of:
        //    target(next(min_slope_he, pmesh), pmesh))
        //    target(next(opposite(min_slope_he, pmesh)
        // in:
        //    PMP::internal::is_outward_oriented
        // is tedious
#if 0
        // find the extremum point
        PID ext_pid = -1;

        for (TID tid : volume_CCs[cc_id]) {
            for (std::size_t j=0; j<3; ++j) {
                PID pid = triangles[tid][j];
                if (ext_pid == PID(-1) ||
                    points[ext_pid].z() < points[pid].z()) {
                    ext_pid = pid;
                }
            }
        }

        CGAL_postcondition(ext_pid != PID(-1));

        // find the edges adjacent to the extrem point
        std::unordered_set<PID> incident_vertices;
        for (TID tid : volume_CCs[cc_id]) {
            for (std::size_t j=0; j<3; ++j) {
                if (triangles[tid][j] == ext_pid) {
                  incident_vertices.insert(triangles[tid][(j+1)%3]);
                  incident_vertices.insert(triangles[tid][(j+2)%3]); // should be sufficient to only do one
                }
            }
        }

        // now we compare slopes, but we need the third points on either side of the edge.......
        ... @todo
#else
        // outward = positive volume, inward = negative volume
        Point3 origin(0,0,0);
        CGAL::FT volume = 0;
        CGAL::internal::Evaluate<CGAL::FT> evaluate;
        for (TID tid : volume_CCs[cc_id]) {
            volume += CGAL::volume(origin,
                                   points[triangles[tid][0]],
                                   points[triangles[tid][1]],
                                   points[triangles[tid][2]]);
            evaluate(volume);
        }

        volume_orientations[cc_id] = (volume > 0) ? Volume_orientation::OUTWARD
                                                  : Volume_orientation::INWARD;

        std::cout << "Volume #" << cc_id
                  << ", volume: " << volume
                  << ", orientation: " << int(volume_orientations[cc_id]) << std::endl;
#endif
    }

    // ----------

    // Purge faces that are not incident to an "interior" volume
    auto removal_predicate = [&](const TID tid) {
        const Volume_orientation first_orientation = volume_orientations[face_volume_IDs[tid][0]];
        const Volume_orientation second_orientation = volume_orientations[face_volume_IDs[tid][1]];

        CGAL_warning(!(volume_orientations[face_volume_IDs[tid][0]] == Volume_orientation::OUTWARD &&
                       volume_orientations[face_volume_IDs[tid][1]] == Volume_orientation::OUTWARD));
        CGAL_warning(!(volume_orientations[face_volume_IDs[tid][0]] == Volume_orientation::INWARD &&
                       volume_orientations[face_volume_IDs[tid][1]] == Volume_orientation::INWARD));

        // discard any face that is not incident to at least one non-rejected volume
        return (first_orientation != Volume_orientation::OUTWARD && second_orientation != Volume_orientation::OUTWARD);
    };
#endif // CGAL_SS3_FILTER_VOLUMES_WITH_ONLY_REACHABLE_FACES

    for (TID tid = 0; tid < triangles.size();) {
        if (removal_predicate(tid)) {
            std::size_t swap_pos = triangles.size() - 1;
            if (tid != swap_pos) {
                std::swap(face_volume_IDs[tid], face_volume_IDs[swap_pos]);
                std::swap(triangles[tid], triangles[swap_pos]);
                std::swap(supporting_planes[tid], supporting_planes[swap_pos]);
                std::swap(speeds[tid], speeds[swap_pos]);
            }
            triangles.pop_back();
            supporting_planes.pop_back();
            speeds.pop_back();
        } else {
            ++tid;
        }
    }

    std::cout << "final has " << triangles.size() << " triangles" << std::endl;

    PMP::remove_isolated_points_in_polygon_soup(points, triangles);

    CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-final_soup.off", points, triangles);

    // ----------

    polyhedron = soup_to_polyhedron(points, triangles, supporting_planes, speeds);
    CGAL_assertion(polyhedron && polyhedron->isConsistent());

    db::_3d::AbstractFile::mergeCoplanarFacets(polyhedron);
    db::_3d::AbstractFile::removeVerticesDegLt3(polyhedron);
    CGAL_postcondition(polyhedron && polyhedron->isConsistent());

    db::_3d::OBJFile::save("results/autoref_event-almost_final.obj", polyhedron, false /*do not triangulate*/);

    // @fixme this is very dangerous because now we will get a different normalization for the planes.
    // It would be better to re-use the planes
    PolyhedronTransformation::harmonizeFacetPlanes(polyhedron);
    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, 0.0);
    CGAL_postcondition(polyhedron && polyhedron->isConsistent());

    polyhedron->initializeAllIDs();

    db::_3d::OBJFile::save("results/autoref_event-final.obj", polyhedron, false /*do not triangulate*/);

#ifndef CGAL_SS3_NO_SKELETON_DS
    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
#endif

    return { polyhedron, event->getOffset() + delta };
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

    bool do_tests = (!isLocked(edge->getVertexSrc()) && !isLocked(edge->getVertexDst()) &&
                     !isLocked(vertices[0]) && !isLocked(vertices[1]) &&
                     !isLocked(vertices[2]) && !isLocked(vertices[3]));

    bool not_flipped_valid = false;
    if(do_tests)
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
            facets_clone[i]->setBasePlaneID(facets[i]->getBasePlaneID()); // @todo useless but just in case for now...
            facets_clone[i]->cachedPlane_ = facets[i]->cachedPlane_;
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

        if (polyhedron_no_flip_offset && !SelfIntersection::hasSelfIntersectingSurface(polyhedron_no_flip_offset)) {
            not_flipped_valid = true;
        }
    }

    std::cout << "not_flipped_valid = " << not_flipped_valid << std::endl;

    bool flipped_valid = false;
#if 0
    if (!do_tests || !not_flipped_valid) { // redundant but for clarity
        flipped_valid = true;
    } else
#endif
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
            facets_clone[i]->setBasePlaneID(facets[i]->getBasePlaneID());
            facets_clone[i]->cachedPlane_ = facets[i]->cachedPlane_; // @todo there should be a constructor/function "Facet::copyPropertiesAndData(otherFacet)"
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

        if (polyhedron_flipped_offset && !SelfIntersection::hasSelfIntersectingSurface(polyhedron_flipped_offset)) {
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
                facets_c[i]->setBasePlaneID(facets[i]->getBasePlaneID());
                facets_c[i]->cachedPlane_ = facets[i]->cachedPlane_;

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
    skel_result_->addEvent(event);
#endif
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

#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
    vertex_data->setNode(event->getNode());
    ArcSPtr arc = createArc(vertex);
    skel_result_->addArc(arc);
    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
#endif

#if 0 // using EPECK for now so we are exact
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
#endif
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

#ifndef CGAL_SS3_NO_SKELETON_DS
    SkelVertexDataSPtr data_offset = SkelVertexData::create(vertex_offset);
    data_offset->setNode(event->getNode());
    ArcSPtr arc = createArc(vertex_offset);
    skel_result_->addArc(arc);
    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
#endif
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

#ifndef CGAL_SS3_NO_SKELETON_DS
    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
#endif
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

#ifndef CGAL_SS3_NO_SKELETON_DS
    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
#endif
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

#ifndef CGAL_SS3_NO_SKELETON_DS
    event->setPolyhedronResult(polyhedron);
    skel_result_->addEvent(event);
#endif
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
    skel_result_->addEvent(event);
#endif
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
    skel_result_->addEvent(event);
#endif
}

void SimpleStraightSkel::handleSurfaceEvent(SurfaceEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());

    std::cout << "########################################" << std::endl;
    std::cout << "#######  Handle Surface Event  #########" << std::endl;
    std::cout << "########################################" << std::endl;

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    std::cout << "Node = " << *(node->getPoint()) << std::endl;
    std::cout << "Edge A = " << event->getEdge1()->toString() << std::endl;
    std::cout << "Edge B = " << event->getEdge2()->toString() << std::endl;

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
    skel_result_->addEvent(event);
#endif
}

void SimpleStraightSkel::handlePolyhedronSplitEvent(PolyhedronSplitEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());

    std::cout << "########################################" << std::endl;
    std::cout << "####  Handle Polyhedron Split Event  ###" << std::endl;
    std::cout << "########################################" << std::endl;

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    std::cout << "Node = " << *(node->getPoint()) << std::endl;
    std::cout << "Edge A = " << event->getEdge1()->toString() << std::endl;
    std::cout << "Edge B = " << event->getEdge2()->toString() << std::endl;

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
    skel_result_->addEvent(event);
#endif
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
    skel_result_->addEvent(event);
#endif

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
    skel_result_->addEvent(event);
#endif
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

#ifdef CGAL_SS3_FILTER_PIERCE_EVENTS_AT_POP_TIME
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
    boost::associative_property_map<std::unordered_map<PCDT_FH, bool> > in_domain(in_domain_map);
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
        std::cout << "Pierce Event rejected" << std::endl;
        event->setPolyhedronResult(polyhedron);
        // skel_result_->addEvent(event);
        return;
    }
#endif // CGAL_SS3_FILTER_PIERCE_EVENTS_AT_POP_TIME

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
    skel_result_->addEvent(event);
#endif
}


StraightSkeletonSPtr SimpleStraightSkel::getResult() const {
#ifndef CGAL_SS3_NO_SKELETON_DS
    std::cerr << "No skeleton to return, it was not built" << std::endl;
    std::exit(1);
#endif
    return this->skel_result_;
}

} }
