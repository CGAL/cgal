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
#include <limits>
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

SimpleStraightSkel::SimpleStraightSkel(PolyhedronSPtr polyhedron, ControllerSPtr controller, const std::list<double>& save_offsets) {
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

SimpleStraightSkelSPtr SimpleStraightSkel::create(PolyhedronSPtr polyhedron, ControllerSPtr controller, const std::list<double>& save_offsets) {
    return SimpleStraightSkelSPtr(new SimpleStraightSkel(polyhedron, controller, save_offsets));
}

void SimpleStraightSkel::initVertexSplitter() {
    use_fast_vertex_splitter_ = true;
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


bool SimpleStraightSkel::isReflex(EdgeSPtr edge) {
    bool result = false;
    VertexSPtr vertex_src = edge->getVertexSrc();
    VertexSPtr vertex_dst = edge->getVertexDst();
    if (vertex_src->getPoint() == vertex_dst->getPoint()) {
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        FacetSPtr facet_src = getFacetSrc(edge);
        FacetSPtr facet_dst = getFacetDst(edge);
        double speed_l = 1.0;
        if (facet_l->hasData()) {
            speed_l = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_l->getData())->getSpeed();
        }
        double speed_r = 1.0;
        if (facet_r->hasData()) {
            speed_r = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_r->getData())->getSpeed();
        }
        double speed_src = 1.0;
        if (facet_src->hasData()) {
            speed_src = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_src->getData())->getSpeed();
        }
        double speed_dst = 1.0;
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
        double speed_l = 1.0;
        if (facet_l->hasData()) {
            speed_l = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_l->getData())->getSpeed();
        }
        double speed_r = 1.0;
        if (facet_r->hasData()) {
            speed_r = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_r->getData())->getSpeed();
        }
        double speed_src = 1.0;
        if (facet_src->hasData()) {
            speed_src = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_src->getData())->getSpeed();
        }
        double speed_dst = 1.0;
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


void SimpleStraightSkel::run() {
    if (controller_) {
        controller_->wait();
        controller_->setDispPolyhedron(polyhedron_);
        controller_->setDispSkel3d(skel_result_);
    }
    DEBUG_PRINT("== Straight Skeleton 3D started ==");
    double t_start = util::Timer::now();
    PolyhedronSPtr polyhedron = polyhedron_->clone();

    // Simple test for weighted straight skeleton.
    //unsigned int j = 0;
    //list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    //while (it_f != polyhedron->facets().end()) {
    //    FacetSPtr facet = *it_f++;
    //    if (j == 0) {
    //        SkelFacetDataSPtr data = SkelFacetData::create(facet);
    //        data->setSpeed(10.0);
    //        break;
    //    }
    //    j++;
    //}

    Point3SPtr p_box_min;
    Point3SPtr p_box_max;
    unsigned int i = 0;
    DEBUG_VAL("Using " << vertex_splitter_->toString()
            << " to initialize polyhedron.");
    if (init(polyhedron)) {
        if (controller_) {
            controller_->wait();
        }
        double offset = 0.0;
        double offset_prev = 0.0;
        AbstractEventSPtr event = nextEvent(polyhedron, offset);
        while (event) {
            DEBUG_VAL("-- Next Event: " << event->toString() << " --");
            if (controller_) {
                controller_->wait();
            }
            offset = event->getOffset();
            polyhedron = shiftFacets(polyhedron, offset - offset_prev);
            if (event->getType() == AbstractEvent::CONST_OFFSET_EVENT) {
                event->setPolyhedronResult(polyhedron);
                skel_result_->addEvent(event);
                bool screenshot_on_const_offset_event =
                        util::Configuration::getInstance()->getBool(
                        "algo_3d_SimpleStraightSkel", "screenshot_on_const_offset_event");
                if (controller_ && screenshot_on_const_offset_event) {
                    controller_->screenshot();
                }
            } else if (event->getType() == AbstractEvent::SAVE_OFFSET_EVENT) {
                event->setPolyhedronResult(polyhedron);
                skel_result_->addEvent(event);
                std::stringstream ss_filename;
                ss_filename << "offset_" << offset << ".obj";
                db::_3d::OBJFile::save(ss_filename.str(), polyhedron);
                save_offsets_.pop_front();
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
            assert(polyhedron->isConsistent());
            assert(skel_result_->isConsistent());
            if (p_box_min && p_box_max) {
                assert(PolyhedronTransformation::isInsideBox(polyhedron,
                        p_box_min, p_box_max));
            } else {
                p_box_min = PolyhedronTransformation::boundingBoxMin(polyhedron);
                p_box_max = PolyhedronTransformation::boundingBoxMax(polyhedron);
            }
            DEBUG_PRINT("-- Finished handling Event --");
            i++;
            DEBUG_VAR(i);
            if (controller_) {
                controller_->wait();
            }
            event = nextEvent(polyhedron, offset);
            offset_prev = offset;
        }
        DEBUG_PRINT("== Straight Skeleton 3D finished ==");
        double time = util::Timer::now() - t_start;
        skel_result_->appendDescription("time=" +
                util::StringFactory::fromDouble(time) + "; ");
        //skel_result_->appendDescription("controller=" +
        //        util::StringFactory::fromBoolean(controller_) + "; ");
        DEBUG_VAR(skel_result_->toString());
    }
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
        if (i >= 3) {
            Vector3SPtr direction;
            Plane3SPtr plane_1 = facets[0]->plane();
            Plane3SPtr plane_2 = facets[1]->plane();
            Plane3SPtr plane_3 = facets[2]->plane();
            double speed_1 = 1.0;
            if (facets[0]->hasData()) {
                speed_1 = std::dynamic_pointer_cast<SkelFacetData>(
                        facets[0]->getData())->getSpeed();
            }
            double speed_2 = 1.0;
            if (facets[1]->hasData()) {
                speed_2 = std::dynamic_pointer_cast<SkelFacetData>(
                        facets[1]->getData())->getSpeed();
            }
            double speed_3 = 1.0;
            if (facets[2]->hasData()) {
                speed_3 = std::dynamic_pointer_cast<SkelFacetData>(
                        facets[2]->getData())->getSpeed();
            }
            Point3SPtr src = KernelWrapper::intersection(plane_1, plane_2, plane_3);
            Plane3SPtr off_1 = KernelWrapper::offsetPlane(plane_1, -speed_1);
            Plane3SPtr off_2 = KernelWrapper::offsetPlane(plane_2, -speed_2);
            Plane3SPtr off_3 = KernelWrapper::offsetPlane(plane_3, -speed_3);
            Point3SPtr dst = KernelWrapper::intersection(off_1, off_2, off_3);
            if (src && dst) {
                direction = KernelFactory::createVector3(*dst - *src);
            }
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
        double speed_l = 1.0;
        if (facet_l->hasData()) {
            SkelFacetDataSPtr data_l = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_l->getData());
            facet_b = data_l->getFacetOrigin();
            speed_l = data_l->getSpeed();
        }
        FacetSPtr facet_f = facet_r;
        double speed_r = 1.0;
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
    it_v = vertices_tosplit.begin();
    while (it_v != vertices_tosplit.end()) {
        VertexSPtr vertex = *it_v++;

        bool equal_speeds = true;
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (facet_wptr.expired()) {
                continue;
            }
            FacetSPtr facet(facet_wptr);
            double speed = 1.0;
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
            AbstractVertexSplitter::splitConvexVertex(vertex);
        } else if (use_fast_vertex_splitter_ && equal_speeds && vertex->isReflex()) {
            AbstractVertexSplitter::splitReflexVertex(vertex);
        } else {
            l.unlock();
            DEBUG_VAR(vertex->toString());
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
    assert(skel_result_->isConsistent());
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

Point3SPtr SimpleStraightSkel::vanishesAt(EdgeSPtr edge) {
    Point3SPtr result = Point3SPtr();
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
            if (KernelWrapper::side(facet->plane(), p_intersect) <= 0) {
                // inside polyhedron
                result = p_intersect;
            }
        }
    }
    return result;
}

Point3SPtr SimpleStraightSkel::crashAt(EdgeSPtr edge_1, EdgeSPtr edge_2) {
    Point3SPtr result = Point3SPtr();
    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(edge_1->getData());
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(edge_2->getData());
    Line3SPtr line_intersection = KernelWrapper::intersection(
            data_1->getSheet()->getPlane(),
            data_2->getSheet()->getPlane());
    if (!line_intersection) {
        // degenerated case
        return result;
    }
    Plane3SPtr plane_1 = edge_1->getFacetL()->plane();
    Plane3SPtr plane_2 = edge_2->getFacetL()->plane();
    Point3SPtr point_1 = KernelWrapper::intersection(plane_1, line_intersection);
    Point3SPtr point_2 = KernelWrapper::intersection(plane_2, line_intersection);
    if (!point_1 || !point_2) {
        // degenerated case
        return result;
    }
    Vector3SPtr direction = KernelFactory::createVector3(*point_2 - *point_1);
    double distance = KernelWrapper::distance(point_2, point_1);
    double facet_speed_1 = 1.0;
    if (edge_1->getFacetL()->hasData()) {
        facet_speed_1 = std::dynamic_pointer_cast<SkelFacetData>(
                edge_1->getFacetL()->getData())->getSpeed();
    }
    double facet_speed_2 = 1.0;
    if (edge_2->getFacetL()->hasData()) {
        facet_speed_2 = std::dynamic_pointer_cast<SkelFacetData>(
                edge_2->getFacetL()->getData())->getSpeed();
    }
    Plane3SPtr offset_plane_1 = KernelWrapper::offsetPlane(plane_1, -facet_speed_1);
    Plane3SPtr offset_plane_2 = KernelWrapper::offsetPlane(plane_2, -facet_speed_2);
    Point3SPtr offset_point_1 = KernelWrapper::intersection(offset_plane_1, line_intersection);
    Point3SPtr offset_point_2 = KernelWrapper::intersection(offset_plane_2, line_intersection);
    Vector3SPtr direction_1 = KernelFactory::createVector3(*offset_point_1 - *point_1);
    Vector3SPtr direction_2 = KernelFactory::createVector3(*offset_point_2 - *point_2);
    double speed_1 = KernelWrapper::distance(offset_point_1, point_1);
    //if (KernelWrapper::angle(direction_1, direction) > M_PI/2.0) {
    if ((*direction_1 * *direction) < 0.0) {
        speed_1 *= -1.0;
    }
    double speed_2 = KernelWrapper::distance(offset_point_2, point_2);
    //if (KernelWrapper::angle(direction_2, direction) > M_PI/2.0) {
    if ((*direction_2 * *direction) < 0.0) {
        speed_2 *= -1.0;
    }
    double dist_1 = (distance * speed_1) / (speed_1 - speed_2);
    Point3SPtr point = KernelWrapper::offsetPoint(point_1, direction, dist_1);

    bool inside_bounds = true;
    FacetSPtr facet_1_src = getFacetSrc(edge_1);
    FacetSPtr facet_1_dst = getFacetDst(edge_1);
    FacetSPtr facet_2_src = getFacetSrc(edge_2);
    FacetSPtr facet_2_dst = getFacetDst(edge_2);
    Vector3SPtr normal_1 = KernelFactory::createVector3(data_1->getSheet()->getPlane());
    Line3SPtr line_normal_1 = KernelFactory::createLine3(point, normal_1);
    if (KernelWrapper::orientation(line(edge_1), line_normal_1) <= 0) {
        inside_bounds = false;
    }
    if (!(facet_1_src == edge_2->getFacetL() ||
            facet_1_src == edge_2->getFacetR())) {
        SkelVertexDataSPtr data_1_src = std::dynamic_pointer_cast<SkelVertexData>(
            edge_1->getVertexSrc()->getData());
        ArcSPtr arc_1_src = data_1_src->getArc();
        if (KernelWrapper::orientation(arc_1_src->line(), line_normal_1) > 0) {
            inside_bounds = false;
        }
    }
    if (!(facet_1_dst == edge_2->getFacetL() ||
            facet_1_dst == edge_2->getFacetR())) {
        SkelVertexDataSPtr data_1_dst = std::dynamic_pointer_cast<SkelVertexData>(
            edge_1->getVertexDst()->getData());
        ArcSPtr arc_1_dst = data_1_dst->getArc();
        if (KernelWrapper::orientation(arc_1_dst->line(), line_normal_1) < 0) {
            inside_bounds = false;
        }
    }
    Vector3SPtr normal_2 = KernelFactory::createVector3(data_2->getSheet()->getPlane());
    Line3SPtr line_normal_2 = KernelFactory::createLine3(point, normal_2);
    if (KernelWrapper::orientation(line(edge_2), line_normal_2) <= 0) {
        inside_bounds = false;
    }
    if (!(facet_2_src == edge_1->getFacetL() ||
            facet_2_src == edge_1->getFacetR())) {
        SkelVertexDataSPtr data_2_src = std::dynamic_pointer_cast<SkelVertexData>(
            edge_2->getVertexSrc()->getData());
        ArcSPtr arc_2_src = data_2_src->getArc();
        if (KernelWrapper::orientation(arc_2_src->line(), line_normal_2) > 0) {
            inside_bounds = false;
        }
    }
    if (!(facet_2_dst == edge_1->getFacetL() ||
            facet_2_dst == edge_1->getFacetR())) {
        SkelVertexDataSPtr data_2_dst = std::dynamic_pointer_cast<SkelVertexData>(
            edge_2->getVertexDst()->getData());
        ArcSPtr arc_2_dst = data_2_dst->getArc();
        if (KernelWrapper::orientation(arc_2_dst->line(), line_normal_2) < 0) {
            inside_bounds = false;
        }
    }
    if (inside_bounds) {
        result = point;
    }
    return result;
}

double SimpleStraightSkel::offsetDist(FacetSPtr facet, Point3SPtr point) {
    Plane3SPtr plane = facet->plane();
    double result = KernelWrapper::distance(plane, point);
    if (KernelWrapper::side(plane, point) < 0) {
        result *= -1.0;
    }
    if (facet->hasData()) {
        double speed = std::dynamic_pointer_cast<SkelFacetData>(
                facet->getData())->getSpeed();
        result /= speed;
    }
    return result;
}


EdgeEventSPtr SimpleStraightSkel::nextEdgeEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    EdgeEventSPtr result = EdgeEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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
        Point3SPtr point = vanishesAt(edge);
        if (!point) {
            continue;
        }
        FacetSPtr facet_src = getFacetSrc(edge);
        FacetSPtr facet_dst = getFacetDst(edge);
        // This does not work when there is more than one edge between both facets.
        // EdgeSPtr edge_2 = facet_src->findEdge(facet_dst);
        bool split_event = false;
        std::list<EdgeSPtr> edges_2 = facet_src->findEdges(facet_dst);
        std::list<EdgeSPtr>::iterator it_e2 = edges_2.begin();
        while (it_e2 != edges_2.end()) {
            EdgeSPtr edge_2 = *it_e2++;
            bool split_event_current = true;
            SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(edge_2->getData());
            Vector3SPtr normal_2 = KernelFactory::createVector3(data_2->getSheet()->getPlane());
            Line3SPtr line_normal_2 = KernelFactory::createLine3(point, normal_2);
            if (KernelWrapper::orientation(line(edge_2), line_normal_2) < 0) {
                // out of bounded area
                split_event_current = false;
            }
            SkelVertexDataSPtr data_2_src = std::dynamic_pointer_cast<SkelVertexData>(
                edge_2->getVertexSrc()->getData());
            ArcSPtr arc_2_src = data_2_src->getArc();
            if (KernelWrapper::orientation(arc_2_src->line(), line_normal_2) > 0) {
                // out of bounded area
                split_event_current = false;
            }
            SkelVertexDataSPtr data_2_dst = std::dynamic_pointer_cast<SkelVertexData>(
                edge_2->getVertexDst()->getData());
            ArcSPtr arc_2_dst = data_2_dst->getArc();
            if (KernelWrapper::orientation(arc_2_dst->line(), line_normal_2) < 0) {
                // out of bounded area
                split_event_current = false;
            }
            if (split_event_current) {
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

        double offset_event = offsetDist(facet_l, point);
        if (offset_event > offset_max) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = EdgeEvent::create();
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(offset + offset_event);
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
            offset_max = offset_event;
        }
    }
    return result;
}

EdgeMergeEventSPtr SimpleStraightSkel::nextEdgeMergeEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    EdgeMergeEventSPtr result = EdgeMergeEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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

        Point3SPtr point = vanishesAt(edge);
        if (!point) {
            continue;
        }
        double offset_event = offsetDist(facet_l, point);
        if (offset_event > offset_max) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = EdgeMergeEvent::create();
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(offset + offset_event);
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
            offset_max = offset_event;
        }
    }
    return result;
}

TriangleEventSPtr SimpleStraightSkel::nextTriangleEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    TriangleEventSPtr result = TriangleEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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

        Point3SPtr point = vanishesAt(edge);
        if (!point) {
            continue;
        }
        if ((KernelWrapper::side(edge->getFacetL()->plane(), point) > 0.0) ||
                KernelWrapper::side(edge->getFacetR()->plane(), point) > 0.0) {
            // triangle may not be a hole
            // after pierce event
            continue;
        }
        double offset_event = offsetDist(facet, point);
        if (offset_event > offset_max) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = TriangleEvent::create();
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(offset + offset_event);
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

            offset_max = offset_event;
        }
    }
    return result;
}

DblEdgeMergeEventSPtr SimpleStraightSkel::nextDblEdgeMergeEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    DblEdgeMergeEventSPtr result = DblEdgeMergeEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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

        Point3SPtr point = vanishesAt(edge);
        if (!point) {
            continue;
        }
        double offset_event = offsetDist(edge->getFacetL(), point);
        if (offset_event > offset_max) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = DblEdgeMergeEvent::create();
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(offset + offset_event);
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
            offset_max = offset_event;
        }
    }
    return result;
}

DblTriangleEventSPtr SimpleStraightSkel::nextDblTriangleEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    DblTriangleEventSPtr result = DblTriangleEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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

        Point3SPtr point = vanishesAt(edge);
        if (!point) {
            continue;
        }
        double offset_event = offsetDist(edge->getFacetL(), point);
        if (offset_event > offset_max) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = DblTriangleEvent::create();
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(offset + offset_event);
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

            offset_max = offset_event;
        }
    }
    return result;
}

TetrahedronEventSPtr SimpleStraightSkel::nextTetrahedronEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    TetrahedronEventSPtr result = TetrahedronEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (isTetrahedron(edge)) {
            FacetSPtr facet = edge->getFacetL();
            Point3SPtr point = vanishesAt(edge);
            if (!point) {
                continue;
            }
            double offset_event = offsetDist(facet, point);
            if (offset_event > offset_max) {
                NodeSPtr node;
                if (!result) {
                    node = Node::create(point);
                    result = TetrahedronEvent::create();
                    result->setNode(node);
                }
                node = result->getNode();
                node->clear();
                node->setOffset(offset + offset_event);
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

                offset_max = offset_event;
            }
        }
    }
    return result;
}

VertexEventSPtr SimpleStraightSkel::nextVertexEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    VertexEventSPtr result = VertexEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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

            Point3SPtr point = crashAt(edge_11, edge_22);
            if (!point) {
                continue;
            }
            double offset_event = offsetDist(edge_11->getFacetL(), point);
            if (offset_event > offset_max) {
                NodeSPtr node;
                if (!result) {
                    node = Node::create(point);
                    result = VertexEvent::create();
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
                node->setOffset(offset + offset_event);
                node->setPoint(point);
                result->setVertex1(vertex_1);
                result->setVertex2(vertex_2);
                result->setFacet1(facet_1);
                result->setFacet2(facet_2);
                offset_max = offset_event;
            }
        }
    }
    return result;
}

FlipVertexEventSPtr SimpleStraightSkel::nextFlipVertexEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    FlipVertexEventSPtr result = FlipVertexEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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

            Point3SPtr point = crashAt(edge_11, edge_22);
            if (!point) {
                continue;
            }
            double offset_event = offsetDist(edge_11->getFacetL(), point);
            if (offset_event > offset_max) {
                NodeSPtr node;
                if (!result) {
                    node = Node::create(point);
                    result = FlipVertexEvent::create();
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
                node->setOffset(offset + offset_event);
                node->setPoint(point);
                result->setVertex1(vertex_1);
                result->setVertex2(vertex_2);
                result->setFacet1(facet_1);
                result->setFacet2(facet_2);
                offset_max = offset_event;
            }
        }
    }
    return result;
}

SurfaceEventSPtr SimpleStraightSkel::nextSurfaceEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    SurfaceEventSPtr result = SurfaceEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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
            Point3SPtr point = crashAt(edge_1, edge_2);
            if (!point) {
                continue;
            }

            // find minimum orthogonal distance
            double offset_event = offsetDist(edge_1->getFacetL(), point);
            if (offset_event > offset_max) {
                NodeSPtr node;
                if (!result) {
                    node = Node::create(point);
                    result = SurfaceEvent::create();
                    result->setNode(node);
                }
                node = result->getNode();
                node->clear();
                node->setOffset(offset + offset_event);
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

                offset_max = offset_event;
            }
        }
    }
    return result;
}

PolyhedronSplitEventSPtr SimpleStraightSkel::nextPolyhedronSplitEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    PolyhedronSplitEventSPtr result = PolyhedronSplitEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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
            Point3SPtr point = crashAt(edge_1, edge_2);
            if (!point) {
                continue;
            }

            // find minimum orthogonal distance
            double offset_event = offsetDist(edge_1->getFacetL(), point);
            if (offset_event > offset_max) {
                NodeSPtr node;
                if (!result) {
                    node = Node::create(point);
                    result = PolyhedronSplitEvent::create();
                    result->setNode(node);
                }
                node = result->getNode();
                node->clear();
                node->setOffset(offset + offset_event);
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

                offset_max = offset_event;
            }
        }
    }
    return result;
}

SplitMergeEventSPtr SimpleStraightSkel::nextSplitMergeEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    SplitMergeEventSPtr result = SplitMergeEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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

            Point3SPtr point = crashAt(edge_11, edge_22);
            if (!point) {
                continue;
            }
            double offset_event = offsetDist(edge_11->getFacetL(), point);
            if (offset_event > offset_max) {
                NodeSPtr node;
                if (!result) {
                    node = Node::create(point);
                    result = SplitMergeEvent::create();
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
                node->setOffset(offset + offset_event);
                node->setPoint(point);
                result->setVertex1(vertex_1);
                result->setVertex2(vertex_2);
                result->setFacet1(facet_1);
                result->setFacet2(facet_2);
                offset_max = offset_event;
            }
        }
    }
    return result;
}

EdgeSplitEventSPtr SimpleStraightSkel::nextEdgeSplitEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    EdgeSplitEventSPtr result = EdgeSplitEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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
            Point3SPtr point = crashAt(edge_1, edge_2);
            if (!point) {
                continue;
            }

            // find minimum orthogonal distance
            double offset_event = offsetDist(edge_1->getFacetL(), point);
            if (offset_event > offset_max) {
                NodeSPtr node;
                if (!result) {
                    node = Node::create(point);
                    result = EdgeSplitEvent::create();
                    result->setNode(node);
                }
                node = result->getNode();
                node->clear();
                node->setOffset(offset + offset_event);
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

                offset_max = offset_event;
            }
        }
    }
    return result;
}

PierceEventSPtr SimpleStraightSkel::nextPierceEvent(PolyhedronSPtr polyhedron, double offset) {
    ReadLock l(polyhedron->mutex());
    PierceEventSPtr result = PierceEventSPtr();
    double offset_max = -std::numeric_limits<double>::max();
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

                if (KernelWrapper::side(facet->plane(), vertex->getPoint()) > 0) {
                    continue;
                }
                if (!IsLineInFacet(facet, arc->line())) {
                    continue;
                }

                FacetSPtr facet_vertex = FacetSPtr(vertex->facets().front());
                double facet_speed_vertex = 1.0;
                if (facet_vertex->hasData()) {
                    facet_speed_vertex = std::dynamic_pointer_cast<SkelFacetData>(
                            facet_vertex->getData())->getSpeed();
                }
                Plane3SPtr plane_vertex_offset = KernelWrapper::offsetPlane(facet_vertex->plane(), -facet_speed_vertex);
                Point3SPtr point_vertex_offset = KernelWrapper::intersection(plane_vertex_offset, arc->line());
                double speed_vertex = KernelWrapper::distance(vertex->getPoint(), point_vertex_offset);

                Point3SPtr point_facet = KernelWrapper::intersection(facet->plane(), arc->line());
                double facet_speed = 1.0;
                if (facet->hasData()) {
                    facet_speed = std::dynamic_pointer_cast<SkelFacetData>(
                            facet->getData())->getSpeed();
                }
                Plane3SPtr plane_facet_offset = KernelWrapper::offsetPlane(facet->plane(), -facet_speed);
                Point3SPtr point_facet_offset = KernelWrapper::intersection(plane_facet_offset, arc->line());
                double speed_facet = KernelWrapper::distance(point_facet, point_facet_offset);

                double distance = KernelWrapper::distance(vertex->getPoint(), point_facet);
                double dist_vertex = (distance * speed_vertex) /
                        (speed_vertex + speed_facet);
                if (KernelWrapper::comparePoints(arc->getDirection(),
                        point_facet, point_facet_offset) < 0) {
                    // for weighted straight skeleton
                    // reflex vertex and facet move into same direction
                    if (speed_facet < speed_vertex) {
                        // facet to slow
                        continue;
                    }
                    dist_vertex = (distance * speed_vertex) /
                            (speed_facet - speed_vertex);
                }
                double offset_event = -dist_vertex / speed_vertex;
                Point3SPtr point = KernelWrapper::offsetPoint(vertex->getPoint(),
                        arc->getDirection(), dist_vertex);
                if (offset_event > offset_max) {
                    NodeSPtr node;
                    if (!result) {
                        node = Node::create(point);
                        result = PierceEvent::create();
                        result->setNode(node);
                    }
                    node = result->getNode();
                    node->clear();
                    node->addArc(arc);
                    node->setOffset(offset + offset_event);
                    node->setPoint(point);
                    result->setFacet(facet);
                    result->setVertex(vertex);
                    offset_max = offset_event;
                }
            }
        }
    }
    return result;
}


AbstractEventSPtr SimpleStraightSkel::nextEvent(PolyhedronSPtr polyhedron, double offset) {
    AbstractEventSPtr result = AbstractEventSPtr();
    if (!polyhedron) {
        return result;
    }
    if (polyhedron->facets().size() == 0) {
        return result;
    }
    AbstractEventSPtr events[15];
    for (unsigned int i = 0; i < 15; i++) {
        events[i] = AbstractEventSPtr();
    }
    double const_offset = util::Configuration::getInstance()->getDouble(
            "algo_3d_SimpleStraightSkel", "const_offset");
    if (const_offset != 0.0) {
        double next_offset = floor(offset/const_offset + 1.0) * const_offset;
        if (next_offset >= offset) {
            next_offset += const_offset;
        }
        events[0] = ConstOffsetEvent::create(next_offset);
    }
    if (!save_offsets_.empty()) {
        events[1] = SaveOffsetEvent::create(save_offsets_.front());
    }
    events[2] = nextEdgeEvent(polyhedron, offset);
    events[3] = nextEdgeMergeEvent(polyhedron, offset);
    events[4] = nextTriangleEvent(polyhedron, offset);
    events[5] = nextDblEdgeMergeEvent(polyhedron, offset);
    events[6] = nextDblTriangleEvent(polyhedron, offset);
    events[7] = nextTetrahedronEvent(polyhedron, offset);
    events[8] = nextVertexEvent(polyhedron, offset);
    events[9] = nextFlipVertexEvent(polyhedron, offset);
    events[10] = nextSurfaceEvent(polyhedron, offset);
    events[11] = nextPolyhedronSplitEvent(polyhedron, offset);
    events[12] = nextSplitMergeEvent(polyhedron, offset);
    events[13] = nextEdgeSplitEvent(polyhedron, offset);
    events[14] = nextPierceEvent(polyhedron, offset);
    for (unsigned int i = 0; i < 15; i++) {
        if (events[i]) {
            if (!result) {
                result = events[i];
            } else if (result->getOffset() < events[i]->getOffset()) {
                result = events[i];
            } else if (result->getOffset() == events[i]->getOffset()) {
                DEBUG_VAL("WARNING: More than one possible next event.");
                DEBUG_VAL(result->toString());
                DEBUG_VAL(events[i]->toString());
                result = events[i];
            }
        }
    }
    result->setHighlight(true);
    return result;
}

PolyhedronSPtr SimpleStraightSkel::shiftFacets(PolyhedronSPtr polyhedron, double offset) {
    PolyhedronSPtr result = Polyhedron::create();

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Plane3SPtr planes[3];
        unsigned int i = 0;
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (i < 3 && it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet = FacetSPtr(facet_wptr);
                double speed = 1.0;
                if (facet->hasData()) {
                    speed = std::dynamic_pointer_cast<SkelFacetData>(
                            facet->getData())->getSpeed();
                }
                planes[i] = KernelWrapper::offsetPlane(facet->plane(), offset*speed);
                i++;
            }
        }
        if (i >= 3) {
            Point3SPtr point = KernelWrapper::intersection(planes[0], planes[1], planes[2]);
            if (!point) {
                result = PolyhedronSPtr();
                DEBUG_SPTR(result);
                return result;
            }
            VertexSPtr offset_vertex = Vertex::create(point);
            // SkelVertexData for each vertex should be created by init
            SkelVertexDataSPtr data;
            if (vertex->hasData()) {
                data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
                SkelVertexDataSPtr offset_data = SkelVertexData::create(offset_vertex);
                offset_data->setArc(data->getArc());
            } else {
                data = SkelVertexData::create(vertex);
            }
            data->setOffsetVertex(offset_vertex);
            result->addVertex(offset_vertex);
        }
    }

    it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        EdgeSPtr edge;
        unsigned int i = 0;
        std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
        while (it_e != vertex->edges().end()) {
            EdgeWPtr edge_wptr = *it_e++;
            if (!edge_wptr.expired()) {
                edge = EdgeSPtr(edge_wptr);
                i++;
            }
        }
        if (i == 1) {
            VertexSPtr vertex_other;
            if (edge->getVertexSrc() == vertex) {
                vertex_other = edge->getVertexDst();
            } else if (edge->getVertexDst() == vertex) {
                vertex_other = edge->getVertexSrc();
            }
            SkelVertexDataSPtr data_other = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_other->getData());
            VertexSPtr offset_vertex_other = data_other->getOffsetVertex();
            Vector3 direction =
                    *(offset_vertex_other->getPoint()) - *(vertex_other->getPoint());
            Point3SPtr point = KernelFactory::createPoint3(
                    *(vertex->getPoint()) + direction);
            VertexSPtr offset_vertex = Vertex::create(point);
            SkelVertexDataSPtr data;
            if (vertex->hasData()) {
                data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
            } else {
                data = SkelVertexData::create(vertex);
            }
            data->setOffsetVertex(offset_vertex);
            result->addVertex(offset_vertex);
        }
    }

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        SkelVertexDataSPtr vertex_src_data = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexSrc()->getData());
        SkelVertexDataSPtr vertex_dst_data = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexDst()->getData());
        if (vertex_src_data && vertex_dst_data) {
            VertexSPtr offset_vertex_src = vertex_src_data->getOffsetVertex();
            VertexSPtr offset_vertex_dst = vertex_dst_data->getOffsetVertex();
            EdgeSPtr offset_edge = Edge::create(offset_vertex_src, offset_vertex_dst);
            SkelEdgeDataSPtr data;
            if (edge->hasData()) {
                data = std::dynamic_pointer_cast<SkelEdgeData>(edge->getData());
                SkelEdgeDataSPtr offset_data = SkelEdgeData::create(offset_edge);
                offset_data->setSheet(data->getSheet());
            } else {
                data = SkelEdgeData::create(edge);
            }
            data->setOffsetEdge(offset_edge);
            result->addEdge(offset_edge);
        }
    }

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        FacetSPtr offset_facet = Facet::create();
        SkelFacetDataSPtr data;
        double speed = 1.0;
        if (facet->hasData()) {
            data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
            speed = data->getSpeed();
            SkelFacetDataSPtr data_offset = SkelFacetData::create(offset_facet);
            data_offset->setFacetOrigin(data->getFacetOrigin());
            data_offset->setSpeed(speed);
        } else {
            data = SkelFacetData::create(facet);
        }
        Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->plane(), offset*speed);
        offset_facet->setPlane(offset_plane);
        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            if (vertex->hasData()) {
                SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex->getData());
                VertexSPtr offset_vertex = vertex_data->getOffsetVertex();
                if (offset_vertex) {
                    offset_facet->addVertex(offset_vertex);
                }
            }
        }
        std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                edge->getData());
            EdgeSPtr offset_edge = edge_data->getOffsetEdge();
            if (facet == edge->getFacetL()) {
                offset_edge->setFacetL(offset_facet);
                offset_facet->addEdge(offset_edge);
            }
            if (facet == edge->getFacetR()) {
                offset_edge->setFacetR(offset_facet);
                offset_facet->addEdge(offset_edge);
            }
        }
        data->setOffsetFacet(offset_facet);
        result->addFacet(offset_facet);
    }

    assert(polyhedron->edges().size() == result->edges().size());
    assert(polyhedron->facets().size() == result->facets().size());

    return result;
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

    bool not_flipped_valid = false;
    {
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
        edge_no_flip->setFacetL(facets_clone[3]);
        edge_no_flip->setFacetR(facets_clone[1]);
        facets_clone[3]->addEdge(edge_no_flip);
        facets_clone[1]->addEdge(edge_no_flip);
        for (unsigned int i = 0; i < 4; i++) {
            edges[i]->setFacetR(facets_clone[(i+3)%4]);
            edges[i]->setFacetL(facets_clone[i]);
            facets_clone[(i+3)%4]->addEdge(edges[i]);
            facets_clone[i]->addEdge(edges[i]);
        }
        PolyhedronSPtr polyhedron_no_flip = Polyhedron::create(4, facets_clone);
        PolyhedronSPtr polyhedron_no_flip_offset = shiftFacets(
                polyhedron_no_flip, -1.0);
        if (!SelfIntersection::hasSelfIntersectingSurface(
                polyhedron_no_flip_offset)) {
            not_flipped_valid = true;
        }
    }
    DEBUG_VAR(not_flipped_valid);

    bool flipped_valid = false;
    {
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
        edge_flipped->setFacetL(facets_clone[0]);
        edge_flipped->setFacetR(facets_clone[2]);
        facets_clone[0]->addEdge(edge_flipped);
        facets_clone[2]->addEdge(edge_flipped);
        for (unsigned int i = 0; i < 4; i++) {
            edges[i]->setFacetR(facets_clone[(i+3)%4]);
            edges[i]->setFacetL(facets_clone[i]);
            facets_clone[(i+3)%4]->addEdge(edges[i]);
            facets_clone[i]->addEdge(edges[i]);
        }
        PolyhedronSPtr polyhedron_flipped = Polyhedron::create(4, facets_clone);
        PolyhedronSPtr polyhedron_flipped_offset = shiftFacets(
                polyhedron_flipped, -1.0);
        if (!SelfIntersection::hasSelfIntersectingSurface(
                polyhedron_flipped_offset)) {
            flipped_valid = true;
        }
    }
    DEBUG_VAR(flipped_valid);

    if (flipped_valid && !not_flipped_valid) {
        flip_edge = true;
    } else if (not_flipped_valid && !flipped_valid) {
        flip_edge = false;
    } else if (flipped_valid && not_flipped_valid) {
        double angle_flipped = KernelWrapper::angle(
                facets[0]->plane(), facets[2]->plane());
        double angle_no_flip = KernelWrapper::angle(
                facets[1]->plane(), facets[3]->plane());
        DEBUG_VAR(angle_no_flip);
        DEBUG_VAR(angle_flipped);
        if (edge_event_ == 0) {
            // convex
            flip_edge = (angle_flipped <= angle_no_flip);
        } else if (edge_event_ == 1) {
            // reflex
            // choose edge that moves faster
            flip_edge = (angle_flipped >= angle_no_flip);
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
                throw std::runtime_error("Not able to handle EdgeEvent.");
            }
        }
    } else {
        throw std::runtime_error("Not able to handle EdgeEvent.");
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

    if (crashAt(edge_1, edge_b) || crashAt(edge_b, edge_1)) {
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
}

void SimpleStraightSkel::handlePolyhedronSplitEvent(PolyhedronSplitEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());

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
}

void SimpleStraightSkel::handleEdgeSplitEvent(EdgeSplitEventSPtr event, PolyhedronSPtr polyhedron) {
    WriteLock l(skel_result_->mutex());

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

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
            event->getVertex()->getData());
    VertexSPtr vertex_offset = vertex_data->getOffsetVertex();
    SkelFacetDataSPtr facet_data = std::dynamic_pointer_cast<SkelFacetData>(
            event->getFacet()->getData());
    FacetSPtr facet_offset = facet_data->getOffsetFacet();
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
