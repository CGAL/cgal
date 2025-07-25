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
 * @file   algo/3d/PolyhedronTransformation.cpp
 * @author Gernot Walzl
 * @date   2012-09-01
 */

#include "algo/3d/PolyhedronTransformation.h"

#include "algo/3d/KernelWrapper.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/skel/SkelVertexData.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "data/3d/skel/SkelFacetData.h"
#include "db/3d/OBJFile.h"

#include "util/StringFactory.h"
#include "util/Configuration.h"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/simplest_rational_in_interval.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/internal/triangle_soup_snap_rounding.h>

#include <boost/property_map/property_map.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <list>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <stack>
#include <unordered_map>
#include <vector>

using namespace data::_3d;

namespace CGAL {
namespace ss_utils {

std::size_t length(const CGAL::FT& n)
{
  std::stringstream ss;
  ss << CGAL::exact(n);
  std::string str = ss.str();
  return str.size();
}

std::size_t length(const Point3& p)
{
  // returns the maximum number of digits between the x, y and z coordinates
  return (std::max)({length(p.x()), length(p.y()), length(p.z())});
}

std::size_t length(const Plane3& plane)
{
  // returns the maximum number of digits between the a, b, c and d coefficients
  return (std::max)({length(plane.a()), length(plane.b()), length(plane.c()), length(plane.d())});
}

} // namespace ss_utils
} // namespace CGAL

namespace algo { namespace _3d {

// @todo use this in OBJ/SM IO functions
auto triangulate_facet_with_CDT2(FacetSPtr facet)
{
    using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
    using PK = CGAL::Projection_traits_3<CGAL::K>;
    using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
    using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
    using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
    using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
    using PCDT_VH = PCDT::Vertex_handle;

    Vector3SPtr n = KernelFactory::createVector3(facet->plane());
    CGAL_precondition(*n != CGAL::NULL_VECTOR);

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
            vh->info() = vertex;
        }
    }

    std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr v0 = edge->src(facet);
        VertexSPtr v1 = edge->dst(facet);
        CGAL_assertion(*(v0->getPoint()) != *(v1->getPoint()));

        PCDT_VH vh0 = face_vhs.at(v0);
        PCDT_VH vh1 = face_vhs.at(v1);

        try {
            pcdt.insert_constraint(vh0, vh1);
        } catch(const typename PCDT::Intersection_of_constraints_exception&) {
            CGAL_SS3_TRANSF_TRACE("Error: Intersection of constraint w/ " << vh0->point() << " " << vh1->point());
            CGAL_SS3_TRANSF_TRACE(facet->toString());
            CGAL_assertion_msg(false, "Intersections are not allowed");
            return PCDT(projection_traits);
        }
    }

    return pcdt;
}

PolyhedronTransformation::PolyhedronTransformation() {
    // intentionally does nothing.
}

PolyhedronTransformation::~PolyhedronTransformation() {
    // intentionally does nothing.
}

Point3SPtr PolyhedronTransformation::boundingBoxMin(PolyhedronSPtr polyhedron) {
    CGAL::FT p_min[3];
    for (unsigned int i = 0; i < 3; i++) {
        p_min[i] = std::numeric_limits<double>::max();
    }
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        for (unsigned int i = 0; i < 3; i++) {
            if ((*p)[i] < p_min[i]) {
                p_min[i] = (*p)[i];
            }
        }
    }
    Point3SPtr result = KernelFactory::createPoint3(p_min[0],p_min[1],p_min[2]);
    return result;
}

Point3SPtr PolyhedronTransformation::boundingBoxMax(PolyhedronSPtr polyhedron) {
    CGAL::FT p_max[3];
    for (unsigned int i = 0; i < 3; i++) {
        p_max[i] = -std::numeric_limits<double>::max();
    }
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        for (unsigned int i = 0; i < 3; i++) {
            if ((*p)[i] > p_max[i]) {
                p_max[i] = (*p)[i];
            }
        }
    }
    Point3SPtr result = KernelFactory::createPoint3(p_max[0],p_max[1],p_max[2]);
    return result;
}

bool PolyhedronTransformation::isInsideBox(PolyhedronSPtr polyhedron,
    Point3SPtr p_box_min, Point3SPtr p_box_max) {
    bool result = true;
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        for (unsigned int i = 0; i < 3; i++) {
            if (!((*p_box_min)[i] <= (*p)[i] &&
                    (*p)[i] <= (*p_box_max)[i])) {
                result = false;
                // CGAL_SS3_TRANSF_TRACE(*p << " is not in the box " << *p_box_min << " " << *p_box_max);
                break;
            }
        }
        if (!result) {
            break;
        }
    }
    return result;
}

void PolyhedronTransformation::translate(PolyhedronSPtr polyhedron, Vector3SPtr v_t) {
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        Point3SPtr p_t = KernelFactory::createPoint3(*p + *v_t);
        vertex->setPoint(p_t);
    }

    polyhedron->initPlanes();

    polyhedron->appendDescription("translate=<" +
            util::StringFactory::fromDouble(CGAL::to_double((*v_t)[0])) + ", " +
            util::StringFactory::fromDouble(CGAL::to_double((*v_t)[1])) + ", " +
            util::StringFactory::fromDouble(CGAL::to_double((*v_t)[2])) + ">; ");
}

void PolyhedronTransformation::scale(PolyhedronSPtr polyhedron, Vector3SPtr v_s) {
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        Point3SPtr p_s = KernelFactory::createPoint3(
                (*p)[0] * (*v_s)[0], (*p)[1] * (*v_s)[1], (*p)[2] * (*v_s)[2]);
        vertex->setPoint(p_s);
    }
    polyhedron->initPlanes();

    polyhedron->appendDescription("scale=<" +
            util::StringFactory::fromDouble(CGAL::to_double((*v_s)[0])) + ", " +
            util::StringFactory::fromDouble(CGAL::to_double((*v_s)[1])) + ", " +
            util::StringFactory::fromDouble(CGAL::to_double((*v_s)[2])) + ">; ");
}

void PolyhedronTransformation::translateNscale(PolyhedronSPtr polyhedron,
    Point3SPtr p_box_min, Point3SPtr p_box_max) {
    Vector3SPtr v_box_min = KernelFactory::createVector3(p_box_min);
    Vector3SPtr v_box_max = KernelFactory::createVector3(p_box_max);
    Vector3SPtr v_size = KernelFactory::createVector3(*v_box_max - *v_box_min);
    Vector3SPtr v_center = KernelFactory::createVector3(
            (*v_box_min + *v_box_max) / 2.0);

    Point3SPtr p_box_min_curr = boundingBoxMin(polyhedron);
    Point3SPtr p_box_max_curr = boundingBoxMax(polyhedron);
    Vector3SPtr v_box_min_curr = KernelFactory::createVector3(p_box_min_curr);
    Vector3SPtr v_box_max_curr = KernelFactory::createVector3(p_box_max_curr);
    Vector3SPtr v_size_curr = KernelFactory::createVector3(
            *v_box_max_curr - *v_box_min_curr);
    Vector3SPtr v_center_curr = KernelFactory::createVector3(
            (*v_box_min_curr + *v_box_max_curr) / 2.0);

    CGAL::FT scale_factor = std::numeric_limits<double>::max(); // do not put FT
    for (unsigned int i = 0; i < 3; i++) {
        CGAL::FT s = (*v_size)[i]/(*v_size_curr)[i];
        if (scale_factor > s) {
            scale_factor = s;
        }
    }
    scale_factor = floor(CGAL::to_double(scale_factor)*1000.0)/1000.0; // @fixme interval
    Vector3SPtr v_s = KernelFactory::createVector3(
            scale_factor, scale_factor, scale_factor);

    Vector3SPtr v_t = KernelFactory::createVector3((*v_center_curr) * -1.0);
    if (v_t->squared_length() > 0.0) {
        translate(polyhedron, v_t);
    }
    if (scale_factor != 1.0) {
        scale(polyhedron, v_s);
    }
    if (v_center->squared_length() > 0.0) {
        translate(polyhedron, v_center);
    }
}

void PolyhedronTransformation::truncatePrecision(PolyhedronSPtr polyhedron)
{
    CGAL_SS3_TRANSF_TRACE("Lower precision of input polyhedron");

    util::ConfigurationSPtr config = util::Configuration::getInstance();
    double range = 1e-10;
    if (config->isLoaded()) {
        range = config->getDouble("main", "truncate_precision");
    }

    CGAL_SS3_TRANSF_TRACE(" truncate precision: " << range);

    double exp = std::ceil(std::log2(1.0 / range));
    double scale = std::pow(2.0, exp);
    CGAL_SS3_TRANSF_TRACE("  scale = " << scale);

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        CGAL_SS3_TRANSF_TRACE("    Truncated from: " << *(vertex->getPoint()));

        double rx = CGAL::Polygon_mesh_processing::autorefine_impl::double_ceil(p->x() * scale) / scale;
        double ry = CGAL::Polygon_mesh_processing::autorefine_impl::double_ceil(p->y() * scale) / scale;
        double rz = CGAL::Polygon_mesh_processing::autorefine_impl::double_ceil(p->z() * scale) / scale;

        vertex->setPoint(KernelFactory::createPoint3(rx, ry, rz));
        CGAL_SS3_TRANSF_TRACE("    Truncated to: " << *(vertex->getPoint()));
    }
}

bool PolyhedronTransformation::resetPoint(VertexSPtr vertex,
                                          const std::array<Plane3SPtr, 3>& planes)
{
    Point3SPtr point = KernelWrapper::intersection(planes[0], planes[1], planes[2]);
    if (!point) {
        CGAL_SS3_TRANSF_TRACE("Warning: triplet of planes does not define a point!");
        Point3SPtr result = Point3SPtr();
        CGAL_SS3_DEBUG_SPTR(result);
        return false;
    }

    vertex->setPoint(point);
    CGAL_SS3_TRANSF_TRACE("  New point = " << *point);

    CGAL_postcondition_code(std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();)
    CGAL_postcondition_code(for (; it_f!=vertex->facets().end(); ++it_f) {)
    CGAL_postcondition_code(    if (FacetSPtr facet = it_f->lock()) {)
    CGAL_postcondition(             facet->getPlane()->has_on(*point));
    CGAL_postcondition_code(    })
    CGAL_postcondition_code(})

    return true;
}

// resets using the first 3 planes, even if there are more
bool PolyhedronTransformation::resetPoint(VertexSPtr vertex)
{
    CGAL_SS3_TRANSF_TRACE("resetPoint() of " << vertex->toString());

    std::array<Plane3SPtr, 3> planes;
    unsigned int i = 0;
    std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (i < 3 && it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (!facet_wptr.expired()) {
            FacetSPtr facet = FacetSPtr(facet_wptr);
            planes[i++] = facet->plane();

            CGAL_SS3_TRANSF_TRACE("  Facet " << facet->getID() << " [" << *(facet->getPlane()) << "]");
        }
    }
    CGAL_postcondition(i == 3);

    return resetPoint(vertex, { planes[0], planes[1], planes[2] });
}

bool PolyhedronTransformation::resetPoints(PolyhedronSPtr polyhedron)
{
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (!resetPoint(vertex)) {
            CGAL_SS3_TRANSF_TRACE("Warning: failed to reset vertex " << vertex->toString());
            return false;
        }
    }
    return true;
}

// @todo this function cannot deal with degree 1 vertices
Point3SPtr PolyhedronTransformation::shiftPoint(VertexSPtr vertex,
                                                const CGAL::FT& offset)
{
    CGAL_SS3_TRANSF_TRACE("shift " << vertex->toString() << "\noffset = " << offset);

    CGAL_precondition(vertex->degree() >= 3);

    Plane3SPtr planes[3];
    unsigned int i = 0;
    std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (i < 3 && it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (!facet_wptr.expired()) {
            FacetSPtr facet = FacetSPtr(facet_wptr);
            Plane3SPtr plane = facet->plane();

            if (vertex->degree() > 3) {
                // planes are _offset_ planes, but it doesn't matter for the tests
                bool independent = true;
                if (i == 1) {
                    independent = !(CGAL::parallel(*(planes[0]), *plane));
                } else if (i == 2) {
                    // @speed avoid recomputing the intersection from scratch later
                    independent = !is_zero(CGAL::determinant(planes[0]->a(), planes[0]->b(), planes[0]->c(),
                                                             planes[1]->a(), planes[1]->b(), planes[1]->c(),
                                                             plane->a(), plane->b(), plane->c()));
                }

                if (!independent) {
                    continue;
                }
            }

            planes[i] = shiftPlane(facet, offset);

            CGAL_SS3_TRANSF_TRACE("facet[" << i << "] = " << facet->getID());
            CGAL_SS3_TRANSF_TRACE("Offset Plane[" << i << "] = " << *(planes[i]));

            ++i;
        }
    }

    if (i < 3) {
        CGAL_SS3_TRANSF_TRACE("Warning: Could not find three independent planes for high-degree vertex" << vertex->toString());
        return { };
    }

    Point3SPtr point = KernelWrapper::intersection(planes[0], planes[1], planes[2]);
    if (!point) {
        CGAL_SS3_TRANSF_TRACE("Error: triplet of planes doesn't define a point!");
        Point3SPtr result = Point3SPtr();
        CGAL_SS3_DEBUG_SPTR(result);
        return result;
    }

    CGAL_SS3_TRANSF_TRACE("  New point = " << *point);

    CGAL_assertion_code(for(Plane3SPtr pi : planes))
    CGAL_assertion(pi->has_on(*point));

    return point;
}

Segment3SPtr PolyhedronTransformation::shiftEdge(EdgeSPtr edge,
                                                 const CGAL::FT& offset)
{
    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    FacetSPtr facet_src = edge->getFacetL()->next(edge->getVertexSrc());
    FacetSPtr facet_dst = edge->getFacetR()->next(edge->getVertexDst());

    CGAL::FT speed_l = 1.0;
    if (facet_l->hasData()) {
        speed_l = std::dynamic_pointer_cast<SkelFacetData>(facet_l->getData())->getSpeed();
    }

    CGAL::FT speed_r = 1.0;
    if (facet_r->hasData()) {
        speed_r = std::dynamic_pointer_cast<SkelFacetData>(facet_r->getData())->getSpeed();
    }

    CGAL::FT speed_src = 1.0;
    if (facet_src->hasData()) {
        speed_src = std::dynamic_pointer_cast<SkelFacetData>(facet_src->getData())->getSpeed();
    }

    CGAL::FT speed_dst = 1.0;
    if (facet_dst->hasData()) {
        speed_dst = std::dynamic_pointer_cast<SkelFacetData>(facet_dst->getData())->getSpeed();
    }

        // Offset the two common planes
    Plane3SPtr offset_plane_l = KernelWrapper::offsetPlane(facet_l->plane(), offset*speed_l);
    Plane3SPtr offset_plane_r = KernelWrapper::offsetPlane(facet_r->plane(), offset*speed_r);
    Plane3SPtr offset_plane_src = KernelWrapper::offsetPlane(facet_src->plane(), offset*speed_src);
    Plane3SPtr offset_plane_dst = KernelWrapper::offsetPlane(facet_dst->plane(), offset*speed_dst);

#if 0
    // leaving it here because it's not that intuitive: factoring the intersection of the
    // two common planes is much slower than computing two 3-plane intersections
    Line3SPtr common_line = KernelWrapper::intersection(offset_plane_l, offset_plane_r);
    Point3SPtr src_point = KernelWrapper::intersection(offset_plane_src, common_line);
    Point3SPtr dst_point = KernelWrapper::intersection(offset_plane_dst, common_line);
#else
    Point3SPtr src_point = KernelWrapper::intersection(offset_plane_src, offset_plane_l, offset_plane_r);
    Point3SPtr dst_point = KernelWrapper::intersection(offset_plane_dst, offset_plane_l, offset_plane_r);
#endif

    CGAL_assertion(bool(src_point) && bool(dst_point));
    return KernelFactory::createSegment3(src_point, dst_point);
}

Plane3SPtr PolyhedronTransformation::shiftPlane(FacetSPtr facet,
                                                const CGAL::FT& offset)
{
    CGAL::FT facet_speed = 1.0;
    if (facet->hasData()) {
        facet_speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
    }

    return KernelWrapper::offsetPlane(facet->plane(), facet_speed*offset);
}

// This function is for the main shift in the event loop.
// @todo it could be adapted to also work for the smaller polyhedra that are used on-the-fly,
// for example in handleEdgeEvent, but then we would need to temporarily store points
// in this degree 3 loop because we need to use old positions in to compute the offset
// of degree 1 vertices.
void PolyhedronTransformation::shiftFacetsInPlace(PolyhedronSPtr polyhedron,
                                                  const CGAL::FT& offset,
                                                  const bool recompute_positions) {
    CGAL_SS3_TRANSF_TRACE("~~~~ Shift polyhedron by " << offset << " [in place]");

    // @speed shift planes once and recompute point and edge positions from this
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;

        // See comment at the top of the function
        CGAL_assertion(vertex->degree() == 3);

        Point3SPtr old_point = vertex->getPoint(), new_point;
        if (offset != 0 || recompute_positions) {
            new_point = shiftPoint(vertex, offset);
            if (!new_point) {
                CGAL_SS3_TRANSF_TRACE("Warning: Failed to shift polyhedron");
                return;
            }
        }

        // the old point position is not used to compute new point position: we need only the planes
        vertex->setPoint(new_point);

        // ugly hack while events still use getOffsetXYZ() to access elements in the shifted polyhedron
        // if getOffsetXYZ is no longer needed, it also means that combinatorics stored
        // in the events can be weak pointers
        SkelVertexDataSPtr data;
        if (vertex->hasData()) {
            data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
        } else {
            data = SkelVertexData::create(vertex);
        }
        data->setOffsetVertex(vertex); // @todo hack
    }

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        SkelEdgeDataSPtr data;
        if (edge->hasData()) {
            data = std::dynamic_pointer_cast<SkelEdgeData>(edge->getData());
        } else {
            data = SkelEdgeData::create(edge);
        }
        data->setOffsetEdge(edge); // @todo hack
    }

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;

        SkelFacetDataSPtr data;
        if (facet->hasData()) {
            data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
        } else {
            data = SkelFacetData::create(facet);
        }
        data->setOffsetFacet(facet); // @todo hack

        CGAL::FT speed = data->getSpeed();
        Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->plane(), offset*speed);
        facet->setPlane(offset_plane);
    }

    CGAL_postcondition(bool(polyhedron) && polyhedron->isConsistent());
}

// @speed plenty of needless recomputations:
// - when we do shiftPoint for adjacent points
// - when we call shiftPoint, and then call shiftPlane later on
// - ...
// @speed get rid of recompute_positions, seems obsolete now
PolyhedronSPtr PolyhedronTransformation::shiftFacets(PolyhedronSPtr polyhedron,
                                                     const CGAL::FT& offset,
                                                     const bool recompute_positions)
{
    CGAL_precondition(offset != 0);

    PolyhedronSPtr result = Polyhedron::create();

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;

        // this could be handled, but they are useless and should have been removed with
        // calls to the function removeVerticesDegLt3()
        CGAL_precondition(vertex->degree() != 2);

        // those are treated in the next loop
        if (vertex->degree() == 1) {
            continue;
        }

        Point3SPtr old_point = vertex->getPoint(), new_point;
        if (offset != 0 || recompute_positions) {
            new_point = shiftPoint(vertex, offset);
            if (!new_point) {
                std::cerr << "Warning: Failed to create shifted point" << std::endl;
                return { };
            }
        }

        VertexSPtr offset_vertex = Vertex::create(new_point);
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

    // Now, deal with degree 1 vertices.
    //
    // Do NOT merge the two vertices loops together: this function assumes that degree 1 vertices
    // are adjacent to degree 3+ vertices that have been offset in the first loop, and offsets
    // degree 1 vertices using the offset of the adjacent degree 3+ vertex.
#if 1
    it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;

        CGAL_precondition(vertex->degree() != 2);

        // those are treated in the first loop
        if (vertex->degree() >= 3) {
            continue;
        }

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

        CGAL_assertion(i == 1);

        VertexSPtr vertex_other = edge->other(vertex);
        CGAL_assertion(vertex_other->degree() >= 3);

        FacetSPtr facet1 = edge->getFacetL();
        FacetSPtr facet2 = edge->getFacetR();
        Plane3SPtr plane1 = facet1->plane();
        Plane3SPtr plane2 = facet2->plane();
        Line3SPtr intersection_line = KernelWrapper::intersection(plane1, plane2);
        CGAL_assertion(bool(intersection_line));

        // Determine the direction using the unshifted positions
        Point3SPtr point_other = vertex_other->getPoint();
        Point3SPtr point = vertex->getPoint();
        CGAL_assertion(*point != *point_other);
        Vector3 direction = intersection_line->to_vector();
        if (CGAL::scalar_product(direction, *point - *point_other) < 0) {
            direction = -direction;
        }

        // Apply the shift to the offset vertex
        SkelVertexDataSPtr data_other = std::dynamic_pointer_cast<SkelVertexData>(vertex_other->getData());
        VertexSPtr offset_vertex_other = data_other->getOffsetVertex();

        // The length of the direction does not matter: when degree 1 are evaluated, the code
        // treats them as infinite rays
        Point3SPtr offset_point = KernelFactory::createPoint3(*(offset_vertex_other->getPoint()) + direction);
        VertexSPtr offset_vertex = Vertex::create(offset_point);
        SkelVertexDataSPtr data;
        if (vertex->hasData()) {
            data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
        } else {
            data = SkelVertexData::create(vertex);
        }
        data->setOffsetVertex(offset_vertex);
        result->addVertex(offset_vertex);
    }
#else
    it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr point = vertex->getPoint();

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
            VertexSPtr vertex_other = (edge->getVertexSrc() == vertex) ? edge->getVertexDst()
                                                                       : edge->getVertexSrc();

            // otherwise vertex_other has not been offset in the previous vertex loop
            CGAL_assertion(vertex_other->degree() >= 3);

            if (offset != 0 || recompute_positions) {
                SkelVertexDataSPtr data_other = std::dynamic_pointer_cast<SkelVertexData>(vertex_other->getData());
                VertexSPtr offset_vertex_other = data_other->getOffsetVertex();
                Vector3 direction = *(offset_vertex_other->getPoint()) - *(vertex_other->getPoint());
                point = KernelFactory::createPoint3(*(vertex->getPoint()) + direction);
            }

            VertexSPtr offset_vertex = Vertex::create(point);
            SkelVertexDataSPtr data;
            if (vertex->hasData()) {
                data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
            } else {
                data = SkelVertexData::create(vertex);
            }
            data->setOffsetVertex(offset_vertex);
            result->addVertex(offset_vertex);
        } else {
            CGAL_assertion_msg(i != 0, "no edge on degree 1 vertex?");
        }
    }
#endif

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;

        SkelVertexDataSPtr vertex_src_data = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexSrc()->getData());
        SkelVertexDataSPtr vertex_dst_data = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexDst()->getData());

        if (vertex_src_data && vertex_dst_data) {
            VertexSPtr vertex_src = edge->getVertexSrc();
            VertexSPtr vertex_dst = edge->getVertexDst();
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
        CGAL::FT speed = 1.0;
        if (facet->hasData()) {
            data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
            speed = data->getSpeed();
            SkelFacetDataSPtr data_offset = SkelFacetData::create(offset_facet);
            data_offset->setFacetOrigin(data->getFacetOrigin());
            data_offset->setSpeed(speed);
        } else {
            data = SkelFacetData::create(facet);
            data->setSpeed(speed);
        }

        Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->plane(), offset*speed);
        offset_facet->setPlane(offset_plane);
        offset_facet->setBasePlaneID(facet->getBasePlaneID());

        // perturbation mechanisms
        offset_facet->cachedSpeed_ = facet->cachedSpeed_;
        offset_facet->cachedPlane_ = facet->cachedPlane_;

        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            if (vertex->hasData()) {
                SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex->getData());
                VertexSPtr offset_vertex = vertex_data->getOffsetVertex();
                if (offset_vertex) {
                    offset_facet->addVertex(offset_vertex);
                    CGAL_assertion(offset_plane->has_on(*(offset_vertex->getPoint())));
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

    result->initializeAllIDs();

    CGAL_postcondition(polyhedron->vertices().size() == result->vertices().size());
    CGAL_postcondition(polyhedron->edges().size() == result->edges().size());
    CGAL_postcondition(polyhedron->facets().size() == result->facets().size());

    return result;
}

std::pair<std::list<VertexSPtr>, std::list<FacetSPtr> >
PolyhedronTransformation::triangulate(FacetSPtr facet,
                                      PolyhedronSPtr polyhedron) {

    CGAL_SS3_TRANSF_TRACE("Triangulate facet " << facet->getID() << " of polyhedron "
                          << polyhedron->getID() << " with " << facet->vertices().size()
                          << " vertices");

    std::list<VertexSPtr> facet_vertices = facet->vertices();

    std::list<FacetSPtr> created_facets;
    CGAL_precondition(facet && polyhedron && facet->vertices().size() >= 3);

    if (facet_vertices.size() == 3) {
        return { facet_vertices, { facet } };
    }

    facet->sortVertices();

    // Prepare containers for triangulation
    auto pcdt = triangulate_facet_with_CDT2(facet);

    using PCDT = decltype(pcdt);
    using PCDT_FH = typename PCDT::Face_handle;

    std::unordered_map<PCDT_FH, bool> in_domain_map;
    boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(pcdt, in_domain);

   // Speed of the subdivded facet is applied to all the subfacets
    CGAL::FT parent_speed = 1.0;
    if (facet->hasData()) {
        SkelFacetDataSPtr parent_data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
        if (parent_data) {
            parent_speed = parent_data->getSpeed();
        }
    }

    polyhedron->removeFacet(facet);

    // Create new facets and edges for each triangle
    for(auto fh : pcdt.finite_face_handles()) {
        if(!get(in_domain, fh)) {
            continue;
        }

        VertexSPtr v0 = fh->vertex(0)->info();
        VertexSPtr v1 = fh->vertex(1)->info();
        VertexSPtr v2 = fh->vertex(2)->info();
        VertexSPtr verts[3] = {v0, v1, v2};
        FacetSPtr new_facet = Facet::create(3, verts);
        Plane3SPtr plane = KernelFactory::createPlane3(v0->getPoint(),
                                                       v1->getPoint(),
                                                       v2->getPoint());
        new_facet->setPlane(plane);
        new_facet->normalizePlaneCoefficients();
        SkelFacetDataSPtr new_data = SkelFacetData::create(new_facet);
        new_data->setSpeed(parent_speed);
        polyhedron->addFacet(new_facet);
        created_facets.push_back(new_facet);
    }

    polyhedron->initializeAllIDs();

    return { facet_vertices , created_facets };
}

// normalize coefficients
void PolyhedronTransformation::normalizeFacetPlanes(PolyhedronSPtr polyhedron)
{
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        facet->normalizePlaneCoefficients();
    }
}

bool PolyhedronTransformation::doAll2PlanesIntersect(PolyhedronSPtr polyhedron)
{
    bool result = true;
    std::list<FacetSPtr>::iterator it_f1 = polyhedron->facets().begin();
    while (it_f1 != polyhedron->facets().end()) {
        FacetSPtr facet1 = *it_f1++;
        std::list<FacetSPtr>::iterator it_f2 = it_f1;
        while (it_f2 != polyhedron->facets().end()) {
            FacetSPtr facet2 = *it_f2++;
            if (!KernelWrapper::intersection(facet1->plane(), facet2->plane())) {
                CGAL_SS3_TRANSF_TRACE("Degenerate facet pair:");
                CGAL_SS3_TRANSF_TRACE("  " << facet1->toString());
                CGAL_SS3_TRANSF_TRACE("  " << facet2->toString());
                result = false;
                break;
            }
        }
        if (!result) {
            break;
        }
    }

    return result;
}

bool PolyhedronTransformation::doAll3PlanesIntersect(PolyhedronSPtr polyhedron) {
    bool result = true;
    std::list<FacetSPtr>::iterator it_f1 = polyhedron->facets().begin();
    while (it_f1 != polyhedron->facets().end()) {
        FacetSPtr facet1 = *it_f1++;
        std::list<FacetSPtr>::iterator it_f2 = it_f1;
        while (it_f2 != polyhedron->facets().end()) {
            FacetSPtr facet2 = *it_f2++;
            std::list<FacetSPtr>::iterator it_f3 = it_f2;
            while (it_f3 != polyhedron->facets().end()) {
                FacetSPtr facet3 = *it_f3++;
                if (!KernelWrapper::intersection(facet1->plane(), facet2->plane(), facet3->plane())) {
                    CGAL_SS3_TRANSF_TRACE("Degenerate facet triplet:");
                    CGAL_SS3_TRANSF_TRACE("  " << facet1->toString());
                    CGAL_SS3_TRANSF_TRACE("  " << facet2->toString());
                    CGAL_SS3_TRANSF_TRACE("  " << facet3->toString());
                    result = false;
                    break;
                }
            }
            if (!result) {
                break;
            }
        }
        if (!result) {
            break;
        }
    }
    return result;
}

bool PolyhedronTransformation::areAllVerticesDegree3(PolyhedronSPtr polyhedron) {
    bool result = true;

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (vertex->degree() != 3) {
            CGAL_SS3_TRANSF_TRACE("High-degree vertex: " << vertex->toString());
            result = false;
            break;
        }
    }

    return result;
}

bool PolyhedronTransformation::isTiltCompatible(PolyhedronSPtr polyhedron) {
    bool result = true;

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        if (facet->numHighDegreeVertices() > 2) {
            CGAL_SS3_TRANSF_TRACE("facet " << facet->getID() << " has too many high-degree vertices "
                                           << "(" << facet->numHighDegreeVertices() << ")");
            result = false;
            break;
        }
    }

    return result;
}

std::array<double, 3> PolyhedronTransformation::randVec(double min, double max) {
    static std::random_device rd;
    unsigned int s = 0; // rd()
    // CGAL_SS3_TRANSF_TRACE("seed = " << s);
    static std::mt19937 gen(s);
    static std::uniform_real_distribution<> rdist(min, max);

    return { rdist(gen), rdist(gen), rdist(gen) };
}

void PolyhedronTransformation::randMovePoints(PolyhedronSPtr polyhedron) {
    double range = 0.001;
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        double value = config->getDouble("main", "rand_move_points_range");
        if (value != 0.0) {
            range = value;
        }
    }

    CGAL_SS3_TRANSF_TRACE("Points will be moved randomly...");
    CGAL_SS3_TRANSF_TRACE("  rand_move_points_range =" << range);

    // Store coefficients for possible un-tilting post-offset
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;

        // If we are doing random point perturbation, the mesh must be a triangle mesh
        // otherwise points will no longer be on the supporting planes of their incident facets
        if (facet->vertices().size() != 3) {
            std::cerr << "Warning: facet " << facet->getID() << " is not a triangle but we are displacing its vertices" << std::endl;
        }

        facet->storePlaneCoefficients();
    }

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;

        // since it's random, move to doubles to get static filters and avoid DAGs
        Point3SPtr p = vertex->getPoint();
        std::array<double, 3> v_r = randVec(-range/2.0, range/2.0);
        Point3SPtr p_t = KernelFactory::createPoint3(CGAL::to_double(p->x()) + v_r[0],
                                                     CGAL::to_double(p->y()) + v_r[1],
                                                     CGAL::to_double(p->z()) + v_r[2]);
        vertex->setPoint(p_t);
    }

    // recompute normalized planes to ensure points are on the supporting planes
    polyhedron->initPlanes();
    normalizeFacetPlanes(polyhedron);
    CGAL_assertion_code(bool success =)
    resetPoints(polyhedron);
    CGAL_assertion(success);
    CGAL_postcondition(polyhedron && polyhedron->isConsistent());

    polyhedron->appendDescription("rand_move_points_range=" +
            util::StringFactory::fromDouble(range) + "; ");
}

void PolyhedronTransformation::randTiltPlanes(PolyhedronSPtr polyhedron) {
    double range = 1e-10;
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        range = config->getDouble("main", "rand_move_points_range");
    }

    // If we only nudged planes with fixed point constraints, we might not ensure generic position,
    // for example if two pairs of constraints are along the same line.
    //
    // @todo could restrict to only high-degree vertices in facets that have 2 high-degree vertices
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        std::array<double, 3> v_r = randVec(-range, range); // @fixme hardcoded values

        double px = CGAL::to_double(p->x()) + v_r[0];
        double py = CGAL::to_double(p->y()) + v_r[1];
        double pz = CGAL::to_double(p->z()) + v_r[2];

        Point3SPtr p_t = KernelFactory::createPoint3(px, py, pz);
        vertex->setPoint(p_t);
    }

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        facet->perturbPlaneCoefficientsHighDegrees(range);
    }
}

void PolyhedronTransformation::randTiltPlanesv3(PolyhedronSPtr polyhedron) {
    CGAL_SS3_TRANSF_TRACE("Random Plane Tilt (v3)");

    double range = 1e-10;
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        range = config->getDouble("main", "rand_move_points_range");
    }

    CGAL_SS3_TRANSF_TRACE("  rand_move_points_range = " << range);

    if (PolyhedronTransformation::isTiltCompatible(polyhedron)) {
        CGAL_SS3_TRANSF_TRACE("Polyhedron can simply be tilted immediately");
        PolyhedronTransformation::randTiltPlanes(polyhedron);
        CGAL_assertion_code(bool success =)
        resetPoints(polyhedron);
        CGAL_assertion(success);
        return;
    }

#ifdef CGAL_SS3_DUMP_FILES
    db::_3d::OBJFile::save("results/tilt-v3_input.obj", polyhedron,
                           false /*do_triangulate*/,
                           true /*convert_to_double*/);
#endif

    CGAL_SS3_TRANSF_TRACE_CODE(unsigned int had_to_triangulate_n = 0;)

    // high-degree vertex --> first 3 incident facets determining the vertex
    std::map<VertexSPtr, std::set<FacetSPtr> > determining_facets;

    // facet --> first 2 determined high-degree vertices
    //
    // The facet becomes fixed at 2 vertices and not 3 vertices despite the vertices being perturbed
    // because if we do 3 random perturbations of vertices, the normal can vary wildly.
    //
    // Ideally, it could be fixed with 3 high-degree vertices and a smarter perturbation (which needs
    // to take into account all incident facets of these 3 fixing vertices...)
    std::map<FacetSPtr, std::set<VertexSPtr> > fixing_vertices;

#ifdef CGAL_SS3_DUMP_FILES
    auto dump_facet = [](std::string filename, FacetSPtr f) {
        auto pcdt = triangulate_facet_with_CDT2(f);

        using PCDT = decltype(pcdt);
        using PCDT_VH = typename PCDT::Vertex_handle;
        using PCDT_FH = typename PCDT::Face_handle;

        std::unordered_map<PCDT_FH, bool> in_domain_map;
        boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);
        CGAL::mark_domain_in_triangulation(pcdt, in_domain);

        std::map<PCDT_VH, std::size_t> point_to_id;
        std::vector<Point3> points;
        std::vector<std::vector<std::size_t> > triangles;
        for (PCDT_VH vh : pcdt.finite_vertex_handles()) {
            point_to_id[vh] = points.size();
            points.push_back(vh->point());
        }

        for (PCDT_FH fh : pcdt.finite_face_handles()) {
            if(!get(in_domain, fh)) {
                continue;
            }

            triangles.push_back({point_to_id[fh->vertex(0)],
                                 point_to_id[fh->vertex(1)],
                                 point_to_id[fh->vertex(2)]});
        }

        CGAL::IO::write_OFF(filename, points, triangles, CGAL::parameters::stream_precision(17));
    };
#endif

    auto is_high_degree = [&](VertexSPtr v) -> bool {
        return (v->degree() > 3);
    };

    auto has_high_degree_vertices = [](FacetSPtr f) -> bool {
        for (VertexSPtr v : f->vertices()) {
            if (v->degree() > 3) {
                return true;
            }
        }
        return false;
    };

    auto is_vertex_determined = [&](VertexSPtr v) -> bool {
        CGAL_SS3_TRANSF_TRACE_V(16, "Checking if V" << v->getID() << " (deg: " << v->degree() << ") is fixed");
        auto it = determining_facets.find(v);
        return (it != determining_facets.end() && it->second.size() == 3);
    };

    auto is_facet_fixed = [&](FacetSPtr f) -> bool {
        CGAL_SS3_TRANSF_TRACE_V(16, "Checking if F" << f->getID() << " (" << f->vertices().size() << " nv) is fixed");
        CGAL_SS3_TRANSF_TRACE_V(16, "  fixing_vertices size: " << fixing_vertices[f].size());
        CGAL_assertion(fixing_vertices[f].size() <= 3);
        return (f->isTriangle() && fixing_vertices[f].size() == 3) ||
               (!f->isTriangle() && fixing_vertices[f].size() == 2);
    };

#ifdef CGAL_SS3_DUMP_FILES
    unsigned int visited_face_id = 0;
    unsigned int nudged_face_id = 0;
#endif

    // Sort by number of high-degree vertices as to avoid triangulating as much as possible
    auto facet_sorter = [&](FacetSPtr a, FacetSPtr b)
    {
        auto hdv_count = [&](FacetSPtr f) -> unsigned int {
            unsigned int hdv_n = 0;
            for (VertexSPtr v : f->vertices()) {
                if (is_high_degree(v)) {
                    ++hdv_n;
                }
            }
            return hdv_n;
        };

        // Give priority to facets with no determined vertices.
        // If both or neither have constrained vertices, give priority to the largest hdv count.
        //
        // The point is to avoid cascading exact number types, even if we have to triangulate a little more
        auto get_determined_count = [&](FacetSPtr f) -> unsigned int
        {
          unsigned int res = 0;
          for (VertexSPtr v : f->vertices()) {
              if (is_vertex_determined(v)) {
                  ++res;
              }
          }
          return res;
        };

        unsigned int a_determined_n = get_determined_count(a);
        unsigned int b_determined_n = get_determined_count(b);

        CGAL_SS3_TRANSF_TRACE_V(64, "F" << a->getID() << " has " << a_determined_n << " determined vertices");
        CGAL_SS3_TRANSF_TRACE_V(64, "F" << b->getID() << " has " << b_determined_n << " determined vertices");

        if (a_determined_n != b_determined_n) {
            // Give priority to the one with the least determined vertices
            return a_determined_n < b_determined_n;
        }

        // same number of determined vertices, give priority to the facet with the most high-degree vertices
        unsigned int a_hdv_n = hdv_count(a);
        unsigned int b_hdv_n = hdv_count(b);

        CGAL_SS3_TRANSF_TRACE_V(64, "F" << a->getID() << " has " << a_hdv_n << " high-degree vertices");
        CGAL_SS3_TRANSF_TRACE_V(64, "F" << b->getID() << " has " << b_hdv_n << " high-degree vertices");

        if (a_hdv_n != b_hdv_n) {
            // Give priority to the one with the most high-degree vertices
            return a_hdv_n > b_hdv_n;
        }

        // same number of determined vertices and high-degree vertices, give priority to the largest facet
        return a->vertices().size() > b->vertices().size();
    };

    // If the facet has no high-degree vertices, we can just tilt it randomly and it will
    // be fine because by definition all of its vertices are degree 3 and will be stable
    // because an unstable configuration results from almost coplanar facets, which have been
    // merged ahead of randomization.
    // UNLESS we have to triangulate a face incident to one vertex of this face without
    // high-degree vertices and then the face now has a high-degree vertex. If that happens,
    // we want the high-degree vertex to constrain to the (up to 2) facets with no high-degree
    // vertices which we are constraining here ahead of the flooding process.
    // Hence, we mark this as fixed with dummy vertices and add 'v' (a non high-degree vertex)
    // to the 'determining_facets' map.
    for (FacetSPtr facet : polyhedron->facets()) {
        if (!has_high_degree_vertices(facet)) {
            CGAL_SS3_TRANSF_TRACE("Nudge and fix F" << facet->getID());
            facet->perturbPlaneCoefficientsNudge(range);

#ifdef CGAL_SS3_DUMP_FILES
            dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_low_degree.OFF", facet);
#endif

            // A low degree facet is a constraining place when nudging a vertex incident to it.
            // Use dummy vertices to get that effect.
            for (VertexSPtr v : facet->vertices()) {
                fixing_vertices[facet].insert(v);
                if (is_facet_fixed(facet)) {
                      break;
                }
            }

            // the point of this is that if the vertex becomes high degree after triangulation,
            // one (or two) facet with low degree vertices will appear in the determining facets
            for (VertexSPtr v : facet->vertices()) {
                determining_facets[v].insert(facet);
                CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " is determined by F" << facet->getID() << " (1)");
            }
        }
    }

    // This is the main list of facets that we will process
    std::list<FacetSPtr> facets_to_process;
    for (FacetSPtr facet : polyhedron->facets()) {
        if (facet->isTriangle() || !has_high_degree_vertices(facet)) {
            continue;
        }

        facets_to_process.push_back(facet);
    }

    // Some preprocessing: if two faces share more than 2 high-degree vertices, we have to triangulate
    // one of them to ensure generic positioning.
    //
    // We triangulate rather than seemingly smarter method of splitting the face because
    // perturbing splitted facets is very difficult because they are coplanar and perturbations
    // are thus unstable unless performed around the splitting edge, which adds a ton of constraints
    // and complexity.
    for(;;) // reset every time we triangulate something to avoid needless subdivisions
    {
        bool did_something = false;

        std::list<FacetSPtr> facets_to_exclude;
        for (FacetSPtr f : facets_to_process) {
            for (EdgeSPtr e : f->edges()) {
                VertexSPtr sv = e->src(f);
                VertexSPtr tv = e->dst(f);
                if (sv->degree() == 3 || tv->degree() == 3) {
                    continue;
                }

                // Find the facets { f' } which appear in both sets of incident facets for the vertices
                std::set<FacetSPtr> sv_facets, tv_facets, common_facets;
                for (FacetWPtr wf : sv->facets()) {
                    if (FacetSPtr fptr = wf.lock()) {
                        sv_facets.insert(fptr);
                    }
                }

                for (FacetWPtr wf : tv->facets()) {
                    if (FacetSPtr fptr = wf.lock()) {
                        tv_facets.insert(fptr);
                    }
                }

                // Find intersection (common facets)
                std::set_intersection(sv_facets.begin(), sv_facets.end(),
                                      tv_facets.begin(), tv_facets.end(),
                                      std::inserter(common_facets, common_facets.begin())
                );

                for (FacetSPtr fprime : common_facets) {
                    if (fprime == f) {
                        continue;
                    }

                    bool has_edge = (tv->next(fprime) == sv);
                    if (!has_edge) {
                        // Mark for triangulation
                        facets_to_exclude.push_back(fprime);

                        CGAL_SS3_TRANSF_TRACE("Facet F" << fprime->getID() << " needs triangulating due to missing high-degree edge between V" << sv->getID() << " and V" << tv->getID());

                        CGAL_SS3_TRANSF_TRACE("Triangulate F" << fprime->getID());
                        CGAL_SS3_TRANSF_TRACE_CODE(++had_to_triangulate_n;)

                        PolyhedronTransformation::triangulate(fprime, polyhedron);

                        did_something = true;
                        break;
                    }
                }
                if (did_something) {
                    break;
                }
            }
            if (did_something) {
                break; // restart
            }
        }

        for (FacetSPtr f : facets_to_exclude) {
            facets_to_process.remove(f);
        }

        if (!did_something) {
            break;
        }
    }

#ifdef CGAL_SS3_DUMP_FILES
    db::_3d::OBJFile::save("results/tilt-v3_preprocessed.obj", polyhedron,
                           false /*do_triangulate*/,
                           true /*convert_to_double*/);
#endif

    // Forward declarations for mutually recursive lambdas
    std::function<void(FacetSPtr, VertexSPtr)> add_fixing_vertex;
    std::function<void(VertexSPtr)> determine_vertex;

    auto nudge_constrained_vertex = [&](VertexSPtr v)
    {
        CGAL_SS3_TRANSF_TRACE("  Nudging V" << v->getID() << " from " << *(v->getPoint()));

        std::vector<Plane3SPtr> constraining_planes;
        for (FacetSPtr df : determining_facets[v]) {
            if (is_facet_fixed(df)) {
                constraining_planes.push_back(df->getPlane());
                CGAL_SS3_TRANSF_TRACE("  F" << df->getID() << " constrains the nudge");
            }
        }

        CGAL_assertion(constraining_planes.size() <= 3);

        const size_t n_fixed = constraining_planes.size();
        if (n_fixed == 3) {
            resetPoint(v, { constraining_planes[0], constraining_planes[1], constraining_planes[2] });
            CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " reset to " << *v->getPoint());
            return;
        }

        Point3SPtr p = v->getPoint();

        // @todo clean of all that stuff + duplicate code in facet.cpp
        static std::random_device rd;
        unsigned int s = 0; // rd()
        // CGAL_SS3_TRANSF_TRACE("seed = " << s);
        static std::mt19937 gen(s);
        static std::uniform_real_distribution<> rdist(-range, range);

        auto nudge = [&](const CGAL::FT& v) {
            // Since we are perturbing, we might as well collapse the DAG of 'v'.
            // the point is also that once 'nv' is a double, its interval will be a singleton,
            // and we will have access to static filters
            double step = rdist(gen);
            double nv = CGAL::to_double(v) + step;
            return nv;
        };

        auto nudge_to_simplest_rational_in_interval = [&](const CGAL::FT& v) {
            double d1 = nudge(v);
            double d2 = nudge(v);
            if (d2 < d1) {
                std::swap(d1, d2);
            }
            CGAL::FT nv = CGAL::simplest_rational_in_interval<CGAL::K::Exact_kernel::FT>(d1, d2);
            return nv;
        };

#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        This does not look quite ready yet: sometimes, the projections are way off
        and with this, sometimes we produce polyhedra in degenerate positions, probably
        because the smallest rational is the same despite the random intervals.
#endif

#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        CGAL::FT x = nudge_to_simplest_rational_in_interval(p->x());
        CGAL::FT y = nudge_to_simplest_rational_in_interval(p->y());
        CGAL::FT z = nudge_to_simplest_rational_in_interval(p->z());
#else
        std::array<double, 3> v_r = randVec(-range/2.0, range/2.0);
        double x = CGAL::to_double(p->x()) + v_r[0];
        double y = CGAL::to_double(p->y()) + v_r[1];
        double z = CGAL::to_double(p->z()) + v_r[2];
#endif
        Point3SPtr p_nudged = KernelFactory::createPoint3(x, y, z);
        CGAL_SS3_TRANSF_TRACE("base nudge: " << x << " " << y << " " << z);

        Point3SPtr p_new;

        if (n_fixed == 0) {
            p_new = p_nudged;
        } else if (n_fixed == 1) {
            Plane3SPtr plane = constraining_planes[0];
#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
            // something similar but a little more subtle:
            // 1. project the point onto the plane
            // 2. express the point as a linear combination of the plane's origin and basis: pp = o + l1 * b1 + l2 * b2
            // 3. nudge l1 and l2 to l1' and l2' with a random interval around l1 and l2, and
            //    simplest_rational_in_interval
            // 4. recompute the point as pp = o + l1' * b1 + l2' * b2
            Point3SPtr pp = KernelWrapper::projection(plane, p_nudged);
            const Point3& o = plane->point();
            const Vector3& b1 = plane->base1();
            const Vector3& b2 = plane->base2();
            CGAL::FT l1 = CGAL::scalar_product(*pp - o, b1);
            CGAL::FT l2 = CGAL::scalar_product(*pp - o, b2);
            CGAL::FT nl1 = nudge_to_simplest_rational_in_interval(l1);
            CGAL::FT nl2 = nudge_to_simplest_rational_in_interval(l2);
            p_new = KernelFactory::createPoint3(o.x() + nl1 * b1.x() + nl2 * b2.x(),
                                                o.y() + nl1 * b1.y() + nl2 * b2.y(),
                                                o.z() + nl1 * b1.z() + nl2 * b2.z());
#else
            p_new = KernelWrapper::projection(plane, p_nudged);
#endif
        } else if (n_fixed == 2) {
            Plane3SPtr plane1 = constraining_planes[0];
            Plane3SPtr plane2 = constraining_planes[1];
            Line3SPtr line = KernelWrapper::intersection(plane1, plane2);
#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
            // something similar but a little more subtle:
            // 1. project the point onto the line
            // 2. express the point as a linear combination of the line's origin and basis: pp = o + l * v
            // 3. nudge l to l' with a random interval around l, and simplest_rational_in_interval
            // 4. recompute the point as pp = o + l' * v
            Point3SPtr pp = KernelWrapper::projection(line, p_nudged);
            const Point3& o = line->point();
            const Vector3& d = line->to_vector();
            CGAL::FT l = CGAL::scalar_product(*pp - o, d);
            CGAL::FT nl = nudge_to_simplest_rational_in_interval(l);
            p_new = KernelFactory::createPoint3(o.x() + nl * d.x(),
                                                o.y() + nl * d.y(),
                                                o.z() + nl * d.z());
#else
            p_new = KernelWrapper::projection(line, p_nudged);
#endif
        }

        CGAL_SS3_TRANSF_TRACE("  Nudged V" << v->getID() << " to " << *p_new);

        v->setPoint(p_new);
    };

    // The face has enough determined vertices to be fixed. Compute its random perturbation,
    // and add the facet ID to its vertices
    add_fixing_vertex = [&](FacetSPtr f, VertexSPtr v)
    {
        CGAL_precondition(fixing_vertices[f].size() <= 3);

        if (is_facet_fixed(f)) {
            CGAL_SS3_TRANSF_TRACE("  F" << f->getID() << " is already fully fixed");
            return;
        }

        CGAL_SS3_TRANSF_TRACE("  Fix F" << f->getID() << " with V" << v->getID());
        fixing_vertices[f].insert(v);

        if (!is_facet_fixed(f)) {
            // nothing to do yet, there are still degrees of freedom in the facet
            return;
        }

        CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
        CGAL_SS3_TRANSF_TRACE_CODE(ss << "F" << f->getID() << " is now fully fixed by");
        CGAL_SS3_TRANSF_TRACE_CODE(for (VertexSPtr fv : fixing_vertices[f]))
        CGAL_SS3_TRANSF_TRACE_CODE(ss << " V" << fv->getID());
        CGAL_SS3_TRANSF_TRACE(ss.str());

        if (f->isTriangle()) {
            CGAL_assertion(fixing_vertices[f].size() == 3); // just to be clear

            // for triangles, all vertices are determined, and there is nothing to nudge
            // (note that vertices were themselves nudged so the facet is nudged).
            f->initPlane();
            f->normalizePlaneCoefficients();

#ifdef CGAL_SS3_DUMP_FILES
            dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_fixed_3.OFF", f);
#endif

            // Here we do not need to add the fixed face to incident determined vertices
            // because all vertices are all fully determined
            return;
        }

        std::vector<Point3SPtr> fixed_points;
        for (VertexSPtr fixed_v : fixing_vertices[f]) {
            fixed_points.push_back(fixed_v->getPoint());
        }

        CGAL_SS3_TRANSF_TRACE("Nudge and fix F" << f->getID() << " with " << fixed_points.size() << " fixed point(s)");
        f->perturbPlaneCoefficientsFixedPoints(range, fixed_points);

#ifdef CGAL_SS3_DUMP_FILES
        dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_fixed_" + std::to_string(fixed_points.size()) + ".OFF", f);
#endif

        CGAL_SS3_TRANSF_TRACE("Newly fixed face determines its high-degree incident vertices...");

        // Need to now tag the vertices of the facet
        for (VertexSPtr v : f->vertices()) {
            std::cout << "incident: " << v->getID() << " (deg=" << v->degree() << "; " << determining_facets[v].size() << " determining facets)" << std::endl;
            if (!is_vertex_determined(v)) {
                determining_facets[v].insert(f);
                CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is determined by F" << f->getID() << " (2)");
            }
        }
    };

    determine_vertex = [&](VertexSPtr v)
    {
        CGAL_precondition(is_high_degree(v));
        CGAL_precondition(is_vertex_determined(v));

        CGAL_SS3_TRANSF_TRACE_CODE(auto it = determining_facets[v].begin();)
        CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " is now fully determined by"
                              << " F" << (*it)->getID() << " depth=" << CGAL::ss_utils::length(*((*it)->getPlane()))
                              << " F" << (*std::next(it))->getID() << " depth=" << CGAL::ss_utils::length(*((*std::next(it))->getPlane()))
                              << " F" << (*std::next(it, 2))->getID() << " depth=" << CGAL::ss_utils::length(*((*std::next(it, 2))->getPlane())));

        // set the nudged position for the vertex: a nudge constrained by already fixed incident facets
        nudge_constrained_vertex(v);

        CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " is now fixed at " << *(v->getPoint()) << " depth=" << CGAL::ss_utils::length(*(v->getPoint())));

        // compute the plane coefficients of any incident facet that becomes fixed
        // by this vertex becoming determined
        for (FacetWPtr wf : v->facets()) {
            if (FacetSPtr f = wf.lock()) {
                add_fixing_vertex(f, v);
            }
        }
    };

    auto is_facet_overconstrained = [&](FacetSPtr f) -> bool
    {
        if (f->isTriangle() || is_facet_fixed(f)) {
            return false;
        }

        // force triangulation if the exact stack is getting too deep
#if 0
        for (VertexSPtr v : f->vertices()) {
            if (CGAL::ss_utils::length(*(v->getPoint())) > 100) {
                std::cout << "stack is too deep, triangulate the face" << std::endl;
                return true;
            }
        }
#endif

        // we cannot handle that facet without triangulation if adding the facet to high-degree
        // vertices would create too many determined vertices (> 2) in any facet incident to
        // the determined high-degree vertices of this facet
        std::map<FacetSPtr, unsigned int> facets_to_test; // facets + number of appearances
        for (VertexSPtr hdv : f->vertices()) {
            if (is_high_degree(hdv)) {
                for (FacetWPtr inc_f : hdv->facets()) {
                    FacetSPtr f = inc_f.lock();
                    if (f) {
                        ++facets_to_test[f];
                    }
                }
            }
        }

        for (const auto& [ft, count] : facets_to_test) {
            // Count the number of high-degree vertices with either:
            // - 3 determining facets
            // - 2 determining facets and incident to 'facet'
            // These are vertices that are determined, or would be determined once we 'add'
            // the facet to its high-degree vertices.
            unsigned int constrain_n = 0;
            for (VertexSPtr v : f->vertices()) {
                if (is_high_degree(v)) {
                    if (is_vertex_determined(v)) {
                        ++constrain_n;
                    } else if (determining_facets[v].size() == 2 && ft->hasVertex(v)) {
                        ++constrain_n;
                    }
                }

                if(constrain_n > 2) {
                    CGAL_SS3_TRANSF_TRACE("F" << ft->getID() << " would be over constrained by fixing of F" << f->getID());
                    return true;
                }
            }
        }

        return false;
    };

    CGAL_SS3_TRANSF_TRACE("== Main facet flood... ==");

    while (!facets_to_process.empty()) {
        facets_to_process.sort(facet_sorter);
        FacetSPtr facet = facets_to_process.front();
        facets_to_process.pop_front();

        CGAL_SS3_TRANSF_TRACE("Pop F" << facet->getID());

        CGAL_assertion(!facet->isTriangle());
        CGAL_assertion(fixing_vertices[facet].size() <= 2);

        CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
        CGAL_SS3_TRANSF_TRACE_CODE(ss << "  Fixing vertices:";)
        CGAL_SS3_TRANSF_TRACE_CODE(for (VertexSPtr fv : fixing_vertices[facet]))
        CGAL_SS3_TRANSF_TRACE_CODE(ss << " V" << fv->getID();)
        CGAL_SS3_TRANSF_TRACE(ss.str());

#ifdef CGAL_SS3_DUMP_FILES
        dump_facet("results/visited_face_" + std::to_string(visited_face_id++) + ".OFF", facet);
#endif

        CGAL_assertion(facet->vertices().size() >= 3);

        if (is_facet_overconstrained(facet)) {
            std::list<FacetSPtr> facets_to_triangulate = { facet };
            std::list<FacetSPtr> final_facets;
            while (!facets_to_triangulate.empty()) {
                FacetSPtr facet_tt = facets_to_triangulate.front();
                facets_to_triangulate.pop_front();

                // the facet is not yet fixed, so no vertex can have it as determining facet
                CGAL_assertion(fixing_vertices[facet_tt].size() < 2);
                CGAL_assertion_code(for (VertexSPtr v : facet_tt->vertices()) {)
                CGAL_assertion(determining_facets[v].size() <= 3);
                CGAL_assertion(determining_facets[v].count(facet_tt) == 0);
                CGAL_assertion_code(})

                CGAL_SS3_TRANSF_TRACE("Triangulate F" << facet_tt->getID());
                CGAL_SS3_TRANSF_TRACE_CODE(++had_to_triangulate_n;)

                auto [local_vertices, new_facets] = PolyhedronTransformation::triangulate(facet_tt, polyhedron);

                // vertices that have become high degree should have their determination information updated
                for (VertexSPtr v : local_vertices) {
                    std::cout << "local vertex " << v->getID() << " (deg=" << v->degree() << "; " << determining_facets[v].size() << " determining facets)" << std::endl;

                    for (FacetWPtr wf : v->facets()) {
                        if (FacetSPtr fptr = wf.lock()) {
                            if (is_facet_fixed(fptr)) {
                                CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is determined by F" << fptr->getID() << " (3)");
                                determining_facets[v].insert(fptr);

                                if (is_vertex_determined(v)) {
                                    break;
                                }
                            }
                        }
                    }
                }

                // existing already-determined vertices are fixed points for the new facets
                for (FacetSPtr nf : new_facets) {
                    CGAL_SS3_TRANSF_TRACE("spawned F" << nf->getID());

                    for (VertexSPtr v : nf->vertices()) {
                        if (is_vertex_determined(v)) {
                            CGAL_SS3_TRANSF_TRACE("newborn F" << nf->getID() << " is constrained by V" << v->getID());
                            fixing_vertices[nf].insert(v);
                        }
                    }

                    CGAL_postcondition(!is_facet_overconstrained(nf));
                }
            }

            continue;
        }

        // Now, adding the facet to the high-degree vertices will not over constrain the facet, so do it:
        for (VertexSPtr v : facet->vertices()) {
            if (!is_vertex_determined(v)) {
                determining_facets[v].insert(facet);
                CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is determined by F" << facet->getID() << " (4)");
                if (is_high_degree(v) && is_vertex_determined(v)) {
                    // When the vertex becomes fixed (its 3 determining facets become known), we need:
                    // - to perturb the position of the vertex
                    // - to update all incident facets to check if they are now fixed and in that case,
                    //   compute their plane coefficients
                    determine_vertex(v);
                }
            }
        }
    }

    // Some facets might have high-degree vertices, but still some freedom of movement after the flooding, fix them
    for (FacetSPtr f : polyhedron->facets()) {
        // if there were 2, it would have been nudged when the 2nd vertex got determined
        if (f->isTriangle() || fixing_vertices[f].size() == 2) {
            continue;
        }

        CGAL_SS3_TRANSF_TRACE("Nudge and fix F" << f->getID() << " [remaining]");

        std::vector<Point3SPtr> fixed_points;
        for (VertexSPtr v : fixing_vertices[f]) {
            CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is a fixing vertex");
            fixed_points.push_back(v->getPoint());
        }

        f->perturbPlaneCoefficientsFixedPoints(range, fixed_points);

        for (VertexSPtr v : f->vertices()) {
            // @todo could avoid these checks with a boost::multi_index_container
            // to have a set that is sorted by insertion order
            if (!is_vertex_determined(v)) {
                determining_facets[v].insert(f);
                CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is determined by F" << f->getID() << " (5)");
            }
        }

        for (VertexSPtr v : f->vertices()) {
            fixing_vertices[f].insert(v);
            if (is_facet_fixed(f)) {
                break;
            }
        }

#ifdef CGAL_SS3_DUMP_FILES
        dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_remaining.OFF", f);
#endif
    }

    // Nudge what can still be nudged, for randomness only
    CGAL_SS3_TRANSF_TRACE("== Nudge unconstrained high-degree vertices... ==");

    // Simply move high-degree vertices that are not yet constrained by 3 facets,
    // and recomputing triangle facets afterwards.
    for (VertexSPtr v : polyhedron->vertices()) {
        if (is_high_degree(v) && determining_facets[v].size() < 3) {
            CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is high degree and not fully determined, nudge it");
            nudge_constrained_vertex(v);

            // determine the vertex
            // since we know only triangle facets are left, we don't need to cascade and check
            // if incident facets become fixed
            for (FacetWPtr wf : v->facets()) {
                if (is_vertex_determined(v)) {
                    CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is now determined");
                    break;
                }
                if (FacetSPtr f = wf.lock()) {
                    determining_facets[v].insert(f);
                }
            }
        }
    }

    // Now handle triangle faces with high degrees
    CGAL_SS3_TRANSF_TRACE("== Deal with remaining triangles... ==");

    for (FacetSPtr f : polyhedron->facets()) {
        if (!f->isTriangle() || !has_high_degree_vertices(f)) {
            continue;
        }

        CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
        CGAL_SS3_TRANSF_TRACE_CODE(ss << "Fix triangle F" << f->getID() << " [");
        CGAL_SS3_TRANSF_TRACE_CODE(auto it = f->vertices().begin();)
        CGAL_SS3_TRANSF_TRACE_CODE(ss << "V" << (*it++)->getID() << " ");
        CGAL_SS3_TRANSF_TRACE_CODE(ss << "V" << (*it++)->getID() << " ");
        CGAL_SS3_TRANSF_TRACE_CODE(ss << "V" << (*it)->getID() << "]");
        CGAL_SS3_TRANSF_TRACE(ss.str());

        // Triangle facets are great because we have freedom to move points and compute planes afterwards
        f->initPlane();
        f->normalizePlaneCoefficients();

#ifdef CGAL_SS3_DUMP_FILES
        dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_triangle.OFF", f);
#endif

        // Here we need to also need to update the determining facets because some neighboring
        // facets could be unfixed high-degree triangles
        for (VertexSPtr v : f->vertices()) {
            if (!is_vertex_determined(v)) {
                determining_facets[v].insert(f);
                CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is determined by F" << f->getID() << " (6)");
                // no need to cascade here, we know only triangles are left
            }
        }
    }

#ifdef CGAL_SS3_DUMP_FILES
    db::_3d::OBJFile::save("results/tilt_v3-pre_reset.obj", polyhedron,
                           false /*do_triangulate*/,
                           true /*convert_to_double*/);
#endif

    CGAL_SS3_TRANSF_TRACE("Reset the position of not-fully-constrained vertices...");

    // Recompute all points which were not fixed (degree 3 vertices)
    for (VertexSPtr v : polyhedron->vertices()) {
        if (!is_high_degree(v)) {
            CGAL_SS3_TRANSF_TRACE("Reset V" << v->getID());
            resetPoint(v);

            // determine (without cascading)
            for (FacetWPtr wf : v->facets()) {
                if (FacetSPtr f = wf.lock()) {
                    determining_facets[v].insert(f);
                }
            }
        } else {
            CGAL_assertion(is_vertex_determined(v)); // high degree vertices have already been determined
        }
    }

    CGAL_SS3_TRANSF_TRACE("All facets processed");

#ifdef CGAL_SS3_DUMP_FILES
    db::_3d::OBJFile::save("results/tilt_v3.obj", polyhedron,
                           false /*do_triangulate*/,
                           true /*convert_to_double*/);
    db::_3d::OBJFile::save("results/tilt_v3-triangulated.obj", polyhedron,
                           true /*do_triangulate*/,
                           true /*convert_to_double*/);
#endif

    CGAL_assertion_code(for (VertexSPtr v : polyhedron->vertices()) {)
    CGAL_assertion(is_vertex_determined(v));
    CGAL_assertion_code(})

    CGAL_assertion_code(for (FacetSPtr f : polyhedron->facets()) {)
    CGAL_assertion(fixing_vertices[f].size() <= 3);
    CGAL_assertion_code(})

    CGAL_assertion_code(for (FacetSPtr facet : polyhedron->facets()) {)
    CGAL_assertion_code(for (VertexSPtr v : facet->vertices()) {)
    CGAL_assertion_code(std::cout << "check V" << v->getID() << " on F" << facet->getID() << std::endl;)
    CGAL_assertion(facet->getPlane()->has_on(*(v->getPoint())));
    CGAL_assertion_code(})
    CGAL_assertion_code(})

    CGAL_SS3_TRANSF_TRACE("Had to triangulate " << had_to_triangulate_n << " facets");

    CGAL_SS3_TRANSF_TRACE_CODE(for (VertexSPtr v : polyhedron->vertices()))
    CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " has depth " << CGAL::depth(*(v->getPoint())));

    CGAL_SS3_TRANSF_TRACE_CODE(for (FacetSPtr f : polyhedron->facets()) )
    CGAL_SS3_TRANSF_TRACE("F" << f->getID() << " has depth " << CGAL::depth(*(f->getPlane())));

    CGAL_SS3_TRANSF_TRACE_CODE(for (VertexSPtr v : polyhedron->vertices()))
    CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " has length " << CGAL::ss_utils::length(*(v->getPoint())));

    CGAL_SS3_TRANSF_TRACE_CODE(for (FacetSPtr f : polyhedron->facets()) )
    CGAL_SS3_TRANSF_TRACE("F" << f->getID() << " has length " << CGAL::ss_utils::length(*(f->getPlane())));
}

} }
