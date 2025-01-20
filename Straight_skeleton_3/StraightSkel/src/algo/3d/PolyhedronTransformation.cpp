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
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/skel/SkelVertexData.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "data/3d/skel/SkelFacetData.h"
#include "db/3d/OBJFile.h"

#include "util/StringFactory.h"
#include "util/Configuration.h"

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <list>
#include <limits>

namespace algo { namespace _3d {

PolyhedronTransformation::PolyhedronTransformation() {
    // intentionally does nothing.
}

PolyhedronTransformation::~PolyhedronTransformation() {
    // intentionally does nothing.
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

// @todo this function cannot deal with degree 1 vertices
Point3SPtr PolyhedronTransformation::shiftPoint(VertexSPtr vertex,
                                                CGAL::FT offset)
{
    // std::cout << "shift vertex " << vertex->toString() << "\noffset = " << offset << std::endl;

    Plane3SPtr planes[3];
    unsigned int i = 0;
    std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (i < 3 && it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (!facet_wptr.expired()) {
            FacetSPtr facet = FacetSPtr(facet_wptr);
            Plane3SPtr plane = facet->plane();
            // std::cout << "Facet #" << facet->toString() << std::endl;

            // @fixme only valid if we are doing perturbations
            if (vertex->degree() != 3) {
                // planes are _offset_ planes, but it doesn't matter for the tests
                bool independent = true;
                if (i == 1) {
                    independent = !(CGAL::parallel(*(planes[0]), *plane));
                } else if (i == 2) {
                    // @todo a little redundant to compute the intersection from scratch later
                    independent = !is_zero(CGAL::determinant(planes[0]->a(), planes[0]->b(), planes[0]->c(),
                                                            planes[1]->a(), planes[1]->b(), planes[1]->c(),
                                                            plane->a(), plane->b(), plane->c()));
                }

                if (!independent) {
                    continue;
                }
            }

            CGAL::FT speed = 1.0;
            if (facet->hasData()) {
                speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
            }

            planes[i] = KernelWrapper::offsetPlane(plane, offset*speed);
            // std::cout << "facet[" << i << "] = " << facet->toString() << std::endl;
            // std::cout << "  Offset Plane[" << i << "] = " << *(planes[i]) << std::endl;

            ++i;
        }
    }

    if (i < 3) {
        std::cerr << "Warning: Couldn't find three independent planes for vertex " << *(vertex->getPoint()) << std::endl;
        return { };
    }

    Point3SPtr point = KernelWrapper::intersection(planes[0], planes[1], planes[2]);
    if (!point) {
        std::cerr << "Error: triplet of planes doesn't define a point!" << std::endl;
        Point3SPtr result = Point3SPtr();
        DEBUG_SPTR(result);
        return result;
    }

    CGAL_assertion_code(for(Plane3SPtr pi : planes))
    CGAL_assertion(pi->has_on(*point));

    return point;
}

Plane3SPtr PolyhedronTransformation::shiftPlane(FacetSPtr facet,
                                                CGAL::FT offset)
{
    CGAL::FT facet_speed = 1.0;
    if (facet->hasData()) {
        facet_speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
    }

    return KernelWrapper::offsetPlane(facet->plane(), facet_speed*offset);
}

// @todo plenty of needless recomputations
PolyhedronSPtr PolyhedronTransformation::shiftFacets(PolyhedronSPtr polyhedron,
                                                     CGAL::FT offset,
                                                     const bool recompute_positions)
{
    std::cout << "~~~~ Shift by " << offset << std::endl;

    std::ofstream shift_out("results/last_shift.polylines.txt");
    shift_out.precision(17);

    PolyhedronSPtr result = Polyhedron::create();

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr point = vertex->getPoint();

        // those are treated in the next loop
        if (vertex->degree() == 1) {
            continue;
        }

        if (offset != 0 || recompute_positions) {
            point = shiftPoint(vertex, offset);
            if (!point) {
                std::cerr << "Warning: Failed to create shifted polyhedron" << std::endl;
                return { };
            }
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

        // std::cout << *(vertex->getPoint()) << " to " << *point << std::endl;
        shift_out << "2 " << *(vertex->getPoint()) << " " << *point << std::endl;
    }

    // Now, deal with degree 1 vertices.
    //
    // Do NOT merge the two vertices loops together: this function assumes that degree 1 vertices
    // are adjacent to degree 3+ vertices that have been offset in the first loop, and offsets
    // degree 1 vertices using the offset of the adjacent degree 3+ vertex.
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

            // std::cout << *(vertex->getPoint()) << " to " << *point << std::endl;
            shift_out << "2 " << *(vertex->getPoint()) << " " << *point << std::endl;
        } else {
            CGAL_assertion_msg(i != 0, "no edge on degree 1 vertex?");
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

// normalize coefficients
void PolyhedronTransformation::normalizeFacetPlanes(PolyhedronSPtr polyhedron)
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

// normalize coefficients and ensure that coplanar facets have the same coefficients
//
// @todo no point doing this with EPECK since we have SQRT errors
// and if the planes were singletons, we would benefit from static filters
void PolyhedronTransformation::harmonizeFacetPlanes(PolyhedronSPtr polyhedron)
{
    std::cout << "\n> harmonizeFacetPlanes()" << std::endl;

    // @todo also harmonize parallel but non equal planes (i.e., opposite normal directions)

    // @todo this needs to be rewritten to use kernel predicates to sort the planes
    // Order all the supporting planes in a global order
    // The order is completely arbitrary, the only thing that we care about
    // is to quickly find which planes are identical
    //
    // for now, abuse unordered containers because it's less code to write than above
    struct FEqual
    {
        bool operator()(FacetSPtr facet_1, FacetSPtr facet_2) const {
            Plane3SPtr plane_1 = facet_1->plane(); // calls initPlane() if needed
            Plane3SPtr plane_2 = facet_2->plane();
            // std::cout << "facet #" << facet_1->getID() << " VS facet #" << facet_2->getID()
            //          << " --> " << (*plane_1 == *plane_2) << std::endl;
            return CGAL::parallel(*plane_1, *plane_2);
        }
    };

    struct FHash
    {
        std::size_t operator()(FacetSPtr) const {
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
# define CGAL_SS3_SINGLETON_PLANE_COEFFICIENTS
# ifdef CGAL_SS3_SINGLETON_PLANE_COEFFICIENTS
            const CGAL::FT a = CGAL::to_double(exact(plane->a()));
            const CGAL::FT b = CGAL::to_double(exact(plane->b()));
            const CGAL::FT c = CGAL::to_double(exact(plane->c()));
            const CGAL::FT d = CGAL::to_double(exact(plane->d()));
            // this should be the only place with unavoidable SQRTs
            const CGAL::FT n = CGAL::approximate_sqrt(CGAL::square(a) + CGAL::square(b) + CGAL::square(c));
# else
            const CGAL::FT a = plane->a();
            const CGAL::FT b = plane->b();
            const CGAL::FT c = plane->c();
            const CGAL::FT d = plane->d();
            // this should be the only place with unavoidable SQRTs
            const CGAL::FT n = CGAL::approximate_sqrt(CGAL::square(a) + CGAL::square(b) + CGAL::square(c));
# endif

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

bool PolyhedronTransformation::hasParallelPlanes(PolyhedronSPtr polyhedron) {
    bool result = false;
    std::list<FacetSPtr>::iterator it_f1 = polyhedron->facets().begin();
    while (it_f1 != polyhedron->facets().end()) {
        FacetSPtr facet1 = *it_f1++;
        std::list<FacetSPtr>::iterator it_f2 = it_f1;
        while (it_f2 != polyhedron->facets().end()) {
            FacetSPtr facet2 = *it_f2++;
            if (!KernelWrapper::intersection(
                    facet1->plane(), facet2->plane())) {
                result = true;
                break;
            }
        }
        if (result) {
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
                if (!KernelWrapper::intersection(
                        facet1->plane(), facet2->plane(), facet3->plane())) {
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

Vector3SPtr PolyhedronTransformation::randVec(double min, double max) {
    double rval[3];
    for (unsigned int i = 0; i < 3; i++) {
        rval[i] = (max-min)*((double)rand()/(double)RAND_MAX) + min;
    }
    Vector3SPtr result = KernelFactory::createVector3(rval[0],rval[1],rval[2]);
    return result;
}

void PolyhedronTransformation::randMovePoints(PolyhedronSPtr polyhedron) {
    std::cout << "Points will be moved randomly..." << std::endl;

    double range = 0.001;
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        double value = config->getDouble("main", "rand_move_points_range");
        if (value != 0.0) {
            range = value;
        }
    }

    std::cout << "rand_move_points_range=" << range << std::endl;

    // Store coefficients for possible un-tilting post-offset
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        facet->storePlaneCoefficients();
    }

    // srand(time(NULL));
    srand(0); // set seed to a const value to reproduce errors

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        Vector3SPtr v_r = randVec(-range/2.0, range/2.0);
        Point3SPtr p_t = KernelFactory::createPoint3(*p + *v_r);
        vertex->setPoint(p_t);
    }
    polyhedron->initPlanes(); // recompute planes

    polyhedron->appendDescription("rand_move_points_range=" +
            util::StringFactory::fromDouble(range) + "; ");

    CGAL_postcondition(polyhedron && polyhedron->isConsistent());
}

PolyhedronSPtr PolyhedronTransformation::perturb(PolyhedronSPtr polyhedron) {
    db::_3d::OBJFile::save("results/perturbation_before.obj", polyhedron, false /*triangulate*/);

    // check if we can just tilt plane coefficients directly without having
    // to keep a triangulated input
    bool canUsePlaneTilts = true;

    // copy the polyhedron because we will merge (almost) coplanar faces and check if the result
    // is a mesh with only degree 3 vertices. Indeed, such a polyhedron can have its faces tilted
    // to create the perturbation.
    PolyhedronSPtr polyhedron_cpy = polyhedron->clone();

    // @todo
    // it would be nice to merge non connected faces too as to give them the same planes.
    // often we have the same "base" plane, but with vertical faces that are going to cut
    // it into separate connected components
    db::_3d::AbstractFile::mergeCoplanarFacets(polyhedron_cpy);

    db::_3d::OBJFile::save("results/pre-tilt_merged.obj", polyhedron_cpy,
                            false /*do_triangulate*/,
                            true /*convert_to_double*/);

    db::_3d::AbstractFile::removeVerticesDegLt3(polyhedron_cpy);
    CGAL_assertion(polyhedron_cpy && polyhedron_cpy->isConsistent());

    std::list<VertexSPtr>::iterator it_v = polyhedron_cpy->vertices().begin();
    while (it_v != polyhedron_cpy->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (vertex->degree() != 3) {
            std::cerr << "Can't use plane tilts because of " << vertex->toString() << std::endl;
            canUsePlaneTilts = false;
            break;
        }
    }

    if (!canUsePlaneTilts) {
        // this is not 'polyhedron_cpy' because the polyhedron must be triangulated
        // for vertices to remain on their incident planes
        randMovePoints(polyhedron);
    } else {
        std::cout << "All vertices are degree 3, tilt the polyhedron's faces" << std::endl;

        // @todo not necessary if we do not try to un-tilt at the end
        harmonizeFacetPlanes(polyhedron_cpy);

        std::list<FacetSPtr>::iterator it_f = polyhedron_cpy->facets().begin();
        while (it_f != polyhedron_cpy->facets().end()) {
            FacetSPtr facet = *it_f++;
            facet->storePlaneCoefficients();
            facet->perturbPlaneCoefficients();
        }

        // recompute vertex positions with the perturbed planes
        polyhedron_cpy = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron_cpy, 0.0);

        if (polyhedron_cpy && polyhedron_cpy->isConsistent()) {
            // @todo try again with another perturbation, or smarter (iterative) perturbation
            std::cerr << "Warning: perturbed polyhedron is no longer valid" << std::endl;
            polyhedron = polyhedron_cpy;
        }
    }

    db::_3d::OBJFile::save("results/perturbation_after.obj", polyhedron, false /*triangulate*/);

    return polyhedron;
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
    DEBUG_VAR(scale_factor);
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
                std::cout << *p << " is not in the box " << *p_box_min << " " << *p_box_max << std::endl;
                break;
            }
        }
        if (!result) {
            break;
        }
    }
    return result;
}

} }
