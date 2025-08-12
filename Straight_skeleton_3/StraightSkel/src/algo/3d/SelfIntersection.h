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
 * @file   algo/3d/SelfIntersection.h
 * @author Gernot Walzl
 * @date   2012-07-18
 */

#ifndef ALGO_3D_SELFINTERSECTION_H
#define ALGO_3D_SELFINTERSECTION_H

#include "debug.h"
#include "data/3d/ptrs.h"

#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"

namespace algo { namespace _3d {

using namespace data::_3d;

class SelfIntersection {
public:
    virtual ~SelfIntersection();

    static bool doEdgesShareAVertex(EdgeSPtr edge1, EdgeSPtr edge2, bool handle_deg1_as_ray);

    // whether the edges intersect in their interior, with edges possibly being rays
    template <typename ProjectionTraits>
    static bool doEdgesIntersect(FacetSPtr,
                                 EdgeSPtr edge1,
                                 EdgeSPtr edge2,
                                 bool handle_degree_1_as_ray, // @todo always true?...
                                 const ProjectionTraits& traits)
    {
        CGAL_precondition(edge1 != edge2);
        CGAL_precondition(edge1->getVertexSrc() != edge1->getVertexDst() &&
                          edge2->getVertexSrc() != edge2->getVertexDst());

        // edges should never be a full line
        CGAL_precondition(edge1->getVertexSrc()->degree() != 1 || edge1->getVertexDst()->degree() != 1);
        CGAL_precondition(edge2->getVertexSrc()->degree() != 1 || edge2->getVertexDst()->degree() != 1);

        // reject any degeneracy as a self-intersection
        if (*(edge1->getVertexSrc()->getPoint()) == *(edge1->getVertexDst()->getPoint()) ||
            *(edge2->getVertexSrc()->getPoint()) == *(edge2->getVertexDst()->getPoint())) {
            return false;
        }

        // we could have something like an overlapping segment and ray, but in that case,
        // we will see the intersection from another pair
        if (doEdgesShareAVertex(edge1, edge2, handle_degree_1_as_ray)) {
            return false;
        }

        auto isRay = [](EdgeSPtr edge) -> bool
        {
            return (edge->getVertexSrc()->degree() == 1 || edge->getVertexDst()->degree() == 1);
        };

        auto treat_o_o = [&](EdgeSPtr, const auto& oa,
                             EdgeSPtr, const auto& ob) -> bool
        {
            return traits.do_intersect_2_object()(oa, ob);
        };

        auto treat_seg_seg = [&](EdgeSPtr a, EdgeSPtr b) -> bool
        {
            CGAL_precondition(!isRay(a) && !isRay(b));
            const Segment3 sa = { *(a->getVertexSrc()->getPoint()),
                                  *(a->getVertexDst()->getPoint()) };
            const Segment3 sb = { *(b->getVertexSrc()->getPoint()),
                                  *(b->getVertexDst()->getPoint()) };
            return treat_o_o(a, sa, b, sb);
        };

        auto treat_seg_ray = [&](EdgeSPtr a, EdgeSPtr b) -> bool
        {
            CGAL_precondition(!isRay(a) && isRay(b));
            const Segment3 s = { *(a->getVertexSrc()->getPoint()),
                                *(a->getVertexDst()->getPoint()) };
            const Ray3 r = (b->getVertexSrc()->degree() == 1) ? Ray3 { *b->getVertexDst()->getPoint(),
                                                                      *b->getVertexSrc()->getPoint() }
                                                              : Ray3 { *b->getVertexSrc()->getPoint(),
                                                                      *b->getVertexDst()->getPoint() };
            return treat_o_o(a, s, b, r);
        };

        auto treat_ray_ray = [&](EdgeSPtr a, EdgeSPtr b) -> bool
        {
            CGAL_precondition(isRay(a) && isRay(b));
            const Ray3 ra = (a->getVertexSrc()->degree() == 1) ? Ray3 { *a->getVertexDst()->getPoint(),
                                                                        *a->getVertexSrc()->getPoint() }
                                                              : Ray3 { *a->getVertexSrc()->getPoint(),
                                                                        *a->getVertexDst()->getPoint() };
            const Ray3 rb = (b->getVertexSrc()->degree() == 1) ? Ray3 { *b->getVertexDst()->getPoint(),
                                                                        *b->getVertexSrc()->getPoint() }
                                                              : Ray3 { *b->getVertexSrc()->getPoint(),
                                                                        *b->getVertexDst()->getPoint() };
            return treat_o_o(a, ra, b, rb);
        };

        if (handle_degree_1_as_ray) {
            if (isRay(edge1)) {
                if (isRay(edge2)) {
                    return treat_ray_ray(edge1, edge2);
                } else {
                    return treat_seg_ray(edge2, edge1);
                }
            } else if (isRay(edge2)) {
                return treat_seg_ray(edge1, edge2);
            }
        }

        return treat_seg_seg(edge1, edge2);
    }

    template <typename ProjectionTraits>
    static bool isSelfIntersectingFacet(FacetSPtr facet,
                                        const ProjectionTraits& traits)
    {
        bool result = false;

        // We can't use is_simple_polygon_2 + projection traits because some edges are rays
        std::list<EdgeSPtr>::iterator it_e1 = facet->edges().begin();
        while (it_e1 != facet->edges().end()) {
            EdgeSPtr edge1 = *it_e1++;
            std::list<EdgeSPtr>::iterator it_e2 = it_e1;
            while (it_e2 != facet->edges().end()) {
                EdgeSPtr edge2 = *it_e2++;
                if (doEdgesIntersect(facet, edge1, edge2, true, traits)) {
                    CGAL_SS3_ALGO_TRACE("edges intersect within the face:");
                    CGAL_SS3_ALGO_TRACE(facet->toString());
                    CGAL_SS3_ALGO_TRACE(edge1->toString());
                    CGAL_SS3_ALGO_TRACE(edge2->toString());
    #ifdef CGAL_SS3_EXIT_ASAP
                    return true;
    #else
                    result = true;
    #endif
                    break;
                }
            }
            if (result) {
                break;
            }
        }
        return result;
    }

    static bool isSelfIntersectingFacet(FacetSPtr facet);
    static unsigned int hasSelfIntersectingFacets(PolyhedronSPtr polyhedron);

    static bool isInsideWithRayShooting(const Point3& point,
                                        FacetSPtr facet,
                                        const bool handle_deg1_as_ray);
    static bool isInsideWithRayShootingV2(Point3SPtr point,
                                          FacetSPtr facet);
    static bool isEdgeInsideFacet(FacetSPtr facet, EdgeSPtr edge, bool handle_deg1_as_ray);
    static bool hasSelfIntersectingSurface(PolyhedronSPtr polyhedron);

protected:
    SelfIntersection();
};

} }

#endif /* ALGO_3D_SELFINTERSECTION_H */
