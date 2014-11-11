// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Jocelyn Meyron
//

#ifndef CGAL_HALFSPACE_INTERSECTION_WITH_CONSTRUCTION_3_H
#define CGAL_HALFSPACE_INTERSECTION_WITH_CONSTRUCTION_3_H

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Origin.h>
#include <CGAL/convex_hull_3.h>

// For interior_polyhedron_3
#include <CGAL/Convex_hull_3/dual/interior_polyhedron_3.h>
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>
#else
#include <CGAL/MP_Float.h>
#endif

namespace CGAL
{
    namespace internal
    {
        template <class Polyhedron>
        class Build_dual_polyhedron :
                public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
        {
            typedef typename Polyhedron::HalfedgeDS HDS;
            typedef typename Polyhedron::Traits::Point_3 Point_3;
            const Polyhedron &_primal;
            Point_3 origin;

            public:
            Build_dual_polyhedron (const Polyhedron & primal,
                                   Point_3 o = Point_3(CGAL::ORIGIN)):
                _primal (primal), origin(o)
            {}

            void operator () (HDS &hds)
            {
                typedef typename Polyhedron::Facet Facet;
                typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
                typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;
                typedef typename Polyhedron::Vertex_const_iterator Vertex_const_iterator;
                typename CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);

                B.begin_surface(_primal.size_of_facets(),
                                _primal.size_of_vertices(),
                                _primal.size_of_halfedges());

                // compute coordinates of extreme vertices in the dual polyhedron
                // from primal faces
                std::map<Facet_const_handle, size_t> extreme_points;
                size_t n = 0;

                for (Facet_const_iterator it = _primal.facets_begin();
                     it != _primal.facets_end(); ++it, ++n)
                {
                    typename Facet::Halfedge_const_handle h = it->halfedge();
                    typename Facet::Plane_3 p ( h->vertex()->point(),
                                                h->next()->vertex()->point(),
                                                h->next()->next()->vertex()->point());
                    // translate extreme vertex
                    Point_3 extreme_p = CGAL::ORIGIN + p.orthogonal_vector () / (-p.d());
                    Point_3 translated_extreme_p(extreme_p.x() + origin.x(),
                                                 extreme_p.y() + origin.y(),
                                                 extreme_p.z() + origin.z());
                    B.add_vertex(translated_extreme_p);
                    extreme_points[it] = n;
                }

                // build faces
                for (Vertex_const_iterator it = _primal.vertices_begin();
                     it != _primal.vertices_end(); ++it)
                {
                    assert (it->is_bivalent() == false);

                    typename Polyhedron::Halfedge_around_vertex_const_circulator
                        h0 = it->vertex_begin(), hf = h0;

                    B.begin_facet();
                    do
                    {
                        B.add_vertex_to_facet(extreme_points[hf->facet()]);
                    } while (++hf != h0);
                    B.end_facet();
                }

                B.end_surface();
            }
        };
    } // namespace internal

    // Compute the intersection of halfspaces by constructing explicitly
    // the dual points with the traits class for convex_hull_3 given
    // as an argument
    template <class PlaneIterator, class Polyhedron, class Traits>
    void halfspace_intersection_with_constructions_3(PlaneIterator pbegin,
                                                     PlaneIterator pend,
                                                     Polyhedron &P,
                                                     boost::optional<typename Polyhedron::Vertex::Point_3> const& origin,
                                                     const Traits & ch_traits) {
            typedef typename Kernel_traits<typename Polyhedron::Vertex::Point_3>::Kernel K;
            typedef typename K::Point_3 Point;
            typedef typename K::Plane_3 Plane;
            typedef typename CGAL::internal::Build_dual_polyhedron<Polyhedron> Builder;

            Point p_origin;
            
            if (origin) {
              p_origin = boost::get(origin);
            } else {
              // choose exact integral type
#ifdef CGAL_USE_GMP
              typedef CGAL::Gmpq ET;
#else
              typedef CGAL::MP_Float ET;
#endif
              // find a point inside the intersection
              typedef Interior_polyhedron_3<K, ET> Interior_polyhedron;
              Interior_polyhedron interior;
              bool res = interior.find(pbegin, pend);
              CGAL_assertion_msg(res, "halfspace_intersection_with_constructions_3: problem when determing a point inside");
              p_origin = interior.inside_point();
            }

            // construct dual points to apply the convex hull
            std::vector<Point> dual_points;
            for (PlaneIterator p = pbegin; p != pend; ++p) {
                // translate plane
                Plane translated_p(p->a(),
                                   p->b(),
                                   p->c(),
                                   p->d() + p_origin.x() * p->a() + p_origin.y() * p->b() + p_origin.z() * p->c());
                dual_points.push_back(CGAL::ORIGIN + translated_p.orthogonal_vector () / (-translated_p.d()));
            }

            Polyhedron ch;
            CGAL::convex_hull_3(dual_points.begin(), dual_points.end(), ch, ch_traits);

            Builder build_dual (ch, p_origin);
            P.delegate(build_dual);
        }

    // Compute the intersection of halfspaces by constructing explicitly
    // the dual points with the traits class for convex_hull_3 given
    // as an argument.
    // An interior point is given.
    template <class PlaneIterator, class Polyhedron, class Traits>
    void halfspace_intersection_with_constructions_3(PlaneIterator pbegin,
                                                     PlaneIterator pend,
                                                     Polyhedron &P,
                                                     typename Polyhedron::Vertex::Point_3 const& origin,
                                                     const Traits & ch_traits) {
      halfspace_intersection_with_constructions_3(pbegin, pend, P,
                                                  boost::optional<typename Polyhedron::Vertex::Point_3>(origin),
                                                  ch_traits);
    }

    // Compute the intersection of halfspaces by constructing explicitly
    // the dual points with the default traits class for convex_hull_3.
    template <class PlaneIterator, class Polyhedron>
    void halfspace_intersection_with_constructions_3 (PlaneIterator pbegin,
                                                      PlaneIterator pend,
                                                      Polyhedron &P,
                                                      boost::optional<typename Polyhedron::Vertex::Point_3> const& origin) {
        typedef typename Kernel_traits<typename Polyhedron::Vertex::Point_3>::Kernel K;
        typedef typename K::Point_3 Point_3;
        typedef typename internal::Convex_hull_3::Default_traits_for_Chull_3<Point_3>::type Traits;

        halfspace_intersection_with_constructions_3(pbegin, pend, P, origin, Traits());
    }

    // Compute the intersection of halfspaces by constructing explicitly
    // the dual points with the default traits class for convex_hull_3.
    // An interior point is given.
    template <class PlaneIterator, class Polyhedron>
    void halfspace_intersection_with_constructions_3 (PlaneIterator pbegin,
                                                      PlaneIterator pend,
                                                      Polyhedron &P,
                                                      typename Polyhedron::Vertex::Point_3 const& origin) {
      halfspace_intersection_with_constructions_3(pbegin, pend, P,
                                                  boost::optional<typename Polyhedron::Vertex::Point_3>(origin));
    }

    // Compute the intersection of halfspaces by constructing explicitly
    // the dual points with the default traits class for convex_hull_3.
    // An interior point is not given.
    template <class PlaneIterator, class Polyhedron>
    void halfspace_intersection_with_constructions_3 (PlaneIterator pbegin,
                                                      PlaneIterator pend,
                                                      Polyhedron &P) {
      halfspace_intersection_with_constructions_3(pbegin, pend, P, boost::none);
    }
} // namespace CGAL

#endif // CGAL_HALFSPACE_INTERSECTION_WITH_CONSTRUCTION_3_H

