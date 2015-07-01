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
// Author(s) : Jocelyn Meyron and Quentin Mérigot
//

#ifndef CGAL_INTERNAL_VCM_VORONOI_COVARIANCE_3_HPP
#define CGAL_INTERNAL_VCM_VORONOI_COVARIANCE_3_HPP

#include <list>
#include <CGAL/array.h>
#include <CGAL/internal/Voronoi_covariance_3/voronoi_covariance_sphere_3.h>
#ifdef CGAL_VORONOI_COVARIANCE_USE_CONSTRUCTIONS
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#else
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#endif

namespace CGAL {
    namespace Voronoi_covariance_3 {
        namespace internal {
            template <class FT>
                inline void
                covariance_matrix_tetrahedron (FT ax, FT ay, FT az,
                                               FT bx, FT by, FT bz,
                                               FT cx, FT cy, FT cz,
                                               cpp11::array<FT,6> &R)
                {
                    const FT det = (ax*cz*by - ax*bz*cy - ay*bx*cz +
                                    ay*cx*bz + az*bx*cy - az*cx*by) / 60.0;

                    R[0] += (ax*ax + ax*bx + ax*cx +
                             bx*bx + bx*cx + cx*cx) * det;
                    R[1] += (ax*ay + ax*by/2.0 + ax*cy/2.0 +
                             bx*ay/2.0 + bx*by + bx*cy/2.0 +
                             cx*ay/2.0 + cx*by/2.0 + cx*cy) * det;
                    R[2] += (ax*az + ax*bz/2.0 + ax*cz/2.0 +
                             bx*az/2.0 + bx*bz + bx*cz/2.0 +
                             cx*az/2.0 + cx*bz/2.0 + cx*cz) * det;

                    R[3] += (ay*ay + ay*by + ay*cy +
                             by*by + by*cy + cy*cy) * det;
                    R[4] += (az*ay + az*by/2.0 + az*cy/2.0 +
                             bz*ay/2.0 + bz*by + bz*cy/2.0 +
                             cz*ay/2.0 + cz*by/2.0 + cz*cy) * det;

                    R[5] += (az*az + az*bz + az*cz +
                             bz*bz + bz*cz + cz*cz) * det;
                }

            template <class FT>
                class Covariance_accumulator_3
                {
                    public:
                        typedef cpp11::array<FT, 6> Result_type;

                    private:
                        Result_type _result;

                    public:
                        Covariance_accumulator_3()
                        {
                            std::fill (_result.begin(), _result.end(), FT(0));
                        }

                        template <class Point>
                            inline void operator () (const Point &a,
                                                     const Point &c,
                                                     const Point &b)
                            {
                                internal::covariance_matrix_tetrahedron (a[0], a[1], a[2],
                                                                         b[0], b[1], b[2],
                                                                         c[0], c[1], c[2],
                                                                         _result);
                            }

                        const Result_type &result() const
                        {
                            return _result;
                        }
                };

            template <class FT>
                class Volume_accumulator_3
                {
                    public:
                        typedef FT Result_type;

                    private:
                        Result_type _result;

                    public:
                        Volume_accumulator_3() : _result(0.0)
                    {}

                        template <class Point>
                            inline void operator () (const Point &a,
                                                     const Point &c,
                                                     const Point &b)
                            {
                                const double  vol = CGAL::volume(a, b, c, Point(CGAL::ORIGIN));
                                //std::cerr << "vol = " << vol << "\n";
                                _result += vol;
                            }

                        const Result_type &result() const
                        {
                            return _result;
                        }
                };


            template <class DT, class Sphere, class F>
                F& tessellate_and_intersect(const DT &dt,
                                            typename DT::Vertex_handle v,
                                            const Sphere &sphere,
                                            F &f)
                {
                    typedef typename DT::Vertex_handle Vertex_handle;
                    typedef typename DT::Geom_traits::Kernel K;
                    typedef typename K::Plane_3 Plane;
                    typedef typename K::Point_3 Point;
                    typedef typename K::Vector_3 Vector;
                    typedef typename CGAL::Convex_hull_traits_3<K> Traits;
                    typedef typename Traits::Polyhedron_3 Polyhedron;

                    std::list<Vertex_handle> vertices;
                    dt.incident_vertices(v,std::back_inserter(vertices));

                    // construct intersection of half-planes using the convex hull function
                    std::list<Plane> planes;
                    for(typename std::list<Vertex_handle>::iterator it = vertices.begin();
                        it != vertices.end(); ++it)
                    {
                        Vector p = ((*it)->point() - v->point())/2;
                        planes.push_back (Plane(CGAL::ORIGIN+p, p));
                    }

                    // add half-planes defining the sphere discretization
                    sphere(std::back_inserter(planes));

                    Polyhedron P;
                    #ifdef CGAL_VORONOI_COVARIANCE_USE_CONSTRUCTIONS
                    halfspace_intersection_with_constructions_3
                    #else
                    halfspace_intersection_3
                    #endif
                      (planes.begin(), planes.end(), P, Point(CGAL::ORIGIN));

                    // apply f to the triangles on the boundary of P
                    for (typename Polyhedron::Facet_iterator it = P.facets_begin();
                         it != P.facets_end(); ++it)
                    {
                        typename Polyhedron::Halfedge_around_facet_circulator
                            h0 = it->facet_begin(), hf = h0--, hs = cpp11::next(hf);

                        while(1)
                        {
                            f (h0->vertex()->point(), hf->vertex()->point(),
                               hs->vertex()->point());
                            if (hs == h0)
                                break;
                            ++hs; ++hf;
                        }
                    }
                    return f;
                }
        } // namespace internal

        template <class DT, class Sphere, class FT>
            void
            voronoi_covariance_3 (const DT &dt,
                                  typename DT::Vertex_handle v,
                                  const Sphere &sphere,
                                  FT covariance[6])
            {
                typename internal::Covariance_accumulator_3<FT> ca;
                internal::tessellate_and_intersect(dt, v, sphere, ca);
                std::copy (ca.result().begin(), ca.result().end(), covariance);
            }

        template <class DT, class Sphere>
            array<typename DT::Geom_traits::FT, 6>
            voronoi_covariance_3 (const DT &dt,
                                  typename DT::Vertex_handle v,
                                  const Sphere &sphere)
            {
                typedef typename DT::Geom_traits::FT FT;
                typename internal::Covariance_accumulator_3<FT> ca;

                return internal::tessellate_and_intersect(dt, v, sphere, ca).result();
            }

        template <class DT, class Sphere>
            typename DT::Geom_traits::FT
            voronoi_volume_3 (const DT &dt,
                              typename DT::Vertex_handle v,
                              const Sphere &sphere)
            {
                typedef typename DT::Geom_traits::FT FT;
                typename internal::Volume_accumulator_3<FT> va;

                return internal::tessellate_and_intersect(dt, v, sphere, va).result();
            }

    } // namespace Voronoi_covariance_3
} // namespace CGAL

#endif // CGAL_INTERNAL_VCM_VORONOI_COVARIANCE_3_HPP

