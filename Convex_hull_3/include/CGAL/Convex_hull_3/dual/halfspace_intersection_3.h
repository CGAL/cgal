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

#ifndef CGAL_HALFSPACE_INTERSECTION_3_H
#define CGAL_HALFSPACE_INTERSECTION_3_H

#include <CGAL/license/Convex_hull_3.h>


#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_3/dual/Convex_hull_traits_dual_3.h>
#include <CGAL/Origin.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/intersections.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Euler_operations.h>
// For interior_polyhedron_3
#include <CGAL/Convex_hull_3/dual/interior_polyhedron_3.h>
#include <CGAL/internal/Exact_type_selector.h>

#include <boost/unordered_map.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/foreach.hpp>

namespace CGAL
{
    namespace Convex_hull_3
    {
        namespace internal
        {
            // Build the primal polyhedron associated to a dual polyhedron
            // We also need the `origin` which represents a point inside the primal polyhedron
          template <class Polyhedron_dual, class Polyhedron, class Point_3>
            void
            build_primal_polyhedron(Polyhedron& primal,
                                    const Polyhedron_dual & _dual,
                                    const Point_3& origin)
            {
              typedef typename Kernel_traits<Point_3>::Kernel Kernel;
              typedef typename Kernel::RT RT;

              // Typedefs for dual
              typedef typename Polyhedron_dual::Facet Facet;
              typedef typename Polyhedron_dual::Facet_const_handle
                Facet_const_handle;
              typedef typename Polyhedron_dual::Facet_const_iterator
                Facet_const_iterator;
              typedef typename Polyhedron_dual::Vertex_const_iterator
                Vertex_const_iterator;

              // typedefs for primal
              typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

              // Typedefs for intersection
              typedef typename Kernel::Plane_3 Plane_3;
              typedef typename Kernel::Line_3 Line_3;
              typedef boost::optional< boost::variant< Point_3,
                                                       Line_3,
                                                       Plane_3 > > result_inter;

                           
              boost::unordered_map <Facet_const_handle, vertex_descriptor> primal_vertices;
              size_t n = 0;

              // First, computing the primal vertices
              for (Facet_const_iterator it = _dual.facets_begin();
                   it != _dual.facets_end(); ++it, ++n) {
                typename Facet::Halfedge_const_handle h = it->halfedge();
                // Build the dual plane corresponding to the current facet
                Plane_3 p1 = h->vertex()->point();
                Plane_3 p2 = h->next()->vertex()->point();
                Plane_3 p3 = h->next()->next()->vertex()->point();

                RT dp1 = p1.d() + origin.x() * p1.a()
                  + origin.y() * p1.b() + origin.z() * p1.c();
                RT dp2 = p2.d() + origin.x() * p2.a()
                  + origin.y() * p2.b() + origin.z() * p2.c();
                RT dp3 = p3.d() + origin.x() * p3.a()
                  + origin.y() * p3.b() + origin.z() * p3.c();

                Plane_3 pp1(p1.a(), p1.b(), p1.c(), dp1);
                Plane_3 pp2(p2.a(), p2.b(), p2.c(), dp2);
                Plane_3 pp3(p3.a(), p3.b(), p3.c(), dp3);

                // Compute the intersection
                result_inter result = CGAL::intersection(pp1, pp2, pp3);
                CGAL_assertion_msg(bool(result),
                                   "halfspace_intersection_3: no intersection");
                CGAL_assertion_msg(boost::get<Point_3>(& *result) != NULL,
                                   "halfspace_intersection_3: intersection is not a point");

                const Point_3* pp = boost::get<Point_3>(& *result);

                // Primal vertex associated to the current dual plane
                Point_3 ppp(origin.x() + pp->x(),
                            origin.y() + pp->y(),
                            origin.z() + pp->z());

                
                primal_vertices[it] = add_vertex(ppp, primal);
              }

              // Then, add facets to the primal polyhedron
              // To do this, for each dual vertex, we circulate around this vertex
              // and we add an edge between each facet we encounter

              for (Vertex_const_iterator it = _dual.vertices_begin();
                   it != _dual.vertices_end(); ++it) {
                std::vector<vertex_descriptor> vertices;
                typename Polyhedron_dual::Halfedge_around_vertex_const_circulator
                  h0 = it->vertex_begin(), hf = h0;
                  do {
                    vertices.push_back(primal_vertices[hf->facet()]);
                  } while (--hf != h0);
                  Euler::add_face(vertices,primal);
              }
            }

            // Test if a point is inside a convex polyhedron
          template <class Polyhedron, class Point>
            bool point_inside_convex_polyhedron (const Polyhedron &P,
                                                 Point const& p) {
            // Compute the equations of the facets of the polyhedron
            typedef typename boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
            typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
            
            typename boost::property_map<Polyhedron, vertex_point_t>::const_type vpmap  = get(CGAL::vertex_point, P);

            BOOST_FOREACH(face_descriptor fd, faces(P))
                {
                  halfedge_descriptor h = halfedge(fd,P), done(h);
                  Point const& p1 = get(vpmap, target(h,P));
                  h = next(h,P);
                  Point const& p2 = get(vpmap, target(h,P));
                  h = next(h,P);
                  Point p3 = get(vpmap, target(h,P));

                  while( h!=done && collinear(p1,p2,p3) )
                  {
                    h = next(h,P);
                    p3 = get(vpmap, target(h,P));
                  }
                  if (h==done) continue; //degenerate facet, skip it

                  if ( orientation (p1, p2, p3, p) != CGAL::NEGATIVE ) return false;
                }

                return true;
            }

            template <class Plane>
            bool collinear_plane (Plane u, Plane v) {
                typedef typename Kernel_traits<Plane>::Kernel Kernel;
                typedef typename Kernel::Vector_3 Vector;

                Vector uu = u.orthogonal_vector();
                Vector vv = v.orthogonal_vector();

                return CGAL::cross_product(uu, vv) == Vector(0, 0, 0);
            }

            template <class Plane>
            bool coplanar_plane (Plane u, Plane v, Plane w) {
                typedef typename Kernel_traits<Plane>::Kernel Kernel;
                typedef typename Kernel::Vector_3 Vector;

                Vector uu = u.orthogonal_vector();
                Vector vv = v.orthogonal_vector();
                Vector ww = w.orthogonal_vector();

                return CGAL::orientation(uu, vv, ww) == CGAL::COPLANAR;
            }

            // Checks if the dimension of intersection of the planes
            // is a polyhedron (dimension == 3)
            template <class InputPlaneIterator>
            bool is_intersection_dim_3 (InputPlaneIterator begin,
                                        InputPlaneIterator end) {
                typedef typename std::iterator_traits<InputPlaneIterator>::value_type Plane;
                typedef typename std::vector<Plane>::iterator PlaneIterator;

                std::vector<Plane> planes(begin, end);
                // Remove same planes
                std::size_t size = planes.size();

                // At least 4 points
                if (size < 4)
                    return false;

                // Look for two non-parallel planes
                PlaneIterator plane1_it = planes.begin();
                PlaneIterator plane2_it = cpp11::next(planes.begin());

                while (plane2_it != planes.end() &&
                       collinear_plane(*plane1_it, *plane2_it)) {
                    ++plane2_it;
                }

                if (plane2_it == planes.end()) return false;

                PlaneIterator plane3_it = cpp11::next(plane2_it);

                // Look for a triple of planes intersecting in a point
                while (plane3_it != planes.end() &&
                       coplanar_plane(*plane1_it, *plane2_it, *plane3_it)) {
                    ++plane3_it;
                }

                if (plane3_it == planes.end()) return false;

                return true;
            }
        } // namespace internal
    } // namespace Convex_hull_3

    // Compute the intersection of halfspaces.
    // If the user gives an origin then the function used it, otherwise, it is
    // computed using a linear program.
    template <class PlaneIterator, class Polyhedron>
    void halfspace_intersection_3 (PlaneIterator begin, PlaneIterator end,
                                   Polyhedron &P,
                                   boost::optional<typename Kernel_traits<typename std::iterator_traits<PlaneIterator>::value_type>::Kernel::Point_3> origin = boost::none) {
        // Checks whether the intersection is a polyhedron
        CGAL_assertion_msg(Convex_hull_3::internal::is_intersection_dim_3(begin, end), "halfspace_intersection_3: intersection not a polyhedron");

        // Types
        typedef typename Kernel_traits<typename std::iterator_traits<PlaneIterator>::value_type>::Kernel K;
        typedef Convex_hull_3::Convex_hull_traits_dual_3<K> Hull_traits_dual_3;
        typedef Polyhedron_3<Hull_traits_dual_3> Polyhedron_dual_3;

        // if a point inside is not provided find one using linear programming
        if (!origin) {
          // choose exact integral type
          typedef typename internal::Exact_field_selector<void*>::Type ET;

          // find a point inside the intersection
          typedef Interior_polyhedron_3<K, ET> Interior_polyhedron;
          Interior_polyhedron interior;
          CGAL_assertion_code(bool interior_point_found = )
          interior.find(begin, end);
          CGAL_assertion_msg(interior_point_found, "halfspace_intersection_3: problem when determing a point inside the intersection");
          origin = boost::make_optional(interior.inside_point());
        }

        // make sure the origin is on the negative side of all the planes
        CGAL_assertion_code(for(PlaneIterator pit=begin;pit!=end;++pit))
          CGAL_assertion(pit->has_on_negative_side(*origin));

        // compute the intersection of the half-space using the dual formulation
        Hull_traits_dual_3 dual_traits(*origin);
        Polyhedron_dual_3 dual_convex_hull;
        CGAL::convex_hull_3(begin, end, dual_convex_hull, dual_traits);
        Convex_hull_3::internal::build_primal_polyhedron(P, dual_convex_hull, *origin);

        // Posterior check if the origin is inside the computed polyhedron
        // The check is done only if the number type is not float or double because in that
        // case we know the construction of dual points is not exact
        CGAL_assertion_msg(
          boost::is_floating_point<typename K::FT>::value ||
          Convex_hull_3::internal::point_inside_convex_polyhedron(P, *origin),
          "halfspace_intersection_3: origin not in the polyhedron"
        );
    }

  #ifndef CGAL_NO_DEPRECATED_CODE
  // Compute the intersection of halfspaces (an interior point is given.)
  template <class PlaneIterator, class Polyhedron>
  void halfspace_intersection_3 (PlaneIterator begin, PlaneIterator end,
                                 Polyhedron &P,
                                 typename Kernel_traits<typename std::iterator_traits<PlaneIterator>::value_type>::Kernel::Point_3 const& origin) {
    halfspace_intersection_3(begin, end , P, boost::make_optional(origin));
  }
  #endif
} // namespace CGAL

#endif // CGAL_HALFSPACE_INTERSECTION_3_H

