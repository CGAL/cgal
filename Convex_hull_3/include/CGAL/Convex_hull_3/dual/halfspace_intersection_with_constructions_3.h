// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jocelyn Meyron
//

#ifndef CGAL_HALFSPACE_INTERSECTION_WITH_CONSTRUCTION_3_H
#define CGAL_HALFSPACE_INTERSECTION_WITH_CONSTRUCTION_3_H

#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Origin.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/assertions.h>

// For interior_polyhedron_3
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_interior_point_3.h>
#include <CGAL/internal/Exact_type_selector.h>

#include <boost/unordered_map.hpp>
#include <list>
#include <vector>

namespace CGAL
{
  namespace Convex_hull_3 {
    namespace internal
    {
      template <class Polyhedron, class Point_3>
      void
      build_dual_polyhedron(const Polyhedron & primal,
                            Polyhedron& dual,
                            Point_3 origin = Point_3(CGAL::ORIGIN))
      {
        typedef typename Kernel_traits<Point_3>::Kernel::Plane_3 Plane_3;

        typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
        typedef typename boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
        typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

        typename boost::property_map<Polyhedron, vertex_point_t>::const_type vpm_primal = get(CGAL::vertex_point, primal);
        typename boost::property_map<Polyhedron, vertex_point_t>::type vpm_dual = get(CGAL::vertex_point, dual);
        // compute coordinates of extreme vertices in the dual polyhedron
        // from primal faces
        boost::unordered_map<face_descriptor, vertex_descriptor> extreme_points;

        for(face_descriptor fd : faces( primal)){
          halfedge_descriptor h = halfedge(fd,primal);
          Plane_3 p (get(vpm_primal, target(h, primal)),
                     get(vpm_primal, target(next(h, primal), primal)),
                     get(vpm_primal, target(next(next(h, primal), primal), primal)));
          // translate extreme vertex
          Point_3 extreme_p = CGAL::ORIGIN + p.orthogonal_vector () / (-p.d());
          Point_3 translated_extreme_p(extreme_p.x() + origin.x(),
                                       extreme_p.y() + origin.y(),
                                       extreme_p.z() + origin.z());
          vertex_descriptor vd = add_vertex(dual);
          extreme_points[fd] = vd;
          put(vpm_dual, vd, translated_extreme_p);
        }

        // build faces
        for(vertex_descriptor vd : vertices(primal)) {
          //CGAL_assertion (it->is_bivalent() == false);

          std::list<vertex_descriptor> vertices;
          for(face_descriptor fd : faces_around_target(halfedge(vd,primal),primal)){
            vertices.push_front(extreme_points[fd]);
          }
        Euler::add_face(vertices,dual);
        }

      }
    } // namespace internal
  } // namespace Convex_hull_3

        // Compute the intersection of halfspaces by constructing explicitly
        // the dual points with the traits class for convex_hull_3 given
        // as an argument
        template <class PlaneIterator, class Polyhedron, class Traits>
          void halfspace_intersection_with_constructions_3(PlaneIterator pbegin,
                                                           PlaneIterator pend,
                                                           Polyhedron &P,
                                                           boost::optional<typename Kernel_traits<typename std::iterator_traits<PlaneIterator>::value_type>::Kernel::Point_3> origin,
                                                           const Traits & ch_traits) {
          typedef typename Kernel_traits<typename std::iterator_traits<PlaneIterator>::value_type>::Kernel K;
          typedef typename K::Point_3 Point;
          typedef typename K::Plane_3 Plane;

          // if a point inside is not provided find one using linear programming
          if (!origin) {
            // find a point inside the intersection
            origin = halfspace_intersection_interior_point_3(pbegin, pend);

            CGAL_assertion_msg(origin!=boost::none, "halfspace_intersection_with_constructions_3: problem when determing a point inside the intersection");
            if (origin==boost::none)
              return;
          }

          const Point p_origin = *origin;

          // construct dual points to apply the convex hull
          std::vector<Point> dual_points;
          for (PlaneIterator p = pbegin; p != pend; ++p) {
            // make sure the origin is on the negative side of all the planes
            CGAL_assertion(p->has_on_negative_side(p_origin));

            // translate plane
            Plane translated_p(p->a(),
                               p->b(),
                               p->c(),
                               p->d() + p_origin.x() * p->a() + p_origin.y() * p->b() + p_origin.z() * p->c());
            dual_points.push_back(CGAL::ORIGIN + translated_p.orthogonal_vector () / (-translated_p.d()));
          }

          Polyhedron ch;
          CGAL::convex_hull_3(dual_points.begin(), dual_points.end(), ch, ch_traits);

          Convex_hull_3::internal::build_dual_polyhedron (ch, P, p_origin);
        }

        // Compute the intersection of halfspaces by constructing explicitly
        // the dual points with the default traits class for convex_hull_3.
        template <class PlaneIterator, class Polyhedron>
          void halfspace_intersection_with_constructions_3 (PlaneIterator pbegin,
                                                            PlaneIterator pend,
                                                            Polyhedron &P,
                                                            boost::optional<typename Kernel_traits<typename std::iterator_traits<PlaneIterator>::value_type>::Kernel::Point_3> const& origin = boost::none) {
          typedef typename Kernel_traits<typename std::iterator_traits<PlaneIterator>::value_type>::Kernel K;
          typedef typename K::Point_3 Point_3;
          typedef typename Convex_hull_3::internal::Default_traits_for_Chull_3<Point_3>::type Traits;

          halfspace_intersection_with_constructions_3(pbegin, pend, P, origin, Traits());
        }


} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_HALFSPACE_INTERSECTION_WITH_CONSTRUCTION_3_H

