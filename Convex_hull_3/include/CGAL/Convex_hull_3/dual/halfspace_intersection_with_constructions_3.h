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

#include <CGAL/license/Convex_hull_3.h>


#include <CGAL/Origin.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/assertions.h>

// For interior_polyhedron_3
#include <CGAL/Convex_hull_3/dual/interior_polyhedron_3.h>
#include <CGAL/internal/Exact_type_selector.h>

#include <boost/foreach.hpp>
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
          
        typename boost::property_map<Polyhedron, vertex_point_t>::const_type vpmap  = get(CGAL::vertex_point, primal);
        // compute coordinates of extreme vertices in the dual polyhedron
        // from primal faces
        boost::unordered_map<face_descriptor, vertex_descriptor> extreme_points;

        BOOST_FOREACH (face_descriptor fd , faces( primal)){
          halfedge_descriptor h = halfedge(fd,primal);
          Plane_3 p (get(vpmap, target(h, primal)),
                     get(vpmap, target(next(h, primal), primal)),
                     get(vpmap, target(next(next(h, primal), primal), primal)));
          // translate extreme vertex
          Point_3 extreme_p = CGAL::ORIGIN + p.orthogonal_vector () / (-p.d());
          Point_3 translated_extreme_p(extreme_p.x() + origin.x(),
                                       extreme_p.y() + origin.y(),
                                       extreme_p.z() + origin.z());
          extreme_points[fd] = add_vertex(translated_extreme_p,dual);
        }
        
        // build faces
        BOOST_FOREACH (vertex_descriptor vd , vertices(primal)) {
          //CGAL_assertion (it->is_bivalent() == false);
          
          std::list<vertex_descriptor> vertices;
          BOOST_FOREACH(face_descriptor fd, faces_around_target(halfedge(vd,primal),primal)){
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
                                                           boost::optional<typename Kernel_traits<typename std::iterator_traits<PlaneIterator>::value_type>::Kernel::Point_3> const& origin,
                                                           const Traits & ch_traits) {
          typedef typename Kernel_traits<typename std::iterator_traits<PlaneIterator>::value_type>::Kernel K;
          typedef typename K::Point_3 Point;
          typedef typename K::Plane_3 Plane;

          Point p_origin;
            
          if (origin) {
            p_origin = boost::get(origin);
          } else {
            // choose exact integral type
            typedef typename internal::Exact_field_selector<void*>::Type ET;

            // find a point inside the intersection
            typedef Interior_polyhedron_3<K, ET> Interior_polyhedron;
            Interior_polyhedron interior;
            CGAL_assertion_code(bool res = )
              interior.find(pbegin, pend);
            CGAL_assertion_msg(res, "halfspace_intersection_with_constructions_3: problem when determing a point inside");
            p_origin = interior.inside_point();
          }

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
          typedef typename internal::Convex_hull_3::Default_traits_for_Chull_3<Point_3>::type Traits;

          halfspace_intersection_with_constructions_3(pbegin, pend, P, origin, Traits());
        }


} // namespace CGAL

#endif // CGAL_HALFSPACE_INTERSECTION_WITH_CONSTRUCTION_3_H

