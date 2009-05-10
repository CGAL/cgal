// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://scm.gforge.inria.fr/svn/cgal/trunk/AABB_tree/include/CGAL/Point_3_triangle_3_point_3_nearest_point_3.h $
// $Id: Point_3_triangle_3_point_3_nearest_point_3.h 48981 2009-04-28 15:14:34Z stayeb $
//
//
// Author(s)     : Camille Wormser, Stephane Tayeb, Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef NEAREST_POINT_SEGMENT_3_H_
#define NEAREST_POINT_SEGMENT_3_H_

#include <CGAL/kernel_basic.h>
#include <CGAL/enum.h>


namespace CGAL {


    namespace CGALi {

        /**
        * @brief returns true if p is inside segment s. If p is not inside s,
        * result is the nearest point of s from p. WARNING: it is assumed that
        * t and p are on the same line.
        * @param query the query point
        * @param s the segment
        * @param closest_point_on_segment if query is not inside s, the nearest point of s from p
        * @param k the kernel
        * @return true if p is inside s
        */
        template <class K>
        inline
            bool
            is_inside_segment_3(const typename K::Point_3& query,
            const typename K::Segment_3 & s,
            typename K::Point_3& closest_point_on_segment,
            const K&)
        {
            typedef typename K::FT FT;
            typedef typename K::Point_3 Point;

            const Point& a = s.source();
            const Point& b = s.target();
            if((b-a)*(query-a) < (FT)0.0)
            {
                closest_point_on_segment = a;
                return false;
            }
            if((a-b)*(query-b) < (FT)0.0)
            {
                closest_point_on_segment = b;
                return false;
            }

            // query is on segment
            return true;
        }

        /**
        * @brief Computes the closest_point from query between bound and
        * any point of segment.
        * @param query the query point
        * @param segment the segment
        * @param bound the farthest point
        * @param k the kernel
        * @return nearest point: bound or a point inside segment
        */
        template <class K>
        typename K::Point_3
            nearest_point_3(const typename K::Point_3& query,
            const typename K::Segment_3& segment,
            const typename K::Point_3& bound,
            const K& k)
        {
            typedef typename K::Point_3 Point_3;
            typedef typename K::Vector_3 Vector_3;
            typedef typename K::FT FT;

            typename K::Compute_squared_distance_3 sq_distance =
                k.compute_squared_distance_3_object();
            typename K::Compare_squared_distance_3 compare_sq_distance =
                k.compare_squared_distance_3_object();
            typename K::Construct_projected_point_3 projection =
                k.construct_projected_point_3_object();

            // Square distance from query to bound
            const FT bound_sq_dist = sq_distance(query, bound);

            // Project query on segment supporting line
            const Point_3 proj = projection(segment.supporting_line(), query);

            // If point is projected outside, return bound
            if ( compare_sq_distance(query, proj, bound_sq_dist) == CGAL::LARGER )
                return bound;

            Point_3 closest_point_on_segment;
            bool inside = is_inside_segment_3(proj,segment,closest_point_on_segment,k);

            // If proj is inside segment, returns it
            if ( inside )
                return proj;

            // Else returns the constructed point (nearest segment' point from proj),
            // if it is closest to query than bound
            if ( compare_sq_distance(query, closest_point_on_segment, bound_sq_dist) == CGAL::LARGER )
                return bound;

            return closest_point_on_segment;
        }

    }  // end namespace CGALi


    template <class K>
    inline
        Point_3<K>
        nearest_point_3(const Point_3<K>& origin,
        const Segment_3<K>& segment,
        const Point_3<K>& bound)
    {
        return CGALi::nearest_point_3(origin, segment, bound, K());
    }

}  // end namespace CGAL


#endif // NEAREST_POINT_SEGMENT_3_H_
