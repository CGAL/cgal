// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Camille Wormser, Stephane Tayeb, Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_NEAREST_POINT_SEGMENT_3_H_
#define CGAL_NEAREST_POINT_SEGMENT_3_H_

#include <CGAL/kernel_basic.h>
#include <CGAL/enum.h>


namespace CGAL {


    namespace internal {

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
            if((b-a)*(query-a) < FT(0))
            {
                closest_point_on_segment = a;
                return false;
            }
            if((a-b)*(query-b) < FT(0))
            {
                closest_point_on_segment = b;
                return false;
            }

            // query is on segment
            return true;
        }

       

    }  // end namespace internal


}  // end namespace CGAL


#endif // CGAL_NEAREST_POINT_SEGMENT_3_H_

// This file uses an indentation width of 4, instead of 2.
// Sets that preference for GNU/Emacs, in a file-local variable.
// 
// Local Variables:
// c-basic-offset: 4
// End:
