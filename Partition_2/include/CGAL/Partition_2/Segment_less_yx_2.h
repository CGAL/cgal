// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_SEGMENT_LESS_YX_2_H
#define CGAL_SEGMENT_LESS_YX_2_H

#include <utility>
#include <CGAL/Partition_2/Turn_reverser.h>

namespace CGAL {

//
// Compares two pairs of points representing two segments. The first is 
// "less than" the second if the second can see some point of the first 
// by looking straight down (i.e., in direction -pi/2). If the first can see
// the second in this way, it is greater than the second.  If neither sees
// the other when looking in the direction -pi/2, the one that is farther
// to the left is less than the one farther to the right.
//
template <class Traits>
class Segment_less_yx_2
{
   typedef typename Traits::Point_2             Point_2;
   typedef std::pair<Point_2, Point_2>          Point_pair;
   typedef typename Traits::Less_xy_2           Less_xy_2;
   typedef typename Traits::Compare_x_2         Compare_x_2;
   typedef typename Traits::Compare_y_2         Compare_y_2;
   typedef typename Traits::Left_turn_2          Left_turn_2;
   typedef Turn_reverser<Point_2, Left_turn_2>   Right_turn_2;

   public:
     Segment_less_yx_2() : 
       _less_xy_2(Traits().less_xy_2_object()),
       _compare_x_2(Traits().compare_x_2_object()),
       _compare_y_2(Traits().compare_y_2_object()),
       _left_turn_2(Traits().left_turn_2_object()),
       _right_turn_2(Right_turn_2(_left_turn_2))
     { }
     

     bool 
     operator()(const Point_pair& p, const Point_pair& q) const
     { 
        Point_2 p_smaller_xy, p_larger_xy;
        Point_2 q_smaller_xy, q_larger_xy;
        // order the point pairs by x value
        if (_less_xy_2(p.first, p.second))
        {
           p_smaller_xy = p.first;
           p_larger_xy = p.second;
        }
        else
        {
           p_smaller_xy = p.second;
           p_larger_xy = p.first;
        }
        if (_less_xy_2(q.first, q.second))
        {
           q_smaller_xy = q.first;
           q_larger_xy = q.second;
        }
        else
        {
           q_smaller_xy = q.second;
           q_larger_xy = q.first;
        }


        // x range of p comes before x range of q
        if (_compare_x_2(p_larger_xy, q_smaller_xy) == SMALLER)
           return true;
        else if (_compare_x_2(p_larger_xy, q_smaller_xy) == EQUAL)
        { 
           // x range of p ends where x range of q starts
           Comparison_result y_comp = _compare_y_2(p_larger_xy, q_smaller_xy);
           if (y_comp == SMALLER)
              return true;
           else if (y_comp == LARGER)
              return false;
           else // y_comp == EQUAL, so p's x range comes before q's
              return true;
        }
        // x range of q comes before x range of p
        else if (_compare_x_2(q_larger_xy, p_smaller_xy) == SMALLER)
           return false;
        else if (_compare_x_2(q_larger_xy, p_smaller_xy) == EQUAL)
        { 
           // x range of p starts where x range of q ends
           Comparison_result y_comp = _compare_y_2(p_smaller_xy, q_larger_xy);
           if (y_comp == SMALLER)
              return true;
           else if (y_comp == LARGER)
              return false;
           else // y_comp == EQUAL, so p's x range comes after q's
              return false;
        }
        // see if one of q's endpoints is contained in p's x range
        else if (_compare_x_2(p_smaller_xy,q_smaller_xy) == SMALLER && 
                 _compare_x_2(q_smaller_xy,p_larger_xy) == SMALLER)
           return _left_turn_2(p_smaller_xy,p_larger_xy,q_smaller_xy);
        else if (_compare_x_2(p_smaller_xy,q_larger_xy) == SMALLER &&
                 _compare_x_2(q_larger_xy,p_larger_xy) == SMALLER)
           return _left_turn_2(p_smaller_xy,p_larger_xy,q_larger_xy);
        //
        // neither of q's endpoints is in p's x-range so see if one of
        // p's endpoints is in q's x-range 
        //
        else if (_compare_x_2(q_smaller_xy,p_smaller_xy) == SMALLER && 
                 _compare_x_2(p_smaller_xy,q_larger_xy) == SMALLER)
           return _right_turn_2(q_smaller_xy,q_larger_xy,p_smaller_xy);
        else if (_compare_x_2(q_smaller_xy,p_larger_xy) == SMALLER &&
                 _compare_x_2(p_larger_xy,q_larger_xy) == SMALLER )
           return _right_turn_2(q_smaller_xy,q_larger_xy,p_larger_xy);
        else // the x ranges are exactly the same
        {
           Comparison_result y_comp = _compare_y_2(p_smaller_xy, q_smaller_xy);
           if (y_comp == SMALLER) 
              return true;
           else if (y_comp == LARGER)
              return false;
           else // look at the other endpoint
           {
              y_comp = _compare_y_2(p_larger_xy, q_larger_xy);
              if (y_comp == SMALLER) 
                 return true;
              else if (y_comp == LARGER)
                 return false;
              else  // point pairs are identical
                 return false;
           }
        }
     }

   private:
      Less_xy_2 _less_xy_2;
      Compare_x_2 _compare_x_2;
      Compare_y_2 _compare_y_2;
      Left_turn_2 _left_turn_2;
      Right_turn_2 _right_turn_2;
};

}

#endif // CGAL_SEGMENT_LESS_YX_2_H
