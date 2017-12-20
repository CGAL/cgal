// Copyright (c) 2015  Università della Svizzera italiana.
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_ORIENTATION_LINF_2_H
#define CGAL_ORIENTATION_LINF_2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>


#include <CGAL/basic.h>
#include <CGAL/enum.h>

namespace CGAL {

typedef Sign OrientationLinf;


template<class K>
class Orientation_Linf_2
{
private:
  typedef typename K::Point_2     Point_2;
  typedef typename K::Compare_x_2 Compare_x_2;
  typedef typename K::Compare_y_2 Compare_y_2;
  typedef typename K::Comparison_result Comparison_result;

  Compare_x_2 compare_x_2;
  Compare_y_2 compare_y_2;


  OrientationLinf predicate(const Point_2& p, const Point_2& q,
	                    const Point_2& r) const
  {
    Comparison_result cmpxpq = compare_x_2(p, q);
    Comparison_result cmpypq = compare_y_2(p, q);
    Comparison_result cmpxpr = compare_x_2(p, r);
    Comparison_result cmpypr = compare_y_2(p, r);
    Comparison_result cmpxqr = compare_x_2(q, r);
    Comparison_result cmpyqr = compare_y_2(q, r);

    //std::cout << "debug Orientation_Linf_2 (p,q,r)= "
    //          << p << ' ' << q << ' ' << r 
    //          << std::endl;

    if ( cmpxpq == EQUAL ) {
      if (cmpypq == EQUAL) {//p and q are same points
        return DEGENERATE;
      } 
      else {//pq forms a vertical line
        if (cmpxpr == EQUAL) {//r lies on the vertical line formed by pq
          return DEGENERATE;
        } 
        else {
          return (cmpypq == cmpxpr) ? RIGHT_TURN : LEFT_TURN ;
        }
      }
    } 
    else if ( cmpypq == EQUAL ) {//p and q forms a horizontal line 
      // here: cmpxpq != EQUAL
      if (cmpypr == EQUAL) {//r lies on the horizontal line formed by pq
        return DEGENERATE;
      } 
      else {
        return (cmpxpq == cmpypr) ? LEFT_TURN : RIGHT_TURN ;
      }
    }
    else {
      // here both cmpxpq != EQUAL and cmpypq != EQUAL
      bool is_monotone = 
          ( ( ( cmpxpr == -cmpxqr ) &&
              ( cmpypr == -cmpyqr )    ) ||
            ( ( cmpxpq == cmpxpr) && ( cmpxpr == cmpxqr ) &&
              ( cmpypq == cmpypr) && ( cmpypr == cmpyqr )    ) ||
            ( (-cmpxpq == cmpxpr) && ( cmpxpr == cmpxqr) &&
              (-cmpypq == cmpypr) && ( cmpypr == cmpyqr)     )    ) ;

      //std::cout << "debug is_monotone=" << is_monotone << std::endl; 

      if (is_monotone) {
        //p, q, r are monotone here
        //std::cout << "debug Orientation_Linf_2 are monotone"
        //      << std::endl;
        return DEGENERATE;
      }
      else { // p, q, r do not form stair case
        //std::cout << "debug Orientation_Linf_2 not monotone"
        //      << std::endl;
        if ( cmpxpq == SMALLER ) {
          if ( cmpypq == SMALLER ) {//CASE-I q lies in the North East of p
            return (cmpyqr == SMALLER ||
                cmpxpr == LARGER  ||
                (cmpxqr == LARGER && cmpypr == SMALLER)) ? 
              LEFT_TURN : RIGHT_TURN;
          }//end of CASE-I
          else {
            // compare_y_2(p, q) == LARGER -- 
            // CASE-II q lies in the South East of p
            return (cmpxqr == SMALLER ||
                cmpypr == SMALLER ||
                (cmpxpr == SMALLER && cmpyqr == SMALLER)) ? 
              LEFT_TURN : RIGHT_TURN;
          }//end of CASE-II
        }
        else {//cmpxpq == LARGER case III and case IV
          //std::cout << "debug Orientation_Linf_2 cmpxpq LARGER"
          //    << std::endl;
          if ( cmpypq == SMALLER ) {//CASE-III q lies in the North West of p
            //std::cout << "debug Orientation_Linf_2 cmpypq SMALLER"
            //    << std::endl;
            return (cmpxpr == SMALLER ||
                cmpyqr == SMALLER ||
                (cmpxqr == SMALLER && cmpypr == SMALLER)) ? 
              RIGHT_TURN : LEFT_TURN;			
          }//end of CASE-III
          else {//compare_y_2(p, q) == LARGER -- CASE-IV q lies in the South West of p
            //std::cout << "debug Orientation_Linf_2 cmpypq LARGER"
            //    << std::endl;
            return (cmpypr == SMALLER ||
                cmpxqr == LARGER ||
                (cmpxpr == LARGER && cmpyqr == SMALLER)) ? 
              RIGHT_TURN : LEFT_TURN;
          }//end of CASE-IV
        }
      }
    }
  }

public:
  typedef OrientationLinf   result_type;
  typedef Point_2           argument_type;

  OrientationLinf operator()(const Point_2& p, const Point_2& q,
	                     const Point_2& r) const
  {
    return predicate(p, q, r);
  }
};

} //namespace CGAL

#endif // CGAL_ORIENTATION_LINF_2_H
