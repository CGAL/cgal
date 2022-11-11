// Copyright (c) 2015  Università della Svizzera italiana.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SIDE_OF_ORIENTED_SQUARE_2_H
#define CGAL_SIDE_OF_ORIENTED_SQUARE_2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>


#include <CGAL/basic.h>
#include <CGAL/Orientation_Linf_2.h>
#include <CGAL/Side_of_bounded_square_2.h>
#include <CGAL/enum.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/basic.h>

namespace CGAL {

    template<class K>
    class Side_of_oriented_square_2
    {
    private:
      typedef typename K::Point_2               Point_2;
      typedef Orientation_Linf_2<K>             Orientation_Linf_2_Type;
      typedef Side_of_bounded_square_2<K>       Side_of_bounded_square_2_Type;
      typedef typename K::Compare_x_2           Compare_x_2;
      typedef typename K::Compare_y_2           Compare_y_2;
      typedef typename K::Comparison_result     Comparison_result;

      Compare_x_2 compare_x_2;
      Compare_y_2 compare_y_2;

      Side_of_bounded_square_2_Type side_of_bounded_square_2;
      Orientation_Linf_2_Type orientation_Linf;

      Oriented_side predicate(const Point_2 &p, const Point_2 &q,
                 const Point_2 &r, const Point_2 &t) const
      {
        CGAL_SDG_DEBUG(std::cout << "debug entering side_of_os (pqrt)= ("
          << p << ") (" << q << ") (" << r << ") (" << t << ")"
          << std::endl;);

        Oriented_side orlpqr = orientation_Linf(p, q, r);

        if (orlpqr == DEGENERATE) {
          // here p,q,r are monotone
          CGAL_SDG_DEBUG(std::cout << "debug side_of_os pqr are monotone" << std::endl;);

          bool is_degenerate_pqt = (orientation_Linf(p,q,t) == DEGENERATE);
          bool is_degenerate_qrt = (orientation_Linf(q,r,t) == DEGENERATE);

          if (is_degenerate_pqt && is_degenerate_qrt) {
            //p,q,r,t are all collinear
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs pqrt all collin" << std::endl;);
            return ON_ORIENTED_BOUNDARY;
          }

          if (! is_degenerate_pqt) {
            //CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs qpt not monotone" << std::endl;);
            return predicate(q,p,t,r);
          }

          if (! is_degenerate_qrt) {
            //CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs rqt not monotone" << std::endl;);
            return predicate(r,q,t,p);
          }

        }
        else { // p,q,r are not monotone
          //CGAL_SDG_DEBUG(std::cout << "side_of_os pqr not monotone" << std::endl;);

          Comparison_result cxtp = compare_x_2(t,p);
          Comparison_result cytp = compare_y_2(t,p);
          Comparison_result cxtq = compare_x_2(t,q);
          Comparison_result cytq = compare_y_2(t,q);
          Comparison_result cxtr = compare_x_2(t,r);
          Comparison_result cytr = compare_y_2(t,r);

          // if t equals any one of p,q,r return zero
          if (((cxtp == EQUAL) && (cytp == EQUAL)) ||
              ((cxtq == EQUAL) && (cytq == EQUAL)) ||
              ((cxtr == EQUAL) && (cytr == EQUAL))   )
          {
            return ON_ORIENTED_BOUNDARY;
          }

          // check if query point is inside any of the following
          // segments: pq, qr, rp
          if (((cxtp == EQUAL) && (cxtq == EQUAL) && (cytp != cytq))
           || ((cxtq == EQUAL) && (cxtr == EQUAL) && (cytq != cytr))
           || ((cxtr == EQUAL) && (cxtp == EQUAL) && (cytr != cytp))
           || ((cytp == EQUAL) && (cytq == EQUAL) && (cxtp != cxtq))
           || ((cytq == EQUAL) && (cytr == EQUAL) && (cxtq != cxtr))
           || ((cytr == EQUAL) && (cytp == EQUAL) && (cxtr != cxtp)))
          {
            CGAL_SDG_DEBUG(std::cout << "debug side_of_os query point in segment"
              << std::endl;);
            return (Oriented_side)
              (( (int) orlpqr ) *
               ( (int) ON_POSITIVE_SIDE ) ) ;
          }

          Bounded_side bspqrt = side_of_bounded_square_2(p,q,r,t);

          CGAL_SDG_DEBUG(std::cout << "debug side_of_os bspqrt="
            << bspqrt << std::endl;);

          if (bspqrt == ON_BOUNDARY) {

            if ((cxtp != EQUAL) && (cytp != EQUAL)
                && (! (orientation_Linf(r,q,t)==DEGENERATE))) {
               // r q t p
               return (Oriented_side)
                 (((int) orientation_Linf(r,q,t)) *
                  ((int) side_of_bounded_square_2(r,q,t,p)) ) ;
            }
            if ((cxtq != EQUAL) && (cytq != EQUAL)
                && (! (orientation_Linf(p,r,t)==DEGENERATE))) {
               // p r t q
               return (Oriented_side)
                 (((int) orientation_Linf(p,r,t)) *
                  ((int) side_of_bounded_square_2(p,r,t,q)) ) ;
            }
            if ((cxtr != EQUAL) && (cytr != EQUAL)
                && (! (orientation_Linf(q,p,t)==DEGENERATE))) {
               // q p t r
               return (Oriented_side)
                 (((int) orientation_Linf(q,p,t)) *
                  ((int) side_of_bounded_square_2(q,p,t,r)) ) ;
            }
            CGAL_SDG_DEBUG(std::cout << "debug side_of_os about to return "
              << "ON_ORIENTED_BOUNDARY" << std::endl;);
            return ON_ORIENTED_BOUNDARY;
          }
          else if ( bspqrt == ON_BOUNDED_SIDE ) {
            return (orlpqr == LEFT_TURN) ?
                    ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE  ;
          }
          else {
            // ( Side_of_bounded_square_2(p,q,r,t) == ON_UNBOUNDED_SIDE )
            CGAL_assertion(bspqrt == ON_UNBOUNDED_SIDE);
            return (orlpqr == LEFT_TURN) ?
                    ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE ;
          }
        }
        CGAL_SDG_DEBUG(std::cout << "should not reach here" << std::endl;);

        CGAL_assertion(false);

        // avoid complaining from certain compilers
        return ON_ORIENTED_BOUNDARY;

      } // end of def of Oriented_side predicate(p, q, r)

    public:

      Oriented_side operator()(const Point_2 &p, const Point_2 &q,
             const Point_2 &r, const Point_2 &t) const
      {
        return predicate(p, q, r, t);
      }
    };


} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_SIDE_OF_ORIENTED_SQUARE_C2_H
