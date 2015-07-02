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
// 
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

//ON_UNBOUNDED_SIDE = -1,
//ON_BOUNDARY = 0,
//ON_BOUNDED_SIDE = 1

#ifndef CGAL_SIDE_OF_BOUNDED_SQUARE_2_H
#define CGAL_SIDE_OF_BOUNDED_SQUARE_2_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/Orientation_Linf_2.h>
#include <CGAL/tuple.h>

namespace CGAL {

    template<class K>
    class Side_of_bounded_square_2
    {
    private:
      typedef typename K::Point_2               Point_2;
      typedef Orientation_Linf_2<K>             Orientation_Linf_2_Type;
      typedef typename K::Compare_x_2           Compare_x_2;
      typedef typename K::Compare_y_2           Compare_y_2;
      typedef typename K::FT                    FT;
      typedef typename K::Comparison_result     Comparison_result;

      typedef CGAL::cpp11::tuple<const bool, const bool,
                 const bool, const bool, const size_t> SmallerEqTuple;

      typedef CGAL::cpp11::tuple<FT const *, FT const *,
                   const bool, const bool, const bool,
                   Point_2 const *, Point_2 const *> MinMaxTuple;

      Compare_x_2 compare_x_2;
      Compare_y_2 compare_y_2;

      Orientation_Linf_2_Type orientation_Linf;

      inline MinMaxTuple
      minmax(const FT &p, const FT &q, const FT &r,
             const Point_2 &pt_p, const Point_2 &pt_q, const Point_2 &pt_r)
      const
      {
        FT const * min_p (NULL);
        FT const * max_p (NULL);
        bool samepq = false;
        bool samepr = false;
        bool sameqr = false;
        const Comparison_result cmppq = compare(p, q);
        switch(cmppq) {
          case SMALLER:
            min_p = &p;
            max_p = &q;
            break;
          case LARGER:
            min_p = &q;
            max_p = &p;
            break;
          default: // EQUAL
            min_p = &p;
            max_p = &q;
            samepq = true;
            break;
        }
        const Comparison_result cmppr = compare(p, r);
        Comparison_result cmpqr;
        if (samepq) {
          cmpqr = cmppr;
          switch(cmppr) {
            case SMALLER:
              max_p = &r;
              break;
            case LARGER:
              min_p = &r;
              break;
            default: // EQUAL is impossible
              CGAL_assertion(false);
              break;
          }
        } else {
          if (min_p == &p) {
            switch(cmppr) {
              case SMALLER:
                cmpqr = compare(q, r);
                switch(cmpqr) {
                  case SMALLER:
                    max_p = &r;
                    break;
                  case LARGER:
                    break;
                  default:
                    // q and r have the same coord
                    sameqr = true;
                    break;
                }
                break;
              case LARGER:
                cmpqr = cmppr;
                min_p = &r;
                break;
              default:
                // p and r have the same coord
                cmpqr = - cmppr;
                samepr = true;
                break;
            }
          } else { // min_p == &q
            switch(cmppr) {
              case SMALLER:
                cmpqr = cmppr;
                max_p = &r;
                break;
              case LARGER:
                cmpqr = compare(q, r);
                switch(cmpqr) {
                  case SMALLER:
                    break;
                  case LARGER:
                    min_p = &r;
                    break;
                  default:
                    // q and r have the same coord
                    sameqr = true;
                    break;
                }
                break;
              default:
                // p and r have the same coord
                cmpqr = - cmppq;
                samepr = true;
                break;
            }
          }
        }
        CGAL_assertion(min_p != NULL);
        CGAL_assertion(max_p != NULL);
        Point_2 const * pt_min_p =
          (min_p == &p)? &pt_p : (min_p == &q)? &pt_q : &pt_r;
        Point_2 const * pt_max_p =
          (max_p == &p)? &pt_p : (max_p == &q)? &pt_q : &pt_r;
        CGAL_SDG_DEBUG(std::cout << "debug minmax cmppq=" << cmppq
            << " cmppr=" << cmppr << " cmpqr=" << cmpqr << std::endl; );
        return CGAL::cpp11::make_tuple(
            min_p, max_p, samepq, samepr, sameqr, pt_min_p, pt_max_p);
      }

      inline SmallerEqTuple analyze_smalleq(
          const Comparison_result & cxtmax,
          const Comparison_result & cytmax,
          const Comparison_result & cxmint,
          const Comparison_result & cymint) const
      {
        CGAL_precondition(cxtmax != LARGER);
        CGAL_precondition(cytmax != LARGER);
        CGAL_precondition(cxmint != LARGER);
        CGAL_precondition(cymint != LARGER);
        size_t count_eq = 0;
        const bool at_rgt (cxtmax == EQUAL);
        bool at_lft (false);
        if (at_rgt) {
          ++count_eq;
          CGAL_assertion( cxmint == SMALLER );
        } else {
          at_lft = cxmint == EQUAL;
          if (at_lft) { ++count_eq; }
        }
        const bool at_top (cytmax == EQUAL);
        bool at_bot (false);
        if (at_top) {
          ++count_eq;
          CGAL_assertion( cymint == SMALLER );
        } else {
          at_bot = cymint == EQUAL;
          if (at_bot) { ++count_eq; }
        }
        return CGAL::cpp11::make_tuple(
            at_rgt, at_top, at_lft, at_bot, count_eq);
      }

      inline Bounded_side predicate(const Point_2 &p, const Point_2 &q,
                  const Point_2 &r, const Point_2 &t) const
      {
        CGAL_SDG_DEBUG(std::cout
            << "debug Side_of_bounded_square_2 entering" << std::endl;);
        CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs (pqrt)= (" << p
            << ") (" << q << ") (" << r << ") (" << t << ")" << std::endl;);

        CGAL_assertion(orientation_Linf(p,q,r) != DEGENERATE);

        //compute the minimum x and maximum x
        const FT px (p.x()), qx (q.x()), rx (r.x());
        MinMaxTuple tupx = minmax(px, qx, rx, p, q, r);
        FT const * lft_p = CGAL::cpp11::get<0>(tupx);
        FT const * rgt_p = CGAL::cpp11::get<1>(tupx);
        const bool samex_pq = CGAL::cpp11::get<2>(tupx);
        const bool samex_pr = CGAL::cpp11::get<3>(tupx);
        const bool samex_qr = CGAL::cpp11::get<4>(tupx);
        CGAL_assertion(lft_p != NULL);
        CGAL_assertion(rgt_p != NULL);

        //compute the minimum y and maximum y
        const FT py (p.y()), qy (q.y()), ry (r.y());
        MinMaxTuple tupy = minmax(py, qy, ry, p, q, r);
        FT const * bot_p = CGAL::cpp11::get<0>(tupy);
        FT const * top_p = CGAL::cpp11::get<1>(tupy);
        const bool samey_pq = CGAL::cpp11::get<2>(tupy);
        const bool samey_pr = CGAL::cpp11::get<3>(tupy);
        const bool samey_qr = CGAL::cpp11::get<4>(tupy);
        CGAL_assertion(bot_p != NULL);
        CGAL_assertion(top_p != NULL);

        CGAL_SDG_DEBUG(std::cout << "debug bs before mirror"
            << "  lft=" << *lft_p << "  rgt=" << *rgt_p
            << "  bot=" << *bot_p << "  top=" << *top_p
            << std::endl; );

        const bool exist_two_with_same_x = samex_pq || samex_pr || samex_qr;
        const bool exist_two_with_same_y = samey_pq || samey_pr || samey_qr;

        bool is_lft_input (true);
        bool is_rgt_input (true);
        bool is_bot_input (true);
        bool is_top_input (true);

        FT dxmirror;
        FT dymirror;
        if (exist_two_with_same_x != exist_two_with_same_y) {
          if (exist_two_with_same_x) {
            FT const *s1 = NULL;
            FT const *s2 = NULL;
            FT const *dx = NULL;
            if (samex_pq) {
              s1 = &py;
              s2 = &qy;
              dx = &ry;
            }
            else if (samex_pr) {
              s1 = &py;
              s2 = &ry;
              dx = &qy;
            }
            else if (samex_qr) {
              s1 = &qy;
              s2 = &ry;
              dx = &py;
            }
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs two same x"
                << std::endl;);
            if ( (bot_p == dx) || (top_p == dx) ) {
              CGAL_assertion (
                 ( ( CGAL::compare(*dx, *s1) == SMALLER ) &&
                   ( CGAL::compare(*dx, *s2) == SMALLER )   ) ||
                 ( ( CGAL::compare(*dx, *s1) == LARGER  ) &&
                   ( CGAL::compare(*dx, *s2) == LARGER  )   )   );
              dxmirror = *s1 + *s2 - *dx;
              if (top_p == dx) {
                CGAL_assertion(CGAL::compare(dxmirror, *bot_p) == SMALLER);
                bot_p = &dxmirror;
                is_bot_input = false;
              } else {
                CGAL_assertion(CGAL::compare(dxmirror, *top_p) == LARGER);
                CGAL_assertion( bot_p == dx );
                top_p = &dxmirror;
                is_top_input = false;
              }
            }
          } else {
            CGAL_assertion( exist_two_with_same_y );
            FT const *s1 = NULL;
            FT const *s2 = NULL;
            FT const *dy = NULL;
            if (samey_pq) {
              s1 = &px;
              s2 = &qx;
              dy = &rx;
            }
            else if (samey_pr) {
              s1 = &px;
              s2 = &rx;
              dy = &qx;
            }
            else if (samey_qr) {
              s1 = &qx;
              s2 = &rx;
              dy = &px;
            }
            CGAL_assertion( dy != NULL );
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs two same y"
                << std::endl;);
            if ( (lft_p == dy) || (rgt_p == dy) ) {
              CGAL_assertion (
                 ( ( CGAL::compare(*dy, *s1) == SMALLER ) &&
                   ( CGAL::compare(*dy, *s2) == SMALLER )   ) ||
                 ( ( CGAL::compare(*dy, *s1) == LARGER  ) &&
                   ( CGAL::compare(*dy, *s2) == LARGER  )   )   );
              dymirror = *s1 + *s2 - *dy;
              if (rgt_p == dy) {
                CGAL_assertion(CGAL::compare(dymirror, *lft_p) == SMALLER);
                lft_p = &dymirror;
                is_lft_input = false;
              } else {
                CGAL_assertion(CGAL::compare(dymirror, *rgt_p) == LARGER);
                CGAL_assertion( lft_p == dy );
                rgt_p = &dymirror;
                is_rgt_input = false;
              }
            }
          } // end of case exist_two_with_same_y
        } // end of if exist_two_with_same_x != exist_two_with_same_y

        const bool is_L_shaped =
          exist_two_with_same_x && exist_two_with_same_y;

        const FT two(2);

        CGAL_SDG_DEBUG( std::cout << "debug bs after mirror"
            << "  lft=" << *lft_p << "  rgt=" << *rgt_p
            << "  bot=" << *bot_p << "  top=" << *top_p
            << std::endl ; );

        Comparison_result cmpsides =
          CGAL::compare(*rgt_p - *lft_p, *top_p - *bot_p);

        FT fix1;
        FT fix2;

        bool are_at_three_corners (false);
        if (cmpsides == EQUAL)
        {
          if (is_L_shaped) {
            are_at_three_corners = true;
          }
        }
        else if (cmpsides == LARGER)
        { //diff x > diff y forms a rectangle
          //need to find the movable side of rectangle
          if (exist_two_with_same_x) {
            // expand rectangle both downwards and upwards
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs move both sides"
                << std::endl;);
            fix1 = (*bot_p + *top_p - *rgt_p + *lft_p)/two;
            is_bot_input = false;
            fix2 = (*top_p + *bot_p + *rgt_p - *lft_p)/two;
            is_top_input = false;
            // update bottom and top
            bot_p = &fix1;
            top_p = &fix2;
          }
          else {
            CGAL_assertion( is_bot_input && is_top_input );
            Point_2 const * const pt_lft_p =
              (is_lft_input)? CGAL::cpp11::get<5>(tupx) : NULL;
            Point_2 const * const pt_rgt_p =
              (is_rgt_input)? CGAL::cpp11::get<6>(tupx) : NULL;
            Point_2 const * const pt_top_p = CGAL::cpp11::get<6>(tupy);
            if ((pt_top_p != pt_lft_p) && (pt_top_p != pt_rgt_p)) {
              // lower the bottom side
              fix1 = *top_p - *rgt_p + *lft_p;
              bot_p = &fix1;
              is_bot_input = false;
            }
            else {
              CGAL_assertion_code(
              Point_2 const * const pt_bot_p = CGAL::cpp11::get<5>(tupy);)
              CGAL_assertion(
                  (pt_bot_p != pt_lft_p) && (pt_bot_p != pt_rgt_p) );
              // augment the top side
              fix2 = *bot_p + *rgt_p - *lft_p;
              top_p = &fix2;
              is_top_input = false;
            }
          } // end of not exist_two_with_same_x case
        } // end of cmpsides == LARGER case
        else
        { // px_max.x() - px_min.x() < py_max.y() - py_min.y())
          // diff x < diff y forms a rectangle
          CGAL_assertion( cmpsides == SMALLER );
          if (exist_two_with_same_y) {
            // change both sides
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs move both sides"
                << std::endl;);

            fix1 = (*lft_p + *rgt_p + *top_p - *bot_p)/two;
            is_rgt_input = false;
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs fatten fix1="
                << fix1 << std::endl;);

            fix2 = (*lft_p + *rgt_p - *top_p + *bot_p)/two;
            is_lft_input = false;
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs fatten fix2="
                << fix2 << std::endl;);

            // update right and left
            rgt_p = &fix1;
            lft_p = &fix2;
          }
          else {
            CGAL_assertion( is_lft_input && is_rgt_input );
            Point_2 const * const pt_lft_p = CGAL::cpp11::get<5>(tupx);
            Point_2 const * const pt_bot_p =
              (is_bot_input) ? CGAL::cpp11::get<5>(tupy) : NULL;
            Point_2 const * const pt_top_p =
              (is_top_input) ? CGAL::cpp11::get<6>(tupy) : NULL;
            // find the movable side or sides of the rectangle
            if ((pt_lft_p != pt_bot_p) && (pt_lft_p != pt_top_p)) {
              // augment the right side
              fix1 = *lft_p + *top_p - *bot_p;
              rgt_p = &fix1;
              is_rgt_input = false;
            } else {
              CGAL_assertion_code(
              Point_2 const * const pt_rgt_p = CGAL::cpp11::get<6>(tupx);)
              CGAL_assertion(
                  (pt_rgt_p != pt_bot_p) && (pt_rgt_p != pt_top_p) );
              // diminish from the left side
              fix2 = *rgt_p - *top_p + *bot_p;
              lft_p = &fix2;
              is_lft_input = false;
            }
          } // end of not exist_two_with_same_y case
        } // end of cmpsides == SMALLER case

        CGAL_SDG_DEBUG( std::cout << "debug bs after side fixing "
            << "lft=" << *lft_p
            << "  rgt=" << *rgt_p << "  bot=" << *bot_p << " "
            << "top=" << *top_p << std::endl ; );

        // comparison of query point t with lrbt
        const Comparison_result cxmint = CGAL::compare(*lft_p, t.x());
        if (cxmint == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        const Comparison_result cxtmax = CGAL::compare(t.x(), *rgt_p);
        if (cxtmax == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        const Comparison_result cymint = CGAL::compare(*bot_p, t.y());
        if (cymint == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        const Comparison_result cytmax = CGAL::compare(t.y(), *top_p);
        if (cytmax == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        // here, each comparison value is either SMALLER or EQUAL
        const SmallerEqTuple tup =
          analyze_smalleq(cxtmax, cytmax, cxmint, cymint);
        const size_t count_eq = CGAL::cpp11::get<4>(tup);
        CGAL_assertion( count_eq <= 2 );
        if (count_eq == 0) {
          CGAL_assertion( cxmint == SMALLER && cxtmax == SMALLER &&
                          cymint == SMALLER && cytmax == SMALLER );
          CGAL_SDG_DEBUG(std::cout
              << "debug Side_of_bs return ON_BOUNDED_SIDE" << std::endl;);
          return ON_BOUNDED_SIDE;
        } else {
          CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs on boundary, "
              << "left=" << cxmint << " right=" << cxtmax
              << " bot=" << cymint << " top  =" << cytmax
              << std::endl; );

          if (count_eq == 2) {
            if (are_at_three_corners) {
              return ON_BOUNDARY;
            } else {
              if (is_L_shaped) {
                return ON_UNBOUNDED_SIDE;
              }
            }
          }

          if (are_at_three_corners) {
            CGAL_assertion( count_eq == 1 );
            return ON_BOUNDED_SIDE;
          }

          const bool at_lft = CGAL::cpp11::get<2>(tup);
          if (is_lft_input && at_lft) {
            CGAL_assertion(cxmint == EQUAL);
            CGAL_SDG_DEBUG(std::cout
                << "debug Side_of_bs t on lft input" << std::endl;);
            const FT lfty = (lft_p == &px) ? py : (lft_p == &qx) ? qy : ry;
            const Comparison_result test =
              test1d(*bot_p, *top_p, lfty, t.y());
            if (test != EQUAL) {
              return (test == SMALLER) ?
                     ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
            }
          }

          const bool at_rgt = CGAL::cpp11::get<0>(tup);
          if (is_rgt_input && at_rgt) {
            CGAL_assertion(cxtmax == EQUAL);
            CGAL_SDG_DEBUG(std::cout
                << "debug Side_of_bs t on rgt input" << std::endl;);
            const FT rgty = (rgt_p == &px) ? py : (rgt_p == &qx) ? qy : ry;
            const Comparison_result test =
              test1d(*bot_p, *top_p, rgty, t.y());
            if (test != EQUAL) {
              return (test == SMALLER) ?
                     ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
            }
          }

          const bool at_bot = CGAL::cpp11::get<3>(tup);
          if (is_bot_input && at_bot) {
            CGAL_assertion(cymint == EQUAL);
            CGAL_SDG_DEBUG(std::cout
                << "debug Side_of_bs t on bot input" << std::endl;);
            const FT botx = (bot_p == &py) ? px : (bot_p == &qy) ? qx : rx;
            const Comparison_result test =
              test1d(*lft_p, *rgt_p, botx, t.x());
            if (test != EQUAL) {
              return (test == SMALLER) ?
                     ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
            }
          }

          const bool at_top = CGAL::cpp11::get<1>(tup);
          if (is_top_input && at_top) {
            CGAL_assertion(cytmax == EQUAL);
            CGAL_SDG_DEBUG(std::cout
                << "debug Side_of_bs t on top input" << std::endl;);
            const FT topx = (top_p == &py) ? px : (top_p == &qy) ? qx : rx;
            const Comparison_result test =
              test1d(*lft_p, *rgt_p, topx, t.x());
            if (test != EQUAL) {
              return (test == SMALLER) ?
                     ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
            }
          }

          CGAL_SDG_DEBUG(std::cout
              << "debug Side_of_bs return ON_BOUNDARY" << std::endl;);
          return ON_BOUNDARY;
        }
      }

      inline Comparison_result test1d(
          const FT& A, const FT& B, const FT&C, const FT& D) const
      {
        const FT two(2);
        CGAL_SDG_DEBUG( std::cout << "debug bs test1d entering with ABCD "
            << A << " " << B << " " << C << " " << D << std::endl; );
        return CGAL::compare(CGAL::abs(A+B-two*D), CGAL::abs(A+B-two*C));
      }

    public:

      inline Bounded_side operator()(const Point_2 &p, const Point_2 &q,
                   const Point_2 &r, const Point_2 &t) const
      {
        return predicate(p, q, r, t);
      }
    };

} //namespace CGAL

#endif // CGAL_SIDE_OF_BOUNDED_SQUARE_2_H
