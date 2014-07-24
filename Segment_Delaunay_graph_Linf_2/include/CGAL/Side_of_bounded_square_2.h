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

      typedef CGAL::cpp11::tuple<bool, bool, bool, bool, size_t>
              SmallerEqTuple;

      Compare_x_2 compare_x_2;
      Compare_y_2 compare_y_2;

      Orientation_Linf_2_Type orientation_Linf;

      template <class Compare, typename E>
      inline void minmax(
          const E &p, const E &q, const E &r,
          Point_2 const * & min_p, Point_2 const * & max_p,
          bool & samepq, bool & samepr, bool & sameqr
          ) const
      {
        Compare cmp;
        samepq = false;
        samepr = false;
        sameqr = false;
        const Comparison_result cmppq = cmp(p, q);
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
        const Comparison_result cmppr = cmp(p, r);
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
                cmpqr = cmp(q, r);
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
                cmpqr = cmp(q, r);
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
        CGAL_SDG_DEBUG(std::cout << "debug minmax cmppq=" << cmppq
            << " cmppr=" << cmppr << " cmpqr=" << cmpqr << std::endl; );
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
        Point_2 const * lft_p (NULL);
        Point_2 const * rgt_p (NULL);
        bool samex_pq (false);
        bool samex_pr (false);
        bool samex_qr (false);
        minmax<Compare_x_2>(p, q, r,
            lft_p, rgt_p, samex_pq, samex_pr, samex_qr);
        CGAL_assertion(lft_p != NULL);
        CGAL_assertion(rgt_p != NULL);

        //compute the minimum y and maximum y
        Point_2 const * bot_p (NULL);
        Point_2 const * top_p (NULL);
        bool samey_pq (false);
        bool samey_pr (false);
        bool samey_qr (false);
        minmax<Compare_y_2>(p, q, r,
            bot_p, top_p, samey_pq, samey_pr, samey_qr);
        CGAL_assertion(bot_p != NULL);
        CGAL_assertion(top_p != NULL);

        CGAL_SDG_DEBUG(std::cout << "debug bs before mirror"
            << "  lft=" << *lft_p << "  rgt=" << *rgt_p
            << "  bot=" << *bot_p << "  top=" << *top_p
            << std::endl; );

        const bool exist_two_with_same_x = samex_pq or samex_pr or samex_qr;
        const bool exist_two_with_same_y = samey_pq or samey_pr or samey_qr;

        bool is_lft_input (true);
        bool is_rgt_input (true);
        bool is_bot_input (true);
        bool is_top_input (true);

        Point_2 dxmirror;
        Point_2 dymirror;
        if (exist_two_with_same_x != exist_two_with_same_y) {
          if (exist_two_with_same_x) {
            Point_2 const *s1 = NULL;
            Point_2 const *s2 = NULL;
            Point_2 const *dx = NULL;
            if (samex_pq) {
              s1 = &p;
              s2 = &q;
              dx = &r;
            }
            else if (samex_pr) {
              s1 = &p;
              s2 = &r;
              dx = &q;
            }
            else if (samex_qr) {
              s1 = &q;
              s2 = &r;
              dx = &p;
            }
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs two same x"
                << std::endl;);
            if ( (bot_p == dx) or (top_p == dx) ) {
              CGAL_assertion (
                 ( ( compare_y_2(*dx, *s1) == SMALLER ) and
                   ( compare_y_2(*dx, *s2) == SMALLER )   ) or
                 ( ( compare_y_2(*dx, *s1) == LARGER  ) and
                   ( compare_y_2(*dx, *s2) == LARGER  )   )   );
              dxmirror = Point_2 (dx->x(), s1->y() + s2->y() - dx->y());
              if (top_p == dx) {
                CGAL_assertion( compare_y_2(dxmirror, *bot_p) == SMALLER );
                bot_p = &dxmirror;
                is_bot_input = false;
              } else {
                CGAL_assertion( compare_y_2(dxmirror, *top_p) == LARGER );
                CGAL_assertion( bot_p == dx );
                top_p = &dxmirror;
                is_top_input = false;
              }
            }
          } else {
            CGAL_assertion( exist_two_with_same_y );
            Point_2 const *s1 = NULL;
            Point_2 const *s2 = NULL;
            Point_2 const *dy = NULL;
            if (samey_pq) {
              s1 = &p;
              s2 = &q;
              dy = &r;
            }
            else if (samey_pr) {
              s1 = &p;
              s2 = &r;
              dy = &q;
            }
            else if (samey_qr) {
              s1 = &q;
              s2 = &r;
              dy = &p;
            }
            CGAL_assertion( dy != NULL );
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs two same y"
                << std::endl;);
            if ( (lft_p == dy) or (rgt_p == dy) ) {
              CGAL_assertion (
                 ( ( compare_x_2(*dy, *s1) == SMALLER ) and
                   ( compare_x_2(*dy, *s2) == SMALLER )   ) or
                 ( ( compare_x_2(*dy, *s1) == LARGER  ) and
                   ( compare_x_2(*dy, *s2) == LARGER  )   )   );
              dymirror = Point_2 (s1->x() + s2->x() - dy->x(), dy->y());
              if (rgt_p == dy) {
                CGAL_assertion( compare_x_2(dymirror, *lft_p) == SMALLER );
                lft_p = &dymirror;
                is_lft_input = false;
              } else {
                CGAL_assertion(compare_x_2(dymirror, *rgt_p) == LARGER);
                CGAL_assertion( lft_p == dy );
                rgt_p = &dymirror;
                is_rgt_input = false;
              }
            }
          } // end of case exist_two_with_same_y
        } // end of if exist_two_with_same_x != exist_two_with_same_y

        const bool is_L_shaped =
          exist_two_with_same_x and exist_two_with_same_y;

        const FT two(2);

        CGAL_SDG_DEBUG( std::cout << "debug bs after mirror"
            << "  lft=" << *lft_p << "  rgt=" << *rgt_p
            << "  bot=" << *bot_p << "  top=" << *top_p
            << std::endl ; );

        Comparison_result cmpsides =
          CGAL::compare(rgt_p->x() - lft_p->x(), top_p->y() - bot_p->y());

        Point_2 fix1;
        Point_2 fix2;

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
          if (compare_x_2(*lft_p, *top_p) == SMALLER and
              compare_x_2(*top_p, *rgt_p) == SMALLER   ) {
            // lower the bottom side
            fix1 = Point_2 (bot_p->x(), top_p->y() - rgt_p->x() + lft_p->x());
            bot_p = &fix1;
            is_bot_input = false;
          }
          else if (compare_x_2(*lft_p, *bot_p) == SMALLER and
                   compare_x_2(*bot_p, *rgt_p) == SMALLER   ){
            // augment the top side
            fix2 = Point_2 (top_p->x(), bot_p->y() + rgt_p->x() - lft_p->x());
            top_p = &fix2;
            is_top_input = false;
          }
          else {
            // expand rectangle both downwards and upwards
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs move both sides"
                << std::endl;);
            fix1 = Point_2 (
                bot_p->x(),
               (bot_p->y() + top_p->y() - rgt_p->x() + lft_p->x())/two);
            is_bot_input = false;
            fix2 = Point_2 (
                top_p->x(),
               (top_p->y() + bot_p->y() + rgt_p->x() - lft_p->x())/two);
            is_top_input = false;

            // update bottom and top
            bot_p = &fix1;
            top_p = &fix2;
          }
        }
        else
        { // px_max.x() - px_min.x() < py_max.y() - py_min.y())
          // diff x < diff y forms a rectangle

          // find the movable side or sides of the rectangle
          if (compare_y_2(*lft_p, *bot_p) == LARGER and
              compare_y_2(*lft_p, *top_p) == SMALLER ) {
            // augment the right side
            fix1 = Point_2 (lft_p->x() + top_p->y() - bot_p->y(), rgt_p->y());
            rgt_p = &fix1;
            is_rgt_input = false;
          }
          else if (compare_y_2(*rgt_p,*bot_p) == LARGER and
                   compare_y_2(*rgt_p,*top_p) == SMALLER ){
            // diminish from the left side
            fix2 = Point_2 (rgt_p->x() - top_p->y() + bot_p->y(), lft_p->y());
            lft_p = &fix2;
            is_lft_input = false;
          }
          else {
            // change both sides
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs move both sides"
                << std::endl;);

            fix1 = Point_2 (
               (lft_p->x() + rgt_p->x() + top_p->y() - bot_p->y())/two,
                rgt_p->y());
            is_rgt_input = false;
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs fatten fix1="
                << fix1 << std::endl;);

            fix2 = Point_2 (
               (lft_p->x() + rgt_p->x() - top_p->y() + bot_p->y())/two,
                lft_p->y());
            is_lft_input = false;
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs fatten fix2="
                << fix2 << std::endl;);

            // update right and left
            rgt_p = &fix1;
            lft_p = &fix2;
          }
        }

        CGAL_SDG_DEBUG( std::cout << "debug bs after side fixing "
            << "lft=" << *lft_p
            << "  rgt=" << *rgt_p << "  bot=" << *bot_p << " "
            << "top=" << *top_p << std::endl ; );

        // comparison of query point t with lrbt
        const Comparison_result cxmint = compare_x_2(*lft_p, t);
        if (cxmint == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        const Comparison_result cxtmax = compare_x_2(t, *rgt_p);
        if (cxtmax == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        const Comparison_result cymint = compare_y_2(*bot_p, t);
        if (cymint == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        const Comparison_result cytmax = compare_y_2(t, *top_p);
        if (cytmax == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        // here, each comparison value is either SMALLER or EQUAL
        const SmallerEqTuple tup =
          analyze_smalleq(cxtmax, cytmax, cxmint, cymint);
        const size_t count_eq = CGAL::cpp11::get<4>(tup);
        CGAL_assertion( count_eq >= 0 );
        CGAL_assertion( count_eq <= 2 );
        if (count_eq == 0) {
          CGAL_assertion( cxmint == SMALLER and cxtmax == SMALLER and
                          cymint == SMALLER and cytmax == SMALLER );
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

          const bool at_lft = CGAL::cpp11::get<2>(tup);
          if (is_lft_input and at_lft) {
            CGAL_assertion(cxmint == EQUAL);
            CGAL_SDG_DEBUG(std::cout
                << "debug Side_of_bs t on lft input" << std::endl;);
            Comparison_result test =
              test1d(bot_p->y(), top_p->y(), lft_p->y(), t.y());
            if (test != EQUAL) {
              return (test == SMALLER) ?
                     ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
            }
          }

          const bool at_rgt = CGAL::cpp11::get<0>(tup);
          if (is_rgt_input and at_rgt) {
            CGAL_assertion(cxtmax == EQUAL);
            CGAL_SDG_DEBUG(std::cout
                << "debug Side_of_bs t on rgt input" << std::endl;);
            Comparison_result test =
              test1d(bot_p->y(), top_p->y(), rgt_p->y(), t.y());
            if (test != EQUAL) {
              return (test == SMALLER) ?
                     ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
            }
          }

          const bool at_bot = CGAL::cpp11::get<3>(tup);
          if (is_bot_input and at_bot) {
            CGAL_assertion(cymint == EQUAL);
            CGAL_SDG_DEBUG(std::cout
                << "debug Side_of_bs t on bot input" << std::endl;);
            Comparison_result test =
              test1d(lft_p->x(), rgt_p->x(), bot_p->x(), t.x());
            if (test != EQUAL) {
              return (test == SMALLER) ?
                     ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
            }
          }

          const bool at_top = CGAL::cpp11::get<1>(tup);
          if (is_top_input and at_top) {
            CGAL_assertion(cytmax == EQUAL);
            CGAL_SDG_DEBUG(std::cout
                << "debug Side_of_bs t on top input" << std::endl;);
            Comparison_result test =
              test1d(lft_p->x(), rgt_p->x(), top_p->x(), t.x());
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
