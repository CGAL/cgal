//ON_UNBOUNDED_SIDE = -1,
//ON_BOUNDARY = 0,
//ON_BOUNDED_SIDE = 1

#ifndef CGAL_SIDE_OF_BOUNDED_SQUARE_2_H
#define CGAL_SIDE_OF_BOUNDED_SQUARE_2_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/Orientation_Linf_2.h>

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

      Compare_x_2 compare_x_2;
      Compare_y_2 compare_y_2;

      Orientation_Linf_2_Type orientation_Linf;


      inline Bounded_side predicate(const Point_2 &p, const Point_2 &q,
                  const Point_2 &r, const Point_2 &t) const
      {
        CGAL_SDG_DEBUG(std::cout
            << "debug Side_of_bounded_square_2 entering" << std::endl;);
        CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs (pqrt)= (" << p
            << ") (" << q << ") (" << r << ") (" << t << ")" << std::endl;);

        CGAL_assertion(orientation_Linf(p,q,r) != DEGENERATE);

        //compute the minimum x and maximum x
        Point_2 const * lft_p;
        Point_2 const * rgt_p;
        bool samex_pq (false);
        bool samex_pr (false);
        bool samex_qr (false);
        Comparison_result cmpxpq = compare_x_2(p, q);
        switch(cmpxpq) {
          case SMALLER:
            lft_p = &p;
            rgt_p = &q;
            break;
          case LARGER:
            lft_p = &q;
            rgt_p = &p;
            break;
          default: // EQUAL
            lft_p = &p;
            rgt_p = &q;
            samex_pq = true;
            break;
        }
        Comparison_result cmpxpr = compare_x_2(p, r);
        Comparison_result cmpxqr;
        if (samex_pq) {
          cmpxqr = cmpxpr;
          switch(cmpxpr) {
            case SMALLER:
              rgt_p = &r;
              break;
            case LARGER:
              lft_p = &r;
              break;
            default: // EQUAL is impossible
              CGAL_assertion(false);
              break;
          }
        } else {
          if (lft_p == &p) {
            switch(cmpxpr) {
              case SMALLER:
                cmpxqr = compare_x_2(q, r);
                switch(cmpxqr) {
                  case SMALLER:
                    rgt_p = &r;
                    break;
                  case LARGER:
                    break;
                  default:
                    // q and r have the same x coord
                    samex_qr = true;
                    break;
                }
                break;
              case LARGER:
                cmpxqr = cmpxpr;
                lft_p = &r;
                break;
              default:
                // p and r have the same x coord
                cmpxqr = - cmpxpr;
                samex_pr = true;
                break;
            }
          } else { // lft_p == &q
            switch(cmpxpr) {
              case SMALLER:
                cmpxqr = cmpxpr;
                rgt_p = &r;
                break;
              case LARGER:
                cmpxqr = compare_x_2(q, r);
                switch(cmpxqr) {
                  case SMALLER:
                    break;
                  case LARGER:
                    lft_p = &r;
                    break;
                  default:
                    // q and r have the same x coord
                    samex_qr = true;
                    break;
                }
                break;
              default:
                // p and r have the same x coord
                cmpxqr = - cmpxpq;
                samex_pr = true;
                break;
            }
          }
        }

        //compute the minimum y and maximum y
        Point_2 const * bot_p;
        Point_2 const * top_p;
        bool samey_pq (false);
        bool samey_pr (false);
        bool samey_qr (false);
        Comparison_result cmpypq = compare_y_2(p, q);
        switch(cmpypq) {
          case SMALLER:
            bot_p = &p;
            top_p = &q;
            break;
          case LARGER:
            bot_p = &q;
            top_p = &p;
            break;
          default: // EQUAL
            bot_p = &p;
            top_p = &q;
            samey_pq = true;
            break;
        }
        Comparison_result cmpypr = compare_y_2(p, r);
        CGAL_SDG_DEBUG( std::cout << "debug bs " << " p=" << p
            << "  r=" << r << "  cmpypr=" << cmpypr << std::endl ; );
        Comparison_result cmpyqr;
        if (samey_pq) {
          cmpyqr = cmpypr;
          switch(cmpypr) {
            case SMALLER:
              top_p = &r;
              break;
            case LARGER:
              bot_p = &r;
              break;
            default: // EQUAL is impossible
              CGAL_assertion(false);
              break;
          }
        } else {
          if (bot_p == &p) {
            switch(cmpypr) {
              case SMALLER:
                cmpyqr = compare_y_2(q, r);
                switch(cmpyqr) {
                  case SMALLER:
                    top_p = &r;
                    break;
                  case LARGER:
                    break;
                  default:
                    // q and r have the same y coord
                    samey_qr = true;
                    break;
                }
                break;
              case LARGER:
                cmpyqr = cmpypr;
                bot_p = &r;
                break;
              default:
                // p and r have the same y coord
                cmpyqr = - cmpypq;
                samey_pr = true;
                break;
            }
          } else { // bot_p == &q
            CGAL_SDG_DEBUG(std::cout << "debug cmpypr=" << cmpypr
                << " and q is lower than p" << std::endl; );
            switch(cmpypr) {
              case SMALLER:
                cmpyqr = cmpypr;
                top_p = &r;
                break;
              case LARGER:
                cmpyqr = compare_y_2(q, r);
                CGAL_SDG_DEBUG( std::cout << "debug bs " << " q=" << q
                    << "  r=" << r << "  cmpyqr=" << cmpyqr << std::endl ; );
                switch(cmpyqr) {
                  case SMALLER:
                    break;
                  case LARGER:
                    bot_p = &r;
                    break;
                  default:
                    // q and r have the same y coord
                    samey_qr = true;
                    break;
                }
                break;
              default:
                // p and r have the same y coord
                cmpyqr = - cmpypq;
                samey_pr = true;
                break;
            }
          }
        }

        bool is_lft_input (true);
        bool is_rgt_input (true);
        bool is_bot_input (true);
        bool is_top_input (true);

        CGAL_SDG_DEBUG(std::cout << "debug bs " << " lft=" << *lft_p <<
            "  rgt=" << *rgt_p << "  bot=" << *bot_p << "  top=" << *top_p
            << std::endl; );
        CGAL_SDG_DEBUG(std::cout << "debug bs cmpypq=" << cmpypq <<
            " cmpypr=" << cmpypr << " cmpyqr=" << cmpyqr
            << std::endl; );

        // check if two points have the same x or y coordinate
        Point_2 const *s1;
        Point_2 const *s2;

        Point_2 const *dx = NULL;
        // check if two points have the same x coordinate
        bool exist_two_with_same_x (false);
        if (samex_pq) {
          exist_two_with_same_x = true;
          s1 = &p;
          s2 = &q;
          dx = &r;
        }
        else if (samex_pr) {
          exist_two_with_same_x = true;
          s1 = &p;
          s2 = &r;
          dx = &q;
        }
        else if (samex_qr) {
          exist_two_with_same_x = true;
          s1 = &q;
          s2 = &r;
          dx = &p;
        }

        Point_2 dxmirror;
        if (exist_two_with_same_x) {
          CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs two same x"
              << std::endl;);
          if ( ( ( compare_y_2(*dx, *s1) == SMALLER ) and
                 ( compare_y_2(*dx, *s2) == SMALLER )   ) or
               ( ( compare_y_2(*dx, *s1) == LARGER  ) and
                 ( compare_y_2(*dx, *s2) == LARGER  )   )   )   {
            dxmirror = Point_2 (dx->x(), s1->y() + s2->y() - dx->y());
            if (compare_y_2(dxmirror, *bot_p) == SMALLER) {
              bot_p = &dxmirror;
              is_bot_input = false;
            }
            if (compare_y_2(dxmirror, *top_p) == LARGER) {
              top_p = &dxmirror;
              is_top_input = false;
            }
          }
        }

        Point_2 const *dy = NULL;
        // check if two points have the same y coordinate
        bool exist_two_with_same_y (false);
        if (samey_pq) {
          exist_two_with_same_y = true;
          s1 = &p;
          s2 = &q;
          dy = &r;
        }
        else if (samey_pr) {
          exist_two_with_same_y = true;
          s1 = &p;
          s2 = &r;
          dy = &q;
        }
        else if (samey_qr) {
          exist_two_with_same_y = true;
          s1 = &q;
          s2 = &r;
          dy = &p;
        }

        Point_2 dymirror;
        if (exist_two_with_same_y) {
          CGAL_assertion( dy != NULL );
          CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs two same y"
              << std::endl;);
          if ( ( ( compare_x_2(*dy, *s1) == SMALLER ) and
                 ( compare_x_2(*dy, *s2) == SMALLER )   ) or
               ( ( compare_x_2(*dy, *s1) == LARGER  ) and
                 ( compare_x_2(*dy, *s2) == LARGER  )   )   )   {
            dymirror = Point_2 (s1->x() + s2->x() - dy->x(), dy->y());
            if (compare_x_2(dymirror, *lft_p) == SMALLER) {
              lft_p = &dymirror;
              is_lft_input = false;
            }
            if (compare_x_2(dymirror, *rgt_p) == LARGER) {
              rgt_p = &dymirror;
              is_rgt_input = false;
            }
          }
        }

        const FT half(0.5);

        CGAL_SDG_DEBUG( std::cout << "debug bs after mirror "
            << "lft=" << *lft_p
            << "  rgt=" << *rgt_p << "  bot=" << *bot_p << " "
            << "top=" << *top_p << std::endl ; );

        Comparison_result cmpsides =
          CGAL::compare(rgt_p->x() - lft_p->x(), top_p->y() - bot_p->y());

        Point_2 fix1;
        Point_2 fix2;

        if (cmpsides == EQUAL)
        {
          // do nothing, lrbt are fine
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
               (bot_p->y() + top_p->y() - rgt_p->x() + lft_p->x())*half);
            is_bot_input = false;
            fix2 = Point_2 (
                top_p->x(),
               (top_p->y() + bot_p->y() + rgt_p->x() - lft_p->x())*half);
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
               (lft_p->x() + rgt_p->x() + top_p->y() - bot_p->y())*half,
                rgt_p->y());
            is_rgt_input = false;
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs fatten fix1="
                << fix1 << std::endl;);

            fix2 = Point_2 (
               (lft_p->x() + rgt_p->x() - top_p->y() + bot_p->y())*half,
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
        Comparison_result cxmint = compare_x_2(*lft_p, t);
        if (cxmint == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        Comparison_result cxtmax = compare_x_2(t, *rgt_p);
        if (cxtmax == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        Comparison_result cymint = compare_y_2(*bot_p, t);
        if (cymint == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        Comparison_result cytmax = compare_y_2(t, *top_p);
        if (cytmax == LARGER) {
          return ON_UNBOUNDED_SIDE;
        }
        // here, each comparison value is either SMALLER or EQUAL
        if( cxmint == SMALLER and
            cxtmax == SMALLER and
            cymint == SMALLER and
            cytmax == SMALLER   ) {
          CGAL_SDG_DEBUG(std::cout
              << "debug Side_of_bs return ON_BOUNDED_SIDE" << std::endl;);
          return ON_BOUNDED_SIDE;
        } else {
          CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs on boundary, "
              << "left=" << cxmint << " right=" << cxtmax
              << " bot=" << cymint << " top  =" << cytmax
              << std::endl; );

          if (is_lft_input and (cxmint == EQUAL)) {
            CGAL_SDG_DEBUG(std::cout
                << "debug Side_of_bs t on lft input" << std::endl;);
            Comparison_result test =
              test1d(bot_p->y(), top_p->y(), lft_p->y(), t.y());
            if (test != EQUAL) {
              return (test == SMALLER) ?
                     ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
            }
          }

          if (is_rgt_input and (cxtmax == EQUAL)) {
            CGAL_SDG_DEBUG(std::cout
                << "debug Side_of_bs t on rgt input" << std::endl;);
            Comparison_result test =
              test1d(bot_p->y(), top_p->y(), rgt_p->y(), t.y());
            if (test != EQUAL) {
              return (test == SMALLER) ?
                     ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
            }
          }

          if (is_bot_input and (cymint == EQUAL)) {
            CGAL_SDG_DEBUG(std::cout
                << "debug Side_of_bs t on bot input" << std::endl;);
            Comparison_result test =
              test1d(lft_p->x(), rgt_p->x(), bot_p->x(), t.x());
            if (test != EQUAL) {
              return (test == SMALLER) ?
                     ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
            }
          }

          if (is_top_input and (cytmax == EQUAL)) {
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
