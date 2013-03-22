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

        Point_2 px_min, px_max, py_min, py_max;

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

        if ( (cmpxpq == LARGER or cmpxpq == EQUAL)
              and (cmpxpr == LARGER or cmpxpr == EQUAL) ){
          px_max = p;
          px_min = (cmpxqr == LARGER or cmpxqr == EQUAL) ? r : q;
          CGAL_assertion(compare_x_2(px_max, *rgt_p) == EQUAL);
          CGAL_assertion(compare_x_2(px_min, *lft_p) == EQUAL);
        }
        else if ((cmpxpq == SMALLER or cmpxpq == EQUAL)
                 and (cmpxqr == LARGER or cmpxqr == EQUAL) ){
          px_max = q;
          px_min = (cmpxpr == LARGER or cmpxpr == EQUAL) ? r : p;
          CGAL_assertion(compare_x_2(px_max, *rgt_p) == EQUAL);
          CGAL_assertion(compare_x_2(px_min, *lft_p) == EQUAL);
        }
        else {
          px_max = r;
          px_min = (cmpxpq == LARGER or cmpxpq == EQUAL) ? q : p;
          CGAL_assertion(compare_x_2(px_max, *rgt_p) == EQUAL);
          CGAL_assertion(compare_x_2(px_min, *lft_p) == EQUAL);
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

        CGAL_SDG_DEBUG(std::cout << "debug bs " << " lft=" << *lft_p <<
            "  rgt=" << *rgt_p << "  bot=" << *bot_p << "  top=" << *top_p
            << std::endl; );
        CGAL_SDG_DEBUG(std::cout << "debug bs cmpypq=" << cmpypq <<
            " cmpypr=" << cmpypr << " cmpyqr=" << cmpyqr
            << std::endl; );

        //compute the minimum y and maximum y
        if ((cmpypq == LARGER or cmpypq == EQUAL)
            and (cmpypr == LARGER or cmpypr == EQUAL) ){
          py_max = p;
          py_min = (cmpyqr == LARGER or cmpyqr == EQUAL) ? r : q;
          CGAL_assertion(compare_y_2(py_max, *top_p) == EQUAL);
          CGAL_assertion(compare_y_2(py_min, *bot_p) == EQUAL);
        }
        else if ((cmpypq == SMALLER or cmpypq == EQUAL)
                 and (cmpyqr == LARGER or cmpyqr == EQUAL) ){
          py_max = q;
          py_min = (cmpypr == LARGER or cmpypr == EQUAL) ? r : p;
          CGAL_assertion(compare_y_2(py_max, *top_p) == EQUAL);
          CGAL_assertion(compare_y_2(py_min, *bot_p) == EQUAL);
        }
        else {
          py_max = r;
          py_min = (cmpypq == LARGER or cmpypq == EQUAL) ? q : p;
          CGAL_assertion(compare_y_2(py_max, *top_p) == EQUAL);
          CGAL_assertion(compare_y_2(py_min, *bot_p) == EQUAL);
        }


        // check if two points have the same x or y coordinate
        Point_2 const *s1;
        Point_2 const *s2;
        Point_2 const *dx;
        Point_2 const *dy;

        bool exist_two_with_same_x (false);
        // check if two points have the same x coordinate
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

        if (exist_two_with_same_x) {
          CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs two same x"
              << std::endl;);
          if ( ( ( compare_y_2(*dx, *s1) == SMALLER ) and
                 ( compare_y_2(*dx, *s2) == SMALLER )   ) or
               ( ( compare_y_2(*dx, *s1) == LARGER  ) and
                 ( compare_y_2(*dx, *s2) == LARGER  )   )   )   {
            Point_2 dxmirror (dx->x(), s1->y() + s2->y() - dx->y());
            if (compare_y_2(dxmirror, py_min) == SMALLER) {
              py_min = dxmirror;
            }
            if (compare_y_2(dxmirror, py_max) == LARGER) {
              py_max = dxmirror;
            }
          }
        }

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

        if (exist_two_with_same_y) {
          CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs two same y"
              << std::endl;);
          if ( ( ( compare_x_2(*dy, *s1) == SMALLER ) and
                 ( compare_x_2(*dy, *s2) == SMALLER )   ) or
               ( ( compare_x_2(*dy, *s1) == LARGER  ) and
                 ( compare_x_2(*dy, *s2) == LARGER  )   )   )   {
            Point_2 dymirror (s1->x() + s2->x() - dy->x(), dy->y());
            if (compare_x_2(dymirror, px_min) == SMALLER) {
              px_min = dymirror;
            }
            if (compare_x_2(dymirror, px_max) == LARGER) {
              px_max = dymirror;
            }
          }
        }

        Point_2 pmin, pmax;
        //pmin and pmax define the bounded square by p,q,r

        const FT half(0.5);

        CGAL_SDG_DEBUG( std::cout << "debug bs " << "px_min=" << px_min
            << "  px_max=" << px_max << "  py_min=" << py_min << " "
            << "py_max=" << py_max << std::endl ; );

        Comparison_result cmpsides =
          CGAL::compare(px_max.x() - px_min.x(), py_max.y() - py_min.y());

        if (cmpsides == EQUAL)
        { //diff x == diff y forms a square
          pmin = Point_2(px_min.x(), py_min.y());
          pmax = Point_2(px_max.x(), py_max.y());
        }
        else if (cmpsides == LARGER)
        { //diff x > diff y forms a rectangle
          //need to find the movable side of rectangle
          if (compare_x_2(px_min,py_max)==SMALLER &&
              compare_x_2(py_max,px_max)==SMALLER   ) {
          //move the lower side
            pmin = Point_2(px_min.x(),py_max.y() - px_max.x() + px_min.x());
            pmax = Point_2(px_max.x(),py_max.y());
          }
          else if (compare_x_2(px_min,py_min)==SMALLER &&
                   compare_x_2(py_min,px_max)==SMALLER   ){
            //move the upper side
            pmin = Point_2(px_min.x(),py_min.y());
            pmax = Point_2(px_max.x(),py_min.y() + px_max.x() - px_min.x());
          }
          else {
            //move both sides
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs move both sides"
                << std::endl;);
            pmin = Point_2(
                px_min.x(),
               (py_max.y() + py_min.y() - px_max.x() + px_min.x())*half);
            pmax = Point_2(
                px_max.x(),
               (py_max.y() + py_min.y() + px_max.x() - px_min.x())*half);
          }
        }
        else
        {//px_max.x() - px_min.x() < py_max.y() - py_min.y())
         // diff x < diff y forms a rectangle
            //need to find the movable side of rectangle
          if (compare_y_2(px_min,py_min)==LARGER &&
              compare_y_2(px_min,py_max)==SMALLER ) {
          //move the right side
            pmin = Point_2(px_min.x(),py_min.y());
            pmax = Point_2(px_min.x() + py_max.y() - py_min.y(),py_max.y());
          }
          else if (compare_y_2(px_max,py_min)==LARGER &&
                   compare_y_2(px_max,py_max)==SMALLER ){
          //move the left side
            pmin = Point_2(px_max.x() - py_max.y() + py_min.y(),py_min.y());
            pmax = Point_2(px_max.x(),py_max.y());
          }
          else {
            //move both sides
            CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs move both sides"
                << std::endl;);
            pmin = Point_2(
                (px_min.x() + px_max.x() - py_max.y() + py_min.y())*half,
                 py_min.y() );
            pmax = Point_2(
                (px_min.x() + px_max.x() + py_max.y() - py_min.y())*half,
                 py_max.y() );
          }
        }

        CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs pmin=" << pmin
                  << " pmax=" << pmax << " t=" << t << std::endl;);

        //Now we answer the predicate bounded side of square
        Comparison_result cxmint = compare_x_2(pmin, t);
        Comparison_result cxtmax = compare_x_2(t, pmax);
        Comparison_result cymint = compare_y_2(pmin, t);
        Comparison_result cytmax = compare_y_2(t, pmax);
        if( cxmint == SMALLER and
            cxtmax == SMALLER and
            cymint == SMALLER and
            cytmax == SMALLER   ) {
          CGAL_SDG_DEBUG(std::cout
              << "debug Side_of_bs return ON_BOUNDED_SIDE" << std::endl;);
          return ON_BOUNDED_SIDE;
        }
        else if (cxmint == LARGER or
                 cxtmax == LARGER or
                 cymint == LARGER or
                 cytmax == LARGER   ) {
          CGAL_SDG_DEBUG(std::cout
              << "debug Side_of_bs return ON_UNBOUNDED_SIDE" << std::endl;);
          return ON_UNBOUNDED_SIDE;
        }
        else {
          CGAL_SDG_DEBUG(std::cout << "debug Side_of_bs on boundary, "
              << "left=" << cxmint << " right=" << cxtmax
              << " bot=" << cymint << " top  =" << cytmax
              << std::endl; );
          CGAL_SDG_DEBUG(std::cout
              << "debug Side_of_bs return ON_BOUNDARY" << std::endl;);
          return ON_BOUNDARY;
        }
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
