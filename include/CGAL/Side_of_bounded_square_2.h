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
      
            
      Bounded_side predicate(const Point_2 &p, const Point_2 &q,
                  const Point_2 &r, const Point_2 &t) const
      {
        std::cout << "debug Side_of_bounded_square_2 entering" 
                  << std::endl;
        std::cout << "debug Side_of_bs (pqrt)= (" << p << ") (" 
                  << q << ") (" << r << ") (" << t << ")" << std::endl;

        Point_2 px_min, px_max, py_min, py_max;

        //should be added as CGAL precondition
        if (orientation_Linf(p,q,r) == DEGENERATE) {
          std::cout << "debug: error: p,q,r are monotone!" << std::endl;
          exit(0);
        }
        
        Comparison_result cmpxpq = compare_x_2(p, q);
        Comparison_result cmpypq = compare_y_2(p, q);
        Comparison_result cmpxpr = compare_x_2(p, r);
        Comparison_result cmpypr = compare_y_2(p, r);
        Comparison_result cmpxqr = compare_x_2(q, r);
        Comparison_result cmpyqr = compare_y_2(q, r);
        
        //compute the minimum x and maximum x
        if ( (cmpxpq == LARGER or cmpxpq == EQUAL) 
              and (cmpxpr == LARGER or cmpxpr == EQUAL) ){
          px_max = p;
          px_min = (cmpxqr == LARGER or cmpxqr == EQUAL) ? r : q;          
        }
        else if ((cmpxpq == SMALLER or cmpxpq == EQUAL) 
                 and (cmpxqr == LARGER or cmpxqr == EQUAL) ){
          px_max = q;
          px_min = (cmpxpr == LARGER or cmpxpr == EQUAL) ? r : p;
        }
        else {
          px_max = r;
          px_min = (cmpxpq == LARGER or cmpxpq == EQUAL) ? q : p; 
        }
        //compute the minimum y and maximum y
        if ((cmpypq == LARGER or cmpypq == EQUAL) 
            and (cmpypr == LARGER or cmpypr == EQUAL) ){
          py_max = p;
          py_min = (cmpyqr == LARGER or cmpyqr == EQUAL) ? r : q;          
        }
        else if ((cmpypq == SMALLER or cmpypq == EQUAL) 
                 and (cmpyqr == LARGER or cmpyqr == EQUAL) ){
          py_max = q;
          py_min = (cmpypr == LARGER or cmpypr == EQUAL) ? r : p;
        }
        else {
          py_max = r;
          py_min = (cmpypq == LARGER or cmpypq == EQUAL) ? q : p; 
        }
        

        // check if two points have the same x or y coordinate 
        bool exist_two_with_same_x;
        bool exist_two_with_same_y;
        Point_2 s1, s2, dx, dy;

        // check if two points have the same x coordinate
        if (cmpxpq == EQUAL) {
          exist_two_with_same_x = true;
          s1 = p;
          s2 = q;
          dx = r;
        } 
        else if (cmpxpr == EQUAL) {
          exist_two_with_same_x = true;
          s1 = p;
          s2 = r;
          dx = q;
        } 
        else if (cmpxqr == EQUAL) {
          exist_two_with_same_x = true;
          s1 = q;
          s2 = r;
          dx = p;
        } 
        else {
          exist_two_with_same_x = false;
        }

        if (exist_two_with_same_x) {
          std::cout << "debug Side_of_bs two same x" << std::endl;
          if ( ( ( compare_y_2(dx, s1) == SMALLER ) and
                 ( compare_y_2(dx, s2) == SMALLER )   ) or
               ( ( compare_y_2(dx, s1) == LARGER  ) and
                 ( compare_y_2(dx, s2) == LARGER  )   )   )   {
            Point_2 dxmirror (dx.x(), s1.y()+s2.y()-dx.y());
            if (compare_y_2(dxmirror, py_min) == SMALLER) {
              py_min = dxmirror;
            }
            if (compare_y_2(dxmirror, py_max) == LARGER) {
              py_max = dxmirror;
            }
          }
        }

        // check if two points have the same y coordinate
        if (cmpypq == EQUAL) {
          exist_two_with_same_y = true;
          s1 = p;
          s2 = q;
          dy = r;
        } 
        else if (cmpypr == EQUAL) {
          exist_two_with_same_y = true;
          s1 = p;
          s2 = r;
          dy = q;
        } 
        else if (cmpyqr == EQUAL) {
          exist_two_with_same_y = true;
          s1 = q;
          s2 = r;
          dy = p;
        } 
        else {
          exist_two_with_same_y = false;
        }

        if (exist_two_with_same_y) {
          std::cout << "debug Side_of_bs two same y" << std::endl;
          if ( ( ( compare_x_2(dy, s1) == SMALLER ) and
                 ( compare_x_2(dy, s2) == SMALLER )   ) or
               ( ( compare_x_2(dy, s1) == LARGER  ) and
                 ( compare_x_2(dy, s2) == LARGER  )   )   )   {
            Point_2 dymirror (s1.x()+s2.x()-dy.x(), dy.y());
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
        
        if (px_max.x() - px_min.x() == py_max.y() - py_min.y()) 
        { //diff x == diff y forms a square
          pmin = Point_2(px_min.x(), py_min.y());
          pmax = Point_2(px_max.x(), py_max.y());
        }
        else if (px_max.x() - px_min.x() > py_max.y() - py_min.y())
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
            std::cout << "debug Side_of_bs move both sides" << std::endl;
            pmin = Point_2(
                px_min.x(), 
               (py_max.y() + py_min.y() - px_max.x() + px_min.x())*FT(0.5));
            pmax = Point_2(
                px_max.x(), 
               (py_max.y() + py_min.y() + px_max.x() - px_min.x())*FT(0.5));
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
            std::cout << "debug Side_of_bs move both sides" << std::endl;
            pmin = Point_2(
                (px_min.x() + px_max.x() - py_max.y() + py_min.y())*FT(0.5),
                 py_min.y() );
            pmax = Point_2(
                (px_min.x() + px_max.x() + py_max.y() - py_min.y())*FT(0.5),
                 py_max.y() );
          }
        }

        std::cout << "debug Side_of_bs pmin=" << pmin
                  << " pmax=" << pmax << " t=" << t << std::endl; 

        //Now we answer the predicate bounded side of square
        if( compare_x_2(pmin, t) == SMALLER && 
            compare_x_2(t, pmax) == SMALLER && 
            compare_y_2(pmin, t) == SMALLER && 
            compare_y_2(t, pmax) == SMALLER   ) {
          std::cout << "debug Side_of_bs return ON_BOUNDED_SIDE" << std::endl;
          return ON_BOUNDED_SIDE;
        }
        else if (compare_x_2(pmin, t) == LARGER || 
                 compare_x_2(t, pmax) == LARGER || 
                 compare_y_2(pmin, t) == LARGER || 
                 compare_y_2(t, pmax) == LARGER   ) {
          std::cout << "debug Side_of_bs return ON_UNBOUNDED_SIDE" << std::endl;
          return ON_UNBOUNDED_SIDE;
        }
        else {
          std::cout << "debug Side_of_bs return ON_BOUNDARY" << std::endl;
          return ON_BOUNDARY;
        }
      }

    public:
      
      Bounded_side operator()(const Point_2 &p, const Point_2 &q,
                   const Point_2 &r, const Point_2 &t) const
      {
        return predicate(p, q, r, t);
      }
    };
    
  
} //namespace CGAL

#endif // CGAL_SIDE_OF_BOUNDED_SQUARE_2_H
