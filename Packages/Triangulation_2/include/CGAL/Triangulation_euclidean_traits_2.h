#ifndef CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
#define CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Ray_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Distance_2.h>

template < class R >
class CGAL_Triangulation_euclidean_traits_2 {
public:
  typedef R Rep;
  typedef CGAL_Point_2<R>  Point;
  typedef CGAL_Segment_2<R> Segment;
  typedef CGAL_Triangle_2<R> Triangle;
  typedef CGAL_Line_2<R> Line;
  typedef CGAL_Direction_2<R> Direction;
  typedef CGAL_Ray_2<R> Ray;



  typedef CGAL_Distance_2<CGAL_Triangulation_euclidean_traits_2<R> > Distance;

    CGAL_Comparison_result compare_x(const Point &p, const Point &q) const
    {
        return CGAL_compare_x(p, q);
    }


    CGAL_Comparison_result compare_y(const Point &p, const Point &q) const
    {
        return CGAL_compare_y(p, q);
    }

  bool compare(const Point &p, const Point &q) const
  {
    return (CGAL_compare_x(p, q)== CGAL_EQUAL &&  
	    (CGAL_compare_x(p, q)== CGAL_EQUAL);
  }

    CGAL_Orientation orientation(const Point &p,
                                 const Point &q,
                                 const Point &r) const
    {
        return CGAL_orientation(p, q, r);
    }


   CGAL_Orientation extremal(const Point &p,
                              const Point &q,
                              const Point &r) const
      {
        if (compare(p,q)) return CGAL_COLLINEAR;
        if (compare(p,r)) return CGAL_COLLINEAR;
        if (compare(r,q)) return CGAL_COLLINEAR;
    
        return CGAL_orientation(p, q, r);
      }
    
    CGAL_Oriented_side side_of_oriented_circle(const Point &p,
                                               const Point &q,
                                               const Point &r,
                                               const Point &s) const
      {
        if (compare(p,s)) return CGAL_ON_ORIENTED_BOUNDARY;
        if (compare(q,s)) return CGAL_ON_ORIENTED_BOUNDARY;
        if (compare(r,s)) return CGAL_ON_ORIENTED_BOUNDARY;
    
        return CGAL_side_of_oriented_circle(p, q, r, s);
      }

     Point circumcenter(const Point &p, const Point &q, const Point &r) const
    {
        return CGAL_circumcenter(p, q, r);
    }

};


#endif // CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
