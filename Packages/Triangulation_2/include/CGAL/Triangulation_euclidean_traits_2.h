#ifndef CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
#define CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Distance_2.h>

template < class R >
class CGAL_Triangulation_euclidean_traits_2 {
public:
  typedef R Rep;
  typedef CGAL_Point_2<R>  Point;
  typedef CGAL_Segment_2<R> Segment;
  typedef CGAL_Triangle_2<R> Triangle;

  typedef CGAL_Distance_2<CGAL_Triangulation_euclidean_traits_2<R> > Distance;

    CGAL_Comparison_result compare_x(const Point &p, const Point &q) const
    {
        return CGAL_compare_x(p, q);
    }


    CGAL_Comparison_result compare_y(const Point &p, const Point &q) const
    {
        return CGAL_compare_y(p, q);
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
        if (p==q) return CGAL_COLLINEAR;
        if (p==r) return CGAL_COLLINEAR;
        if (r==q) return CGAL_COLLINEAR;

        return CGAL_orientation(p, q, r);
    }

    CGAL_Oriented_side side_of_oriented_circle(const Point &p,
                                               const Point &q,
                                               const Point &r,
                                               const Point &s) const
    {
        if (p==s) return CGAL_ON_ORIENTED_BOUNDARY;
        if (q==s) return CGAL_ON_ORIENTED_BOUNDARY;
        if (r==s) return CGAL_ON_ORIENTED_BOUNDARY;

        return CGAL_side_of_oriented_circle(p, q, r, s);
    }
};


#endif // CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
