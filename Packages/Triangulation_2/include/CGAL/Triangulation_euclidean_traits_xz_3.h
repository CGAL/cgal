// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-12 $
// release_date  : $CGAL_Date: 1999/04/28 $
//
// file          : include/CGAL/Triangulation_euclidean_traits_xz_3.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================


#ifndef CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_XZ_3_H
#define CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_XZ_3_H

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Distance_2.h>


#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template < class R >
class Triangulation_euclidean_traits_xz_3 {
public:
    typedef Triangulation_euclidean_traits_xz_3<R> Traits;
    typedef R Rp;
    typedef Point_3<R>  Point;
    typedef Segment_3<R> Segment;
    typedef Triangle_3<R> Triangle;
    typedef Line_3<R>   Line;
    typedef Ray_3<R>    Ray;
    typedef Direction_3<R> Direction;
    
    
  Triangulation_euclidean_traits_xz_3(){}
  Triangulation_euclidean_traits_xz_3(
				const Triangulation_euclidean_traits_xz_3& et){}
  Triangulation_euclidean_traits_xz_3 &operator=(
		       const Triangulation_euclidean_traits_xz_3&  et){return *this;}

  typename Rp::FT x(const Point &p) const { return p.x(); }
  typename Rp::FT y(const Point &p) const { return p.z(); }

    Comparison_result compare_x(const Point &p, const Point &q) const
      {
        return CGAL::compare(x(p), x(q));
      }
    Comparison_result compare_y(const Point &p, const Point &q) const
      {
        return CGAL::compare(y(p), y(q));
      }
    bool compare(const Point &p, const Point &q) const
      {
        return (compare_x(p, q)== EQUAL &&  
		compare_y(p, q)== EQUAL);
      }
    
    Orientation orientation(const Point &p,
                                 const Point &q,
                                 const Point &r) const
      {
        return orientationC2(x(p), y(p), x(q), y(q), x(r), y(r));
      }
    

    Oriented_side side_of_oriented_circle(const Point &p,
					  const Point &q,
					  const Point &r,
					  const Point &s) const
    {
      return side_of_oriented_circleC2(x(p), y(p),
				       x(q), y(q),
				       x(r), y(r),
				       x(s), y(s));
    }
        
    
    class Distance : public Distance_2<Traits> 
    {
      typedef typename Distance_2<Traits>::Point Point;

    public:
        Distance(const Point& p0,
                 const Traits* traits = NULL)
            : Distance_2<Traits>(p0, traits) { }
    
    
        Distance(const Point& p0,
                 const Point& p1,
                 const Traits* traits = NULL)
            : Distance_2<Traits>(p0,p1,traits) { }
    
        Distance(const Point& p0,
                 const Point& p1,
                 const Point& p2,
                 const Traits* traits = NULL)
            : Distance_2<Traits>(p0,p1,p2,traits) { }
    
        Comparison_result
        compare() const
        {
          Point p0 = get_point(0);
          Point p1 = get_point(1);
          Point p2 = get_point(2);
          return cmp_dist_to_pointC2(x(p0),y(p0),x(p1),y(p1),x(p2),y(p2));
        }
    };
    
};

CGAL_END_NAMESPACE


#endif // CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_XZ_3_H
