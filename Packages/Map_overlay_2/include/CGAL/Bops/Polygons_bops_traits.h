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
// release       : $CGAL_Revision: CGAL-2.5-I-11 $
// release_date  : $CGAL_Date: 2002/08/04 $
//
// file          : include/CGAL/Polygons_bops_traits.h
// package       : Map_overlay (1.12)
// maintainer    : Efi Fogel <efif@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra          <estere@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_POLYGONS_BOPS_TRAITS_2_H
#define CGAL_POLYGONS_BOPS_TRAITS_2_H

#include <CGAL/Polygon_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Ray_2_Segment_2_intersection.h>

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Polygons_bops_traits_2 : public Arr_segment_traits_2<Kernel_>
{
  typedef Kernel_                       Kernel;
  typedef typename Kernel::FT           FT;

public:
  typedef Arr_segment_traits_2<Kernel>  Base;
  
  typedef typename Base::Point_2        Point_2;
  typedef typename Base::X_curve_2      X_curve_2;
  typedef X_curve_2                     Curve_2;
  
  typedef typename Kernel::Ray_2        Ray_2;   // a new type for this class.
  typedef typename Kernel::Direction_2  Direction_2;
  
  typedef typename Base::Curve_point_status Curve_point_status;
  
  // Obsolete, for backward compatibility
  typedef Point_2                       Point;
  typedef X_curve_2                     X_curve;
  typedef Curve_2                       Curve;

  //protected:
  //typedef typename Base::Curve_status   Curve_status;
  
public:
  
  Polygons_bops_traits_2() : Base() 
  {}

  // Calculate the intersection point of a vertical ray enamating from p and
  // cv. We consider the ray as a line since, it is given that the ray
  // hits cv, hence we don't have to calculate the whether p is
  // above/below cv.
  Point_2 calc_hitting_point(const X_curve_2& cv, 
                             const Point_2& p)
  {
    Curve_point_status status = curve_get_point_status(cv,p); 
    
    // It is given that p is on cv range.
    CGAL_assertion(status != CURVE_NOT_IN_RANGE);
    
    if (status == ON_CURVE)
          return p;

    if (curve_is_vertical(cv))
      {
        // Else below or above curve.
        //    |            *
        //    |     OR     |
        //    *            | 
        
        if (status == UNDER_CURVE)
          return lowest(cv.source(),cv.target());
        else // status == ABOVE_CURVE
          return highest(cv.source(),cv.target());
      }
    
    Direction_2 dir = (status == UNDER_CURVE)? 
      Direction_2(FT(0),FT(1)) : Direction_2(FT(0),FT(-1));
    Ray_2 ray(p,dir);
    
    //double alpha = (status == UNDER_CURVE)? 90 : -90;
    //Ray_2 ray(p,alpha);
    Point_2 hitting_p;
    
    CGAL_assertion(do_intersect(cv,ray));
    
    Object obj=CGAL::intersection(cv, ray);
    assign(hitting_p, obj);

    return hitting_p;
  }
private:
  bool is_left(const Point_2 &p1, const Point_2 &p2) const 
  { return (compare_x(p1, p2) == SMALLER); }
  bool is_right(const Point_2 &p1, const Point_2 &p2) const 
  { return (compare_x(p1, p2) == LARGER); }
  bool is_same_x(const Point_2 &p1, const Point_2 &p2) const 
  { return (compare_x(p1, p2) == EQUAL); }
  bool is_lower(const Point_2 &p1, const Point_2 &p2) const 
  { return (compare_y(p1, p2) == SMALLER); }
  bool is_higher(const Point_2 &p1, const Point_2 &p2) const 
  { return (compare_y(p1, p2) == LARGER); }
  bool is_same_y(const Point_2 &p1, const Point_2 &p2) const 
  { return (compare_y(p1, p2) == EQUAL); }
  bool is_same(const Point_2 &p1, const Point_2 &p2) const
  {
    return (compare_x(p1, p2) == EQUAL) &&
      (compare_y(p1, p2) == EQUAL);
  }
  const Point_2& leftmost(const Point_2 &p1, const Point_2 &p2) const
  { return (is_left(p1, p2) ? p1 : p2); }

  const Point_2& rightmost(const Point_2 &p1, const Point_2 &p2) const
  { return (is_right(p1, p2) ? p1 : p2); }
  
  const Point_2& lowest(const Point_2 &p1, const Point_2 &p2) const
  { return (is_lower(p1, p2) ? p1 : p2); }
  
  const Point_2& highest(const Point_2 &p1, const Point_2 &p2) const
  { return (is_higher(p1, p2) ? p1 : p2); }
};

CGAL_END_NAMESPACE

#endif




