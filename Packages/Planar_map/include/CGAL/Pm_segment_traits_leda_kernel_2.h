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
// release       : $$
// release_date  : $$
//
// file          : include/CGAL/Pm_segment_traits_leda_kernel.h
// package       : Planar_map (5.87)
// maintainer    : Efi Fogel         <efif@post.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Shai Hirsch       <shaihi@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PM_SEGMENT_TRAITS_LEDA_KERNEL
#define CGAL_PM_SEGMENT_TRAITS_LEDA_KERNEL

#include <CGAL/rat_leda_in_CGAL_2.h>
#include <CGAL/Planar_map_2/Pm_segment_utilities_2.h>

CGAL_BEGIN_NAMESPACE

/*!
 */
template< class FT_ >
class my_rat_direction : public leda_rat_point {
public:
public:
  typedef FT_                                   RT;
  typedef FT_                                   FT;

  my_rat_direction(const leda_rat_point & p) : leda_rat_point(p) {}
  my_rat_direction(const FT x, const FT y) : leda_rat_point(x,y) {}

  /*!
   */
  Sign sign_of_determinant2x2(const FT & a00, const FT & a01,
                              const FT & a10, const FT & a11) const
  {
    return
      static_cast<Sign>(static_cast<int>(CGAL_NTS compare( a00*a11, a10*a01)));
  }
    
  /*!
   */
  Comparison_result
  compare_angle_with_x_axis_2(const FT & dx1, const FT & dy1,
                              const FT & dx2, const FT & dy2) const
  {
    // angles are in [-pi,pi], and the angle between Ox and d1 is compared
    // with the angle between Ox and d2
    int quadrant_1 = (dx1 >= FT(0)) ? ((dy1 >= FT(0))?1:4)
      : ((dy1 >= FT(0))?2:3);
    int quadrant_2 = (dx2 >= FT(0)) ? ((dy2 >= FT(0))?1:4)
      : ((dy2 >= FT(0))?2:3);
    // We can't use CGAL_NTS compare(quadrant_1,quadrant_2) because in case
    // of tie, we need additional computation
    if (quadrant_1 > quadrant_2) return LARGER;
    else if (quadrant_1 < quadrant_2) return SMALLER;
    return Comparison_result(-sign_of_determinant2x2(dx1,dy1,dx2,dy2));
  }
      
  /*!
   */
  Comparison_result
  compare_angle_with_x_axis(const my_rat_direction & d1,
                            const my_rat_direction & d2) const
  { return compare_angle_with_x_axis_2(d1.xcoord(), d1.ycoord(),
                                       d2.xcoord(), d2.ycoord()); }

  /*!
   */
  bool operator<(const my_rat_direction & d) const
  { return compare_angle_with_x_axis(*this, d) == SMALLER; }

  /*!
   */
  bool operator<=(const my_rat_direction & d) const
  { return compare_angle_with_x_axis(*this, d) != LARGER; }
      
  /*!
   */
  bool counterclockwise_in_between(const my_rat_direction & d1,
                                   const my_rat_direction & d2) const
  {
    if (d1 <*this) return (*this < d2 || d2 <= d1);
    return (*this < d2 && d2 <= d1);
  }
};

template< class FT_ >
class Pm_segment_traits_leda_kernel_2 {
private:
  typedef Pm_segment_traits_leda_kernel_2<FT_>  Self;
    
public:
  typedef FT_                                   RT;
  typedef FT_                                   FT;
  typedef leda_rat_point                        Point_2;
  typedef leda_rat_segment                      Segment_2;
  typedef my_rat_direction<FT>                  Direction_2;

  /*! Functor
   */
  class Is_vertical_2 {
  public:
    /*! \todo Can cv be a point, or should it be prohobited by a precondition.
     */
    bool operator()(const Segment_2 & cv) const
    {
      if (cv.is_trivial()) return true;
      return cv.is_vertical();
    }
  };

  /*! Functor
   */
  class Equal_2 {
  public:
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    { return (p1 == p2); }

    bool operator()(const Segment_2 & c1, const Segment_2 & c2) const
    { return (c1 == c2); }
  };

  /*! Functor
   */
  class Compare_x_2 {
  public:
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      int res = Self::Point_2::cmp_x(p1, p2);
      return ((res < 0) ? SMALLER : ((res > 0) ? LARGER : EQUAL));
    }
  };

  /*! Functor
   */
  class Compare_y_2 {
  public:
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      int res = Self::Point_2::cmp_y(p1, p2);
      return ((res < 0) ? SMALLER : ((res > 0) ? LARGER : EQUAL));
    }
  };

  /*! Functor
   */
  class Less_x_2 {
  public:
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    { return (Self::Point_2::cmp_x(p1, p2) < 0); }
  };
    
  /*! Functor
   */
  class Less_y_2 {
  public:
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    { return (Self::Point_2::cmp_y(p1, p2) < 0); }
  };
    
  /*! Functor
   */
  class Construct_vertex_2 {
  public:
    Point_2 operator()(const Segment_2 & cv, int id) const
    { return ((id == 0) ? cv.source() : cv.target()); }
  };

  /*! Functor
   */
  class Compare_y_at_x_2 {
  private:
    int cmp_x(const Point_2 & p1, const Point_2 & p2) const
    { return Self::Point_2::cmp_x(p1, p2); }

    int cmp_y(const Point_2 & p1, const Point_2 & p2) const
    { return Self::Point_2::cmp_y(p1, p2); }
      
    bool curve_is_in_x_range(const Segment_2 & cv, const Point_2 & p) const
    {
      return
        !(((cmp_x(p, cv.source()) < 0) && (cmp_x(p, cv.target()) < 0)) ||
         ((cmp_x(p, cv.source()) > 0) && (cmp_x(p, cv.target()) > 0)));
    }
	
    bool curve_is_in_y_range(const Segment_2 & cv, const Point_2 & p) const
    { 
      return
        !(((cmp_y(p, cv.source()) < 0) && (cmp_y(p, cv.target()) < 0)) ||
         ((cmp_y(p, cv.source()) > 0) && (cmp_y(p, cv.target()) > 0)));
    }

    Orientation orientation(const Point_2 & p, const Point_2 & q,
                            const Point_2 & r) const
    { return CGAL::orientation(p, q, r); }
	
  public:
    Comparison_result operator()(const Point_2 & q,
                                 const Segment_2 & cv1, const Segment_2 & cv2)
      const
    {
      if ((!curve_is_in_x_range(cv1, q)) || (!curve_is_in_x_range(cv2, q)))
        return EQUAL;
		
      // bug ??? in LEDA - 
      // cmp_segments_at_xcoord returns wrong answer if
      // cv1 (or cv2) are from right to left
      // cv1_ and cv2_ are the same as cv1 and cv2 - 
      //   oriented from left to right
      Segment_2 cv1_ = cv1;
      Segment_2 cv2_ = cv2;
      if (lexicographically_xy_larger(cv1.source(), cv1.target()))
        cv1_ = cv1.reversal();
      if (lexicographically_xy_larger(cv2.source(), cv2.target()))
        cv2_ = cv2.reversal();
  		
      // checking verical curves.
      if (cv1_.is_vertical()) {
        if (cv2_.is_vertical()) {
          // both cv1 and cv2 are vertical
          int res = cmp_y(cv1_.target(), cv2_.source());
          return ((res < 0) ? SMALLER : ((res > 0) ? LARGER : EQUAL));
        }

        // only cv1 is vertical.
        if (orientation(cv2_.source(), cv2_.target(), cv1_.source()) > 0)
          return LARGER;
                      
        if (orientation(cv2_.source(), cv2_.target(), cv1_.target()) < 0)
          return SMALLER;
  
        return EQUAL;
      }
                  
      if (cv2_.is_vertical()) {
        // only cv2 is vertical:
        if (orientation(cv1_.source(), cv1_.target(), cv2_.source()) > 0 )
          return SMALLER;
                      
        if (orientation(cv1_.source(), cv1_.target(), cv2_.target()) < 0)
          return LARGER;
  
        return EQUAL;  
      }
                    
      // Non of the curves are vertical:
      int res = cmp_segments_at_xcoord(cv1_, cv2_, q);
      return ((res < 0) ? SMALLER : ((res > 0) ? LARGER : EQUAL));
    }

    Comparison_result operator()(const Point_2 & p, const Segment_2 & cv)
    {
      if (!curve_is_in_x_range(cv, p)) return EQUAL;

      if (cv.is_vertical()) {
        if ((cmp_y(p, cv.source()) < 0) && (cmp_y(p, cv.target()) < 0))
          return SMALLER;
        if ((cmp_y(p, cv.source()) > 0) && (cmp_y(p, cv.target()) > 0))
          return LARGER;
        return EQUAL;
      }

      Orientation o = (cmp_x(cv.source(), cv.target()) < 0) ?
        orientation(cv.source(), cv.target(), p) :
        orientation(cv.target(), cv.source(), p);
			
      return ((o < 0) ? SMALLER : ((o > 0) ? LARGER : EQUAL));
    }
  };

  /*!
   */
  class Compare_slope_2 {
  public:
    Comparison_result operator()(const Segment_2 & cv1, const Segment_2 & cv2)
      const
    {
      int res = cmp_slopes(cv1, cv2);
      return ((res < 0) ? SMALLER : ((res > 0) ? LARGER : EQUAL));
    }
  };

  /*!
   */
  class Counterclockwise_in_between_2 {
  public:
    bool operator()(const Direction_2 & d,
                    const Direction_2 & d1, const Direction_2 & d2) const
    {
      return d.counterclockwise_in_between(d1, d2);
    }
  };

  /*!
   */
  class Construct_direction_2 {
  public:
    Direction_2 operator()(const Segment_2 & cv) const
    { return Direction_2(cv.target() - cv.source()); }
  };

  /*!
   */
  class Construct_opposite_direction_2 {
  public:
    Direction_2 operator()(const Direction_2 & d) const
    { return Direction_2(-d.xcoord(), -d.ycoord()); }
  };

  // creators:
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(); }
  Equal_2 equal_2_object() const { return Equal_2(); }
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(); }
  Compare_y_2 compare_y_2_object() const { return Compare_y_2(); }
  Construct_vertex_2 construct_vertex_2_object() const
    { return Construct_vertex_2(); }
  Compare_y_at_x_2 compare_y_at_x_2_object() const
    { return Compare_y_at_x_2(); }
  Less_x_2 less_x_2_object() const { return Less_x_2(); }
  Compare_slope_2 compare_slope_2_object() const { return Compare_slope_2(); }
  Counterclockwise_in_between_2 counterclockwise_in_between_2_object() const
    { return Counterclockwise_in_between_2(); }
  Construct_direction_2 construct_direction_2_object() const
    { return Construct_direction_2(); }
  Construct_opposite_direction_2 construct_opposite_direction_2_object() const
    { return Construct_opposite_direction_2(); }
};

CGAL_END_NAMESPACE

#endif
