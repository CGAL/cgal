#ifndef CGAL_PM_SEGMENT_TRAITS_LEDA_KERNEL
#define CGAL_PM_SEGMENT_TRAITS_LEDA_KERNEL

#include <CGAL/rat_leda_in_CGAL_2.h>
#include <CGAL/Planar_map_2/Pm_segment_utilities_2.h>

CGAL_BEGIN_NAMESPACE

template< class FT_ >
class Pm_segment_traits_leda_kernel_2 {
private:
  typedef Pm_segment_traits_leda_kernel_2       Self;
    
public:
  typedef FT_                                   RT;
  typedef FT_                                   FT;
  typedef leda_rat_point                        Point_2;
  typedef leda_rat_segment                      Segment_2;

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
    
  // creators:
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(); }
  Equal_2 equal_2_object() const { return Equal_2(); }
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(); }
  Compare_y_2 compare_y_2_object() const { return Compare_y_2(); }
  Construct_vertex_2 construct_vertex_2_object() const { return Construct_vertex_2(); }
  Compare_y_at_x_2 compare_y_at_x_2_object() const { return Compare_y_at_x_2(); }
  Less_x_2 less_x_2_object() const { return Less_x_2(); }
  Compare_slope_2 compare_slope_2_object() const { return Compare_slope_2(); }
};

CGAL_END_NAMESPACE

#endif
