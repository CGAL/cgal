#ifdef CGAL_PM_SEGMENT_TRAITS_LEDA_TRAITS
#define CGAL_PM_SEGMENT_TRAITS_LEDA_TRAITS

#include <CGAL/rat_leda_in_CGAL_2.h>
#include <CGAL/Planar_map_2/Pm_segment_utilities_2.h>
#include <CGAL/Planar_map_2/Pm_point_utilities_2.h>

CGAL_BEGIN_NAMESPACE

class Pm_segment_traits_leda_kernel_2
{
private:
  typedef Pm_segment_traits_leda_kernel_2 Self

public:
  typedef leda_rat_point      Point_2
  typedef leda_rat_segment    Segment_2

  /*!
   */
  struct Is_vertical_2_ {
    /*! \todo Can cv be a point, or should it be prohobited by a precondition.
     */
    bool operator()(const X_curve_2 & cv) const
    {
      if (cv.is_trivial())
        return true;
      return cv.is_vertical();
    }
  };

  /*!
   */
  struct Equal_2_ {
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    { return is_same(p1, p2); }

    bool operator()(const X_curve_2 & c1, const X_curve_2 & c2) const
    { return is_same(c1, c2); }
  };

  /*!
   */
  Comparison_result Comparison_result_from_int(int res) const 
  {
    if (res < 0) return SMALLER;
    if (res > 0) return LARGER;
    return EQUAL;
  }
    
  /*!
   */
  struct Compare_x_2_ {
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      return Comparison_result_from_int(
#if (__LEDA__ >= 380)
        Self::Point_2::cmp_x(p1,p2)
#else // backward compatability to LEDA
        compare(p1.xcoord(),p2.xcoord())
#endif
      );
    }
  };

  /*!
   */
  struct Compare_y_2_ {
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      return Comparison_result_from_int(

#if (__LEDA__ >= 380)
        Self::Point_2::cmp_y(p1,p2)
#else // backward compatability to LEDA   
        compare(p1.ycoord(),p2.ycoord()) 
#endif
      );
    }
  };

  /*!
   */
  struct Construct_vertex_2_ {
    Point_2 operator()() const
    {
    }
  };

  /*!
   */
  struct Compare_y_at_x_2_ {
  private:
    bool curve_is_in_x_range(const X_curve_2 & cv, const Point_2 & q) const
    { 
      return !( is_right(q, rightmost(cv.source(), cv.target())) ||
                is_left(q, leftmost(cv.source(), cv.target()))	 );
    }
	
    bool curve_is_in_y_range(const X_curve_2 &cv, const Point_2 & q) const
    { 
      bool r = !( is_lower(q, lowest(cv.source(), cv.target())) ||
                  is_higher(q, highest(cv.source(), cv.target())) );
      return r;
    }

    Orientation orientation(const Point_2 & p, const Point_2 & q,
                            const Point_2 & r) const
    {
      return CGAL::orientation(p, q, r);
    }
	
  public:
    Comparison_result operator()(const X_curve_2 & cv1, const X_curve_2 & cv2, 
                                 const Point_2 & q) const
    {
      if ((!curve_is_in_x_range(cv1, q)) || 
          (!curve_is_in_x_range(cv2, q)))
        return EQUAL;
		
      int res;
      // bug ??? in LEDA - 
      // cmp_segments_at_xcoord returns wrong answer if
      // cv1 (or cv2) are from right to left
      // cv1_ and cv2_ are the same as cv1 and cv2 - 
      //   oriented from left to right
      X_curve_2 cv1_ = cv1;
      X_curve_2 cv2_ = cv2;
      if (lexicographically_xy_larger(cv1.source(), cv1.target()))
        cv1_ = cv1.reversal();
      if (lexicographically_xy_larger(cv2.source(), cv2.target()))
        cv2_ = cv2.reversal();
  		
                  
      // checking verical curves.
      struct Is_vertical_2_ is_vertical;
  
      if (is_vertical(cv1_)) {
                    
        if (is_vertical(cv2_)) {
          // both cv1 and cv2 are vertical
          if ( is_lower(cv1_.target(), cv2_.source()) )
            return SMALLER;
          if ( is_higher(cv1_.source(), cv2_.target()) )
            return LARGER;
          return EQUAL; // overlapping. 
        } // end  both cv1 and cv2 are vertical.
        else { // only cv1 is vertical.
          if (orientation(cv2_.source(), 
                          cv2_.target(), 
                          cv1_.source()) > 0 )
            return LARGER;
                      
          if (orientation(cv2_.source(), 
                          cv2_.target(), 
                          cv1_.target()) < 0)
            return SMALLER;
  
          return EQUAL;
        }
      }
                  
      if (is_vertical(cv2_)) { // only cv2 is vertical.
        if (orientation(cv1_.source(), 
                        cv1_.target(), 
                        cv2_.source()) > 0 )
          return SMALLER;
                      
        if (orientation(cv1_.source(), 
                        cv1_.target(), 
                        cv2_.target()) < 0)
          return LARGER;
  
        return EQUAL;  
      }
                    
      // end checking verical curves.
  
      res = cmp_segments_at_xcoord(cv1_, cv2_, q);
  		
      if (res < 0) 
        return SMALLER;
      if (res > 0) 
        return LARGER;
      return EQUAL;
    }
  };

  typedef Is_vertical_2_        Is_vertical_2;
  typedef Equal_2_              Equal_2;
  typedef Compare_x_2_          Compare_x_2;
  typedef Compare_y_2_          Compare_y_2;
  typedef Construct_vertex_2_   Construct_vertex_2;
  typedef Compare_y_at_x_2_     Compare_y_at_x_2;

  Is_vertical_2 is_vertical_2_object() { return Is_vertical_2(); }
  Equal_2_object() { return Equal_2(); }
  Compare_x_2_object() { return Compare_x_2(); }
  Compare_y_2_object() { return Compare_y_2(); }
  Construct_vertex_2_object() { return Construct_vertex_2(); }
  Compare_y_at_x_2_object() { return Compare_y_at_x_2(); }
    
CGAL_END_NAMESPACE

#endif
