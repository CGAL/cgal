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
// file          : include/CGAL/Arr_debug_traits_2.h
// package       : Arrangement (5.87)
// maintainer    : Efi Fogel         <efif@post.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Ron Wein          <wein@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <danha@post.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_ARR_DEBUG_TRAITS_2_H
#define CGAL_ARR_DEBUG_TRAITS_2_H

#include <CGAL/tags.h>
 
CGAL_BEGIN_NAMESPACE

/*!
 * A traits class for debugging a generic traits class.
 */
template <class Base_traits_>
class Arr_debug_traits_2
{
private:

  Base_traits_                 base;
  mutable std::ostream&        os;

public:

  typedef Base_traits_                             Base_traits;
  typedef typename Base_traits::Has_left_category  Has_left_category;
  typedef CGAL::Tag_false                          Has_reflect_category;
  typedef typename Base_traits::Point_2            Point_2;
  typedef typename Base_traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Base_traits::Curve_2            Curve_2;

  // For backward compatability:    
  typedef Point_2                         Point;
  typedef X_monotone_curve_2              X_curve;
  typedef X_monotone_curve_2              Curve;

  /*!
   * Default constructor.
   */
  Arr_debug_traits_2() :
    base(),
    os(std::cout)
  {}
  
  /*!
   * Constructor with an output file.
   */
  Arr_debug_traits_2 (std::ostream& _os) :
    base(),
    os(_os)
  {}
    
  /*!
   * Compare the x-coordinates of two given points.
   * \param p1 The first point.
   * \param p2 The second point.
   * \return LARGER if x(p1) > x(p2); SMALLER if x(p1) < x(p2); or else EQUAL.
   */
  Comparison_result compare_x (const Point_2& p1, 
			       const Point_2& p2) const
  {
    Comparison_result res = base.compare_x (p1, p2);

    os << "compare_x (" << p1 << "   ," 
       << '\t' << p2 << ")" << std::endl;
    os << "-> " << res << std::endl;

    return (res);
  }

  /*! 
   * Compares lexigoraphically the two points: by x, then by y.
   * \param p1 Te first point.
   * \param p2 The second point.
   * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2); 
   *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
   *         or else EQUAL.
   */
  Comparison_result compare_xy (const Point_2& p1, 
				const Point_2& p2) const
  {
    Comparison_result res = base.compare_xy (p1, p2);

    os << "compare_xy (" << p1 << "   ," 
       << '\t' << p2 << ")" << std::endl;
    os << "-> " << res << std::endl;

    return (res);
  }

  /*!
   * Check whether the given curve is a vertical segment.
   * \param cv The curve.
   * \return (true) if the curve is vertical.
   */
  bool curve_is_vertical(const X_monotone_curve_2& cv) const 
  {
    bool        res = base.curve_is_vertical (cv);

    os << "curve_is_vertical (" << cv << std::endl;
    os << "-> " << res << std::endl;

    return (res);
  } 

  /*!
   * Check whether the given point is in the x-range of the given curve.
   * \param cv The curve.
   * \param q The point.
   * \return (true) if q is in the x-range of cv.
   */
  bool point_in_x_range (const X_monotone_curve_2& cv, 
			 const Point_2& q) const
  {
    bool        res = base.point_in_x_range (cv, q);

    os << "point_in_x_range (" << cv << "   ,"
       << '\t' << q << ")" << std::endl;
    os << "-> " << res << std::endl;

    return (res);
  }

  /*!
   * Get the relative status of two curves at the x-coordinate of a given 
   * point.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param q The point.
   * \pre The point q is in the x-range of the two curves.
   * \return LARGER if cv1(x(q)) > cv2(x(q)); SMALLER if cv1(x(q)) < cv2(x(q));
   *  or else EQUAL.
   */
  Comparison_result curves_compare_y_at_x (const X_monotone_curve_2& cv1, 
					   const X_monotone_curve_2& cv2, 
					   const Point_2& q) const
  {
    Comparison_result res = base.curves_compare_y_at_x (cv1, cv2, q);

    os << "curves_compare_y_at_x (" << cv1 << "   ,"
       << '\t' << cv2 << "   ,"
       << '\t' << q << ")" << std::endl;
    os << "-> " << res << std::endl;

    return (res);
  }

  /*!
   * Compares the y value of two curves immediately to the left of the 
   * x-value of their intersection point.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param q The point.
   * \pre The point q is in the x range of the two curves, and both of them 
   * must be also be defined to its left. Furthermore, cv1(x(q) == cv2(x(q)).
   * \return The relative position of cv1 with respect to cv2 to the left of
   * x(q): LARGER, SMALLER or EQUAL.
   */
  Comparison_result curves_compare_y_at_x_left (const X_monotone_curve_2& cv1,
						const X_monotone_curve_2& cv2, 
						const Point_2& q) const 
  {
    Comparison_result res = base.curves_compare_y_at_x_left (cv1, cv2, q);

    os << "curves_compare_y_at_x_left (" << cv1 << "   ,"
       << '\t' << cv2 << "   ,"
       << '\t' << q << ")" << std::endl;
    os << "-> " << res << std::endl;

    return (res);
  }

  /*!
   * Compares the y value of two curves in an epsilon environment to the right
   * of the x-value of their intersection point.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param q The point.
   * \pre The point q is in the x range of the two curves, and both of them 
   * must be also be defined to its right. Furthermore, cv1(x(q) == cv2(x(q)).
   * \return The relative position of cv1 with respect to cv2 to the right of
   * x(q): LARGER, SMALLER or EQUAL.
   */
  Comparison_result curves_compare_y_at_x_right (const X_monotone_curve_2& cv1,
						 const X_monotone_curve_2& cv2,
						 const Point_2& q) const
  {
    Comparison_result res = base.curves_compare_y_at_x_right (cv1, cv2, q);

    os << "curves_compare_y_at_x_right (" << cv1 << "   ,"
       << '\t' << cv2 << "   ,"
       << '\t' << q << ")" << std::endl;
    os << "-> " << res << std::endl;

    return (res);
  }
    
  /*! 
   * Return the location of the given point with respect to the input curve.
   * \param cv The curve.
   * \param p The point.
   * \pre p is in the x-range of cv.
   * \return SMALLER if y(p) < cv(x(p)), that is the point is below the curve;
   *         LARGER if y(p) > cv(x(p)), that is the point is above the curve;
   *         or else (if p is on the curve) EQUAL.
   */
  Comparison_result curve_compare_y_at_x (const Point_2& p,
					  const X_monotone_curve_2& cv) const
  {
    Comparison_result res = base.curve_compare_y_at_x (p, cv);

    os << "curve_compare_y_at_x (" << p << "   ,"
       << '\t' << cv << ")" << std::endl;
    os << "-> " << res << std::endl;

    return (res);
  }

  /*! 
   * Check if the two curves are the same (have the same graph).
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \return (true) if the two curves are the same.
   */
  bool curve_equal (const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
  {
    bool        res = base.curve_equal (cv1, cv2);

    os << "curve_equal (" << cv1 << "   ,"
       << '\t' << cv2 << ")" << std::endl;
    os << "-> " << res << std::endl;

    return (res);
  }

  /*!
   * Check if the two points are the same.
   * \param p1 The first point.
   * \param p2 The second point.
   * \return (true) if p1 == p2.
   */
  bool point_equal (const Point_2& p1,
		    const Point_2& p2) const
  {
    bool        res = base.point_equal (p1, p2);

    os << "point_equal (" << p1 << "   ,"
       << '\t' << p2 << ")" << std::endl;
    os << "-> " << res << std::endl;

    return (res);
  }
  
  /*!
   * Get the curve source.
   * \param cv The curve.
   * \return The source point.
   */
  Point_2 curve_source (const X_monotone_curve_2& cv) const 
  { 
    Point_2        p = base.curve_source (cv);

    os << "curve_source (" << cv << ")" << std::endl;
    os << "-> " << p << std::endl;

    return (p);
  }

  /*!
   * Get the curve target.
   * \param cv The curve.
   * \return The target point.
   */
  Point_2 curve_target(const X_monotone_curve_2& cv) const 
  { 
    Point_2        p = base.curve_target (cv);

    os << "curve_target (" << cv << ")" << std::endl;
    os << "-> " << p << std::endl;

    return (p);
  }

  /*! 
   * Cut the given curve into x-monotone subcurves and insert them to the
   * given output iterator.
   * \param cv The curve.
   * \param o The output iterator
   * \return The past-the-end iterator
   */
  template<class OutputIterator>
  OutputIterator curve_make_x_monotone (const Curve_2& cv,
                                        OutputIterator oi) const
  {
    std::list<X_monotone_curve_2>                          xcvs;
    typename std::list<X_monotone_curve_2>::const_iterator it;

    base.curve_make_x_monotone (cv, std::back_inserter(xcvs));
    
    os << "curve_make_x_monotone (" << cv << ")" << std::endl;
    for (it = xcvs.begin(); it != xcvs.end(); ++it)
    {
      os << (it == xcvs.begin() ? "-> " : "   ") << *it << std::endl;
      *oi++ = *it;
    }
    
    return (oi);
  } 

  /*!
   * Split a given curve at a given split point into two sub-curves.
   * \param cv The curve to split
   * \param c1 The output first part of the split curve. 
   *           Its source will be the source of the original curve.
   * \param c2 The output second part of the split curve.
   *           Its target will be the target of the original curve.
   * \param p The split point.
   * \pre p lies on cv but is not one of its endpoints.
   */
  void curve_split (const X_monotone_curve_2& cv, 
		    X_monotone_curve_2& c1, X_monotone_curve_2& c2, 
		    const Point_2& p) const
  {
    base.curve_split (cv, c1, c2, p);

    os << "curve_split (" << cv << "   ,"
       << '\t' << p << ")" << std::endl;
    os << "-> " << c1 << std::endl
       << "   " << c2 << std::endl;

    return;
  }

  /*!
   * Find the nearest intersection point (or points) of two given curves to
   * the right lexicographically of a given reference point.
   */
  CGAL::Object nearest_intersection_to_right (const X_monotone_curve_2& cv1,
					      const X_monotone_curve_2& cv2,
					      const Point_2& p) const
  {
    CGAL::Object       res = base.nearest_intersection_to_right (cv1, cv2, p);
    Point_2            ip;
    X_monotone_curve_2 icv;

    os << "nearest_intersection_to_right (" << cv1 << "   ,"
       << '\t' << cv2 << "   ,"
       << '\t' << p << ")" << std::endl;
    
    if (res.is_empty())
    {
      os << "-> NONE!" << std::endl;
    }
    else if (CGAL::assign (ip, res))
    {
      os << "-> " << ip << std::endl;
    }
    else if (CGAL::assign (icv, res))
    {
      os << "-> " << icv << std::endl;
    }
    else
    {
      os << "-> ERROR!" << std::endl;
    }

    return (res);
  }

  /*!
   * Find the nearest intersection point (or points) of two given curves to
   * the left lexicographically of a given reference point.
   */
  CGAL::Object nearest_intersection_to_left (const X_monotone_curve_2& cv1,
					     const X_monotone_curve_2& cv2,
					     const Point_2& p) const
  {
    CGAL::Object       res = base.nearest_intersection_to_left (cv1, cv2, p);
    Point_2            ip;
    X_monotone_curve_2 icv;

    os << "nearest_intersection_to_left (" << cv1 << "   ,"
       << '\t' << cv2 << "   ,"
       << '\t' << p << ")" << std::endl;
    
    if (res.is_empty())
    {
      os << "-> NONE!" << std::endl;
    }
    else if (CGAL::assign (ip, res))
    {
      os << "-> " << ip << std::endl;
    }
    else if (CGAL::assign (icv, res))
    {
      os << "-> " << icv << std::endl;
    }
    else
    {
      os << "-> ERROR!" << std::endl;
    }

    return (res);
  }

  /*!
   * Check whether the two given curves overlap.
   * \patam cv1 The first curve.
   * \patam cv2 The second curve.
   * \return (true) if the two curves overlap in a one-dimensional subcurve
   * (i.e., not in a finite number of points). Otherwise, if they have a finite
   * number of intersection points, or none at all, return (false).
   */
  bool curves_overlap (const X_monotone_curve_2& cv1,
                       const X_monotone_curve_2& cv2) const
  {
    bool       res = base.curves_overlap (cv1, cv2);

    os << "curves_overlap (" << cv1 << "   ,"
       << '\t' << cv2 << std::endl;
    os << "-> " << res << std::endl;

    return (res);
  }

};

CGAL_END_NAMESPACE

#endif
