// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
#ifndef CGAL_CURVE_DATA_TRAITS_2_H
#define CGAL_CURVE_DATA_TRAITS_2_H

CGAL_BEGIN_NAMESPACE

/*! A generic traits class for maintaining an arrangement of curves that have
 * an extra data field. This traits class is templated with a Data class an
 * an ordinary traits class which is also used as a based traits class to
 * inherit from. It extracts the original Curve_2 and X_monotone_curve_2 types
 * from the ordinary traits class, and redefines them to have Data as an extra
 * field.
 *
 * The Data field is updated when the curves are converted from Curve_2 to
 * X_monotone_curve_2, and when the X_monotone_curve_2 curves are split. All
 * other methods are inherited from the base ordinary traits class.
 */

template <class Traits_, class Data_>
class Arr_curve_data_traits_2 : public Traits_ {
public:
  typedef Traits_                               Traits;
  typedef Data_                                 Data;
  typedef typename Traits::Curve_2              Org_curve_2;
  typedef typename Traits::X_monotone_curve_2   Org_x_monotone_curve_2;
  typedef typename Traits::Point_2              Point_2;

  
  class Curve_2 : public Org_curve_2 {
  private:
    Data m_data;

  public:
    Curve_2() {}
    
    Curve_2(const Org_curve_2 & cv, const Data & data) :
      Org_curve_2(cv),
      m_data(data)
    {}

    const Data & get_data() const { return m_data; }

    void set_data(const Data & data) { m_data = data; }
  };
  
  class X_monotone_curve_2 : public Org_x_monotone_curve_2 {
  private:
    Data m_data;

  public:

    X_monotone_curve_2() {}
    
    X_monotone_curve_2(const Org_x_monotone_curve_2 & cv, const Data & data) :
      Org_x_monotone_curve_2(cv),
      m_data(data)
    {}

    const Data & get_data() const { return m_data; }

    void set_data(const Data & data) { m_data = data; }
  };

  // For backward compatibility:
  typedef Point_2                               Point;
  typedef X_monotone_curve_2                    X_curve;
  typedef Curve_2                               Curve;
  
public:
  Arr_curve_data_traits_2 () : Traits() {}
  
  /*! Cut the given curve into x-monotone subcurves and insert them to the
   * given output iterator. While segments are x_monotone, still need to pass
   * them out.
   * \param cv The curve.
   * \param o The output iterator
   * \return The past-the-end iterator
   *
   * why this function dose not compile as const???
   */
  template<class OutputIterator>
  OutputIterator curve_make_x_monotone(const Curve_2 & cv,
                                       OutputIterator o)// const
  {
    std::list<Org_x_monotone_curve_2>  org_x_curves;
    
    Traits::curve_make_x_monotone (cv, std::back_inserter(org_x_curves));

    typename std::list<Org_x_monotone_curve_2>::const_iterator it;
    for (it = org_x_curves.begin(); it != org_x_curves.end(); it++) {
      *o++ = X_monotone_curve_2 (*it, cv.get_data());
    }
    
    return o;
  } 

  /*! Split a given curve at a given split point into two sub-curves.
   * \param cv the curve to split
   * \param c1 the output first part of the split curve. Its source is the
   * source of the original curve.
   * \param c2 the output second part of the split curve. Its target is the
   * target of the original curve.
   * \param p the split point.
   * \pre p lies on cv but is not one of its end-points.
   */
  void curve_split(const X_monotone_curve_2 & cv, 
		   X_monotone_curve_2 & c1, X_monotone_curve_2 & c2, 
                   const Point_2 & p) const
  {
    Org_x_monotone_curve_2 org_c1, org_c2;
    Traits::curve_split (cv, org_c1, org_c2, p);
    c1 = X_monotone_curve_2 (org_c1, cv.get_data());
    c2 = X_monotone_curve_2 (org_c2, cv.get_data());
  }

  /*!
   * Find the nearest intersection of the two given curves to the right of 
   * a given reference point.
   * Nearest is defined as the lexicographically nearest point, not including 
   * the point reference point itself.
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
   * there is an overlapping subcurve, that contains the reference point in
   * its x-range, the function should return an X_monotone_curve_2 whose 
   * interior is strictly to the right of the reference point (that is, whose
   * left endpoint is the projection of the reference point onto the 
   * overlapping subcurve).
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \return An empty object if there is no intersection to the right of p.
   *         An object wrapping a Point_2 in case of a simple intersection.
   *         An object wrapping an X_monotone_curve_2 in case of an overlap.
   */
  Object nearest_intersection_to_right (const X_monotone_curve_2 & cv1,
				        const X_monotone_curve_2 & cv2,
					const Point_2 & p) const
  {
    // Use the base method and resturn the result as is in case there is
    // no intersection or if the intersection is a single point.
    CGAL::Object res = Traits::nearest_intersection_to_right (cv1, cv2, p);
    Point_2      ip;

    if (res.is_empty() || CGAL::assign (ip, res))
      return (res);

    // Attach data to the overlapping curve and return it.
    Org_x_monotone_curve_2  icv;

    CGAL::assign (icv, res);
    return (CGAL::make_object (X_monotone_curve_2 (icv, cv1.get_data())));
  }

  /*!
   * Find the nearest intersection of the two given curves to the left of 
   * a given reference point.
   * Nearest is defined as the lexicographically nearest point, not including 
   * the point reference point itself.
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
   * there is an overlapping subcurve, that contains the reference point in
   * its x-range, the function should return an X_monotone_curve_2 whose 
   * interior is strictly to the left of the reference point (that is, whose
   * right endpoint is the projection of the reference point onto the 
   * overlapping subcurve).
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \return An empty object if there is no intersection to the left of p.
   *         An object wrapping a Point_2 in case of a simple intersection.
   *         An object wrapping an X_monotone_curve_2 in case of an overlap.
   */
  Object nearest_intersection_to_left (const X_monotone_curve_2 & cv1,
				       const X_monotone_curve_2 & cv2,
				       const Point_2 & p) const
  {
    // Use the base method and resturn the result as is in case there is
    // no intersection or if the intersection is a single point.
    CGAL::Object res = Traits::nearest_intersection_to_left (cv1, cv2, p);
    Point_2      ip;

    if (res.is_empty() || CGAL::assign (ip, res))
      return (res);

    // Attach data to the overlapping curve and return it.
    Org_x_monotone_curve_2  icv;

    CGAL::assign (icv, res);
    return (CGAL::make_object (X_monotone_curve_2 (icv, cv1.get_data())));
  }
};

CGAL_END_NAMESPACE

#endif
