// Copyright (c) 2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POLYCURVE_2_H
#define CGAL_ARR_POLYCURVE_2_H

/*! \file
 * Header file for the polycurve classes used by the Arr_polycurve_traits_2
 * class.
 */

#include <list>
#include <iostream>
#include <vector>
#include <iterator>

#include <CGAL/Bbox_2.h>

namespace CGAL {


/*! \class
 * Representation of a polycurve.
 */
template <class CurveTraits>
class Polycurve_2
{
public:
  typedef CurveTraits                                   Curve_traits_2;
  typedef typename Curve_traits_2::Point_2              Point_2;
  typedef typename Curve_traits_2::X_monotone_curve_2   Curve_2;
  typedef typename Curve_traits_2::Curve_2              Curve_2;

  typedef std::vector<Curve_2>::const_iterator          const_iterator;
  typedef std::vector<Curve_2>::const_reverse_iterator  const_reverse_iterator;
  
protected:
  // The segments that comprise the polycurve:
  std::vector<Curve_2> curves;

public:
  /*! Default constructor. */
  Polycurve_2() {}

  /*!
   * Constructor from a range of curves.
   * \param begin An iterator pointing to the first curve in the range.
   * \param end An iterator pointing after the last curve in the range.
   * \pre the target of the i-th curve is identical to the source of the
   * i+1-th curve for i = 0,1,...,n-1
   */
  template <class InputIterator>
  Polycurve_2(InputIterator begin, InputIterator end)
  {
    std::copy(begin, end, std::back_inserter(curves));
  }

  /*! Append a curve to the polycurve.
   * \param curve The new curve
   * \pre the target of the current last curve must be identical to the
   * source of curve
   */
  void push_back(const Curve & curve)
  {
    curves.push_back (curve);
  }

  /*!
   * Create a bounding-box for the polycurve.
   * \return The bounding-box.
   */
  Bbox_2 bbox() const
  {
    // Compute the union of the bounding boxes of all curves.
    unsigned int n = size();
    Bbox_2 bbox;
    unsigned int  i;

    const_iterator ci = begin();
    bbox = ci->bbox(*ci);
    for (++ci; ci != end(); ++ci)
      bbox = bbox + ci->bbox(*ci);
    return bbox;
  }

  /*! Obtain an iterator for the polycurve points. */
  const_iterator begin() const { return curves.begin(); }

  /*! Obtain a past-the-end iterator for the polycurve points. */
  const_iterator end() const { return curves.end(); }

  /*! Obtain an reverse iterator for the polycurve points. */
  const_reverse_iterator rbegin() const { return curves.rbegin(); }

  /*! Obtain a reverse past-the-end iterator for the polycurve points. */
  const_reverse_iterator rend() const { return return curves.rend(); }

  /*!
   * Obtain the number of curves that comprise the polycurve.
   * \return The number of curves.
   */
  inline unsigned int size() const { return curves.size(); }

  /*!
   * Obtain the i'th segment of the polycurve.
   * \param i The segment index(from 0 to size()-1).
   * \return A const reference to the segment.
   */
  inline const Curve_2 & operator[](const unsigned int i) const
  { return curves[i]; }

  /*! Clear the polycurve. */
  inline void clear() { curves.clear(); }
};

/*! \class
 * Representation of an x-monotone polycurve.
 * An x-monotone polycurve is always directed from left to right.
 */
template <class CurveTraits>
class X_monotone_polycurve_2 : public Polycurve_2<CurveTraits>
{
public:
  typedef CurveTraits                           Curve_traits_2;
  typedef Polycurve_2<Curve_traits_2>           Base;
  typedef typename Curve_traits_2::Point_2      Point_2;
  typedef typename Curve_traits_2::Curve_2      Curve_2;

  /*! Default constructor. */
  X_monotone_polycurve_2() : Base() {}

  /*!
   * Constructor from a range of points, defining the endpoints of the
   * polycurve curves.
   */
  template <class InputIterator>
  X_monotone_polycurve_2(InputIterator begin, InputIterator end) :
    Base(begin, end)
  {
    // Make sure the range of points contains at least two points.
    Curve_traits_2 seg_traits;
    InputIterator ps = begin;
    CGAL_precondition(ps != end);
    InputIterator pt = ps;
    ++pt;
    CGAL_precondition(pt != end);

    CGAL_precondition_code(
      typename Curve_traits_2::Compare_x_2 compare_x =
        seg_traits.compare_x_2_object();
      );
    CGAL_precondition_code(
      typename Curve_traits_2::Compare_xy_2 compare_xy =
        seg_traits.compare_xy_2_object();
      );
    
    
    // Make sure there is no change of directions as we traverse the polycurve.
    CGAL_precondition_code(
      const Comparison_result cmp_x_res = compare_x(*ps, *pt);
    );
    const Comparison_result cmp_xy_res = compare_xy(*ps, *pt);
    CGAL_precondition(cmp_xy_res != EQUAL);
    ++ps; ++pt;
    while (pt != end) {
      CGAL_precondition (compare_xy(*ps, *pt) == cmp_xy_res);
      CGAL_precondition (compare_x(*ps, *pt) == cmp_x_res);
      ++ps; ++pt;
    }

    // Reverse the polycurve so it always directed from left to right.
    if (cmp_xy_res == LARGER)
      _reverse();
  }

  /*!
   * Append a segment to the polycurve.
   * \param seg The new segment to be appended to the polycurve.
   * \pre If the polycurve is not empty, the segment source must be the
   *      same as the target point of the last segment in the polycurve
   *      (thus it must extend it to the right).
   */
  inline void push_back (const Curve_2& seg)
  {
    CGAL_precondition_code (Curve_traits_2   seg_tr);
    CGAL_precondition_code (const unsigned int n = this->size());
    CGAL_precondition (seg_tr.compare_xy_2_object() (seg.source(),
						     seg.target()) == SMALLER);
    CGAL_precondition (n == 0 ||
		       seg_tr.equal_2_object() (this->curves[n - 1].target(),
						seg.source()));

    this->curves.push_back (seg);
  }

private:

  /*! Reverse the polycurve. */
  void _reverse()
  {
    typename Base::const_reverse_iterator  ps = this->rbegin();
    typename Base::const_reverse_iterator  pt = ps;
    ++pt;

    std::vector<Curve_2>  rev_segs (this->size());
    unsigned int            i = 0;

    while (pt != this->rend())
    {
      rev_segs[i] = Curve_2 (*ps, *pt);
      ++ps; ++pt;
      i++;
    }

    this->curves = rev_segs;
    return;
  }

};

/*! Output operator for a polycurve. */
template <class CurveTraits>
std::ostream& operator<< (std::ostream & os,
			  const Polycurve_2<CurveTraits>& cv)
{
  typename Polycurve_2<CurveTraits>::const_iterator  iter = cv.begin();

  // Print the number of points:
  os << cv.points();

  while (iter != cv.end())
  {
    os << "  " << *iter;
    ++iter;
  }
  return (os);
}

/*! Input operator for a polycurve. */
template <class CurveTraits>
std::istream & operator>> (std::istream & is,
                           Polycurve_2<CurveTraits> & pl)
{
  std::cerr << "Not implemented yet!" << std::endl;
  return is;
}


} //namespace CGAL

#endif
