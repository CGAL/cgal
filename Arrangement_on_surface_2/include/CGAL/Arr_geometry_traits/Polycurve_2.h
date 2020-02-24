// Copyright (c) 2009,2010,2011,2015 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ron Wein  <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//                 Dror Atariah <dror.atariah@fu-berlin.de>

#ifndef CGAL_POLYCURVE_2_H
#define CGAL_POLYCURVE_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Header file for the polyline classes used by the
 * Arr_polycurve_basic_traits_2, Arr_polycurve_traits_2, and
 * Arr_polyline_traits_2 classes.
 */

#include <list>
#include <vector>
#include <iterator>
#include <limits>

#include <CGAL/Bbox_2.h>

namespace CGAL {

namespace internal {

/*! \class
 * Representation of a polycurve.
 */
template <typename SubcurveType_2, typename PointType_2>
class Polycurve_2 {
public:
  typedef SubcurveType_2                                Subcurve_type_2;
  typedef PointType_2                                   Point_type_2;

protected:
  // The subcurves that comprise the poyline:
  typedef typename std::vector<Subcurve_type_2>         Subcurves_container;

  Subcurves_container m_subcurves;

public:
  typedef typename Subcurves_container::size_type       Size;
  typedef typename Subcurves_container::size_type       size_type;

  /*! Construct default. */
  Polycurve_2() : m_subcurves() {}

  /*! Construct from a subcurve. */
  Polycurve_2(const Subcurve_type_2& seg) : m_subcurves()
  { m_subcurves.push_back(seg); }

  /* Roadmap:
   * - Currently construction from points is marked as deprecated.
   *   For the sake of backwards compatibility, it is not removed yet.
   * - In the next release, construction from points will be removed.
   * - Then, the remaining construction from subcurves will become private
   *   in this class as construction is provided by the traits.
   *   Furthermore, the polycurve traits class will become a friend of
   *   this class.
   */

  /*! Construct from a range. The range can be either:
   *
   * For the sake of backwards compatibility we have to keep the possibility
   * of construction from a range of points. Therefore, we have to test the
   * type of the elements of the range again.
   * Since ultimately the construction from points will be deprecated, the
   * actual constructor that we implemented which handles a range of point
   * is already deprecated and SHOULD NOT BE USED.
   *
   * If you want to construct a polycurve from a range of points, use the
   * construction functors from the traits class.
   *
   * - Range of points, and the polycurve is defined by the order of the
   *   points.
   * - Range of linear object. The polycurve is the sequence of linear
   *   objects.
   * \param begin An iterator pointing to the first point in the range.
   * \param end An iterator pointing after the last point in the range.
   * \pre Depends on the range's content. See the implementations for
   *      further details.
   */
  template <typename InputIterator>
  Polycurve_2(InputIterator begin, InputIterator end) :
    m_subcurves()
  {
    typedef typename std::iterator_traits<InputIterator>::value_type VT;
    typedef typename boost::is_same<VT, Point_type_2>::type Is_point;
    construct_polycurve(begin, end, Is_point());
  }

  /*! Construct a polycurve from a range of subcurves.
   * \param begin An iterator pointing to the first subcurve in the range.
   * \param end An iterator pointing after the past-the-end subcurve
   *        in the range.
   * \pre The end of subcurve n should be the beginning of subcurve n+1.
   */
  template <typename InputIterator>
  void construct_polycurve(InputIterator begin, InputIterator end,
                          boost::false_type)
  { m_subcurves.assign(begin, end); }

  /*! Construct a polycurve from a range of points.
   * \param begin An iterator pointing to the first point in the range.
   * \param end An iterator pointing after the last point in the range.
   * \pre There are at least 2 points in the range.
   *      In other cases, an empty polycurve will be created.
   */
  template <typename InputIterator>
  CGAL_DEPRECATED void construct_polycurve(InputIterator begin,
                                          InputIterator end,
                                          boost::true_type)
  {
    // Check if there are no points in the range:
    InputIterator ps = begin;
    if (ps == end) return;

    InputIterator pt = ps;
    ++pt;

    // The range contains only one point. A degenerated polycurve is
    // constructed.
    // With one degenerated subcurve, where source=target.
    if (pt == end) {
      m_subcurves.push_back(Subcurve_type_2(*ps, *ps));
      return;
    }

    // Construct a subcurve from each to adjacent points.
    // The container has to contain at least two points.
    while (pt != end) {
      m_subcurves.push_back(Subcurve_type_2(*ps, *pt));
      ++ps;
      ++pt;
    }
  }

  /* Roadmap: Make this private in the next version (after the traits
   *          becomes friendly...)
   */

  /*! Append a subcurve to the (x-monotone) polycurve.
   * Warning: This is a risky function! Don't use it! Prefer the
   *          provided implementation in the traits class.
   * \param seg The new subcurve to be appended to the polycurve.
   * \pre If the polycurve is not empty, seg source must be the
   *      same as the target point of the last subcurve in the polycurve.
   */
  inline void push_back(const Subcurve_type_2& seg)
  { this->m_subcurves.push_back(seg); }

  /*! Append a subcurve to the (x-monotone) polycurve.
   * Warning: This is a risky function! Don't use it! Prefer the
   *          provided implementation in the traits class.
   * \param seg The new subcurve to be appended to the polycurve.
   * \pre If the polycurve is not empty, seg source must be the
   *      same as the target point of the last subcurve in the polycurve.
   */
  inline void push_front(const Subcurve_type_2& seg)
  { this->m_subcurves.insert(this->m_subcurves.begin(), seg); }

  /*! Append a point to the polycurve.
   * To properly implemented this function the traits class is needed,
   * thus it is deprecated.
   */
  CGAL_DEPRECATED void push_back(const Point_type_2& p)
  {
    CGAL_assertion(!m_subcurves.empty());
    Point_type_2 pt = p;
    Point_type_2 ps = m_subcurves.back().target();
    m_subcurves.push_back(Subcurve_type_2(ps, pt));
  }

  /*! Compute the bounding box of the polycurve.
   * \return The bounding-box.
   */
  Bbox_2 bbox() const
  {
    // Compute the union of the bounding boxes of all subcurves.
    size_type n = this->number_of_subcurves();
    Bbox_2 bbox;
    for (std::size_t i = 0; i < n; ++i)
      bbox = (i > 0) ? (bbox + (*this)[i].bbox()) : (*this)[i].bbox();
    return bbox;
  }

  class Point_const_iterator;
  friend class Point_const_iterator;
  typedef std::reverse_iterator<Point_const_iterator>
    Point_const_reverse_iterator;

  /*! An iterator for the polycurve points. */
  class Point_const_iterator {
  public:
    // Type definitions:
    typedef std::bidirectional_iterator_tag     iterator_category;
    typedef Point_type_2                        value_type;
    typedef std::ptrdiff_t                      difference_type;
    typedef size_t                              size_type;
    typedef const value_type&                   reference;
    typedef const value_type*                   pointer;

  private:
    // The polycurve curve.
    const Polycurve_2<Subcurve_type_2, Point_type_2>* m_cvP;
    size_type m_num_pts;                        // Its number of points.
    size_type m_index;                          // The current point.

    /*! Private constructor.
     * \param cv The scanned curve.
     * \param index The index of the subcurve.
     */
    Point_const_iterator(const Polycurve_2<Subcurve_type_2, Point_type_2>* cvP,
                         size_type index) :
      m_cvP(cvP),
      m_index(index)
    {
      m_num_pts = (m_cvP == nullptr) ? 0 :
        ((m_cvP->number_of_subcurves() == 0) ?
         0 : (m_cvP->number_of_subcurves() + 1));
    }

    bool is_index_valid() const
    { return m_index != std::numeric_limits<size_type>::max BOOST_PREVENT_MACRO_SUBSTITUTION (); }

  public:
    /*! Default constructor. */
    Point_const_iterator() :
      m_cvP(nullptr),
      m_num_pts(0),
      m_index(std::numeric_limits<size_type>::max BOOST_PREVENT_MACRO_SUBSTITUTION ())
    {}

    /*! Dereference operator.
     * \return The current point.
     */
    const Point_type_2& operator*() const
    {
      CGAL_assertion(m_cvP != nullptr);
      CGAL_assertion((is_index_valid()) && (m_index < m_num_pts));

      // First point is the source of the first subcurve.
      // Else return the target of the(i-1)'st subcurve.
      if (m_index == 0) return ((*m_cvP)[0]).source();
      else return ((*m_cvP)[m_index - 1]).target();
    }

    /*! Arrow operator.
     * \return A pointer to the current point.
     */
    const Point_type_2* operator->() const { return (&(this->operator* ())); }

    /*! Increment operators. */
    Point_const_iterator& operator++()
    {
      if ((m_cvP != nullptr) && (m_index < m_num_pts)) ++m_index;
      return (*this);
    }

    Point_const_iterator operator++(int)
    {
      Point_const_iterator temp = *this;
      if ((m_cvP != nullptr) && (m_index < m_num_pts)) ++m_index;
      return temp;
    }

    /*! Decrement operators. */
    Point_const_iterator& operator--()
    {
      if ((m_cvP != nullptr) && (is_index_valid())) --m_index;
      return (*this);
    }

    Point_const_iterator operator--(int)
    {
      Point_const_iterator temp = *this;
      if ((m_cvP != nullptr) && (is_index_valid())) --m_index;
      return temp;
    }

    /*! Equality operators. */
    bool operator==(const Point_const_iterator& it) const
    { return ((m_cvP == it.m_cvP) && (m_index == it.m_index)); }

    bool operator!=(const Point_const_iterator& it) const
    { return ((m_cvP != it.m_cvP) || (m_index != it.m_index)); }

    friend class Polycurve_2<Subcurve_type_2, Point_type_2>;
  };

  // Backward compatibility
  typedef Point_const_iterator                  const_iterator;
  typedef Point_const_reverse_iterator          const_reverse_iterator;

  /*! Obtain an iterator for the polycurve points.*/
  Point_const_iterator points_begin() const
  {
    if (number_of_subcurves() == 0) return (Point_const_iterator());
    else return (Point_const_iterator(this, 0));
  }

  /*! Obtain an iterator for the polycurve points.*/
  CGAL_DEPRECATED Point_const_iterator begin() const { return points_begin(); }

  /*! Obtain a past-the-end iterator for the polycurve points.*/
  Point_const_iterator points_end() const
  {
    if (number_of_subcurves() == 0) return (Point_const_iterator());
    else return (Point_const_iterator(this, number_of_subcurves() + 1));
  }

  /*! Obtain a past-the-end iterator for the polycurve points.*/
  CGAL_DEPRECATED Point_const_iterator end() const { return points_end(); }

  /*! Obtain a reverse iterator for the polycurve points. */
  Point_const_reverse_iterator points_rbegin() const
  { return Point_const_reverse_iterator(points_end()); }

  /*! Obtain a reverse iterator for the polycurve points. */
  CGAL_DEPRECATED Point_const_reverse_iterator rbegin() const
  { return points_rbegin(); }

  /*! Obtain a reverse past-the-end iterator for the polycurve points. */
  Point_const_reverse_iterator points_rend() const
  { return Point_const_reverse_iterator(points_begin()); }

  /*! Obtain a reverse past-the-end iterator for the polycurve points. */
  CGAL_DEPRECATED Point_const_reverse_iterator rend() const
  { return points_rend(); }

  typedef typename Subcurves_container::iterator Subcurve_iterator;
  typedef typename Subcurves_container::const_iterator
                                                 Subcurve_const_iterator;
  typedef typename std::reverse_iterator<Subcurve_const_iterator>
                                                 Subcurve_const_reverse_iterator;

  /*! Obtain an iterator for the polycurve subcurves. */
  Subcurve_const_iterator subcurves_begin() const
  { return m_subcurves.begin(); }

  /*! Deprecated! */
  CGAL_DEPRECATED Subcurve_const_iterator begin_subcurves() const
  { return subcurves_begin(); }

  /*! Obtain a past-the-end iterator for the polycurve subcurves. */
  Subcurve_const_iterator subcurves_end() const
  { return m_subcurves.end(); }

  /*! Deprecated! */
  CGAL_DEPRECATED Subcurve_const_iterator end_subcurves() const
  { return subcurves_end(); }

  /*! Obtain a reverse iterator for the polycurve subcurves. */
  Subcurve_const_reverse_iterator subcurves_rbegin() const
  { return (Subcurve_const_reverse_iterator(subcurves_end())); }

  /*! Obtain a reverse past-the-end iterator for the polycurve points. */
  Subcurve_const_reverse_iterator subcurves_rend() const
  { return (Subcurve_const_reverse_iterator(subcurves_begin())); }

  /*! Deprecated!
   * Obtain the number of points contained in the polycurve.
   * In general (for example if the polycurve is not bounded), then the
   * number of vertices cannot be read-off from the number of subcurves, and
   * the traits class is needed.
   * \return The number of points.
   */
  CGAL_DEPRECATED std::size_t points() const
  { return (number_of_subcurves() == 0) ? 0 : number_of_subcurves() + 1; }

  /*! Obtain the number of subcurves that comprise the poyline.
   * \return The number of subcurves.
   */
  size_type number_of_subcurves() const
  { return m_subcurves.size(); }

  /*! Deprecated! Replaced by number_of_subcurves(). */
  CGAL_DEPRECATED size_type size() const { return number_of_subcurves(); }

  /*! Obtain the ith subcurve of the polycurve.
   * \param[in] i The subcurve index(from 0 to size()-1).
   * \return A const reference to the subcurve.
   */
  inline const Subcurve_type_2& operator[](const std::size_t i) const
  {
    CGAL_assertion(i < number_of_subcurves());
    return (m_subcurves[i]);
  }

  /*! Clear the polycurve. */
  inline void clear() { m_subcurves.clear(); }
};

/*! \class
 * Representation of an x-monotone polycurve.
 * An x-monotone polycurve is always directed from left to right.
 */
template <typename SubcurveType_2, typename PointType_2>
class X_monotone_polycurve_2 :
    public Polycurve_2<SubcurveType_2, PointType_2>
{
public:
  typedef SubcurveType_2                             Subcurve_type_2;
  typedef PointType_2                                Point_type_2;

  typedef Polycurve_2<Subcurve_type_2, Point_type_2> Base;

  /*! Construct default. */
  X_monotone_polycurve_2() : Base() {}

  /*! Construct from a subcurve. */
  X_monotone_polycurve_2(Subcurve_type_2 seg) : Base(seg) {}

  /*! Construct from a range.
   * Similar to the constructor of a general polycurve.
   * Like in the case of general polycurve, for the sake of backwards
   * compatibility we have to keep an implementation of construction
   * from a range of points. DO NOT USE THIS CONSTRUCTION.
   */
  template <typename InputIterator>
  X_monotone_polycurve_2(InputIterator begin, InputIterator end) :
    Base(begin, end)
  {
    typedef typename std::iterator_traits<InputIterator>::value_type VT;
    typedef typename boost::is_same<VT, Point_type_2>::type Is_point;
    construct_x_monotone_polycurve(begin, end, Is_point());
  }

  /*! Construct from a range of subcurves.
   * This constructor is expected to be called only from the
   * traits class, after the input was verified there.
   * \pre The range of subcurves form an x-monotone polycurve.
   */
  template <typename InputIterator>
  void construct_x_monotone_polycurve(InputIterator, InputIterator,
                                      boost::false_type)
  {}

  /*! Construct from a range of points, defining the endpoints of the
   * polycurve subcurves.
   */
  template <typename InputIterator>
  CGAL_DEPRECATED void construct_x_monotone_polycurve(InputIterator,
                                                     InputIterator,
                                                     boost::true_type)
  {}
};

} // namespace polycurve
} //namespace CGAL

#endif
