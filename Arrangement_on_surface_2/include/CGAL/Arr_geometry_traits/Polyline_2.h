// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein  <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POLYLINE_2_H
#define CGAL_ARR_POLYLINE_2_H

/*! \file
 * Header file for the polyline classes used by the Arr_polyline_traits_2
 * class.
 */

#include <list>
#include <iostream>
#include <vector>
#include <iterator>

#include <CGAL/Bbox_2.h>

namespace CGAL {


/*! \class
 * Representation of a polyline.
 */
template <typename SegmentTraits_>
class _Polyline_2
{
public:

  typedef SegmentTraits_                        Segment_traits_2;
  typedef typename Segment_traits_2::Point_2    Point_2;
  typedef typename Segment_traits_2::Curve_2    Segment_2;

protected:
  // The segments that comprise the poyline:
  typedef typename std::vector<Segment_2>        Segments_container;
  typedef typename Segments_container::size_type Seg_container_size;
  Segments_container                             m_segments;

public:
  /*! Default constructor. */
  _Polyline_2() : m_segments() {}

  /*!
   * Constructor from a range. The range can be either:
   * - Range of points, and the polyline is defined by the order of the points
   * - Range of linear object. The polyline is the sequence of linear objects.
   *   Note, that the end of object n should be the beginning of object n+1.
   *   See furtehr details in the corresponding implementation.
   * \param begin An iterator pointing to the first point in the range.
   * \param end An iterator pointing after the last point in the range.
   * \pre Depends on the range's content. See the implementations
   *      for furtehr details.
   */
  template <typename InputIterator>
  _Polyline_2(InputIterator begin, InputIterator end) :
    m_segments()
  {
    construct_polyline(begin, end, *begin);
  }

  /*!
   * Construct a polyline from a range of segments.
   * \param begin An iterator pointing to the first segment in the range.
   * \param end An iterator pointing after the past-the-end segment in the range.
   */
  template <typename InputIterator>
  void construct_polyline(InputIterator begin, InputIterator end, 
                          const Segment_2& /* */)
  {
    m_segments.assign(begin, end);
  }
  
  /*!
   * Construct a polyline from a range of points.
   * \param begin An iterator pointing to the first point in the range.
   * \param end An iterator pointing after the last point in the range.
   * \pre There are at least 2 points in the range.
   *      In other cases, an empty polyline will be created.
   */
  template <typename InputIterator>
  void construct_polyline(InputIterator begin, InputIterator end, 
                          const Point_2& /* */)
  {
    // Check if there are no points in the range:
    InputIterator  ps = begin;

    if (ps == end)
      return;

    // Construct a segment from each to adjacent points.
    InputIterator pt = ps;
    ++pt;

    while (pt != end) {
      m_segments.push_back(Segment_2(*ps, *pt));
      ++ps;
      ++pt;
    }
  }

  /*! Append a point to the polyline. */
  void push_back (const Point_2 & p)
  {
    Point_2 pt = p;
    Point_2 ps = m_segments.back().target();
    m_segments.push_back (Segment_2 (ps, pt));
  }

  /*!
   * Create a bounding-box for the polyline.
   * \return The bounding-box.
   */
  Bbox_2 bbox() const
  {
    // Compute the union of the bounding boxes of all segments.
    unsigned int  n = this->size();
    Bbox_2        bbox;
    unsigned int  i;
    
    for (i = 0; i < n; ++i)
    {
      if (i > 0)
	bbox = bbox +(*this)[i].bbox();
      else
	bbox = (*this)[i].bbox();
    }

    return (bbox);
  }

  class const_iterator;
  friend class const_iterator;

  /*! An iterator for the polyline points. */
  /*
   * TODO: How to iterate over points when the polyline
   * has an open end(s)?
   */
  class const_iterator
  {
  public:

    // Type definitions:
    typedef std::bidirectional_iterator_tag     iterator_category;
    typedef Point_2                             value_type;
    typedef std::ptrdiff_t                      difference_type;
    typedef size_t                              size_type;
    typedef const value_type&                   reference;
    typedef const value_type*                   pointer;

  private:
    
    const _Polyline_2<SegmentTraits_> * m_cvP;  // The polyline curve.
    int   m_num_pts;                            // Its number of points.
    int   m_index;                              // The current point.

    /*!
     * Private constructor.
     * \param cv The scanned curve.
     * \param index The index of the segment.
     */
    const_iterator (const _Polyline_2<SegmentTraits_>* cvP, int index) :
      m_cvP(cvP),
      m_index(index)
    {
      if (m_cvP == NULL)
	m_num_pts = 0;
      else
	m_num_pts =
	  (m_cvP->size() == 0) ? 0 : static_cast<int>(m_cvP->size() + 1);
    }

  public:

    /*! Default constructor. */
    const_iterator() :
      m_cvP(NULL),
      m_num_pts(0),
      m_index(-1)
    {}

    /*! 
     * Dereference operator.
     * \return The current point.
     */
    const Point_2& operator*() const
    {
      CGAL_assertion(m_cvP != NULL);
      CGAL_assertion(m_index >= 0 && m_index < m_num_pts);
      
      if (m_index == 0)
      {
        // First point is the source of the first segment.
        return ((*m_cvP)[0]).source();
      }
      else
      {
        // Return the target of the(i-1)'st segment.
	return ((*m_cvP)[m_index - 1]).target();
      }
    }

    /*! 
     * Arrow operator.
     * \return A pointer to the current point.
     */
    const Point_2* operator-> () const
    {
      return (&(this->operator* ()));
    }

    /*! Increment operators. */
    const_iterator& operator++() 
    {
      if (m_cvP != NULL && m_index < m_num_pts)
	++m_index;
      return (*this);
    }

    const_iterator operator++ (int)
    {
      const_iterator  temp = *this;
      if (m_cvP != NULL && m_index < m_num_pts)
	++m_index;
      return (temp);
    }

    /*! Decrement operators. */
    const_iterator& operator-- ()
    {
      if (m_cvP != NULL && m_index >= 0)
	--m_index;
      return (*this);
    }

    const_iterator operator--(int)
    {
      const_iterator  temp = *this;
      if (m_cvP != NULL && m_index >= 0)
	--m_index;
      return (temp);
    }

    /*! Equality operators. */
    bool operator==(const const_iterator& it) const
    {
      return (m_cvP == it.m_cvP && m_index == it.m_index);
    }

    bool operator!=(const const_iterator& it) const
    {
      return (m_cvP != it.m_cvP || m_index != it.m_index);
    }

    friend class _Polyline_2<SegmentTraits_>;
  };

  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  /*! Get an iterator for the polyline points. */
  const_iterator begin() const
  {
    if (size() == 0)
      return (const_iterator (NULL, -1));
    else
      return (const_iterator (this, 0));
  }

  /*! Get a past-the-end iterator for the polyline points. */
  const_iterator end() const
  {
    if (size() == 0)
      return (const_iterator (NULL, -1));
    else
      return (const_iterator (this, size() + 1));
  }

  /*! Get a reverse iterator for the polyline points. */
  const_reverse_iterator rbegin() const
  {
    return (const_reverse_iterator (end()));
  }

  /*! Get a reverse past-the-end iterator for the polyline points. */
  const_reverse_iterator rend() const
  {
    return (const_reverse_iterator (begin()));
  }

  class const_segments_iterator;
  friend class const_segments_iterator;

  /*! An iterator for the polyline points. */
  class const_segments_iterator
  {
  public:

    // Type definitions:
    typedef std::bidirectional_iterator_tag     iterator_category;
    typedef Segment_2                           value_type;
    typedef std::ptrdiff_t                      difference_type;
    // TODO: Shouldn't size_t be replaced by Seg_container_size?
    typedef size_t                              size_type;
    typedef const value_type&                   reference;
    typedef const value_type*                   pointer;

  private:
    
    const _Polyline_2<SegmentTraits_>* m_cvP;  // The polyline curve.
    // TODO: Use size() directly.
    int   m_num_segs;                           // Its number of segments
    int   m_index;                              // The current segment.

    /*!
     * Private constructor.
     * \param cv The scanned curve.
     * \param index The index of the segment.
     */
    const_segments_iterator (const _Polyline_2<SegmentTraits_>* cvP, 
			     int index) :
      m_cvP(cvP),
      m_index(index)
    {
      if (m_cvP == NULL)
	{
	  m_num_segs = 0;
	}
      else
	{
	  m_num_segs = m_cvP->size();
	}
    }

  public:

    /*! Default constructor. */
    const_segments_iterator() :
      m_cvP(NULL),
      m_num_segs(0),
      m_index(-1)
    { }

    /*! 
     * Dereference operator.
     * \return The current point.
     */

    const Segment_2& operator*() const
    {
      CGAL_assertion(m_cvP != NULL);
      CGAL_assertion(m_index >= 0 && m_index < m_num_segs);
      return (*m_cvP)[m_index];
    }

    /*! 
     * Arrow operator.
     * \return A pointer to the current sgement.
     */
    const Segment_2* operator-> () const
    {
      return (&(this->operator* ()));
    }

    /*! Increment operators. */
    const_segments_iterator& operator++() 
    {
      CGAL_assertion(m_index >= 0 && m_index < m_num_segs);
      if (m_cvP != NULL && m_index < m_num_segs)
	++m_index;
      return (*this);
    }

    const_segments_iterator operator++ (int)
    {
      const_segments_iterator  temp = *this;
      if (m_cvP != NULL && m_index < m_num_segs)
	++m_index;
      return (temp);
    }

    /*! Decrement operators. */
    const_segments_iterator& operator-- ()
    {
      // TODO: Shouldn't m_index > 0??
      if (m_cvP != NULL && m_index >= 0)
	--m_index;
      return (*this);
    }

    const_segments_iterator operator--(int)
    {
      const_iterator  temp = *this;
      // TODO: Shouldn't m_index > 0??
      if (m_cvP != NULL && m_index >= 0)
	--m_index;
      return (temp);
    }

    /*! Equality operators. */
    bool operator==(const const_segments_iterator& it) const
    {
      return (m_cvP == it.m_cvP && m_index == it.m_index);
    }

    bool operator!=(const const_segments_iterator& it) const
    {
      return (m_cvP != it.m_cvP || m_index != it.m_index);
    }

    friend class _Polyline_2<SegmentTraits_>;
  };

  typedef std::reverse_iterator<const_segments_iterator>
    const_reverse_segments_iterator;

/*! Get an iterator for the polyline's segments. */
const_segments_iterator begin_segments() const
  {
    if (size() == 0)
      return (const_segments_iterator (NULL, -1));
    else
      return (const_segments_iterator (this, 0));
  }

  /*! Get a past-the-end iterator for the polyline's segments. */
  const_segments_iterator end_segments() const
  {
    if (size() == 0)
      return (const_segments_iterator (NULL, -1));
    else
      //      return (const_segments_iterator (this, size() + 1));
      return (const_segments_iterator (this, size()));
  }

  /*! Get a reverse iterator for the polyline's segments. */
  const_reverse_segments_iterator rbegin_segments() const
  {
    return (const_reverse_segments_iterator (end_segments()));
  }

  /*! Get a reverse past-the-end iterator for the polyline points. */
  const_reverse_segments_iterator rend_segments() const
  {
    return (const_reverse_segments_iterator (begin_segments()));
  }

  /*!
   * Get the number of points contained in the polyline.
   * TODO: Has to be updated once the polyline will be
   * unbounded.
   * \return The number of points.
   */
  unsigned int points() const
  {
    return (size() == 0) ? 0 : size() + 1;
  }

  /*!
   * Get the number of segments that comprise the poyline.
   * \return The number of segments.
   */
  Seg_container_size size() const
  {
    return Seg_container_size(m_segments.size());
  }

  /*!
   * Get the ith segment of the polyline.
   * \param i The segment index(from 0 to size()-1).
   * \return A const reference to the segment.
   */
  inline const Segment_2& operator[] (const unsigned int i) const
  {
    CGAL_assertion (i < size());
    return (m_segments[i]);
  }

  /*! Clear the polyline. */
  inline void clear ()
  {
    m_segments.clear();
  }
};

/*! \class
 * Representation of an x-monotone polyline.
 * An x-monotone polyline is always directed from left to right.
 */
template <typename SegmentTraits_>
class _X_monotone_polyline_2 : public _Polyline_2<SegmentTraits_>
{
public:
  typedef SegmentTraits_                        Segment_traits_2;
  typedef _Polyline_2<SegmentTraits_>           Base;
  typedef typename Segment_traits_2::Point_2    Point_2;
  typedef typename Segment_traits_2::Curve_2    Segment_2;

  /*! Default constructor. */
  _X_monotone_polyline_2() : Base () {}

  /*! Constructor */
  template <typename InputIterator>
  _X_monotone_polyline_2(InputIterator begin, InputIterator end) :
    Base(begin, end)
  {
    construct_x_monotone_polyline(begin, end, *begin);
  }
  
  /*!
   * Constructs from a range of segments.
   */
  template <typename InputIterator>
  void construct_x_monotone_polyline(InputIterator begin, InputIterator end,
                                     const Segment_2& /* */)
  {
    for (InputIterator it = begin ; it != end; ++it )
      {
	std::cout << "\033[1;32mP1: Seg is: " << *it << "\033[0m\n";
      }
  }
  
  /*!
   * Constructs from a range of points, defining the endpoints of the
   * polyline segments.
   */
  template <typename InputIterator>
  void construct_x_monotone_polyline(InputIterator begin, InputIterator end,
                                     const Point_2& /* */)
  {
    // Make sure the range of points contains at least two points.
    Segment_traits_2 seg_traits;
    InputIterator ps = begin;
    CGAL_precondition (ps != end);
    InputIterator pt = ps;
    ++pt;
    CGAL_precondition (pt != end);

    CGAL_precondition_code(
      typename Segment_traits_2::Compare_x_2 compare_x =
        seg_traits.compare_x_2_object();
      );
    CGAL_precondition_code(
      typename Segment_traits_2::Compare_xy_2 compare_xy =
        seg_traits.compare_xy_2_object();
      );
    
    
    // Make sure there is no change of directions as we traverse the polyline.
    CGAL_precondition_code (
      const Comparison_result cmp_x_res = compare_x(*ps, *pt);
    );
    const Comparison_result cmp_xy_res = compare_xy(*ps, *pt);
    CGAL_precondition (cmp_xy_res != EQUAL);
    ++ps; ++pt;
    while (pt != end) {
      CGAL_precondition (compare_xy(*ps, *pt) == cmp_xy_res);
      CGAL_precondition (compare_x(*ps, *pt) == cmp_x_res);
      ++ps; ++pt;
    }

    // Reverse the polyline so it always directed from left to right.
    if (cmp_xy_res == LARGER)
      _reverse();
  }

  /*!
   * Append a segment to the polyline.
   * \param seg The new segment to be appended to the polyline.
   * \pre If the polyline is not empty, the segment source must be the
   *      same as the target point of the last segment in the polyline
   *      (thus it must extend it to the right).
   */
  inline void push_back (const Segment_2& seg)
  {
    CGAL_precondition_code 
      (
       Segment_traits_2   seg_tr;
       const unsigned int n = this->size()
       );
    // TODO: This precondition seems to cause problems. Probably it is due to
    //       the usage of source/target.
    // CGAL_precondition (seg_tr.compare_xy_2_object() (seg.source(),
    // 						     seg.target()) == SMALLER);
    CGAL_precondition (n == 0 ||
		       seg_tr.equal_2_object() (this->m_segments[n - 1].target(),
						seg.source()));

    this->m_segments.push_back (seg);
  }

private:
  /*! Reverse the polyline. */
  void _reverse()
  {
    typename Base::const_reverse_iterator  ps = this->rbegin();
    typename Base::const_reverse_iterator  pt = ps;
    ++pt;

    std::vector<Segment_2>  rev_segs (this->size());
    unsigned int            i = 0;

    while (pt != this->rend())
    {
      rev_segs[i] = Segment_2 (*ps, *pt);
      ++ps; ++pt;
      i++;
    }

    this->m_segments = rev_segs;
    return;
  }

};

/*! Output operator for a polyline. */
template <typename SegmentTraits>
std::ostream& operator<< (std::ostream & os,
			  const _Polyline_2<SegmentTraits>& cv)
{
  // std::cout << "\033[1;31m\nP2: << operator in Polyline_2.h\033[0m\n";
  typename _Polyline_2<SegmentTraits>::const_iterator  iter = cv.begin();

  // Print the number of points:
  os << cv.points();

  while (iter != cv.end())
  {
    os << "  " << *iter;
    ++iter;
  }
  return (os);
}

/*! Input operator for a polyline. */
template <typename SegmentTraits>
std::istream& operator>> (std::istream& is,
			 _Polyline_2<SegmentTraits>& pl)
{
  typedef _Polyline_2<SegmentTraits>  Curve_2;
  typedef typename Curve_2::Point_2   Point_2;

  // Read the number of input points.
  unsigned int        n_pts;

  is >> n_pts;

  // Read m_num_pts points to a list.
  Point_2             p;
  std::list<Point_2>  pts;
  unsigned int        i;

  for (i = 0; i < n_pts; ++i)
  {
    is >> p;
    pts.push_back(p);
  }

  // Create the polyline curve.
  pl = Curve_2(pts.begin(), pts.end());

  return (is);
}


} //namespace CGAL

#endif
