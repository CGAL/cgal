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
//                 Dror Atariah <dror.atariah@fu-berlin.de>

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

  namespace polyline{

    /*! \class
     * Representation of a polyline.
     */
    template <typename Segment_type_T, typename Point_type_T>
    class Polyline_2 {
    public:
      typedef Segment_type_T                         Segment_type_2;
      typedef Point_type_T                           Point_type_2;

    protected:
      // The segments that comprise the poyline:
      typedef typename std::vector<Segment_type_2>   Segments_container;
      Segments_container                             m_segments;

    public:
      typedef typename Segments_container::size_type Segments_size_type;

      /*! Default constructor. */
      Polyline_2() : m_segments() {}

      Polyline_2(const Segment_type_2& seg) : m_segments()
      {
        m_segments.push_back(seg);
      }

      /*
       * Roadmap:
       * - Currently construction from points is marked as deprecated.
       *   For the sake of backwards compatibility, it is not removed yet.
       * - In the next release, construction from points will be removed.
       * - Then, the remaining construction from segments will become private
       *   in this class as construction is provided by the traits.
       *   Furthermore, the polyline traits class will become a friend of
       *   this class.
       */
      /*!
       * Constructor from a range. The range can be either:
       *
       * For the sake of backwards compatibility we have to keep the possibility
       * of construction from a range of points. Therefore, we have to test the
       * type of the elements of the range again.
       * Since ultimately the construction from points will be deprecated, the
       * actual constructor that we implemented which handles a range of point
       * is already deprecated and SHOULD NOT BE USED.
       *
       * If you want to construct a polyline from a range of points, use the
       * construction functors from the traits class.
       *
       * - Range of points, and the polyline is defined by the order of the
       *   points.
       * - Range of linear object. The polyline is the sequence of linear
       *   objects.
       * \param begin An iterator pointing to the first point in the range.
       * \param end An iterator pointing after the last point in the range.
       * \pre Depends on the range's content. See the implementations for
       *      further details.
       */
      template <typename InputIterator>
      Polyline_2(InputIterator begin, InputIterator end) :
        m_segments()
      {
        typedef typename std::iterator_traits<InputIterator>::value_type VT;
        typedef typename boost::is_same<VT, Point_type_2>::type Is_point;
        construct_polyline(begin, end, Is_point());
      }

      /*!
       * Construct a polyline from a range of segments.
       * \param begin An iterator pointing to the first segment in the range.
       * \param end An iterator pointing after the past-the-end segment
       *        in the range.
       * \pre The end of segment n should be the beginning of segment n+1.
       */
      template <typename InputIterator>
      void construct_polyline(InputIterator begin, InputIterator end,
                              boost::false_type)
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
      CGAL_DEPRECATED void construct_polyline(InputIterator begin,
                                              InputIterator end,
                                              boost::true_type)
      {
        // Check if there are no points in the range:
        InputIterator  ps = begin;

        if (ps == end)
          return;

        InputIterator pt = ps;
        ++pt;

        // The range contains only one point. A degenerated polyline is
        // constructed.
        // With one degenerated segment, where source=target.
        if (pt == end) {
          m_segments.push_back(Segment_type_2(*ps, *ps));
          return;
        }

        // Construct a segment from each to adjacent points.
        // The container has to contain at least two points.
        while (pt != end) {
          m_segments.push_back(Segment_type_2(*ps, *pt));
          ++ps;
          ++pt;
        }
      }

      /*
       * Roadmap: Make this private in the next version (after the traits
       *          becomes friendly...)
       */
      /*!
       * Append a segment to the (x-monotone) polyline.
       * Warning: This is a risky function! Don't use it! Prefer the
       *          provided implementation in the traits class.
       * \param seg The new segment to be appended to the polyline.
       * \pre If the polyline is not empty, seg source must be the
       *      same as the target point of the last segment in the polyline.
       */
      inline void push_back(const Segment_type_2& seg)
      {
        this->m_segments.push_back(seg);
      }

      /*!
       * Append a segment to the (x-monotone) polyline.
       * Warning: This is a risky function! Don't use it! Prefer the
       *          provided implementation in the traits class.
       * \param seg The new segment to be appended to the polyline.
       * \pre If the polyline is not empty, seg source must be the
       *      same as the target point of the last segment in the polyline.
       */
      inline void push_front(const Segment_type_2& seg)
      {
        this->m_segments.insert(this->m_segments.begin(), seg);
      }

      /*!
       * Append a point to the polyline.
       * To properly implemented this function the traits class is needed,
       * thus it is deprecated.
       */
      CGAL_DEPRECATED void push_back(const Point_type_2& p)
      {
        CGAL_assertion(!m_segments.empty());
        Point_type_2 pt = p;
        Point_type_2 ps = m_segments.back().target();
        m_segments.push_back(Segment_type_2(ps, pt));
      }

      /*!
       * \return The bounding-box.
       */
      Bbox_2 bbox() const
      {
        // Compute the union of the bounding boxes of all segments.
        Segments_size_type n = this->number_of_segments();
        Bbox_2 bbox;
        for (std::size_t i = 0; i < n; ++i) {
          bbox = (i > 0) ? (bbox + (*this)[i].bbox()) : (*this)[i].bbox();
        }
        return bbox;
      }

      class const_iterator;
      friend class const_iterator;
      typedef std::reverse_iterator<const_iterator>
      const_reverse_iterator;

      /*! An iterator for the polyline points. */
      class CGAL_DEPRECATED const_iterator {
      public:
        // Type definitions:
        typedef std::bidirectional_iterator_tag     iterator_category;
        typedef Point_type_2                        value_type;
        typedef std::ptrdiff_t                      difference_type;
        typedef size_t                              size_type;
        typedef const value_type&                   reference;
        typedef const value_type*                   pointer;

      private:
        // The polyline curve.
        const Polyline_2<Segment_type_2, Point_type_2>* m_cvP;
        int   m_num_pts;                            // Its number of points.
        int   m_index;                              // The current point.

        /*!
         * Private constructor.
         * \param cv The scanned curve.
         * \param index The index of the segment.
         */
        const_iterator(const Polyline_2<Segment_type_2, Point_type_2>* cvP,
                       int index) :
          m_cvP(cvP),
          m_index(index)
        {
          if (m_cvP == NULL)
            m_num_pts = 0;
          else
            m_num_pts = (m_cvP->number_of_segments() == 0) ?
              0 : (m_cvP->number_of_segments() + 1);
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
        const Point_type_2& operator*() const
        {
          CGAL_assertion(m_cvP != NULL);
          CGAL_assertion((m_index >= 0) && (m_index < m_num_pts));

          if (m_index == 0) {
            // First point is the source of the first segment.
            return ((*m_cvP)[0]).source();
          }
          else {
            // Return the target of the(i-1)'st segment.
            return ((*m_cvP)[m_index - 1]).target();
          }
        }

        /*!
         * Arrow operator.
         * \return A pointer to the current point.
         */
        const Point_type_2* operator->() const
        {
          return (&(this->operator* ()));
        }

        /*! Increment operators. */
        const_iterator& operator++()
        {
          if ((m_cvP != NULL) && (m_index < m_num_pts))
            ++m_index;
          return (*this);
        }

        const_iterator operator++(int)
        {
          const_iterator  temp = *this;
          if ((m_cvP != NULL) && (m_index < m_num_pts))
            ++m_index;
          return temp;
        }

        /*! Decrement operators. */
        const_iterator& operator-- ()
        {
          if ((m_cvP != NULL) && (m_index >= 0))
            --m_index;
          return (*this);
        }

        const_iterator operator--(int)
        {
          const_iterator  temp = *this;
          if ((m_cvP != NULL) && (m_index >= 0))
            --m_index;
          return temp;
        }

        /*! Equality operators. */
        bool operator==(const const_iterator& it) const
        {
          return ((m_cvP == it.m_cvP) && (m_index == it.m_index));
        }

        bool operator!=(const const_iterator& it) const
        {
          return ((m_cvP != it.m_cvP) || (m_index != it.m_index));
        }

        friend class Polyline_2<Segment_type_2, Point_type_2>;
      };

      /* ! Get an iterator for the polyline points.*/
      CGAL_DEPRECATED const_iterator begin() const
      {
        if (number_of_segments() == 0)
          return (const_iterator(NULL, -1));
        else
          return (const_iterator(this, 0));
      }

      /*! Get a past-the-end iterator for the polyline points.*/
      CGAL_DEPRECATED const_iterator end() const
      {
        if (number_of_segments() == 0)
          return (const_iterator(NULL, -1));
        else
          return (const_iterator(this, number_of_segments() + 1));
      }

      /*! Get a reverse iterator for the polyline points. */
      CGAL_DEPRECATED const_reverse_iterator rbegin() const
      {
        return (const_reverse_iterator(end()));
      }

      /*! Get a reverse past-the-end iterator for the polyline points. */
      CGAL_DEPRECATED const_reverse_iterator rend() const
      {
        return (const_reverse_iterator(begin()));
      }

      typedef typename Segments_container::const_iterator
      Segment_const_iterator;
      typedef typename std::reverse_iterator<Segment_const_iterator>
      Segment_const_reverse_iterator;

      /*! Get an iterator for the polyline's segments. */
      Segment_const_iterator begin_segments() const
      { return m_segments.begin(); }

      /*! Get a past-the-end iterator for the polyline's segments. */
      Segment_const_iterator end_segments() const
      { return m_segments.end(); }

      /*! Get a reverse iterator for the polyline's segments. */
      Segment_const_reverse_iterator rbegin_segments() const
      { return (Segment_const_reverse_iterator(end_segments())); }

      /*! Get a reverse past-the-end iterator for the polyline points. */
      Segment_const_reverse_iterator rend_segments() const
      { return (Segment_const_reverse_iterator(begin_segments())); }

      /*! Deprecated!
       * Get the number of points contained in the polyline.
       * In general (for example if the polyline is not bounded), then the
       * number of vertices cannot be read-off from the number of segments, and
       * the traits class is needed.
       * \return The number of points.
       */
      CGAL_DEPRECATED std::size_t points() const
      {
        return (number_of_segments() == 0) ? 0 : number_of_segments() + 1;
      }

      /*! Deprecated! Replaced by number_of_segments()
       * Get the number of segments that comprise the poyline.
       * \return The number of segments.
       */
      CGAL_DEPRECATED Segments_size_type size() const
      { return m_segments.size(); }

      /*!
       * Get the number of segments that comprise the poyline.
       * \return The number of segments.
       */
      Segments_size_type number_of_segments() const
      { return m_segments.size(); }

      /*!
       * Get the ith segment of the polyline.
       * \param i The segment index(from 0 to size()-1).
       * \return A const reference to the segment.
       */
      inline const Segment_type_2& operator[](const std::size_t i) const
      {
        CGAL_assertion(i < number_of_segments());
        return (m_segments[i]);
      }

      /*! Clear the polyline. */
      inline void clear()
      {
        m_segments.clear();
      }
    };

    /*! \class
     * Representation of an x-monotone polyline.
     * An x-monotone polyline is always directed from left to right.
     */
    template <typename Segment_type_2_T, typename Point_type_2_T>
    class X_monotone_polyline_2 :
      public Polyline_2<Segment_type_2_T, Point_type_2_T>
    {
    public:
      typedef Segment_type_2_T                          Segment_type_2;
      typedef Point_type_2_T                            Point_type_2;

      typedef Polyline_2<Segment_type_2, Point_type_2>  Base;

      /*! Default constructor. */
      X_monotone_polyline_2() : Base() {}

      X_monotone_polyline_2(Segment_type_2 seg) : Base(seg){ }

      /*! Constructor
       * Similar to the constructor of a general polyline.
       * Like in the case of general polyline, for the sake of backwards
       * compatibility we have to keep an implementation of construction
       * from a range of points. DO NOT USE THIS CONSTRUCTION.
       */
      template <typename InputIterator>
      X_monotone_polyline_2(InputIterator begin, InputIterator end) :
        Base(begin, end)
      {
        typedef typename std::iterator_traits<InputIterator>::value_type VT;
        typedef typename boost::is_same<VT, Point_type_2>::type Is_point;
        construct_x_monotone_polyline(begin, end, Is_point());
      }

      /*!
       * Constructs from a range of segments.
       * This constructor is expected to be called only from the
       * traits class, after the input was verified there.
       * \pre The range of segments form an x-monotone polyline.
       */
      template <typename InputIterator>
      void construct_x_monotone_polyline(InputIterator, InputIterator,
                                         boost::false_type)
      {}

      /*!
       * Constructs from a range of points, defining the endpoints of the
       * polyline segments.
       */
      template <typename InputIterator>
      CGAL_DEPRECATED void construct_x_monotone_polyline(InputIterator,
                                                         InputIterator,
                                                         boost::true_type)
      {}
    };

    /*! Output operator for a polyline. */
    template <typename Segment_type_2_T, typename Point_type_2_T>
    std::ostream&
    operator<<(std::ostream & os,
               const Polyline_2<Segment_type_2_T, Point_type_2_T>& cv)
    {
      typedef Segment_type_2_T                          Segment_type_2;
      typedef Point_type_2_T                            Point_type_2;
      typedef Polyline_2<Segment_type_2, Point_type_2>  Curve_2;

      typename Curve_2::Segment_const_iterator iter = cv.begin_segments();
      while (iter != cv.end_segments()) {
        if (iter == cv.begin_segments()) {
          os << " " << *iter;
          ++iter;
        }
        else {
          os << " <-> " << *iter;
          ++iter;
        }
      }
      return (os);
    }


    /*! Input operator for a polyline. */
    template <typename Segment_type_2_T, typename Point_type_2_T>
    std::istream&
    operator>>(std::istream& is,
               Polyline_2<Segment_type_2_T, Point_type_2_T>& pl)
    {
      typedef Segment_type_2_T                          Segment_type_2;
      typedef Point_type_2_T                            Point_type_2;
      typedef Polyline_2<Segment_type_2, Point_type_2>  Curve_2;

      // Read the number of input points.
      unsigned int        n_pts;

      is >> n_pts;

      CGAL_precondition_msg(n_pts > 1,
                            "Input must contain at least two points");

      // Read m_num_pts points to a list.
      Point_type_2 p;
      std::list<Point_type_2> pts;
      for (unsigned int i = 0; i < n_pts; ++i) {
        is >> p;
        pts.push_back(p);
      }

      std::list<Segment_type_2> segments;
      typename std::list<Point_type_2>::iterator curr = pts.begin();
      typename std::list<Point_type_2>::iterator next = pts.begin();
      ++next;
      while (next != pts.end()) {
        segments.push_back(Segment_type_2(*curr, *next));
        ++curr;
        ++next;
      }

      // Create the polyline curve.
      pl = Curve_2(segments.begin(), segments.end());

      return (is);
    }
  } // namespace polyline
} //namespace CGAL
#endif
