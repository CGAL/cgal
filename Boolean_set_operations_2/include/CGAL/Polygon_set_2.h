// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_POLYGON_SET_2_H
#define CGAL_POLYGON_SET_2_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_default_dcel.h>

#include <vector>

namespace CGAL {

template <class Kernel,
          typename Containter = std::vector<typename Kernel::Point_2>,
          class Dcel_ =
            Gps_default_dcel<Gps_segment_traits_2<Kernel, Containter> > >
class Polygon_set_2 :
  public General_polygon_set_2<Gps_segment_traits_2<Kernel, Containter>, Dcel_>
{  
private:
  typedef General_polygon_set_2<Gps_segment_traits_2<Kernel, Containter>, Dcel_>
                                                          Base;
  typedef typename Base::Base                             Base_of_base;
  typedef Polygon_set_2<Kernel, Containter, Dcel_>        Self;

public:
  typedef  typename Base::Traits_2                        Traits_2;    
  typedef  typename Base::Polygon_2                       Polygon_2;
  typedef  typename Base::Polygon_with_holes_2            Polygon_with_holes_2;
  typedef  typename Base::Arrangement_2                   Arrangement_2;
  typedef  typename Base::Size                            Size;

  /*! Default constructor. */
  Polygon_set_2 () :
    Base()
  {}

  /*! Constructor from the base class. */
  Polygon_set_2 (const Base& base) :
    Base (base)
  {}

  /*! Constructor with traits object. */
  Polygon_set_2 (Traits_2& tr) :
    Base(tr)
  {}

  /*! Constructor from a polygon. */
  explicit Polygon_set_2 (const Polygon_2& pgn) :
    Base (pgn)
  {}

  /*! Constructor from a polygon with holes. */
  explicit Polygon_set_2 (const Polygon_with_holes_2& pwh):
    Base (pwh)
  {}

  
  inline void intersection(const Polygon_2& pgn)
  {
    Base::intersection(pgn);
  }

  inline void intersection(const Polygon_with_holes_2& pgn)
  {
    Base::intersection(pgn);
  }

  inline void intersection(const Self& other)
  {
    Base::intersection(other);
  }

  inline void intersection(const Self& ps1, const Self& ps2)
  {
    Base::intersection(base(ps1), base(ps2));
  }

  template <class InputIterator>
  inline void intersection(InputIterator begin, InputIterator end)
  {
    Base::intersection(begin, end);
  }

  template <class InputIterator1, class InputIterator2>
  inline void intersection(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2)
  {
    Base::intersection(begin1, end1, begin2, end2);
  }

  inline void join(const Polygon_2& pgn)
  {
    Base::join(pgn);
  }

  inline void join(const Polygon_with_holes_2& pgn)
  {
    Base::join(pgn);
  }

  inline void join(const Self& other)
  {
    Base::join(other);
  }

  inline void join(const Self& ps1, const Self& ps2)
  {
    Base::join(base(ps1), base(ps2));
  }

  template <class InputIterator>
  inline void join(InputIterator begin, InputIterator end)
  {
    Base::join(begin, end);
  }

  template <class InputIterator1, class InputIterator2>
  inline void join(InputIterator1 begin1, InputIterator1 end1,
                   InputIterator2 begin2, InputIterator2 end2)
  {
    Base::join(begin1, end1, begin2, end2);
  }

  inline void difference(const Polygon_2& pgn)
  {
    Base::difference(pgn);
  }

  inline void difference(const Polygon_with_holes_2& pgn)
  {
    Base::difference(pgn);
  }

  inline void difference(const Self& other)
  {
    Base::difference(other);
  }

  inline void difference(const Self& ps1, const Self& ps2)
  {
    Base::difference(base(ps1), base(ps2));
  }

  template <class InputIterator>
  inline void difference(InputIterator begin, InputIterator end)
  {
    Base::difference(begin, end);
  }

  template <class InputIterator1, class InputIterator2>
  inline void difference(InputIterator1 begin1, InputIterator1 end1,
                         InputIterator2 begin2, InputIterator2 end2)
  {
    Base::difference(begin1, end1, begin2, end2);
  }

  inline void symmetric_difference(const Polygon_2& pgn)
  {
    Base::symmetric_difference(pgn);
  }

  inline void symmetric_difference(const Polygon_with_holes_2& pgn)
  {
    Base::symmetric_difference(pgn);
  }

  inline void symmetric_difference(const Self& other)
  {
    Base::symmetric_difference(other);
  }

  inline void symmetric_difference(const Self& ps1, const Self& ps2)
  {
    Base::symmetric_difference(base(ps1), base(ps2));
  }

  template <class InputIterator>
  inline void symmetric_difference(InputIterator begin, InputIterator end)
  {
    Base::symmetric_difference(begin, end);
  }

  template <class InputIterator1, class InputIterator2>
  inline void symmetric_difference(InputIterator1 begin1, InputIterator1 end1,
                                   InputIterator2 begin2, InputIterator2 end2)
  {
    Base::symmetric_difference(begin1, end1, begin2, end2);
  }

  inline bool do_intersect(const Polygon_2& pgn)
  {
    return (Base::do_intersect(pgn));
  }

  inline bool do_intersect(const Polygon_with_holes_2& pgn)
  {
    return (Base::do_intersect(pgn));
  }

  inline bool do_intersect(const Self& other)
  {
    return (Base::do_intersect(other));
  }

  template <class InputIterator>
  inline bool do_intersect(InputIterator begin, InputIterator end)
  {
    return (Base::do_intersect(begin, end));
  }

  template <class InputIterator1, class InputIterator2>
  inline bool do_intersect(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2)
  {
    return (Base::do_intersect(begin1, end1, begin2, end2));
  }

private:

  inline const Base_of_base& base(const Self& other) const
  {
    return (static_cast<const Base_of_base&>(other));
  }

};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
