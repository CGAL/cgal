// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_CELL_ATTRIBUTE_WITH_POINT_H
#define CGAL_CELL_ATTRIBUTE_WITH_POINT_H 1

#include <CGAL/Cell_attribute.h>

namespace CGAL {

  /** @file Cell_attribute_with_point.h
   * Definition of cell attribute with point, with or without info.
   */

  /// Point associated with a cell.
  template < class Point >
  class Point_for_cell
  {
  public:
    /// Contructor without parameter.
    Point_for_cell()
    {}

    /// Contructor with a point in parameter.
    Point_for_cell(const Point& apoint) : mpoint(apoint)
    {}

    /// Get the point associated with the cell.
    Point& point()
    { return mpoint; }

    /// Get the point associated with the cell.
    const Point& point() const
    { return mpoint; }

  protected:
    /// The point associated with the cell.
    Point mpoint;
  };

  /// Attribute associated with a point and an info.
  template < class LCC, class Info_=void, class Tag=Tag_true,
             class Functor_on_merge_=Null_functor,
             class Functor_on_split_=Null_functor,
             class WithID=Tag_false >
  class Cell_attribute_with_point :
    public Cell_attribute<LCC, Info_, Tag,
                          Functor_on_merge_, Functor_on_split_, WithID>,
    public Point_for_cell<typename LCC::Point>
  {
    template <class, class, class, class>
    friend class Compact_container;

    template <class, class>
    friend class Concurrent_compact_container;

  public:
    typedef Cell_attribute_with_point<LCC, Info_, Tag, Functor_on_merge_,
                                      Functor_on_split_> Self;

    typedef Cell_attribute<LCC, Info_, Tag,
                           Functor_on_merge_, Functor_on_split_> Base1;
    typedef Point_for_cell<typename LCC::Point> Base2;

    typedef typename LCC::Point             Point;
    typedef typename LCC::Dart_handle       Dart_handle;
    typedef typename LCC::Dart_const_handle Dart_const_handle;

    typedef Info_                Info;
    typedef Functor_on_merge_    Functor_on_merge;
    typedef Functor_on_split_    Functor_on_split;

    using Base1::info;

    bool operator==(const Self& other) const
    { return Base1::operator==(other) && this->point()==other.point(); }

    bool operator!=(const Self& other) const
    { return !operator==(other); }

  protected:
    /// Default contructor.
    Cell_attribute_with_point()
    {}

    /// Contructor with a point in parameter.
    Cell_attribute_with_point(const Point& apoint) : Base2(apoint)
    {}

    /// Contructor with a point and an attribute in parameters.
    Cell_attribute_with_point(const Point& apoint, const Info& ainfo) :
      Base1(ainfo),
      Base2(apoint)
    {}
  };

  /// Attribute associated with a point and without info.
  template < class LCC, class Tag,
             class Functor_on_merge_,
             class Functor_on_split_,
             class WithID>
  class Cell_attribute_with_point<LCC, void, Tag,
                                  Functor_on_merge_, Functor_on_split_,
                                  WithID>:
    public Cell_attribute<LCC, void,
                          Tag, Functor_on_merge_, Functor_on_split_, WithID>,
    public Point_for_cell<typename LCC::Point>
  {
    template <class, class, class, class>
    friend class Compact_container;

    template <class, class>
    friend class Concurrent_compact_container;

  public:
    typedef Cell_attribute<LCC, void, Tag,
                           Functor_on_merge_, Functor_on_split_, WithID> Base1;
    typedef Point_for_cell<typename LCC::Point> Base2;

    typedef void                            Info;
    typedef typename LCC::Point             Point;
    typedef typename LCC::Dart_handle       Dart_handle;
    typedef typename LCC::Dart_const_handle Dart_const_handle;

    typedef Functor_on_merge_ Functor_on_merge;
    typedef Functor_on_split_ Functor_on_split;

    bool operator==(const Cell_attribute_with_point& other) const
    { return Base1::operator==(other) && this->point()==other.point(); }

    bool operator!=(const Cell_attribute_with_point& other) const
    { return !operator==(other); }

    template<typename Cellattr>
    bool operator==(const Cellattr&) const
    { return false; }

  protected:
    /// Default contructor.
    Cell_attribute_with_point()
    {}

    /// Contructor with a point in parameter.
    Cell_attribute_with_point(const Point& apoint) : Base2(apoint)
    {}
  };

} // namespace CGAL

#endif // CGAL_CELL_ATTRIBUTE_WITH_POINT_H //
// EOF //
