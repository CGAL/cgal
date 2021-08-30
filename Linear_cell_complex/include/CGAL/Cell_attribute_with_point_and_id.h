// Copyright (c) 2017 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_CELL_ATTRIBUTE_WITH_POINT_AND_ID_H
#define CGAL_CELL_ATTRIBUTE_WITH_POINT_AND_ID_H 1

#include <CGAL/Cell_attribute_with_point.h>

namespace CGAL {

  // A cell attribute with point and id, when Info_!=void
  template <class Refs, class Info_=void, class Tag_=Tag_true,
            class OnMerge=Null_functor,
            class OnSplit=Null_functor>
  class Cell_attribute_with_point_and_id: public
      Cell_attribute_with_point<Refs, Info_, Tag_, OnMerge, OnSplit, Tag_true>
  {
    typedef Cell_attribute_with_point
          <Refs, Info_, Tag_, OnMerge, OnSplit, Tag_true> Base;

    template <class, class, class, class>
    friend class Compact_container;

    template <class, class>
    friend class Concurrent_compact_container;

  public:
    typedef typename Base::Point Point;

  protected:
    /// Default contructor.
    Cell_attribute_with_point_and_id()
    {}

    /// Contructor with an info in parameter.
    Cell_attribute_with_point_and_id(const Point& apoint) : Base(apoint)
    {}

    /// Contructor with a point and an attribute in parameters.
    Cell_attribute_with_point_and_id(const Point& apoint, const Info_& ainfo) :
      Base(apoint, ainfo)
    {}
  };

  /// Specialization when Info==void.
  template <class Refs, class Tag_, class OnMerge, class OnSplit>
  class Cell_attribute_with_point_and_id<Refs, void, Tag_, OnMerge, OnSplit>:
      public Cell_attribute_with_point<Refs, void, Tag_, OnMerge, OnSplit, Tag_true>
  {
    typedef Cell_attribute_with_point
         <Refs, void, Tag_, OnMerge, OnSplit, Tag_true> Base;

    template <class, class, class, class>
    friend class Compact_container;

    template <class, class>
    friend class Concurrent_compact_container;

  public:
    typedef typename Base::Point Point;

  protected:
    /// Default contructor.
    Cell_attribute_with_point_and_id()
    {}

    /// Contructor with a point in parameter.
    Cell_attribute_with_point_and_id(const Point& apoint) : Base(apoint)
    {}
  };
} // namespace CGAL

#endif // CGAL_CELL_ATTRIBUTE_WITH_POINT_AND_ID_H //
// EOF //
