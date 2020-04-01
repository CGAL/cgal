// Copyright (c) 2017 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_COMPARE_HANDLES_WITH_OR_WITHOUT_TIMESTAMPS_H
#define CGAL_COMPARE_HANDLES_WITH_OR_WITHOUT_TIMESTAMPS_H

#include <CGAL/Has_timestamp.h>
#include <CGAL/tags.h>
#include <iterator>

namespace CGAL {

struct Compare_handles_with_or_without_timestamps
{
  template<typename Handle>
  bool operator()(const Handle h1, const Handle h2) const
  {
    typedef typename std::iterator_traits<Handle>::value_type Type;
    return less(h1, h2,
                Boolean_tag<CGAL::internal::Has_timestamp<Type>::value>());
  }

  template<typename Handle>
  bool less(const Handle h1, const Handle h2, Tag_false) const
  {
    return &*h1 < &*h2;
  }

  template<typename Handle>
  bool less(const Handle h1, const Handle h2, Tag_true) const
  {
    if(h1 == Handle())      return (h2 != Handle());
    else if(h2 == Handle()) return false;
    else                    return h1->time_stamp() < h2->time_stamp();
  }
};

} // end namespace CGAL

#endif // CGAL_COMPARE_HANDLES_WITH_OR_WITHOUT_TIMESTAMPS_H
