// Copyright (c) 2023  GeometryFactory Sarl (France).
// All rights reserved.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_SIMPLEX_BASE_WITH_TIME_STAMP_H
#define CGAL_SIMPLEX_BASE_WITH_TIME_STAMP_H

#include <CGAL/tags.h> // for Tag_true
#include <cstdint>     // for std::size_t
#include <utility>     // for std::forward

namespace CGAL {

/// @brief A base class with a time stamp.
///
/// This class is a base class that can be used to add a time stamp to any class.
/// \cgalModels{BaseWithTimeStamp}
template <typename BaseWithTSBase>
class Triangulation_simplex_base_with_time_stamp
  : public BaseWithTSBase
 {
  std::size_t time_stamp_ = std::size_t(-2);

public:
  using BaseWithTSBase::BaseWithTSBase;

  using Has_timestamp = CGAL::Tag_true;

  std::size_t time_stamp() const {
    return time_stamp_;
  }
  void set_time_stamp(const std::size_t& ts) {
    time_stamp_ = ts;
  }

  template < class TDS >
  struct Rebind_TDS {
    typedef typename BaseWithTSBase::template Rebind_TDS<TDS>::Other Base2;
    typedef Triangulation_simplex_base_with_time_stamp<Base2> Other;
  };
};

} // namespace CGAL

#endif // CGAL_SIMPLEX_BASE_WITH_TIME_STAMP_H
