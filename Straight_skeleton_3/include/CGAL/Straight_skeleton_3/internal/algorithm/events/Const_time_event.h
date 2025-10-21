// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   data/3d/skel/ConstTimeEvent.h
 * author Gernot Walzl
 * date   2012-03-27
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_CONST_TIME_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_CONST_TIME_EVENT_H

#include <CGAL/Straight_skeleton_3/internal/algorithm/events/Abstract_event.h>

#include <memory>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class Const_time_event
  : public Abstract_event<Traits>
{
  using Base = Abstract_event<Traits>;
  using Const_time_event_sptr = std::shared_ptr<Const_time_event<Traits> >;

public:
  Const_time_event()
    : Base(Base::CONST_TIME_EVENT)
  { }

  static Const_time_event_sptr create()
  {
    return std::make_shared<Const_time_event>();
  }

  bool operator==(const Const_time_event& other) const
  {
    return (Base::time() == other.time());
  }
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_CONST_OFFSET_EVENT_H */
