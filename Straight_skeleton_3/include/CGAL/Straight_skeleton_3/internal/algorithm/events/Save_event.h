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
 * file   data/3d/skel/SaveEvent.h
 * author Gernot Walzl
 * date   2013-12-28
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SAVE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SAVE_EVENT_H

#include <CGAL/Straight_skeleton_3/internal/algorithm/events/Abstract_event.h>

#include <memory>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class SaveEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using SaveEventSPtr = std::shared_ptr<SaveEvent<Traits> >;

public:
  SaveEvent()
    : Base(AbstractEvent<Traits>::SAVE_EVENT)
  { }

  static SaveEventSPtr create()
  {
    return std::make_shared<SaveEvent>();
  }

  bool operator==(const SaveEvent& other) const
  {
    return (Base::getTime() == other.getTime());
  }
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SAVE_EVENT_H */
