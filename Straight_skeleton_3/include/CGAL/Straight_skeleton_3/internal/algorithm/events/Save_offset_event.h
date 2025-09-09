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
 * @file   data/3d/skel/SaveOffsetEvent.h
 * @author Gernot Walzl
 * @date   2013-12-28
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SAVE_OFFSET_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SAVE_OFFSET_EVENT_H

#include <CGAL/Straight_skeleton_3/internal/algorithm/events/Abstract_event.h>

#include <memory>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class SaveOffsetEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using SaveOffsetEventSPtr = std::shared_ptr<SaveOffsetEvent<Traits> >;

private:
  using FT = typename Traits::FT;

public:
  SaveOffsetEvent()
    : Base(AbstractEvent<Traits>::SAVE_OFFSET_EVENT),
      offset_(-1)
  { }

  SaveOffsetEvent(const FT& offset)
    : Base(AbstractEvent<Traits>::SAVE_OFFSET_EVENT),
      offset_(offset)
  { }

  static SaveOffsetEventSPtr create()
  {
    return std::make_shared<SaveOffsetEvent>();
  }

  static SaveOffsetEventSPtr create(const FT& offset)
  {
    return std::make_shared<SaveOffsetEvent>(offset);
  }

  const FT& getOffset() const
  {
    return this->offset_;
  }

  void setOffset(const FT& offset)
  {
    this->offset_ = offset;
  }

  bool operator==(const SaveOffsetEvent& other) const
  {
    return (this->offset_ == other.offset_);
  }

protected:
  FT offset_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SAVE_OFFSET_EVENT_H */
