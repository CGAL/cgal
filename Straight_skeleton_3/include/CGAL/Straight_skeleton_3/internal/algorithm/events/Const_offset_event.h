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
 * @file   data/3d/skel/ConstOffsetEvent.h
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_CONST_OFFSET_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_CONST_OFFSET_EVENT_H

#include <CGAL/Straight_skeleton_3/internal/algorithm/events/Abstract_event.h>

#include <memory>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class ConstOffsetEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using ConstOffsetEventSPtr = std::shared_ptr<ConstOffsetEvent<Traits> >;

private:
  using FT = typename Traits::FT;

public:
  ConstOffsetEvent()
    : Base(Base::CONST_OFFSET_EVENT),
      offset_(-1)
  { }

  ConstOffsetEvent(const FT& offset = FT{-1})
    : Base(Base::CONST_OFFSET_EVENT),
      offset_(offset)
  { }

  static ConstOffsetEventSPtr create()
  {
    return std::make_shared<ConstOffsetEvent>();
  }

  static ConstOffsetEventSPtr create(const FT& offset)
  {
    return std::make_shared<ConstOffsetEvent>(offset);
  }

  const FT& getTime() const
  {
    return this->offset_;
  }

  void setTime(const FT& offset)
  {
    this->offset_ = offset;
  }

  bool operator==(const ConstOffsetEvent& other) const
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

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_CONST_OFFSET_EVENT_H */
