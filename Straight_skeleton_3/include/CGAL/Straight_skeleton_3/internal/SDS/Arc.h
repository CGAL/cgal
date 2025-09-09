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
 * @file   data/3d/skel/Arc.h
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_ARC_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_ARC_H

#include <CGAL/Straight_skeleton_3/IO/String_factory.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>

#include <list>
#include <memory>
#include <sstream>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace SDS {

template <typename Traits>
class Node;

template <typename Traits>
class Sheet;

template <typename Traits>
class StraightSkeleton;

template <typename Traits>
class Arc
{
  using Vector_3 = typename Traits::Vector_3;
  using Line_3 = typename Traits::Line_3;

  using Vector3SPtr = std::shared_ptr<Vector_3>;
  using Line3SPtr = std::shared_ptr<Line_3>;

private:
  using NodeSPtr = std::shared_ptr<Node<Traits> >;
  using ArcWPtr = std::weak_ptr<Arc<Traits> >;
  using ArcSPtr = std::shared_ptr<Arc<Traits> >;
  using SheetWPtr = std::weak_ptr<Sheet<Traits> >;
  using SheetSPtr = std::shared_ptr<Sheet<Traits> >;
  using StraightSkeletonWPtr = std::weak_ptr<StraightSkeleton<Traits> >;
  using StraightSkeletonSPtr = std::shared_ptr<StraightSkeleton<Traits> >;

  using KernelFactory = kernel::KernelFactory<Traits>;

public:
  Arc(NodeSPtr node_src, Vector3SPtr direction)
  {
    node_src_ = node_src;
    direction_ = direction;
    id_ = -1;
  }

  Arc(NodeSPtr node_src, NodeSPtr node_dst)
  {
    node_src_ = node_src;
    node_dst_ = node_dst;
    id_ = -1;
  }

  ~Arc() {
    node_src_.reset();
    node_dst_.reset();
    direction_.reset();
    sheets_.clear();
  }

  static ArcSPtr create(NodeSPtr node_src, Vector3SPtr direction)
  {
    ArcSPtr result = std::make_shared<Arc>(node_src, direction);
    node_src->addArc(result);
    return result;
  }

  static ArcSPtr create(NodeSPtr node_src, NodeSPtr node_dst)
  {
    ArcSPtr result = std::make_shared<Arc>(node_src, node_dst);
    node_src->addArc(result);
    node_dst->addArc(result);
    return result;
  }

  NodeSPtr getNodeSrc() const
  {
    CGAL_SS3_DEBUG_SPTR(node_src_);
    return node_src_;
  }

  void setNodeSrc(NodeSPtr node_src)
  {
    this->node_src_ = node_src;
  }

  typename std::list<ArcWPtr>::iterator getNodeSrcListIt() const
  {
    return this->node_src_list_it_;
  }

  void setNodeSrcListIt(typename std::list<ArcWPtr>::iterator node_src_list_it)
  {
    this->node_src_list_it_ = node_src_list_it;
  }

  NodeSPtr getNodeDst() const
  {
    CGAL_SS3_DEBUG_SPTR(node_dst_);
    return node_dst_;
  }

  void setNodeDst(NodeSPtr node_dst)
  {
    this->node_dst_ = node_dst;
  }

  typename std::list<ArcWPtr>::iterator getNodeDstListIt() const
  {
    return this->node_dst_list_it_;
  }

  void setNodeDstListIt(typename std::list<ArcWPtr>::iterator node_dst_list_it)
  {
    this->node_dst_list_it_ = node_dst_list_it;
  }

  Vector3SPtr getDirection() const
  {
    CGAL_SS3_DEBUG_SPTR(direction_);
    return this->direction_;
  }

  void setDirection(Vector3SPtr direction)
  {
    this->direction_ = direction;
  }

  StraightSkeletonSPtr getSkel() const
  {
    return this->skel_.lock();
  }

  void setSkel(StraightSkeletonSPtr skel)
  {
    this->skel_ = skel;
  }

  typename std::list<ArcSPtr>::iterator getListIt() const
  {
    return this->list_it_;
  }

  void setListIt(typename std::list<ArcSPtr>::iterator list_it)
  {
    this->list_it_ = list_it;
  }

  int getID() const
  {
    return this->id_;
  }

  void setID(int id)
  {
    this->id_ = id;
  }

  void addSheet(SheetSPtr sheet)
  {
    sheets_.insert(sheets_.end(), SheetWPtr(sheet));
  }

  bool removeSheet(SheetSPtr sheet)
  {
    bool result = false;
    SheetWPtr sheet_wptr;
    typename std::list<SheetWPtr>::iterator it = sheets_.begin();
    while (it != sheets_.end()) {
      sheet_wptr = *it;
      if (sheet_wptr.lock() == sheet) {
        sheets_.erase(it);
        result = true;
        break;
      }
      it++;
    }
    return result;
  }

  std::list<SheetWPtr>& sheets()
  {
    return this->sheets_;
  }

  Line3SPtr line() const
  {
    Line3SPtr result = Line3SPtr();
    if (node_dst_) {
      result = KernelFactory::createLine3(node_src_->getPoint(), node_dst_->getPoint());
    } else if (direction_) {
      result = KernelFactory::createLine3(node_src_->getPoint(), direction_);
    }
    return result;
  }

  bool hasNodeDst() const
  {
    bool result = false;
    if (node_dst_) {
      result = true;
    }
    return result;
  }

  std::string toString() const
  {
    std::string result("Arc(");
    if (id_ != -1) {
      result += "id=" + IO::StringFactory::fromInteger(id_) + ", ";
    } else {
      result += IO::StringFactory::fromPointer(this) + ", ";
    }
    result += "src=" + node_src_->toString() + ", ";
    if (node_dst_) {
      result += "dst=" + node_dst_->toString();
    } else {
      result += "dir=<" + IO::StringFactory::fromDouble(CGAL::to_double((*direction_)[0])) + ", " +
                          IO::StringFactory::fromDouble(CGAL::to_double((*direction_)[1])) + ", " +
                          IO::StringFactory::fromDouble(CGAL::to_double((*direction_)[2])) + ">";
    }
    result += ")";
    return result;
  }

protected:
  NodeSPtr node_src_;
  typename std::list<ArcWPtr>::iterator node_src_list_it_;
  NodeSPtr node_dst_;
  typename std::list<ArcWPtr>::iterator node_dst_list_it_;
  Vector3SPtr direction_;
  typename std::list<SheetWPtr> sheets_; // every arc has 3 sheets
  StraightSkeletonWPtr skel_;

  typename std::list<ArcSPtr>::iterator list_it_;

  int id_;
};

} // namespace SDS
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_ARC_H */

