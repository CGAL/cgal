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
 * file   data/3d/skel/Arc.h
 * author Gernot Walzl
 * date   2012-03-27
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
  : public std::enable_shared_from_this<Arc<Traits> >
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
  Arc(const NodeSPtr& node_src, const Vector3SPtr& direction)
  {
    node_src_ = node_src;
    direction_ = direction;
    id_ = next_id_++;
  }

  Arc(const NodeSPtr& node_src, const NodeSPtr& node_dst)
  {
    node_src_ = node_src;
    node_dst_ = node_dst;
    id_ = next_id_++;
  }

  ~Arc() {
    node_src_.reset();
    node_dst_.reset();
    direction_.reset();
    sheets_.clear();
  }

  static ArcSPtr create(const NodeSPtr& node_src, const Vector3SPtr& direction)
  {
    ArcSPtr result = std::make_shared<Arc>(node_src, direction);
    node_src->addArc(result);
    return result;
  }

  static ArcSPtr create(const NodeSPtr& node_src, const NodeSPtr& node_dst)
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

  void setNodeSrc(const NodeSPtr& node_src)
  {
    CGAL_SS3_DEBUG_SPTR(node_src);
    CGAL_precondition(!hasNodeSrc() && !hasNodeDst());
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

  bool hasNodeSrc() const
  {
    return bool(node_src_);
  }

  NodeSPtr getNodeDst() const
  {
    CGAL_SS3_DEBUG_SPTR(node_dst_);
    return node_dst_;
  }

  void setNodeDst(const NodeSPtr& node_dst)
  {
    CGAL_SS3_DEBUG_SPTR(node_dst);
    CGAL_precondition(hasNodeSrc() && !hasNodeDst());
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

  bool hasNodeDst() const
  {
    return bool(node_dst_);
  }

  NodeSPtr other(const NodeSPtr& node) const
  {
    CGAL_precondition(node == node_src_ || node == node_dst_);
    if (node == node_src_) {
      return node_dst_;
    } else {
      return node_src_;
    }
  }

  void closeArc(const NodeSPtr& node_dst)
  {
    CGAL_SS3_DEBUG_SPTR(node_dst);
    CGAL_precondition(hasNodeSrc() && !hasNodeDst());
    setNodeDst(node_dst);
    node_dst->addArc(this->shared_from_this());
  }

  Vector3SPtr getDirection() const
  {
    CGAL_SS3_DEBUG_SPTR(direction_);
    return this->direction_;
  }

  void setDirection(const Vector3SPtr& direction)
  {
    this->direction_ = direction;
  }

  StraightSkeletonSPtr getSkel() const
  {
    return this->skel_.lock();
  }

  void setSkel(const StraightSkeletonSPtr& skel)
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

  void setID(const int id)
  {
    this->id_ = id;
  }

  void addSheet(const SheetSPtr& sheet)
  {
    CGAL_SS3_DEBUG_SPTR(sheet);
    CGAL_precondition(!containsSheet(sheet));
    sheets_.insert(sheets_.end(), SheetWPtr(sheet));
  }

  bool removeSheet(const SheetSPtr& sheet)
  {
    CGAL_SS3_DEBUG_SPTR(sheet);
    CGAL_precondition(containsSheet(sheet));
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

  bool containsSheet(const SheetSPtr& sheet) const
  {
    CGAL_SS3_DEBUG_SPTR(sheet);
    SheetWPtr sheet_wptr = SheetWPtr(sheet);
    bool result = (sheets_.end() != STL_Extension::internal::weak_find(sheets_.begin(), sheets_.end(), sheet_wptr));
    return result;
  }

  std::list<SheetWPtr>& sheets()
  {
    return this->sheets_;
  }

  ArcSPtr next(const SheetSPtr& sheet, const NodeSPtr& node) const
  {
    CGAL_precondition(containsSheet(sheet));
    CGAL_precondition(node == node_src_ || node == node_dst_);
    ArcSPtr result = ArcSPtr();
    // check the arcs incident to the node
    std::list<ArcWPtr> warcs = node->arcs();
    typename std::list<ArcWPtr>::const_iterator it = warcs.begin();
    while (it != warcs.end()) {
      ArcWPtr arc_wptr = *it++;
      CGAL_SS3_DEBUG_WPTR(arc_wptr);
      if (ArcSPtr arc = arc_wptr.lock()) {
        if (arc == this->shared_from_this()) {
          continue;
        }
        if (arc->containsSheet(sheet)) {
          result = arc;
          break;
        }
      }
    }
    return result;
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

  std::string toString() const
  {
    std::string result("Arc(");
    // Arc ID
    if (id_ != -1) {
      result += "id=" + IO::StringFactory::fromInteger(id_) + ", ";
    } else {
      result += IO::StringFactory::fromPointer(this) + ", ";
    }
    // Node IDs
    result += "src=";
    if (node_src_) {
      result += IO::StringFactory::fromInteger(node_src_->getID());
    } else {
      result += "-1";
    }
    result += ", dst=";
    if (node_dst_) {
      result += IO::StringFactory::fromInteger(node_dst_->getID());
    } else {
      result += "-1";
    }
    // Incident sheet IDs
    result += ", sheets={";
    bool first = true;
    for (const SheetWPtr& sheet_wptr : sheets_) {
      if (SheetSPtr sheet = sheet_wptr.lock()) {
        if (!first) result += ", ";
        result += IO::StringFactory::fromInteger(sheet->getID());
        first = false;
      }
    }
    result += "}";

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

private:
  static int next_id_;
};

template <typename Traits>
int Arc<Traits>::next_id_ = 0;

} // namespace SDS
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_ARC_H */

