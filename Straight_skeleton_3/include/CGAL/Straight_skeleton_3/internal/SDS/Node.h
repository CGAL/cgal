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
 * @file   data/3d/skel/Node.h
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_NODE_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_NODE_H

#include <CGAL/Straight_skeleton_3/internal/weak_find.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>

#include <list>
#include <memory>
#include <sstream>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace SDS {

template <typename Traits>
class Arc;

template <typename Traits>
class Sheet;

template <typename Traits>
class StraightSkeleton;

template <typename Traits>
class Node
  : public std::enable_shared_from_this<Node<Traits> >
{
  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;

  using Point3SPtr = std::shared_ptr<Point_3>;

private:
  using NodeSPtr = std::shared_ptr<Node<Traits> >;
  using ArcWPtr = std::weak_ptr<Arc<Traits> >;
  using ArcSPtr = std::shared_ptr<Arc<Traits> >;
  using SheetWPtr = std::weak_ptr<Sheet<Traits> >;
  using SheetSPtr = std::shared_ptr<Sheet<Traits> >;
  using StraightSkeletonWPtr = std::weak_ptr<StraightSkeleton<Traits> >;
  using StraightSkeletonSPtr = std::shared_ptr<StraightSkeleton<Traits> >;

public:
  Node()
    : offset_(0), id_(-1)
  { }

  Node(Point3SPtr point)
    : point_(point), offset_(0), id_(-1)
  { }

  static NodeSPtr create()
  {
    return std::make_shared<Node>();
  }

  static NodeSPtr create(Point3SPtr point)
  {
    return std::make_shared<Node>(point);
  }

  Point3SPtr getPoint() const
  {
    CGAL_SS3_DEBUG_SPTR(point_);
    return point_;
  }

  void setPoint(Point3SPtr point)
  {
    this->point_ = point;
  }

  const FT& getTime() const
  {
    return offset_;
  }

  void setTime(const FT& offset)
  {
    this->offset_ = offset;
  }

  StraightSkeletonSPtr getSkel() const
  {
    return this->skel_.lock();
  }

  void setSkel(StraightSkeletonSPtr skel)
  {
    this->skel_ = skel;
  }

  typename std::list<NodeSPtr>::iterator getListIt() const
  {
    return this->list_it_;
  }

  void setListIt(typename std::list<NodeSPtr>::iterator list_it)
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

  void addArc(ArcSPtr arc)
  {
    typename std::list<ArcWPtr>::iterator it = arcs_.insert(arcs_.end(), ArcWPtr(arc));
    if (arc->getNodeSrc() == this->shared_from_this()) {
      arc->setNodeSrcListIt(it);
    } else if (arc->hasNodeDst()) {
      if (arc->getNodeDst() == this->shared_from_this()) {
        arc->setNodeDstListIt(it);
      }
    }
  }

  bool removeArc(ArcSPtr arc)
  {
    bool result = false;
    if (arc->getNodeSrc() == this->shared_from_this()) {
      arcs_.erase(arc->getNodeSrcListIt());
      arc->setNodeSrcListIt(typename std::list<ArcWPtr>::iterator());
      result = true;
    } else if (arc->hasNodeDst()) {
      if (arc->getNodeDst() == this->shared_from_this()) {
        arcs_.erase(arc->getNodeDstListIt());
        arc->setNodeDstListIt(typename std::list<ArcWPtr>::iterator());
        result = true;
      }
    }
    return result;
  }

  void addSheet(SheetSPtr sheet)
  {
    sheets_.insert(sheets_.end(), SheetWPtr(sheet));
  }

  bool removeSheet(SheetSPtr sheet)
  {
    bool result = false;
    typename std::list<SheetWPtr>::iterator it = sheets_.begin();
    while (it != sheets_.end()) {
      typename std::list<SheetWPtr>::iterator it_current = it;
      SheetWPtr sheet_wptr = *it++;
      if (sheet_wptr.lock() == sheet) {
        sheets_.erase(it_current);
        result = true;
        break;
      }
    }
    return result;
  }

  bool containsArc(ArcSPtr arc) const
  {
    ArcWPtr arc_wptr = ArcWPtr(arc);
    bool result = (arcs_.end() != STL_Extension::internal::weak_find(arcs_.begin(), arcs_.end(), arc_wptr));
    return result;
  }

  bool containsSheet(SheetSPtr sheet) const
  {
    SheetWPtr sheet_wptr = SheetWPtr(sheet);
    bool result = (sheets_.end() != STL_Extension::internal::weak_find(sheets_.begin(), sheets_.end(), sheet_wptr));
    return result;
  }

  void clear()
  {
    typename std::list<SheetWPtr>::iterator it_s = sheets_.begin();
    while (it_s != sheets_.end()) {
      SheetWPtr sheet_wptr = *it_s++;
      if (SheetSPtr sheet = sheet_wptr.lock()) {
        removeSheet(sheet);
      }
    }
    sheets_.clear();
    typename std::list<ArcWPtr>::iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
      ArcWPtr arc_wptr = *it_a++;
      if (ArcSPtr arc = arc_wptr.lock()) {
        removeArc(arc);
      }
    }
    arcs_.clear();
  }

  std::list<ArcWPtr>& arcs()
  {
    return this->arcs_;
  }

  std::list<SheetWPtr>& sheets()
  {
    return this->sheets_;
  }

  unsigned int degree() const
  {
    unsigned int result = 0;
    typename std::list<ArcWPtr>::const_iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
      ArcWPtr arc_wptr = *it_a++;
      if (!arc_wptr.expired()) {
        ++result;
      }
    }
    return result;
  }

  std::string toString() const
  {
    std::string result("Node(");
    if (id_ != -1) {
      result += "id=" + IO::StringFactory::fromInteger(id_) + ", ";
    } else {
      result += IO::StringFactory::fromPointer(this) + ", ";
    }
    result += "<" + IO::StringFactory::fromDouble(CGAL::to_double(getPoint()->x())) + " ";
    result += IO::StringFactory::fromDouble(CGAL::to_double(getPoint()->y())) + " ";
    result += IO::StringFactory::fromDouble(CGAL::to_double(getPoint()->z())) + ">)";
    return result;
  }

protected:
  Point3SPtr point_;
  FT offset_;

  std::list<ArcWPtr> arcs_;
  std::list<SheetWPtr> sheets_;
  StraightSkeletonWPtr skel_;

  typename std::list<NodeSPtr>::iterator list_it_;

  int id_;
};

} // namespace SDS
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_NODE_H */
