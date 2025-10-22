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
 * file   data/3d/skel/Straight_skeleton_3.h
 * author Gernot Walzl
 * date   2011-11-26
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_SDS_STRAIGHT_SKELETON_3_H
#define CGAL_STRAIGHT_SKELETON_3_SDS_STRAIGHT_SKELETON_3_H

#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/events/Abstract_event.h>
#include <CGAL/Straight_skeleton_3/internal/weak_find.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>

#include <list>
#include <memory>
#include <string>
#include <unordered_set>

namespace CGAL {

/*!
 * \ingroup PkgStraightSkeleton3Classes
 *
 * This is the 3D straight skeleton.
 */
template <typename Traits_>
class Straight_skeleton_3
  : public std::enable_shared_from_this<Straight_skeleton_3<Traits_> >
{
public:
  using Traits = Traits_;
  using Polyhedron = Straight_skeletons_3::internal::HDS::Polyhedron<Traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using StraightSkeletonWPtr = std::weak_ptr<Straight_skeleton_3<Traits> >;
  using StraightSkeletonSPtr = std::shared_ptr<Straight_skeleton_3<Traits> >;

private:
  using String_factory = Straight_skeletons_3::IO::String_factory;

public:
  class Node;
  class Arc;
  class Sheet;

public:
  /*!
  * \ingroup PkgStraightSkeleton3Classes
  *
  * This is the 3D straight skeleton node class.
  */
  class Node
    : public std::enable_shared_from_this<Node>
  {
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Point3SPtr = std::shared_ptr<Point_3>;

  private:
    using NodeSPtr = std::shared_ptr<Node>;
    using ArcWPtr = std::weak_ptr<Arc>;
    using ArcSPtr = std::shared_ptr<Arc>;
    using SheetWPtr = std::weak_ptr<Sheet>;
    using SheetSPtr = std::shared_ptr<Sheet>;

  public:
    Node() {
      id_ = next_id_++;
    }

    static NodeSPtr create()
    {
      return std::make_shared<Node>();
    }

    Point3SPtr point() const
    {
      CGAL_SS3_DEBUG_SPTR(point_);
      return point_;
    }

    void set_point(const Point3SPtr& point)
    {
      CGAL_SS3_DEBUG_SPTR(point);
      this->point_ = point;
    }

    const FT& time() const
    {
      return time_;
    }

    void set_time(const FT& time)
    {
      this->time_ = time;
    }

    StraightSkeletonSPtr get_skeleton() const
    {
      return this->skel_.lock();
    }

    void set_skeleton(const StraightSkeletonSPtr& skel)
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

    int get_ID() const
    {
      return this->id_;
    }

    void set_ID(const int id)
    {
      this->id_ = id;
    }

    void add_arc(const ArcSPtr& arc)
    {
      CGAL_SS3_DEBUG_SPTR(arc);
      CGAL_precondition(!has_incident_arc(arc));
      typename std::list<ArcWPtr>::iterator it = arcs_.insert(arcs_.end(), ArcWPtr(arc));
      if (arc->getNodeSrc() == this->shared_from_this()) {
        arc->setNodeSrcListIt(it);
      } else if (arc->hasNodeDst()) {
        if (arc->getNodeDst() == this->shared_from_this()) {
          arc->setNodeDstListIt(it);
        }
      }
    }

    bool remove_arc(const ArcSPtr& arc)
    {
      CGAL_SS3_DEBUG_SPTR(arc);
      CGAL_precondition(has_incident_arc(arc));
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

    void add_sheet(const SheetSPtr& sheet)
    {
      CGAL_SS3_DEBUG_SPTR(sheet);
      CGAL_precondition(!has_incident_sheet(sheet));
      sheets_.insert(sheets_.end(), SheetWPtr(sheet));
    }

    bool remove_sheet(const SheetSPtr& sheet)
    {
      CGAL_SS3_DEBUG_SPTR(sheet);
      CGAL_precondition(has_incident_sheet(sheet));
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

    bool has_incident_arc(const ArcSPtr& arc) const
    {
      CGAL_SS3_DEBUG_SPTR(arc);
      ArcWPtr arc_wptr = ArcWPtr(arc);
      bool result = (arcs_.end() != STL_Extension::internal::weak_find(arcs_.begin(), arcs_.end(), arc_wptr));
      return result;
    }

    bool has_incident_sheet(const SheetSPtr& sheet) const
    {
      CGAL_SS3_DEBUG_SPTR(sheet);
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
          remove_sheet(sheet);
        }
      }
      sheets_.clear();
      typename std::list<ArcWPtr>::iterator it_a = arcs_.begin();
      while (it_a != arcs_.end()) {
        ArcWPtr arc_wptr = *it_a++;
        if (ArcSPtr arc = arc_wptr.lock()) {
          remove_arc(arc);
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

    std::string to_string() const
    {
      std::string result("Node(");
      if (id_ != -1) {
        result += "id=" + String_factory::fromInteger(id_) + ", ";
      } else {
        result += String_factory::fromPointer(this) + ", ";
      }
      result += "<" + String_factory::fromDouble(CGAL::to_double(point()->x())) + " ";
      result += String_factory::fromDouble(CGAL::to_double(point()->y())) + " ";
      result += String_factory::fromDouble(CGAL::to_double(point()->z())) + ">, ";
      result += ", arcs={";
      bool first = true;
      for (const ArcWPtr& arc_wptr : arcs_) {
        if (ArcSPtr arc = arc_wptr.lock()) {
          if (!first) result += ", ";
          result += String_factory::fromInteger(arc->get_ID());
          first = false;
        }
      }
      result += "}";
      result += ", sheets={";
      first = true;
      for (const SheetWPtr& sheet_wptr : sheets_) {
        if (SheetSPtr sheet = sheet_wptr.lock()) {
          if (!first) result += ", ";
          result += String_factory::fromInteger(sheet->get_ID());
          first = false;
        }
      }
      result += "}";

      result += ")";
      return result;
    }

  protected:
    Point3SPtr point_;
    FT time_;
    int id_;

    std::list<ArcWPtr> arcs_;
    std::list<SheetWPtr> sheets_;
    StraightSkeletonWPtr skel_;

    typename std::list<NodeSPtr>::iterator list_it_;

  private:
    static int next_id_;
  };

public:
  class Arc
    : public std::enable_shared_from_this<Arc>
  {
    using Vector_3 = typename Traits::Vector_3;
    using Line_3 = typename Traits::Line_3;

    using Vector3SPtr = std::shared_ptr<Vector_3>;
    using Line3SPtr = std::shared_ptr<Line_3>;

  private:
    using NodeSPtr = std::shared_ptr<Node>;
    using ArcWPtr = std::weak_ptr<Arc>;
    using ArcSPtr = std::shared_ptr<Arc>;
    using SheetWPtr = std::weak_ptr<Sheet>;
    using SheetSPtr = std::shared_ptr<Sheet>;

    using Kernel_factory = Straight_skeletons_3::internal::kernel::Kernel_factory<Traits>;

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
      node_src->add_arc(result);
      return result;
    }

    static ArcSPtr create(const NodeSPtr& node_src, const NodeSPtr& node_dst)
    {
      ArcSPtr result = std::make_shared<Arc>(node_src, node_dst);
      node_src->add_arc(result);
      node_dst->add_arc(result);
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
      node_dst->add_arc(this->shared_from_this());
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

    StraightSkeletonSPtr get_skeleton() const
    {
      return this->skel_.lock();
    }

    void set_skeleton(const StraightSkeletonSPtr& skel)
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

    int get_ID() const
    {
      return this->id_;
    }

    void set_ID(const int id)
    {
      this->id_ = id;
    }

    void add_sheet(const SheetSPtr& sheet)
    {
      CGAL_SS3_DEBUG_SPTR(sheet);
      CGAL_precondition(!has_incident_sheet(sheet));
      sheets_.insert(sheets_.end(), SheetWPtr(sheet));
    }

    bool remove_sheet(const SheetSPtr& sheet)
    {
      CGAL_SS3_DEBUG_SPTR(sheet);
      CGAL_precondition(has_incident_sheet(sheet));
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

    bool has_incident_sheet(const SheetSPtr& sheet) const
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
      CGAL_precondition(has_incident_sheet(sheet));
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
          if (arc->has_incident_sheet(sheet)) {
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
        result = Kernel_factory::createLine3(node_src_->point(), node_dst_->point());
      } else if (direction_) {
        result = Kernel_factory::createLine3(node_src_->point(), direction_);
      }
      return result;
    }

    std::string to_string() const
    {
      std::string result("Arc(");
      // Arc ID
      if (id_ != -1) {
        result += "id=" + String_factory::fromInteger(id_) + ", ";
      } else {
        result += String_factory::fromPointer(this) + ", ";
      }
      // Node IDs
      result += "src=";
      if (node_src_) {
        result += String_factory::fromInteger(node_src_->get_ID());
      } else {
        result += "-1";
      }
      result += ", dst=";
      if (node_dst_) {
        result += String_factory::fromInteger(node_dst_->get_ID());
      } else {
        result += "-1";
      }
      // Incident sheet IDs
      result += ", sheets={";
      bool first = true;
      for (const SheetWPtr& sheet_wptr : sheets_) {
        if (SheetSPtr sheet = sheet_wptr.lock()) {
          if (!first) result += ", ";
          result += String_factory::fromInteger(sheet->get_ID());
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

public:
  class Sheet
    : public std::enable_shared_from_this<Sheet>
  {
    using Plane_3 = typename Traits::Plane_3;

    using Plane3SPtr = std::shared_ptr<Plane_3>;

  private:
    using NodeSPtr = std::shared_ptr<Node>;
    using ArcSPtr = std::shared_ptr<Arc>;
    using SheetSPtr = std::shared_ptr<Sheet>;

    using EdgeWPtr = typename Polyhedron::EdgeWPtr;
    using FacetSPtr = typename Polyhedron::FacetSPtr;

  public:
    Sheet() {
      id_ = next_id_++;
    }

    ~Sheet()
    {
      facet_b_.reset();
      facet_f_.reset();
    }

    static SheetSPtr create()
    {
      return std::make_shared<Sheet>();
    }

    FacetSPtr getFacetB() const
    {
      CGAL_SS3_DEBUG_SPTR(facet_b_);
      return facet_b_;
    }

    void setFacetB(const FacetSPtr& facet_b)
    {
      CGAL_SS3_DEBUG_SPTR(facet_b);
      facet_b_ = facet_b;
    }

    FacetSPtr getFacetF() const
    {
      CGAL_SS3_DEBUG_SPTR(facet_f_);
      return facet_f_;
    }

    void setFacetF(const FacetSPtr& facet_f)
    {
      CGAL_SS3_DEBUG_SPTR(facet_f);
      facet_f_ = facet_f;
    }

    void add_edge(const EdgeWPtr& edge)
    {
      CGAL_precondition(STL_Extension::internal::weak_find(edges_.begin(), edges_.end(), edge) == edges_.end());
      edges_.push_back(edge);
    }

    const std::list<EdgeWPtr>& edges() const
    {
      return edges_;
    }

    StraightSkeletonSPtr get_skeleton() const
    {
      return this->skel_.lock();
    }

    void set_skeleton(const StraightSkeletonSPtr& skel)
    {
      this->skel_ = skel;
    }

    typename std::list<SheetSPtr>::iterator getListIt() const
    {
      return this->list_it_;
    }

    void setListIt(typename std::list<SheetSPtr>::iterator list_it)
    {
      this->list_it_ = list_it;
    }

    int get_ID() const
    {
      return this->id_;
    }

    void set_ID(const int id)
    {
      this->id_ = id;
    }

    Plane3SPtr get_plane() const
    {
      CGAL_SS3_DEBUG_SPTR(this->plane_);
      return this->plane_;
    }

    void set_plane(const Plane3SPtr& plane)
    {
      CGAL_SS3_DEBUG_SPTR(plane);
      this->plane_ = plane;
    }

    void add_node(const NodeSPtr& node)
    {
      CGAL_SS3_DEBUG_SPTR(node);
      CGAL_precondition((std::find(nodes_.begin(), nodes_.end(), node) == nodes_.end()));
      /*std::list<NodeSPtr>::iterator it =*/ nodes_.insert(nodes_.end(), node);
      if (!node->has_incident_sheet(this->shared_from_this())) {
        node->add_sheet(this->shared_from_this());
      }
    }

    bool remove_node(const NodeSPtr& node)
    {
      CGAL_SS3_DEBUG_SPTR(node);
      bool result = false;
      nodes_.remove(node);
      result = node->remove_sheet(this->shared_from_this());
      return result;
    }

    void add_arc(const ArcSPtr& arc)
    {
      CGAL_SS3_DEBUG_SPTR(arc);
      /*std::list<ArcSPtr>::iterator it =*/ arcs_.insert(arcs_.end(), arc);
      arc->add_sheet(this->shared_from_this());
    }

    bool remove_arc(const ArcSPtr& arc)
    {
      CGAL_SS3_DEBUG_SPTR(arc);
      bool result = false;
      arcs_.remove(arc);
      result = arc->remove_sheet(this->shared_from_this());
      return result;
    }

    void addContour(const ArcSPtr& contour)
    {
      CGAL_SS3_DEBUG_SPTR(contour);
      /*std::list<ArcSPtr>::iterator it =*/ contours_.insert(contours_.end(), contour);
      contour->add_sheet(this->shared_from_this());
    }

    /**
    * merge 'sheet' into this sheet.
    */
    void merge(const SheetSPtr& sheet)
    {
      CGAL_SS3_DEBUG_SPTR(sheet);

      CGAL_precondition((this->facet_b_ == sheet->getFacetB() && this->facet_f_ == sheet->getFacetF()) ||
                        (this->facet_b_ == sheet->getFacetF() && this->facet_f_ == sheet->getFacetB()));

      // copy all arcs that do not yet exist in the sheet
      for (ArcSPtr arc : sheet->arcs()) {
        if (std::find(arcs_.begin(), arcs_.end(), arc) == arcs_.end()) {
          add_arc(arc);
        }
      }

      // copy all contours that do not yet exist in the sheet
      for (ArcSPtr contour : sheet->contours()) {
        if (std::find(contours_.begin(), contours_.end(), contour) == contours_.end()) {
          addContour(contour);
        }
      }

      // copy all nodes that do not yet exist in the sheet
      for (NodeSPtr node : sheet->nodes()) {
        if (std::find(nodes_.begin(), nodes_.end(), node) == nodes_.end()) {
          add_node(node);
        }
      }
    }

    std::list<ArcSPtr>& arcs()
    {
      return this->arcs_;
    }

    std::list<ArcSPtr>& contours()
    {
      return this->contours_;
    }

    std::list<NodeSPtr>& nodes()
    {
      return this->nodes_;
    }

    std::string to_string() const
    {
      std::stringstream result;
      result << "Sheet(";
      if (id_ != -1) {
        result << "id=" << String_factory::fromInteger(id_) << ", ";
      } else {
        result << String_factory::fromPointer(this) << ", ";
      }

      // node IDs
      result << "nodes={";
      bool first = true;
      for (const NodeSPtr& node : nodes_) {
        if (!first) result << ", ";
        result << String_factory::fromInteger(node->get_ID());
        first = false;
      }
      result << "}, ";

      // contour IDs
      result << "contours={";
      first = true;
      for (const ArcSPtr& contour : contours_) {
        if (!first) result << ", ";
        result << String_factory::fromInteger(contour->get_ID());
        first = false;
      }
      result << "}, ";

      // arc IDs
      result << "arcs={";
      first = true;
      for (const ArcSPtr& arc : arcs_) {
        if (!first) result << ", ";
        result << String_factory::fromInteger(arc->get_ID());
        first = false;
      }
      result << "}";

      result << ")";
      return result.str();
    }

  protected:
    // facets of the static polyhedron
    FacetSPtr facet_b_;
    FacetSPtr facet_f_;

    std::list<ArcSPtr> arcs_;
    // Squatting the arc type for contours, which are edges of the static polyhedron.
    // Contours are needed to generate a closed polygon (with holes) when drawing sheets.
    std::list<ArcSPtr> contours_;
    std::list<NodeSPtr> nodes_;
    StraightSkeletonWPtr skel_;

    typename std::list<SheetSPtr>::iterator list_it_;

    // these are edges of the *mutable* polyhedron, whose active movement traces the sheet
    std::list<EdgeWPtr> edges_;

    int id_;
    Plane3SPtr plane_;

  private:
    static int next_id_;
  };

public:
  using NodeWPtr = std::weak_ptr<Node>;
  using NodeSPtr = std::shared_ptr<Node>;
  using ArcWPtr = std::weak_ptr<Arc>;
  using ArcSPtr = std::shared_ptr<Arc>;
  using SheetWPtr = std::weak_ptr<Sheet>;
  using SheetSPtr = std::shared_ptr<Sheet>;

private:
  using Abstract_event = Straight_skeletons_3::internal::algorithm::Abstract_event<Traits>;
  using Abstract_event_sptr = std::shared_ptr<Abstract_event>;

public:
  ~Straight_skeleton_3()
  {
    polyhedron_.reset();
    sheets_.clear();
    arcs_.clear();
    nodes_.clear();
  }

  static StraightSkeletonSPtr create()
  {
    return std::make_shared<Straight_skeleton_3>();
  }

  void add_node(NodeSPtr node)
  {
    typename std::list<NodeSPtr>::iterator it = nodes_.insert(nodes_.end(), node);
    node->set_skeleton(this->shared_from_this());
    node->setListIt(it);
  }

  bool remove_node(NodeSPtr node)
  {
    bool result = false;
    if (node->get_skeleton() == this->shared_from_this()) {
      nodes_.erase(node->getListIt());
      node->set_skeleton(StraightSkeletonSPtr());
      node->setListIt(std::list<NodeSPtr>::iterator());
      result = true;
    }
    return result;
  }

  void add_arc(ArcSPtr arc)
  {
    typename std::list<ArcSPtr>::iterator it = arcs_.insert(arcs_.end(), arc);
    arc->set_skeleton(this->shared_from_this());
    arc->setListIt(it);
  }

  bool remove_arc(ArcSPtr arc)
  {
    bool result = false;
    if (arc->get_skeleton() == this->shared_from_this()) {
      arcs_.erase(arc->getListIt());
      arc->set_skeleton(StraightSkeletonSPtr());
      arc->setListIt(typename std::list<ArcSPtr>::iterator());
      result = true;
    }
    return result;
  }

  void add_sheet(SheetSPtr sheet)
  {
    typename std::list<SheetSPtr>::iterator it = sheets_.insert(sheets_.end(), sheet);
    sheet->set_skeleton(this->shared_from_this());
    sheet->setListIt(it);
  }

  bool remove_sheet(SheetSPtr sheet)
  {
    bool result = false;
    if (sheet->get_skeleton() == this->shared_from_this()) {
      sheets_.erase(sheet->getListIt());
      sheet->set_skeleton(StraightSkeletonSPtr());
      sheet->setListIt(typename std::list<SheetSPtr>::iterator());
      result = true;
    }
    return result;
  }

  void merge_sheets(SheetSPtr sheet_into,
                    SheetSPtr sheet_from)
  {
    CGAL_SS3_DEBUG_SPTR(sheet_into);
    CGAL_SS3_DEBUG_SPTR(sheet_from);

    sheet_into->merge(sheet_from);

    for (NodeSPtr node : sheet_from->nodes()) {
      node->remove_sheet(sheet_from);
    }
    for (ArcSPtr arc : sheet_from->arcs()) {
      arc->remove_sheet(sheet_from);
    }
    remove_sheet(sheet_from);
  }

  PolyhedronSPtr get_polyhedron() const {
      return this->polyhedron_;
  }

  void set_polyhedron(PolyhedronSPtr polyhedron) {
      this->polyhedron_ = polyhedron;
  }

public:
  /*!
   * returns the sheets of the straight skeleton.
   */
  std::list<SheetSPtr>& sheets()
  {
    return this->sheets_;
  }

  /*!
   * returns the arcs of the straight skeleton.
   */
  std::list<ArcSPtr>& arcs()
  {
    return this->arcs_;
  }

  /*!
   * returns the nodes of the straight skeleton.
   */
  std::list<NodeSPtr>& nodes()
  {
    return this->nodes_;
  }

public:
  void reset_all_IDs()
  {
    typename std::list<SheetSPtr>::iterator it_s = sheets_.begin();
    while (it_s != sheets_.end()) {
      SheetSPtr sheet = *it_s++;
      sheet->set_ID(-1);
    }
    typename std::list<ArcSPtr>::iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
      ArcSPtr arc = *it_a++;
      arc->set_ID(-1);
    }
    typename std::list<NodeSPtr>::iterator it_n = nodes_.begin();
    while (it_n != nodes_.end()) {
      NodeSPtr node = *it_n++;
      node->set_ID(-1);
    }
  }

  bool is_consistent(bool is_partial = true) const
  {
    bool result = true;

    // Check for duplicate elements
    {
      std::unordered_set<NodeSPtr> node_set;
      for (const NodeSPtr& node : nodes_) {
        if (!node_set.insert(node).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate node in nodes_ list");
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          result = false;
          break;
        }
      }
      std::unordered_set<ArcSPtr> arc_set;
      for (const ArcSPtr& arc : arcs_) {
        if (!arc_set.insert(arc).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate arc in arcs_ list");
          CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
          result = false;
          break;
        }
      }
      std::unordered_set<SheetSPtr> sheet_set;
      for (const SheetSPtr& sheet : sheets_) {
        if (!sheet_set.insert(sheet).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate sheet in sheets_ list");
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          result = false;
          break;
        }
      }
    }

    // Check uniqueness of Node IDs
    {
      std::unordered_set<int> node_ids;
      for (const NodeSPtr& node : nodes_) {
        int id = node->get_ID();
        if (!node_ids.insert(id).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate Node ID in nodes_ list: " << id);
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          result = false;
          break;
        }
      }
    }
    // Check uniqueness of Arc IDs
    {
      std::unordered_set<int> arc_ids;
      for (const ArcSPtr& arc : arcs_) {
        int id = arc->get_ID();
        if (!arc_ids.insert(id).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate Arc ID in arcs_ list: " << id);
          CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
          result = false;
          break;
        }
      }
    }
    // Check uniqueness of Sheet IDs
    {
      std::unordered_set<int> sheet_ids;
      for (const SheetSPtr& sheet : sheets_) {
        int id = sheet->get_ID();
        if (!sheet_ids.insert(id).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate Sheet ID in sheets_ list: " << id);
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          result = false;
          break;
        }
      }
    }

    typename std::list<NodeSPtr>::const_iterator it_n = nodes_.begin();
    while (it_n != nodes_.end()) {
      NodeSPtr node = *it_n++;
      CGAL_SS3_DEBUG_SPTR(node);
      if (node->get_skeleton() != this->shared_from_this()) {
        CGAL_SS3_SKEL_DS_TRACE(node->to_string());
        result = false;
        break;
      }
      typename std::list<ArcWPtr>::const_iterator it_a = node->arcs().begin();
      while (it_a != node->arcs().end()) {
        ArcWPtr arc_wptr = *it_a++;
        CGAL_SS3_DEBUG_WPTR(arc_wptr);
        ArcSPtr arc = arc_wptr.lock();
        if (arc_wptr.expired() || !arc) {
          CGAL_SS3_SKEL_DS_TRACE("Error: invalid ArcWPtr (expired or cannot lock) in node's arc list");
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          result = false;
          break;
        }

        if (std::find(arcs_.begin(), arcs_.end(), arc) == arcs_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: node references an arc not present in the skeleton's arc list");
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
          result = false;
          break;
        }

        if (node != arc->getNodeSrc() && node != arc->getNodeDst()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: incident arc does not have node as endpoint");
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
          result = false;
          break;
        }
      }
      typename std::list<SheetWPtr>::const_iterator it_s = node->sheets().begin();
      while (it_s != node->sheets().end()) {
        SheetWPtr sheet_wptr = *it_s++;
        CGAL_SS3_DEBUG_WPTR(sheet_wptr);
        SheetSPtr sheet = sheet_wptr.lock();
        if (sheet_wptr.expired() || !sheet) {
          CGAL_SS3_SKEL_DS_TRACE("Error: invalid SheetWPtr (expired or cannot lock) in node's sheet list");
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          result = false;
          break;
        }

        if (std::find(sheets_.begin(), sheets_.end(), sheet) == sheets_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: node references a sheet not present in the skeleton's sheet list");
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          result = false;
          break;
        }

        std::list<NodeSPtr> nodes = sheet->nodes();
        if (nodes.end() == std::find(nodes.begin(), nodes.end(), node)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: incident sheet does not have node");
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          result = false;
          break;
        }
      }
    }

    typename std::list<ArcSPtr>::const_iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
      ArcSPtr arc = *it_a++;
      CGAL_SS3_DEBUG_SPTR(arc);
      ArcWPtr arc_wptr(arc);

      if (arc->get_skeleton() != this->shared_from_this()) {
        CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
        result = false;
        break;
      }

      if (!arc->getNodeSrc()) {
        CGAL_SS3_SKEL_DS_TRACE("Error: arc does not have source");
        result = false;
        break;
      }

      if (std::find(nodes_.begin(), nodes_.end(), arc->getNodeSrc()) == nodes_.end()) {
        CGAL_SS3_SKEL_DS_TRACE("Error: arc source is not a node present in the skeleton's node list");
        CGAL_SS3_SKEL_DS_TRACE(arc->getNodeSrc()->to_string());
        result = false;
        break;
      }

      std::list<ArcWPtr> warcs = arc->getNodeSrc()->arcs();
      if (warcs.end() == STL_Extension::internal::weak_find(warcs.begin(), warcs.end(), arc_wptr)) {
        CGAL_SS3_SKEL_DS_TRACE("Error: arc's source node does not have arc");
        CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
        CGAL_SS3_SKEL_DS_TRACE(arc->getNodeSrc()->to_string());
        result = false;
        break;
      }

      if (arc->hasNodeDst()) {
        if (std::find(nodes_.begin(), nodes_.end(), arc->getNodeDst()) == nodes_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc destination is not a node present in the skeleton's node list");
          CGAL_SS3_SKEL_DS_TRACE(arc->getNodeDst()->to_string());
          result = false;
          break;
        }

        warcs = arc->getNodeDst()->arcs();
        if (warcs.end() == STL_Extension::internal::weak_find(warcs.begin(), warcs.end(), arc_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc's target node does not have arc");
          CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
          CGAL_SS3_SKEL_DS_TRACE(arc->getNodeDst()->to_string());
          result = false;
          break;
        }
      } else {
        if (is_partial) {
          if (!arc->getDirection()) {
            CGAL_SS3_SKEL_DS_TRACE("Error: arc does not have destination nor direction (partial skeleton)");
            result = false;
            break;
          }
        } else {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc does not have destination (full skeleton)");
          result = false;
          break;
        }
      }
      unsigned int num_sheets = 0;
      typename std::list<SheetWPtr>::const_iterator it_s = arc->sheets().begin();
      while (it_s != arc->sheets().end()) {
        SheetWPtr sheet_wptr = *it_s++;
        CGAL_SS3_DEBUG_WPTR(sheet_wptr);
        SheetSPtr sheet = sheet_wptr.lock();
        if (sheet_wptr.expired() || !sheet) {
          CGAL_SS3_SKEL_DS_TRACE("Error: invalid SheetWPtr (expired or cannot lock) in arc's sheet list");
          CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
          result = false;
          break;
        }

        if (std::find(sheets_.begin(), sheets_.end(), sheet) == sheets_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc references a sheet that is not present in the skeleton's sheet list");
          CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          result = false;
          break;
        }

        ++num_sheets;

        std::list<ArcSPtr> arcs = sheet->arcs();
        if (arcs.end() == std::find(arcs.begin(), arcs.end(), arc)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc's sheet does not have arc");
          CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          result = false;
          break;
        }
      }

      if (num_sheets != 3) {
        CGAL_SS3_SKEL_DS_TRACE("Error: Arc does not have 3 sheets.");
        CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
        CGAL_SS3_SKEL_DS_TRACE("num_sheets = " << num_sheets);
        result = false;
        break;
      }
    }
    typename std::list<SheetSPtr>::const_iterator it_s = sheets_.begin();
    while (it_s != sheets_.end()) {
      SheetSPtr sheet = *it_s++;
      CGAL_SS3_DEBUG_SPTR(sheet);
      SheetWPtr sheet_wptr(sheet);
      if (sheet->get_skeleton() != this->shared_from_this()) {
        CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
        result = false;
        break;
      }
      if (!polyhedron_->has_facet(sheet->getFacetB())) {
        CGAL_SS3_SKEL_DS_TRACE("Error: sheet's facet B is not in the polyhedron");
        CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
        result = false;
      }
      if (!polyhedron_->has_facet(sheet->getFacetF())) {
        CGAL_SS3_SKEL_DS_TRACE("Error: sheet's facet F is not in the polyhedron");
        CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
        result = false;
      }
      typename std::list<NodeSPtr>::const_iterator it_n = sheet->nodes().begin();
      while (it_n != sheet->nodes().end()) {
        NodeSPtr node = *it_n++;
        CGAL_SS3_DEBUG_SPTR(node);

        if (std::find(nodes_.begin(), nodes_.end(), node) == nodes_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: sheet references a node that is not present in the skeleton's node list");
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          result = false;
          break;
        }

        std::list<SheetWPtr> wsheets = node->sheets();
        if (wsheets.end() == STL_Extension::internal::weak_find(wsheets.begin(), wsheets.end(), sheet_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: sheet's node does not have sheet");
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          result = false;
          break;
        }
      }

      // check that a node only has two incident arcs at most in the same sheet
      std::unordered_map<NodeSPtr, int> arc_nodes;

      typename std::list<ArcSPtr>::const_iterator it_a = sheet->arcs().begin();
      while (it_a != sheet->arcs().end()) {
        ArcSPtr arc = *it_a++;
        CGAL_SS3_DEBUG_SPTR(arc);

        if (std::find(arcs_.begin(), arcs_.end(), arc) == arcs_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: sheet references an arc that is not present in the skeleton's arc list");
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
          result = false;
          break;
        }

        std::list<SheetWPtr> wsheets = arc->sheets();
        if (wsheets.end() == STL_Extension::internal::weak_find(wsheets.begin(), wsheets.end(), sheet_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: sheet's arc does not have sheet");
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
          result = false;
          break;
        }

        arc_nodes[arc->getNodeSrc()]++;
        if (arc->hasNodeDst()) {
          arc_nodes[arc->getNodeDst()]++;
        }
        for (auto [node, count] : arc_nodes) {
          if (count > 2) {
            CGAL_SS3_SKEL_DS_TRACE("Error: sheet has a node with more than two incident arcs in the sheet");
            CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
            CGAL_SS3_SKEL_DS_TRACE(node->to_string());
            result = false;
            break;
          }
          CGAL_SS3_DEBUG_SPTR(node);
          std::list<NodeSPtr> nodes = sheet->nodes();
          if (nodes.end() == std::find(nodes.begin(), nodes.end(), node)) {
            CGAL_SS3_SKEL_DS_TRACE("Error: Nodes of arc of sheet should be nodes of sheet");
            CGAL_SS3_SKEL_DS_TRACE(node->to_string());
            CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
            CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
            result = false;
            break;
          }
        }
      }

      // Contour checks
      for (const ArcSPtr& contour : sheet->contours()) {
        // Check that contour has exactly one incident sheet (the current sheet)
        int sheet_count = 0;
        for (const SheetWPtr& contour_sheet_wptr : contour->sheets()) {
          SheetSPtr contour_sheet = contour_sheet_wptr.lock();
          if (contour_sheet == sheet) {
            ++sheet_count;
          }
        }
        if (sheet_count != 1) {
          CGAL_SS3_SKEL_DS_TRACE("Error: contour arc does not have exactly one incident sheet (should be the current sheet)");
          CGAL_SS3_SKEL_DS_TRACE(contour->to_string());
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          result = false;
          break;
        }

        // Check that contour has two valid, distinct incident nodes
        NodeSPtr node_src = contour->getNodeSrc();
        NodeSPtr node_dst = contour->getNodeDst();
        if (!node_src || !node_dst) {
          CGAL_SS3_SKEL_DS_TRACE("Error: contour does not have two valid incident nodes");
          CGAL_SS3_SKEL_DS_TRACE(contour->to_string());
          result = false;
          break;
        }

        // Check that both nodes appear in the sheet's node list
        const auto& sheet_nodes = sheet->nodes();
        if (std::find(sheet_nodes.begin(), sheet_nodes.end(), node_src) == sheet_nodes.end() ||
            std::find(sheet_nodes.begin(), sheet_nodes.end(), node_dst) == sheet_nodes.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: contour's nodes do not appear in sheet's node list");
          CGAL_SS3_SKEL_DS_TRACE(contour->to_string());
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          result = false;
          break;
        }
      }
    }

    // @todo check contours

    return result;
  }

  std::string to_string() const
  {
    std::stringstream sstr;
    sstr << "Straight_skeleton_3\n";
    sstr << "  Nodes:  " << nodes_.size() << std::endl;
    for (const NodeSPtr& node : nodes_) {
      sstr << "    " << node->to_string() << "\n";
    }
    sstr << "  Arcs:   " << arcs_.size() << std::endl;
    for (const ArcSPtr& arc : arcs_) {
      sstr << "    " << arc->to_string() << "\n";
    }
    sstr << "  Sheets: " << sheets_.size() << std::endl;
    for (const SheetSPtr& sheet : sheets_) {
      sstr << "    " << sheet->to_string() << "\n";
    }
    return sstr.str();
  }

protected:
  PolyhedronSPtr polyhedron_;

  std::list<NodeSPtr> nodes_;
  std::list<ArcSPtr> arcs_;
  std::list<SheetSPtr> sheets_;
};

} // namespace CGAL

template <typename Traits>
int CGAL::Straight_skeleton_3<Traits>::Node::next_id_ = 0;

template <typename Traits>
int CGAL::Straight_skeleton_3<Traits>::Arc::next_id_ = 0;

template <typename Traits>
int CGAL::Straight_skeleton_3<Traits>::Sheet::next_id_ = 0;

#endif /* CGAL_STRAIGHT_SKELETON_3_SDS_STRAIGHT_SKELETON_3_H */
