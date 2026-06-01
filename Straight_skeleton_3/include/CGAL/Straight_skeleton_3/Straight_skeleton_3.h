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

#include <CGAL/license/Straight_skeleton_3.h>

#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Geom_utils.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/HDS_utils.h>
#include <CGAL/Straight_skeleton_3/internal/weak_find.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>

#include <list>
#include <memory>
#include <optional>
#include <string>
#include <unordered_set>

namespace CGAL {

namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename GeomTraits>
class Hds_utils;

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3

/*!
 * \ingroup PkgStraightSkeleton3Classes
 *
 * The class `Straight_skeleton_3` is a data structure providing access
 * to the combinatorial and geometrical information of a 3D straight skeleton.
 *
 * A 3D straight skeleton is composed of nodes, arcs, and sheets.
 *
 * - Nodes are points at which topological changes have occurred during the propagation of the faces
 *   of the polyhedron.
 * - Arcs connect two nodes, and represent the trace of a vertex of the polyhedron during
 *   the propagation of the faces of the polyhedron. The vertex spawned at the source of the arc,
 *   and disappears at the target of the arc.
 * - Sheets are 2D surfaces connecting arcs and nodes, representing the trace of the edges
 *   of the polyhedron during the propagation of the faces.
 *
 * If a maximal propagation time was set during construction of the skeleton (this is necessarily
 * the case for outward offsetting), the straight skeleton may be partial, i.e., some arcs
 * and sheets may be unbounded.
 *
 * \tparam GeomTraits must be a model of `Kernel` and be the same type as the one used to
 *         construct the straight skeleton.
 */
template <typename GeomTraits>
class Straight_skeleton_3
  : public std::enable_shared_from_this<Straight_skeleton_3<GeomTraits> >
{
  /*!
  * Convenience typedef
  */
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Plane_3 = typename Geom_traits::Plane_3;

  using Polyhedron = Straight_skeletons_3::internal::HDS::Polyhedron<Geom_traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

public:
  using StraightSkeletonWPtr = std::weak_ptr<Straight_skeleton_3<Geom_traits> >;
  using StraightSkeletonSPtr = std::shared_ptr<Straight_skeleton_3<Geom_traits> >;

private:
  using String_factory = Straight_skeletons_3::IO::String_factory;
  using Kernel_wrapper = Straight_skeletons_3::internal::kernel::Kernel_wrapper<Geom_traits>;
  using Geom_utils = Straight_skeletons_3::internal::algorithm::Geom_utils<Geom_traits>;
  using Hds_utils = Straight_skeletons_3::internal::algorithm::Hds_utils<Geom_traits>;

public:
  class Node;
  class Arc;
  class Sheet;

public:
  /*!
  * \ingroup PkgStraightSkeleton3Classes
  *
  * The class `Node` represents a node of the straight skeleton, i.e., a combinatorial entity
  * representing a point in 3D space at which an event has occurred during the propagation
  * of the faces of the polyhedron.
  *
  * A node stores its 3D position, the time at which it was created, and the incident arcs and sheets.
  *
  * \sa `CGAL::Straight_skeleton_3`
  * \sa `CGAL::Straight_skeleton_3::Arc`
  * \sa `CGAL::Straight_skeleton_3::Sheet`
  */
  class Node
    : public std::enable_shared_from_this<Node>
  {
  public:
    /*!
    * the field type used for coordinates and time values.
    */
    using FT = typename Geom_traits::FT;

    /*!
    * the point type used for the position of the node.
    */
    using Point_3 = typename Geom_traits::Point_3;

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

    /*!
    * returns the position of the node.
    */
    const Point_3& point() const
    {
      return point_;
    }

    void set_point(const Point_3& point)
    {
      this->point_ = point;
    }

    /*!
    * returns the time at which the node was created.
    */
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

    typename std::list<NodeSPtr>::iterator get_list_it() const
    {
      return this->list_it_;
    }

    void set_list_it(typename std::list<NodeSPtr>::iterator list_it)
    {
      this->list_it_ = list_it;
    }

    /*!
    * returns a unique ID of the node.
    */
    int id() const
    {
      return this->id_;
    }

    void set_id(const int id)
    {
      this->id_ = id;
    }

    void add_arc(const ArcSPtr& arc)
    {
      CGAL_SS3_DEBUG_SPTR(arc);
      CGAL_precondition(!has_incident_arc(arc));
      typename std::list<ArcWPtr>::iterator it = arcs_.insert(arcs_.end(), ArcWPtr(arc));
      if (arc->source() == this->shared_from_this()) {
        arc->set_source_list_it(it);
      } else if (arc->has_target()) {
        if (arc->target() == this->shared_from_this()) {
          arc->set_target_list_it(it);
        }
      }
    }

    bool remove_arc(const ArcSPtr& arc)
    {
      CGAL_SS3_DEBUG_SPTR(arc);
      CGAL_precondition(has_incident_arc(arc));
      bool result = false;
      if (arc->source() == this->shared_from_this()) {
        arcs_.erase(arc->get_source_list_it());
        arc->set_source_list_it(typename std::list<ArcWPtr>::iterator());
        result = true;
      } else if (arc->has_target()) {
        if (arc->target() == this->shared_from_this()) {
          arcs_.erase(arc->get_target_list_it());
          arc->set_target_list_it(typename std::list<ArcWPtr>::iterator());
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

    /*!
    * return the list of incident arcs.
    */
    std::list<std::weak_ptr<Arc> >& arcs()
    {
      return this->arcs_;
    }

    /*!
    * return the list of incident sheets.
    */
    std::list<std::weak_ptr<Sheet> >& sheets()
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

    /*!
    * returns combinatorial and geometrical information about the node as a string.
    */
    std::string to_string() const
    {
      std::string result("Node(");
      if (id_ != -1) {
        result += "id=" + String_factory::fromInteger(id_) + ", ";
      } else {
        result += String_factory::fromPointer(this) + ", ";
      }
      result += "<" + String_factory::fromDouble(CGAL::to_double(point().x())) + " ";
      result += String_factory::fromDouble(CGAL::to_double(point().y())) + " ";
      result += String_factory::fromDouble(CGAL::to_double(point().z())) + ">, ";
      result += ", arcs={";
      bool first = true;
      for (const ArcWPtr& arc_wptr : arcs_) {
        if (ArcSPtr arc = arc_wptr.lock()) {
          if (!first) result += ", ";
          result += String_factory::fromInteger(arc->id());
          first = false;
        }
      }
      result += "}";
      result += ", sheets={";
      first = true;
      for (const SheetWPtr& sheet_wptr : sheets_) {
        if (SheetSPtr sheet = sheet_wptr.lock()) {
          if (!first) result += ", ";
          result += String_factory::fromInteger(sheet->id());
          first = false;
        }
      }
      result += "}";

      result += ")";
      return result;
    }

  protected:
    Point_3 point_;
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
  /*!
  * \ingroup PkgStraightSkeleton3Classes
  *
  * The class `Arc` represents an arc of the straight skeleton, i.e., a connection between two nodes.
  *
  * If the skeleton is partial because the propagation was interrupted, some arcs may be unbounded
  * and have only a source node. In this case, the arc stores a direction vector instead
  * of a target node.
  *
  * \sa `CGAL::Straight_skeleton_3`
  * \sa `CGAL::Straight_skeleton_3::Node`
  * \sa `CGAL::Straight_skeleton_3::Sheet`
  */
  class Arc
    : public std::enable_shared_from_this<Arc>
  {
  public:
    /*!
    * the vector type used for the direction of the arc.
    */
    using Vector_3 = typename Geom_traits::Vector_3;
    /*!
    * the line type used for the supporting line of the arc.
    */
    using Line_3 = typename Geom_traits::Line_3;

    using NodeSPtr = std::shared_ptr<Node>;
    using ArcWPtr = std::weak_ptr<Arc>;
    using ArcSPtr = std::shared_ptr<Arc>;
    using SheetWPtr = std::weak_ptr<Sheet>;
    using SheetSPtr = std::shared_ptr<Sheet>;

  public:
    Arc(const NodeSPtr& node_src, const Vector_3& direction)
    {
      node_src_ = node_src;
      direction_ = direction;
      id_ = next_id_++;
    }

    Arc(const NodeSPtr& node_src, const NodeSPtr& node_tgt)
    {
      node_src_ = node_src;
      node_tgt_ = node_tgt;
      id_ = next_id_++;
    }

    ~Arc() {
      node_src_.reset();
      node_tgt_.reset();
      direction_.reset();
      sheets_.clear();
    }

    static ArcSPtr create(const NodeSPtr& node_src, const Vector_3& direction)
    {
      ArcSPtr result = std::make_shared<Arc>(node_src, direction);
      node_src->add_arc(result);
      return result;
    }

    static ArcSPtr create(const NodeSPtr& node_src, const NodeSPtr& node_tgt)
    {
      ArcSPtr result = std::make_shared<Arc>(node_src, node_tgt);
      node_src->add_arc(result);
      node_tgt->add_arc(result);
      return result;
    }

    /*!
    * returns the source node of the arc.
    */
    std::shared_ptr<Node> source() const
    {
      CGAL_SS3_DEBUG_SPTR(node_src_);
      return node_src_;
    }

    void set_source(const NodeSPtr& node_src)
    {
      CGAL_SS3_DEBUG_SPTR(node_src);
      CGAL_precondition(!has_source() && !has_target());
      this->node_src_ = node_src;
    }

    typename std::list<ArcWPtr>::iterator get_source_list_it() const
    {
      return this->node_src_list_it_;
    }

    void set_source_list_it(typename std::list<ArcWPtr>::iterator node_src_list_it)
    {
      this->node_src_list_it_ = node_src_list_it;
    }

    bool has_source() const
    {
      return bool(node_src_);
    }

    /*!
    * returns the target node of the arc.
    * Note that for unbounded arcs, this function returns a null pointer.
    * \sa `has_target()`
    */
    std::shared_ptr<Node> target() const
    {
      CGAL_SS3_DEBUG_SPTR(node_tgt_);
      return node_tgt_;
    }

    void set_target(const NodeSPtr& node_tgt)
    {
      CGAL_SS3_DEBUG_SPTR(node_tgt);
      CGAL_precondition(has_source() && !has_target());
      this->node_tgt_ = node_tgt;
    }

    typename std::list<ArcWPtr>::iterator get_target_list_it() const
    {
      return this->tgt_list_it_;
    }

    void set_target_list_it(typename std::list<ArcWPtr>::iterator node_tgt_list_it)
    {
      this->tgt_list_it_ = node_tgt_list_it;
    }

    /*!
    * returns whether the arc has a target node.
    */
    bool has_target() const
    {
      return bool(node_tgt_);
    }

    NodeSPtr other(const NodeSPtr& node) const
    {
      CGAL_precondition(node == node_src_ || node == node_tgt_);
      if (node == node_src_) {
        return node_tgt_;
      } else {
        return node_src_;
      }
    }

    void close_arc(const NodeSPtr& node_tgt)
    {
      CGAL_SS3_DEBUG_SPTR(node_tgt);
      CGAL_precondition(has_source() && !has_target());
      set_target(node_tgt);
      node_tgt->add_arc(this->shared_from_this());
    }

    bool has_direction() const
    {
      return this->direction_.has_value();
    }

    /*!
    * returns the direction of the arc, i.e., the direction in which the vertex of the polyhedron
    * is moving.
    */
    const Vector_3& direction() const
    {
      CGAL_SS3_DEBUG_SPTR(this->direction_->has_value());
      return *(this->direction_);
    }

    void set_direction(const Vector_3& direction)
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

    typename std::list<ArcSPtr>::iterator get_list_it() const
    {
      return this->list_it_;
    }

    void set_list_it(typename std::list<ArcSPtr>::iterator list_it)
    {
      this->list_it_ = list_it;
    }

    int id() const
    {
      return this->id_;
    }

    void set_id(const int id)
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

    /*!
    * returns the list of incident sheets.
    */
    std::list<std::weak_ptr<Sheet> >& sheets()
    {
      return this->sheets_;
    }

    /*!
    * returns the other arc incident to `node` in the given `sheet`.
    * \pre `sheet` contains this arc, and `node` is either the source or target of this arc.
    */
    std::shared_ptr<Arc> next(const std::shared_ptr<Sheet>& sheet, const std::shared_ptr<Node>& node) const
    {
      CGAL_precondition(has_incident_sheet(sheet));
      CGAL_precondition(node == node_src_ || node == node_tgt_);
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

    /*!
    * returns the 3D line supporting the arc.
    */
    Line_3 line() const
    {
      Line_3 result;
      if (node_tgt_) {
        result = { node_src_->point(), node_tgt_->point() };
      } else if (direction_.has_value()) {
        result = { node_src_->point(), *direction_ };
      }
      return result;
    }

    /*!
    * returns combinatorial and geometrical information about the arc as a string.
    */
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
        result += String_factory::fromInteger(node_src_->id());
      } else {
        result += "-1";
      }
      result += ", tgt=";
      if (node_tgt_) {
        result += String_factory::fromInteger(node_tgt_->id());
      } else {
        result += "-1";
      }
      // Incident sheet IDs
      result += ", sheets={";
      bool first = true;
      for (const SheetWPtr& sheet_wptr : sheets_) {
        if (SheetSPtr sheet = sheet_wptr.lock()) {
          if (!first) result += ", ";
          result += String_factory::fromInteger(sheet->id());
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
    NodeSPtr node_tgt_;
    typename std::list<ArcWPtr>::iterator tgt_list_it_;
    std::optional<Vector_3> direction_;
    typename std::list<SheetWPtr> sheets_; // every arc has 3 sheets
    StraightSkeletonWPtr skel_;

    typename std::list<ArcSPtr>::iterator list_it_;

    int id_;

  private:
    static int next_id_;
  };

public:
  /*!
  * \ingroup PkgStraightSkeleton3Classes
  *
  * The class `Sheet` represents a sheet of the straight skeleton, i.e., a 2D surface
  * bounded by arcs and nodes.
  *
  * A special type of arcs, called *contours*, are also stored in the sheet. Contours are not
  * skeleton arcs per se, but rather correspond to the edges of the input polyhedron,
  * which are needed for the union of arcs to form a closed polygon (when the sheet is bounded).
  *
  * In the case of a partial straight skeleton, the sheet may be unbounded.
  *
  * \sa `CGAL::Straight_skeleton_3`
  * \sa `CGAL::Straight_skeleton_3::Node`
  * \sa `CGAL::Straight_skeleton_3::Sheet`
  */
  class Sheet
    : public std::enable_shared_from_this<Sheet>
  {
  public:
    /*!
    * the plane type used for the supporting plane of the sheet.
    */
    using Plane_3 = typename Geom_traits::Plane_3;

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

    FacetSPtr get_facet_B() const
    {
      CGAL_SS3_DEBUG_SPTR(facet_b_);
      return facet_b_;
    }

    void set_facet_B(const FacetSPtr& facet_b)
    {
      CGAL_SS3_DEBUG_SPTR(facet_b);
      facet_b_ = facet_b;
    }

    FacetSPtr get_facet_F() const
    {
      CGAL_SS3_DEBUG_SPTR(facet_f_);
      return facet_f_;
    }

    void set_facet_F(const FacetSPtr& facet_f)
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

    typename std::list<SheetSPtr>::iterator get_list_it() const
    {
      return this->list_it_;
    }

    void set_list_it(typename std::list<SheetSPtr>::iterator list_it)
    {
      this->list_it_ = list_it;
    }

    int id() const
    {
      return this->id_;
    }

    void set_id(const int id)
    {
      this->id_ = id;
    }

    /*!
    * returns the supporting plane of the sheet.
    */
    const Plane_3& get_plane() const
    {
      return this->plane_;
    }

    void set_plane(const Plane_3& plane)
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

    void add_contour(const ArcSPtr& contour)
    {
      CGAL_SS3_DEBUG_SPTR(contour);
      /*std::list<ArcSPtr>::iterator it =*/ contours_.insert(contours_.end(), contour);
      contour->add_sheet(this->shared_from_this());
    }

    // merge 'sheet' into this sheet.
    void merge(const SheetSPtr& sheet)
    {
      CGAL_SS3_DEBUG_SPTR(sheet);

      CGAL_precondition((this->facet_b_ == sheet->get_facet_B() && this->facet_f_ == sheet->get_facet_F()) ||
                        (this->facet_b_ == sheet->get_facet_F() && this->facet_f_ == sheet->get_facet_B()));

      // copy all arcs that do not yet exist in the sheet
      for (ArcSPtr arc : sheet->arcs()) {
        if (std::find(arcs_.begin(), arcs_.end(), arc) == arcs_.end()) {
          add_arc(arc);
        }
      }

      // copy all contours that do not yet exist in the sheet
      for (ArcSPtr contour : sheet->contours()) {
        if (std::find(contours_.begin(), contours_.end(), contour) == contours_.end()) {
          add_contour(contour);
        }
      }

      // copy all nodes that do not yet exist in the sheet
      for (NodeSPtr node : sheet->nodes()) {
        if (std::find(nodes_.begin(), nodes_.end(), node) == nodes_.end()) {
          add_node(node);
        }
      }
    }

    /*!
    * returns the list of arcs of the sheet.
    */
    std::list<std::shared_ptr<Arc> >& arcs()
    {
      return this->arcs_;
    }

    /*!
    * returns the list of contours of the sheet.
    * Note that this returns a list of arcs, as a sheet might be incident to multiple input edges.
    */
    std::list<std::shared_ptr<Arc> >& contours()
    {
      return this->contours_;
    }

    /*!
    * returns the list of nodes of the sheet.
    */
    std::list<std::shared_ptr<Node> >& nodes()
    {
      return this->nodes_;
    }

    /*!
    * returns combinatorial and geometrical information about the sheet as a string.
    */
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
        result << String_factory::fromInteger(node->id());
        first = false;
      }
      result << "}, ";

      // contour IDs
      result << "contours={";
      first = true;
      for (const ArcSPtr& contour : contours_) {
        if (!first) result << ", ";
        result << String_factory::fromInteger(contour->id());
        first = false;
      }
      result << "}, ";

      // arc IDs
      result << "arcs={";
      first = true;
      for (const ArcSPtr& arc : arcs_) {
        if (!first) result << ", ";
        result << String_factory::fromInteger(arc->id());
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
    Plane_3 plane_;

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
    node->set_list_it(it);
  }

  bool remove_node(NodeSPtr node)
  {
    bool result = false;
    if (node->get_skeleton() == this->shared_from_this()) {
      nodes_.erase(node->get_list_it());
      node->set_skeleton(StraightSkeletonSPtr());
      node->set_list_it(std::list<NodeSPtr>::iterator());
      result = true;
    }
    return result;
  }

  void add_arc(ArcSPtr arc)
  {
    typename std::list<ArcSPtr>::iterator it = arcs_.insert(arcs_.end(), arc);
    arc->set_skeleton(this->shared_from_this());
    arc->set_list_it(it);
  }

  bool remove_arc(ArcSPtr arc)
  {
    bool result = false;
    if (arc->get_skeleton() == this->shared_from_this()) {
      arcs_.erase(arc->get_list_it());
      arc->set_skeleton(StraightSkeletonSPtr());
      arc->set_list_it(typename std::list<ArcSPtr>::iterator());
      result = true;
    }
    return result;
  }

  void add_sheet(SheetSPtr sheet)
  {
    typename std::list<SheetSPtr>::iterator it = sheets_.insert(sheets_.end(), sheet);
    sheet->set_skeleton(this->shared_from_this());
    sheet->set_list_it(it);
  }

  bool remove_sheet(SheetSPtr sheet)
  {
    bool result = false;
    if (sheet->get_skeleton() == this->shared_from_this()) {
      sheets_.erase(sheet->get_list_it());
      sheet->set_skeleton(StraightSkeletonSPtr());
      sheet->set_list_it(typename std::list<SheetSPtr>::iterator());
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
  std::list<std::shared_ptr<Sheet> >& sheets()
  {
    return this->sheets_;
  }

  /*!
   * returns the arcs of the straight skeleton.
   */
  std::list<std::shared_ptr<Arc> >& arcs()
  {
    return this->arcs_;
  }

  /*!
   * returns the nodes of the straight skeleton.
   */
  std::list<std::shared_ptr<Node> >& nodes()
  {
    return this->nodes_;
  }

public:
  void reset_all_IDs()
  {
    typename std::list<SheetSPtr>::iterator it_s = sheets_.begin();
    while (it_s != sheets_.end()) {
      SheetSPtr sheet = *it_s++;
      sheet->set_id(-1);
    }
    typename std::list<ArcSPtr>::iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
      ArcSPtr arc = *it_a++;
      arc->set_id(-1);
    }
    typename std::list<NodeSPtr>::iterator it_n = nodes_.begin();
    while (it_n != nodes_.end()) {
      NodeSPtr node = *it_n++;
      node->set_id(-1);
    }
  }

  bool is_consistent(bool no_unbounded_elements = false) const
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
        int id = node->id();
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
        int id = arc->id();
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
        int id = sheet->id();
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
        CGAL_SS3_SKEL_DS_TRACE("Error: node's skeleton pointer is incorrect");
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

        if (node != arc->source() && node != arc->target()) {
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
        CGAL_SS3_SKEL_DS_TRACE("Error: arc's skeleton pointer is incorrect");
        CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
        result = false;
        break;
      }

      if (!arc->source()) {
        CGAL_SS3_SKEL_DS_TRACE("Error: arc does not have source");
        result = false;
        break;
      }

      if (std::find(nodes_.begin(), nodes_.end(), arc->source()) == nodes_.end()) {
        CGAL_SS3_SKEL_DS_TRACE("Error: arc source is not a node present in the skeleton's node list");
        CGAL_SS3_SKEL_DS_TRACE(arc->source()->to_string());
        result = false;
        break;
      }

      std::list<ArcWPtr> warcs = arc->source()->arcs();
      if (warcs.end() == STL_Extension::internal::weak_find(warcs.begin(), warcs.end(), arc_wptr)) {
        CGAL_SS3_SKEL_DS_TRACE("Error: arc's source node does not have arc");
        CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
        CGAL_SS3_SKEL_DS_TRACE(arc->source()->to_string());
        result = false;
        break;
      }

      if (arc->has_target()) {
        if (std::find(nodes_.begin(), nodes_.end(), arc->target()) == nodes_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc target is not a node present in the skeleton's node list");
          CGAL_SS3_SKEL_DS_TRACE(arc->target()->to_string());
          result = false;
          break;
        }

        warcs = arc->target()->arcs();
        if (warcs.end() == STL_Extension::internal::weak_find(warcs.begin(), warcs.end(), arc_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc's target node does not have arc");
          CGAL_SS3_SKEL_DS_TRACE(arc->to_string());
          CGAL_SS3_SKEL_DS_TRACE(arc->target()->to_string());
          result = false;
          break;
        }
      } else {
        if (no_unbounded_elements) {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc does not have target (full skeleton)");
          result = false;
          break;
        } else {
          if (!arc->has_direction()) {
            CGAL_SS3_SKEL_DS_TRACE("Error: arc has neither target nor direction");
            result = false;
            break;
          }
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
        CGAL_SS3_SKEL_DS_TRACE("Error: sheet's skeleton pointer is incorrect");
        CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
        result = false;
        break;
      }
      if (!polyhedron_->has_facet(sheet->get_facet_B())) {
        CGAL_SS3_SKEL_DS_TRACE("Error: sheet's facet B is not in the polyhedron");
        CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
        result = false;
      }
      if (!polyhedron_->has_facet(sheet->get_facet_F())) {
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
      CGAL::unordered_flat_map<NodeSPtr, int> arc_nodes;

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

        arc_nodes[arc->source()]++;
        if (arc->has_target()) {
          arc_nodes[arc->target()]++;
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
        NodeSPtr node_src = contour->source();
        NodeSPtr node_tgt = contour->target();
        if (!node_src || !node_tgt) {
          CGAL_SS3_SKEL_DS_TRACE("Error: contour does not have two valid incident nodes");
          CGAL_SS3_SKEL_DS_TRACE(contour->to_string());
          result = false;
          break;
        }

        // Check that both nodes appear in the sheet's node list
        const auto& sheet_nodes = sheet->nodes();
        if (std::find(sheet_nodes.begin(), sheet_nodes.end(), node_src) == sheet_nodes.end() ||
            std::find(sheet_nodes.begin(), sheet_nodes.end(), node_tgt) == sheet_nodes.end()) {
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

  // Dump skeleton nodes in an .xyz format
  void dump_nodes(const std::string& filename) const
  {
    std::ofstream nodes_out(filename);
    nodes_out.precision(17);
    for (NodeSPtr node : nodes_) {
      nodes_out << node->point() << "\n";
    }
    nodes_out.close();
  }

  // Dump skeleton arcs in a CGAL polylines format
  void dump_arcs(const std::string& filename,
                 const FT& time = FT(0)) const
  {
    std::ofstream arcs_out(filename);
    arcs_out.precision(17);
    for (ArcSPtr arc : arcs_) {
      arcs_out << "2 ";
      arcs_out << arc->source()->point() << " ";
      if (arc->has_target()) {
        arcs_out << arc->target()->point() << "\n";
      } else {
        std::set<FacetSPtr> incident_faces;
        CGAL_assertion(arc->sheets().size() == 3);
        for (SheetWPtr sheet_wptr : arc->sheets()) {
          if (SheetSPtr sheet = sheet_wptr.lock()) {
            incident_faces.insert(sheet->get_facet_B());
            incident_faces.insert(sheet->get_facet_F());
          }
        }

        CGAL_assertion(incident_faces.size() == 3);

        std::array<Plane_3, 3> offset_planes;
        unsigned int i = 0;
        for (FacetSPtr inc_f : incident_faces) {
          offset_planes[i++] = Geom_utils::offset_plane(inc_f->get_plane(), // static polyhedron's
                                                        Hds_utils::get_speed(inc_f) * time);
        }
        CGAL_postcondition(i == 3);

        std::optional<Point_3> res = Kernel_wrapper::intersection(offset_planes[0], offset_planes[1], offset_planes[2]);
        CGAL_assertion(res.has_value());

        arcs_out << *res << "\n";
      }
    }
    arcs_out.close();
  }

protected:
  PolyhedronSPtr polyhedron_;

  std::list<NodeSPtr> nodes_;
  std::list<ArcSPtr> arcs_;
  std::list<SheetSPtr> sheets_;
};

} // namespace CGAL

template <typename GeomTraits>
int CGAL::Straight_skeleton_3<GeomTraits>::Node::next_id_ = 0;

template <typename GeomTraits>
int CGAL::Straight_skeleton_3<GeomTraits>::Arc::next_id_ = 0;

template <typename GeomTraits>
int CGAL::Straight_skeleton_3<GeomTraits>::Sheet::next_id_ = 0;

#endif /* CGAL_STRAIGHT_SKELETON_3_SDS_STRAIGHT_SKELETON_3_H */
