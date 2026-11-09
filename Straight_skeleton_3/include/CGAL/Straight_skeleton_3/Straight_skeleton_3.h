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
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>
#include <CGAL/Straight_skeleton_3/internal/weak_find.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>

#include <CGAL/unordered_flat_map.h>

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

template <typename GeomTraits>
class Polyhedron_transformation;

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3

/*!
 * \ingroup PkgStraightSkeleton3Classes
 *
 * The class `Straight_skeleton_3` is a data structure providing access
 * to the combinatorial and geometrical information of a 3D straight skeleton.
 *
 * A 3D straight skeleton is composed of nodes, spokes, and sheets.
 *
 * - Nodes are points at which topological changes have occurred during the propagation of the faces
 *   of the polyhedron.
 * - Spokes connect two nodes, and represent the trace of a vertex of the polyhedron during
 *   the propagation of the faces of the polyhedron. The vertex spawned at the source of the spoke,
 *   and disappears at the target of the spoke.
 * - Sheets are 2D surfaces connecting spokes and nodes, representing the trace of the edges
 *   of the polyhedron during the propagation of the faces.
 *
 * If a maximal propagation time was set during construction of the skeleton (this is necessarily
 * the case for outward offsetting), the straight skeleton may be partial, i.e., some spokes
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
  using Transformation = Straight_skeletons_3::internal::algorithm::Polyhedron_transformation<GeomTraits>;

public:
  class Node;
  class Spoke;
  class Sheet;

public:
  /*!
  * \ingroup PkgStraightSkeleton3Classes
  *
  * The class `Node` represents a node of the straight skeleton, i.e., a combinatorial entity
  * representing a point in 3D space at which an event has occurred during the propagation
  * of the faces of the polyhedron.
  *
  * A node stores its 3D position, the time at which it was created, and the incident spokes and sheets.
  *
  * \sa `CGAL::Straight_skeleton_3`
  * \sa `CGAL::Straight_skeleton_3::Spoke`
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
    using SpokeWPtr = std::weak_ptr<Spoke>;
    using SpokeSPtr = std::shared_ptr<Spoke>;
    using SheetWPtr = std::weak_ptr<Sheet>;
    using SheetSPtr = std::shared_ptr<Sheet>;

  public:
    Node() {
      id_ = next_node_id_++;
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

    void add_spoke(const SpokeSPtr& spoke)
    {
      CGAL_SS3_DEBUG_SPTR(spoke);
      CGAL_precondition(!has_incident_spoke(spoke));
      typename std::list<SpokeWPtr>::iterator it = spokes_.insert(spokes_.end(), SpokeWPtr(spoke));
      if (spoke->source() == this->shared_from_this()) {
        spoke->set_source_list_it(it);
      } else if (spoke->has_target()) {
        if (spoke->target() == this->shared_from_this()) {
          spoke->set_target_list_it(it);
        }
      }
    }

    bool remove_spoke(const SpokeSPtr& spoke)
    {
      CGAL_SS3_DEBUG_SPTR(spoke);
      CGAL_precondition(has_incident_spoke(spoke));
      bool result = false;
      if (spoke->source() == this->shared_from_this()) {
        spokes_.erase(spoke->get_source_list_it());
        spoke->set_source_list_it(typename std::list<SpokeWPtr>::iterator());
        result = true;
      } else if (spoke->has_target()) {
        if (spoke->target() == this->shared_from_this()) {
          spokes_.erase(spoke->get_target_list_it());
          spoke->set_target_list_it(typename std::list<SpokeWPtr>::iterator());
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

    bool has_incident_spoke(const SpokeSPtr& spoke) const
    {
      CGAL_SS3_DEBUG_SPTR(spoke);
      SpokeWPtr spoke_wptr = SpokeWPtr(spoke);
      bool result = (spokes_.end() != STL_Extension::internal::weak_find(spokes_.begin(), spokes_.end(), spoke_wptr));
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
      typename std::list<SpokeWPtr>::iterator it_a = spokes_.begin();
      while (it_a != spokes_.end()) {
        SpokeWPtr spoke_wptr = *it_a++;
        if (SpokeSPtr spoke = spoke_wptr.lock()) {
          remove_spoke(spoke);
        }
      }
      spokes_.clear();
    }

    /*!
    * return the list of incident spokes.
    */
    std::list<std::weak_ptr<Spoke> >& spokes()
    {
      return this->spokes_;
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
      typename std::list<SpokeWPtr>::const_iterator it_a = spokes_.begin();
      while (it_a != spokes_.end()) {
        SpokeWPtr spoke_wptr = *it_a++;
        if (!spoke_wptr.expired()) {
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
      result += ", spokes={";
      bool first = true;
      for (const SpokeWPtr& spoke_wptr : spokes_) {
        if (SpokeSPtr spoke = spoke_wptr.lock()) {
          if (!first) result += ", ";
          result += String_factory::fromInteger(spoke->id());
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

    std::list<SpokeWPtr> spokes_;
    std::list<SheetWPtr> sheets_;
    StraightSkeletonWPtr skel_;

    typename std::list<NodeSPtr>::iterator list_it_;

  private:
    static int next_node_id_;
  };

public:
  /*!
  * \ingroup PkgStraightSkeleton3Classes
  *
  * The class `Spoke` represents a spoke of the straight skeleton, i.e., a connection between two nodes.
  *
  * If the skeleton is partial because the propagation was interrupted, some spokes may be unbounded
  * and have only a source node. In this case, the spoke stores a direction vector instead
  * of a target node.
  *
  * \sa `CGAL::Straight_skeleton_3`
  * \sa `CGAL::Straight_skeleton_3::Node`
  * \sa `CGAL::Straight_skeleton_3::Sheet`
  */
  class Spoke
    : public std::enable_shared_from_this<Spoke>
  {
  public:
    /*!
    * the vector type used for the direction of the spoke.
    */
    using Vector_3 = typename Geom_traits::Vector_3;
    /*!
    * the line type used for the supporting line of the spoke.
    */
    using Line_3 = typename Geom_traits::Line_3;

    using NodeSPtr = std::shared_ptr<Node>;
    using SpokeWPtr = std::weak_ptr<Spoke>;
    using SpokeSPtr = std::shared_ptr<Spoke>;
    using SheetWPtr = std::weak_ptr<Sheet>;
    using SheetSPtr = std::shared_ptr<Sheet>;

  public:
    Spoke(const NodeSPtr& node_src, const Vector_3& direction)
    {
      node_src_ = node_src;
      direction_ = direction;
      id_ = next_spoke_id_++;
    }

    Spoke(const NodeSPtr& node_src, const NodeSPtr& node_tgt)
    {
      node_src_ = node_src;
      node_tgt_ = node_tgt;
      id_ = next_spoke_id_++;
    }

    ~Spoke() {
      node_src_.reset();
      node_tgt_.reset();
      direction_.reset();
      sheets_.clear();
    }

    static SpokeSPtr create(const NodeSPtr& node_src, const Vector_3& direction)
    {
      SpokeSPtr result = std::make_shared<Spoke>(node_src, direction);
      node_src->add_spoke(result);
      return result;
    }

    static SpokeSPtr create(const NodeSPtr& node_src, const NodeSPtr& node_tgt)
    {
      SpokeSPtr result = std::make_shared<Spoke>(node_src, node_tgt);
      node_src->add_spoke(result);
      node_tgt->add_spoke(result);
      return result;
    }

    /*!
    * returns the source node of the spoke.
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

    typename std::list<SpokeWPtr>::iterator get_source_list_it() const
    {
      return this->node_src_list_it_;
    }

    void set_source_list_it(typename std::list<SpokeWPtr>::iterator node_src_list_it)
    {
      this->node_src_list_it_ = node_src_list_it;
    }

    bool has_source() const
    {
      return bool(node_src_);
    }

    /*!
    * returns the target node of the spoke.
    * Note that for unbounded spokes, this function returns a null pointer.
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

    typename std::list<SpokeWPtr>::iterator get_target_list_it() const
    {
      return this->tgt_list_it_;
    }

    void set_target_list_it(typename std::list<SpokeWPtr>::iterator node_tgt_list_it)
    {
      this->tgt_list_it_ = node_tgt_list_it;
    }

    /*!
    * returns whether the spoke has a target node.
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

    void close_spoke(const NodeSPtr& node_tgt)
    {
      CGAL_SS3_DEBUG_SPTR(node_tgt);
      CGAL_precondition(has_source() && !has_target());
      set_target(node_tgt);
      node_tgt->add_spoke(this->shared_from_this());
    }

    bool has_direction() const
    {
      return this->direction_.has_value();
    }

    /*!
    * returns the direction of the spoke, i.e., the direction in which the vertex of the polyhedron
    * is moving.
    */
    const Vector_3& direction() const
    {
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

    typename std::list<SpokeSPtr>::iterator get_list_it() const
    {
      return this->list_it_;
    }

    void set_list_it(typename std::list<SpokeSPtr>::iterator list_it)
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
    * returns the other spoke incident to `node` in the given `sheet`.
    * \pre `sheet` contains this spoke, and `node` is either the source or target of this spoke.
    */
    std::shared_ptr<Spoke> next(const std::shared_ptr<Sheet>& sheet, const std::shared_ptr<Node>& node) const
    {
      CGAL_precondition(has_incident_sheet(sheet));
      CGAL_precondition(node == node_src_ || node == node_tgt_);
      SpokeSPtr result = SpokeSPtr();
      // check the spokes incident to the node
      std::list<SpokeWPtr> wspokes = node->spokes();
      typename std::list<SpokeWPtr>::const_iterator it = wspokes.begin();
      while (it != wspokes.end()) {
        SpokeWPtr spoke_wptr = *it++;
        CGAL_SS3_DEBUG_WPTR(spoke_wptr);
        if (SpokeSPtr spoke = spoke_wptr.lock()) {
          if (spoke == this->shared_from_this()) {
            continue;
          }
          if (spoke->has_incident_sheet(sheet)) {
            result = spoke;
            break;
          }
        }
      }
      return result;
    }

    /*!
    * returns the 3D line supporting the spoke.
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
    * returns combinatorial and geometrical information about the spoke as a string.
    */
    std::string to_string() const
    {
      std::string result("Spoke(");
      // Spoke ID
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
    typename std::list<SpokeWPtr>::iterator node_src_list_it_;
    NodeSPtr node_tgt_;
    typename std::list<SpokeWPtr>::iterator tgt_list_it_;
    std::optional<Vector_3> direction_;
    typename std::list<SheetWPtr> sheets_; // every spoke has 3 sheets
    StraightSkeletonWPtr skel_;

    typename std::list<SpokeSPtr>::iterator list_it_;

    int id_;

  private:
    static int next_spoke_id_;
  };

public:
  /*!
  * \ingroup PkgStraightSkeleton3Classes
  *
  * The class `Sheet` represents a sheet of the straight skeleton, i.e., a 2D surface
  * bounded by spokes and nodes.
  *
  * A special type of spokes, called *contours*, are also stored in the sheet. Contours are not
  * skeleton spokes per se, but rather correspond to the edges of the input polyhedron,
  * which are needed for the union of spokes to form a closed polygon (when the sheet is bounded).
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
    using SpokeSPtr = std::shared_ptr<Spoke>;
    using SheetSPtr = std::shared_ptr<Sheet>;

    using EdgeWPtr = typename Polyhedron::EdgeWPtr;
    using FacetSPtr = typename Polyhedron::FacetSPtr;

  public:
    Sheet() {
      id_ = next_sheet_id_++;
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

    void add_spoke(const SpokeSPtr& spoke)
    {
      CGAL_SS3_DEBUG_SPTR(spoke);
      /*std::list<SpokeSPtr>::iterator it =*/ spokes_.insert(spokes_.end(), spoke);
      spoke->add_sheet(this->shared_from_this());
    }

    bool remove_spoke(const SpokeSPtr& spoke)
    {
      CGAL_SS3_DEBUG_SPTR(spoke);
      bool result = false;
      spokes_.remove(spoke);
      result = spoke->remove_sheet(this->shared_from_this());
      return result;
    }

    void add_contour(const SpokeSPtr& contour)
    {
      CGAL_SS3_DEBUG_SPTR(contour);
      /*std::list<SpokeSPtr>::iterator it =*/ contours_.insert(contours_.end(), contour);
      contour->add_sheet(this->shared_from_this());
    }

    // merge 'sheet' into this sheet.
    void merge(const SheetSPtr& sheet)
    {
      CGAL_SS3_DEBUG_SPTR(sheet);

      CGAL_precondition((this->facet_b_ == sheet->get_facet_B() && this->facet_f_ == sheet->get_facet_F()) ||
                        (this->facet_b_ == sheet->get_facet_F() && this->facet_f_ == sheet->get_facet_B()));

      // copy all spokes that do not yet exist in the sheet
      for (SpokeSPtr spoke : sheet->spokes()) {
        if (std::find(spokes_.begin(), spokes_.end(), spoke) == spokes_.end()) {
          add_spoke(spoke);
        }
      }

      // copy all contours that do not yet exist in the sheet
      for (SpokeSPtr contour : sheet->contours()) {
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
    * returns the list of spokes of the sheet.
    */
    std::list<std::shared_ptr<Spoke> >& spokes()
    {
      return this->spokes_;
    }

    /*!
    * returns the list of contours of the sheet.
    * Note that this returns a list of spokes, as a sheet might be incident to multiple input edges.
    */
    std::list<std::shared_ptr<Spoke> >& contours()
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
      for (const SpokeSPtr& contour : contours_) {
        if (!first) result << ", ";
        result << String_factory::fromInteger(contour->id());
        first = false;
      }
      result << "}, ";

      // spoke IDs
      result << "spokes={";
      first = true;
      for (const SpokeSPtr& spoke : spokes_) {
        if (!first) result << ", ";
        result << String_factory::fromInteger(spoke->id());
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

    std::list<SpokeSPtr> spokes_;
    // Squatting the spoke type for contours, which are edges of the static polyhedron.
    // Contours are needed to generate a closed polygon (with holes) when drawing sheets.
    std::list<SpokeSPtr> contours_;
    std::list<NodeSPtr> nodes_;
    StraightSkeletonWPtr skel_;

    typename std::list<SheetSPtr>::iterator list_it_;

    // these are edges of the *mutable* polyhedron, whose active movement traces the sheet
    std::list<EdgeWPtr> edges_;

    int id_;
    Plane_3 plane_;

  private:
    static int next_sheet_id_;
  };

public:
  using NodeWPtr = std::weak_ptr<Node>;
  using NodeSPtr = std::shared_ptr<Node>;
  using SpokeWPtr = std::weak_ptr<Spoke>;
  using SpokeSPtr = std::shared_ptr<Spoke>;
  using SheetWPtr = std::weak_ptr<Sheet>;
  using SheetSPtr = std::shared_ptr<Sheet>;

public:
  ~Straight_skeleton_3()
  {
    polyhedron_.reset();
    sheets_.clear();
    spokes_.clear();
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

  void add_spoke(SpokeSPtr spoke)
  {
    typename std::list<SpokeSPtr>::iterator it = spokes_.insert(spokes_.end(), spoke);
    spoke->set_skeleton(this->shared_from_this());
    spoke->set_list_it(it);
  }

  bool remove_spoke(SpokeSPtr spoke)
  {
    bool result = false;
    if (spoke->get_skeleton() == this->shared_from_this()) {
      spokes_.erase(spoke->get_list_it());
      spoke->set_skeleton(StraightSkeletonSPtr());
      spoke->set_list_it(typename std::list<SpokeSPtr>::iterator());
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
    for (SpokeSPtr spoke : sheet_from->spokes()) {
      spoke->remove_sheet(sheet_from);
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
   * returns the spokes of the straight skeleton.
   */
  std::list<std::shared_ptr<Spoke> >& spokes()
  {
    return this->spokes_;
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
    typename std::list<SpokeSPtr>::iterator it_a = spokes_.begin();
    while (it_a != spokes_.end()) {
      SpokeSPtr spoke = *it_a++;
      spoke->set_id(-1);
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
      CGAL::unordered_flat_set<NodeSPtr> node_set;
      for (const NodeSPtr& node : nodes_) {
        if (!node_set.insert(node).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate node in nodes_ list");
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          result = false;
          break;
        }
      }
      CGAL::unordered_flat_set<SpokeSPtr> spoke_set;
      for (const SpokeSPtr& spoke : spokes_) {
        if (!spoke_set.insert(spoke).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate spoke in spokes_ list");
          CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
          result = false;
          break;
        }
      }
      CGAL::unordered_flat_set<SheetSPtr> sheet_set;
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
      CGAL::unordered_flat_set<int> node_ids;
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
    // Check uniqueness of Spoke IDs
    {
      CGAL::unordered_flat_set<int> spoke_ids;
      for (const SpokeSPtr& spoke : spokes_) {
        int id = spoke->id();
        if (!spoke_ids.insert(id).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate Spoke ID in spokes_ list: " << id);
          CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
          result = false;
          break;
        }
      }
    }
    // Check uniqueness of Sheet IDs
    {
      CGAL::unordered_flat_set<int> sheet_ids;
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
      typename std::list<SpokeWPtr>::const_iterator it_a = node->spokes().begin();
      while (it_a != node->spokes().end()) {
        SpokeWPtr spoke_wptr = *it_a++;
        CGAL_SS3_DEBUG_WPTR(spoke_wptr);
        SpokeSPtr spoke = spoke_wptr.lock();
        if (spoke_wptr.expired() || !spoke) {
          CGAL_SS3_SKEL_DS_TRACE("Error: invalid SpokeWPtr (expired or cannot lock) in node's spoke list");
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          result = false;
          break;
        }

        if (std::find(spokes_.begin(), spokes_.end(), spoke) == spokes_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: node references a spoke not present in the skeleton's spoke list");
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
          result = false;
          break;
        }

        if (node != spoke->source() && node != spoke->target()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: incident spoke does not have node as endpoint");
          CGAL_SS3_SKEL_DS_TRACE(node->to_string());
          CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
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

    typename std::list<SpokeSPtr>::const_iterator it_a = spokes_.begin();
    while (it_a != spokes_.end()) {
      SpokeSPtr spoke = *it_a++;
      CGAL_SS3_DEBUG_SPTR(spoke);
      SpokeWPtr spoke_wptr(spoke);

      if (spoke->get_skeleton() != this->shared_from_this()) {
        CGAL_SS3_SKEL_DS_TRACE("Error: spoke's skeleton pointer is incorrect");
        CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
        result = false;
        break;
      }

      if (!spoke->source()) {
        CGAL_SS3_SKEL_DS_TRACE("Error: spoke does not have source");
        result = false;
        break;
      }

      if (std::find(nodes_.begin(), nodes_.end(), spoke->source()) == nodes_.end()) {
        CGAL_SS3_SKEL_DS_TRACE("Error: spoke source is not a node present in the skeleton's node list");
        CGAL_SS3_SKEL_DS_TRACE(spoke->source()->to_string());
        result = false;
        break;
      }

      std::list<SpokeWPtr> wspokes = spoke->source()->spokes();
      if (wspokes.end() == STL_Extension::internal::weak_find(wspokes.begin(), wspokes.end(), spoke_wptr)) {
        CGAL_SS3_SKEL_DS_TRACE("Error: spoke's source node does not have spoke");
        CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
        CGAL_SS3_SKEL_DS_TRACE(spoke->source()->to_string());
        result = false;
        break;
      }

      if (spoke->has_target()) {
        if (std::find(nodes_.begin(), nodes_.end(), spoke->target()) == nodes_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: spoke target is not a node present in the skeleton's node list");
          CGAL_SS3_SKEL_DS_TRACE(spoke->target()->to_string());
          result = false;
          break;
        }

        wspokes = spoke->target()->spokes();
        if (wspokes.end() == STL_Extension::internal::weak_find(wspokes.begin(), wspokes.end(), spoke_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: spoke's target node does not have spoke");
          CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
          CGAL_SS3_SKEL_DS_TRACE(spoke->target()->to_string());
          result = false;
          break;
        }
      } else {
        if (no_unbounded_elements) {
          CGAL_SS3_SKEL_DS_TRACE("Error: spoke does not have target (full skeleton)");
          result = false;
          break;
        } else {
          if (!spoke->has_direction()) {
            CGAL_SS3_SKEL_DS_TRACE("Error: spoke has neither target nor direction");
            result = false;
            break;
          }
        }
      }
      unsigned int num_sheets = 0;
      typename std::list<SheetWPtr>::const_iterator it_s = spoke->sheets().begin();
      while (it_s != spoke->sheets().end()) {
        SheetWPtr sheet_wptr = *it_s++;
        CGAL_SS3_DEBUG_WPTR(sheet_wptr);
        SheetSPtr sheet = sheet_wptr.lock();
        if (sheet_wptr.expired() || !sheet) {
          CGAL_SS3_SKEL_DS_TRACE("Error: invalid SheetWPtr (expired or cannot lock) in spoke's sheet list");
          CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
          result = false;
          break;
        }

        if (std::find(sheets_.begin(), sheets_.end(), sheet) == sheets_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: spoke references a sheet that is not present in the skeleton's sheet list");
          CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          result = false;
          break;
        }

        ++num_sheets;

        std::list<SpokeSPtr> spokes = sheet->spokes();
        if (spokes.end() == std::find(spokes.begin(), spokes.end(), spoke)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: spoke's sheet does not have spoke");
          CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          result = false;
          break;
        }
      }

      if (num_sheets != 3) {
        CGAL_SS3_SKEL_DS_TRACE("Error: Spoke does not have 3 sheets.");
        CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
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

      // check that a node only has two incident spokes at most in the same sheet
      CGAL::unordered_flat_map<NodeSPtr, int> spoke_nodes;

      typename std::list<SpokeSPtr>::const_iterator it_a = sheet->spokes().begin();
      while (it_a != sheet->spokes().end()) {
        SpokeSPtr spoke = *it_a++;
        CGAL_SS3_DEBUG_SPTR(spoke);

        if (std::find(spokes_.begin(), spokes_.end(), spoke) == spokes_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: sheet references a spoke that is not present in the skeleton's spoke list");
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
          result = false;
          break;
        }

        std::list<SheetWPtr> wsheets = spoke->sheets();
        if (wsheets.end() == STL_Extension::internal::weak_find(wsheets.begin(), wsheets.end(), sheet_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: sheet's spoke does not have sheet");
          CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
          CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
          result = false;
          break;
        }

        spoke_nodes[spoke->source()]++;
        if (spoke->has_target()) {
          spoke_nodes[spoke->target()]++;
        }
        for (auto [node, count] : spoke_nodes) {
          if (count > 2) {
            CGAL_SS3_SKEL_DS_TRACE("Error: sheet has a node with more than two incident spokes in the sheet");
            CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
            CGAL_SS3_SKEL_DS_TRACE(node->to_string());
            result = false;
            break;
          }
          CGAL_SS3_DEBUG_SPTR(node);
          std::list<NodeSPtr> nodes = sheet->nodes();
          if (nodes.end() == std::find(nodes.begin(), nodes.end(), node)) {
            CGAL_SS3_SKEL_DS_TRACE("Error: Nodes of spoke of sheet should be nodes of sheet");
            CGAL_SS3_SKEL_DS_TRACE(node->to_string());
            CGAL_SS3_SKEL_DS_TRACE(spoke->to_string());
            CGAL_SS3_SKEL_DS_TRACE(sheet->to_string());
            result = false;
            break;
          }
        }
      }

      // Contour checks
      for (const SpokeSPtr& contour : sheet->contours()) {
        // Check that contour has exactly one incident sheet (the current sheet)
        int sheet_count = 0;
        for (const SheetWPtr& contour_sheet_wptr : contour->sheets()) {
          SheetSPtr contour_sheet = contour_sheet_wptr.lock();
          if (contour_sheet == sheet) {
            ++sheet_count;
          }
        }
        if (sheet_count != 1) {
          CGAL_SS3_SKEL_DS_TRACE("Error: contour spoke does not have exactly one incident sheet (should be the current sheet)");
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
    sstr << "  Spokes:   " << spokes_.size() << std::endl;
    for (const SpokeSPtr& spoke : spokes_) {
      sstr << "    " << spoke->to_string() << "\n";
    }
    sstr << "  Sheets: " << sheets_.size() << std::endl;
    for (const SheetSPtr& sheet : sheets_) {
      sstr << "    " << sheet->to_string() << "\n";
    }
    return sstr.str();
  }

  // dumps skeleton nodes in an .xyz format
  void dump_nodes(const std::string& filename,
                  const FT& time) const
  {
    std::ofstream nodes_out(filename);
    nodes_out.precision(17);
    for (NodeSPtr node : nodes_) {
      if (node->time() >= time) {
        nodes_out << node->point() << "\n";
      }
    }
    nodes_out.close();
  }

  // dumps skeleton spokes in a CGAL polylines format
  // If the skeleton isn't complete, partial spokes are drawn using the current direction
  // of the subtending trisector, and, if `time` is chosen too large, might overshoot
  // the real spoke of the complete skeleton.
  void dump_spokes(const std::string& filename,
                 const FT& time) const
  {
    std::ofstream spokes_out(filename);
    spokes_out.precision(17);
    for (const SpokeSPtr& spoke : spokes_) {
      CGAL_assertion(spoke->has_source());
      if (spoke->source()->time() <= time) {
        continue;
      }

      spokes_out << "2 ";
      spokes_out << spoke->source()->point() << " ";
      if (spoke->has_target() && spoke->target()->time() >= time) {
        spokes_out << spoke->target()->point() << "\n";
      } else {
        CGAL::unordered_flat_set<FacetSPtr> incident_faces;
        CGAL_assertion(spoke->sheets().size() == 3);
        for (SheetWPtr sheet_wptr : spoke->sheets()) {
          if (SheetSPtr sheet = sheet_wptr.lock()) {
            incident_faces.insert(sheet->get_facet_B());
            incident_faces.insert(sheet->get_facet_F());
          }
        }

        CGAL_assertion(incident_faces.size() == 3);

        std::array<Plane_3, 3> offset_planes;
        unsigned int i = 0;
        for (const FacetSPtr& inc_f : incident_faces) {
          // 'inc_f' is a face in the base polyhedron, it lives at t=0
          offset_planes[i++] = Transformation::shift_plane(inc_f, time);
        }

        std::optional<Point_3> res = Kernel_wrapper::intersection(offset_planes[0], offset_planes[1], offset_planes[2]);
        CGAL_assertion(res.has_value());

        spokes_out << *res << "\n";
      }
    }
    spokes_out.close();
  }

protected:
  PolyhedronSPtr polyhedron_;

  std::list<NodeSPtr> nodes_;
  std::list<SpokeSPtr> spokes_;
  std::list<SheetSPtr> sheets_;
};

} // namespace CGAL

template <typename GeomTraits>
int CGAL::Straight_skeleton_3<GeomTraits>::Node::next_node_id_ = 0;

template <typename GeomTraits>
int CGAL::Straight_skeleton_3<GeomTraits>::Spoke::next_spoke_id_ = 0;

template <typename GeomTraits>
int CGAL::Straight_skeleton_3<GeomTraits>::Sheet::next_sheet_id_ = 0;

#endif /* CGAL_STRAIGHT_SKELETON_3_SDS_STRAIGHT_SKELETON_3_H */
