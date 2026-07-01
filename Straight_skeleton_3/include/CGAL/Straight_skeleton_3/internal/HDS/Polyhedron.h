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
 * file   data/3d/Polyhedron.h
 * author Gernot Walzl
 * date   2011-11-26
 */

#ifndef CGAL_STRAIGHT_SKELETON_INTERNAL_HDS_POLYHEDRON_H
#define CGAL_STRAIGHT_SKELETON_INTERNAL_HDS_POLYHEDRON_H

#include <CGAL/license/Straight_skeleton_3.h>

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/internal/weak_find.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>

#include <CGAL/number_utils.h>
#include <CGAL/Kernel/global_functions.h>

#include <algorithm>
#include <fstream>
#include <list>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace CGAL {

template <typename GeomTraits>
class Straight_skeleton_3;

namespace Straight_skeletons_3 {
namespace internal {
namespace HDS {

template <typename GeomTraits>
class Polyhedron
  : public std::enable_shared_from_this<Polyhedron<GeomTraits> >
{
public:
  using Geom_traits = GeomTraits;

  using PolyhedronSPtr = std::shared_ptr<Polyhedron>;
  using PolyhedronWPtr = std::weak_ptr<Polyhedron>;

private:
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Segment_3 = typename Geom_traits::Segment_3;
  using Vector_3 = typename Geom_traits::Vector_3;
  using Line_3 = typename Geom_traits::Line_3;
  using Plane_3 = typename Geom_traits::Plane_3;

public:
  class Vertex;
  class Edge;
  class Facet;

public:
  class Vertex
    : public std::enable_shared_from_this<Vertex>
  {
  public:
    class Vertex_data
    {
      using VertexWPtr = std::weak_ptr<Vertex>;
      using VertexSPtr = std::shared_ptr<Vertex>;

      using VertexDataSPtr = std::shared_ptr<Vertex_data>;

    public:
      Vertex_data() { /*intentionally does nothing*/ }
      virtual ~Vertex_data() { /*intentionally does nothing*/ }

      static VertexDataSPtr create(const VertexSPtr& vertex)
      {
        CGAL_SS3_DEBUG_SPTR(vertex);
        VertexDataSPtr result = std::make_shared<Vertex_data>();
        vertex->set_data(result);
        return result;
      }

      virtual VertexDataSPtr clone() const {
        return std::make_shared<Vertex_data>(*this);
      }
    };

    class Skeleton_vertex_data
      : public Vertex_data
    {
      using VertexWPtr = std::weak_ptr<Vertex>;
      using VertexSPtr = std::shared_ptr<Vertex>;

      using NodeWPtr = typename CGAL::Straight_skeleton_3<Geom_traits>::NodeWPtr;
      using NodeSPtr = typename CGAL::Straight_skeleton_3<Geom_traits>::NodeSPtr;
      using ArcWPtr = typename CGAL::Straight_skeleton_3<Geom_traits>::ArcWPtr;
      using ArcSPtr = typename CGAL::Straight_skeleton_3<Geom_traits>::ArcSPtr;

      using VertexDataSPtr = std::shared_ptr<Vertex_data>;
      using SkelVertexDataSPtr = std::shared_ptr<Skeleton_vertex_data>;

    public:
      Skeleton_vertex_data() { /*intentionally does nothing*/ }
      virtual ~Skeleton_vertex_data() { /*intentionally does nothing*/ }

      static SkelVertexDataSPtr create(const VertexSPtr& vertex)
      {
        CGAL_SS3_DEBUG_SPTR(vertex);
        SkelVertexDataSPtr result = std::make_shared<Skeleton_vertex_data>();
        vertex->set_data(result);
        return result;
      }

      virtual VertexDataSPtr clone() const override {
        return std::make_shared<Skeleton_vertex_data>(*this);
      }

      ArcSPtr get_arc() const
      {
        CGAL_SS3_DEBUG_WPTR(arc_);
        return this->arc_.lock();
      }

      void set_arc(const ArcSPtr& arc)
      {
        CGAL_SS3_DEBUG_SPTR(arc);
        this->arc_ = arc;
      }

      NodeWPtr get_wnode() const
      {
        return this->node_;
      }

      NodeSPtr get_node() const
      {
        CGAL_SS3_DEBUG_WPTR(node_);
        return this->node_.lock();
      }

      void set_wnode(NodeWPtr node)
      {
        this->node_ = node;
      }

      void set_node(const NodeSPtr& node)
      {
        CGAL_SS3_DEBUG_SPTR(node);
        this->node_ = node;
      }

      bool has_final_point() const
      {
        return final_point_.has_value();
      }

      const Point_3& get_final_point() const
      {
        return *(final_point_);
      }

      void set_final_point(const std::optional<Point_3>& point)
      {
        final_point_ = point;
      }

    protected:
      ArcWPtr arc_;
      NodeWPtr node_;
      std::optional<Point_3> final_point_;
    };

  private:
    using VertexSPtr = std::shared_ptr<Vertex>;
    using EdgeWPtr = std::weak_ptr<Edge>;
    using EdgeSPtr = std::shared_ptr<Edge>;
    using FacetWPtr = std::weak_ptr<Facet>;
    using FacetSPtr = std::shared_ptr<Facet>;

    using VertexDataSPtr = std::shared_ptr<Vertex_data>;
    using SkelVertexDataSPtr = std::shared_ptr<Skeleton_vertex_data>;

  public:
    Vertex(const Point_3& point)
      : point_(point), id_(-1)
    { }

    virtual ~Vertex()
    {
      facets_.clear();
      edges_.clear();
    }

    static VertexSPtr create(const Point_3& point)
    {
      CGAL_SS3_DEBUG_SPTR(point);
      return std::make_shared<Vertex>(point);
    }

    VertexSPtr clone() const
    {
      VertexSPtr result = std::make_shared<Vertex>(point());
      CGAL_SS3_DEBUG_SPTR(result);
      result->set_id(id());
      if (has_data()) {
        result->set_data(get_data()->clone());
      }
      return result;
    }

    const Point_3& point() const
    {
      return point_;
    }

    void set_point(const Point_3& point)
    {
      point_ = point;
    }

    std::list<EdgeWPtr>& edges()
    {
      return edges_;
    }

    std::list<FacetWPtr>& facets()
    {
      return facets_;
    }

    PolyhedronSPtr get_polyhedron() const
    {
      return polyhedron_.lock();
    }

    void set_polyhedron(const PolyhedronSPtr& polyhedron)
    {
      polyhedron_ = polyhedron;
    }

    typename std::list<VertexSPtr>::iterator getPolyhedronListIt() const
    {
      return polyhedron_list_it_;
    }

    void setPolyhedronListIt(typename std::list<VertexSPtr>::iterator list_it)
    {
      polyhedron_list_it_ = list_it;
    }

    int id() const
    {
      return id_;
    }

    void set_id(const int id)
    {
      id_ = id;
    }

    VertexDataSPtr get_data() const
    {
      CGAL_SS3_DEBUG_SPTR(this->data_);
      return this->data_;
    }

    void set_data(const VertexDataSPtr& data)
    {
      this->data_ = data;
    }

    bool has_data() const
    {
      bool result = false;
      if (data_) {
        result = true;
      }
      return result;
    }

    unsigned int degree() const
    {
      unsigned int result = 0;
      for (EdgeWPtr edge_wptr : edges_) {
        if (!edge_wptr.expired()) {
          ++result;
        }
      }
      return result;
    }

    void add_edge(const EdgeSPtr& edge)
    {
      EdgeWPtr edge_wptr(edge);
      typename std::list<EdgeWPtr>::iterator it = edges_.insert(edges_.end(), edge_wptr);
      VertexSPtr vertex_src = edge->source();
      VertexSPtr vertex_tgt = edge->target();
      if (vertex_src == this->shared_from_this() && vertex_tgt == this->shared_from_this()) {
        typename std::list<EdgeWPtr>::iterator it_e =
          STL_Extension::internal::weak_find(edges_.begin(), edges_.end(), edge_wptr);
        if (it_e == edge->getVertexSrcListIt()) {
          edge->setVertexTgtListIt(it);
        } else {
          edge->setVertexSrcListIt(it);
        }
      } else if (vertex_src == this->shared_from_this()) {
        edge->setVertexSrcListIt(it);
      } else if (vertex_tgt == this->shared_from_this()) {
        edge->setVertexTgtListIt(it);
      }
    }

    bool remove_edge(const EdgeSPtr& edge)
    {
      bool result = false;
      if (edge->source() == this->shared_from_this()) {
        edges_.erase(edge->getVertexSrcListIt());
        edge->setVertexSrcListIt(typename std::list<EdgeWPtr>::iterator());
        result = true;
      } else if (edge->target() == this->shared_from_this()) {
        edges_.erase(edge->getVertexTgtListIt());
        edge->setVertexTgtListIt(typename std::list<EdgeWPtr>::iterator());
        result = true;
      }
      return result;
    }

    EdgeSPtr first_edge() const
    {
      EdgeSPtr result = EdgeSPtr();
      for (EdgeWPtr edge_wptr : edges_) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          result = edge;
          break;
        }
      }
      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    EdgeSPtr edge(unsigned int i)
    {
      EdgeSPtr edge = first_edge();
      for (unsigned int j = 0; j < i; ++j) {
        edge = edge->next(this->shared_from_this());
      }
      return edge;
    }

    /**
      * Searches for an edge to the given target.
      * The orientation of the edge is ignored.
      */
    EdgeSPtr find_edge(const VertexSPtr& tgt) const
    {
      EdgeSPtr result = EdgeSPtr();
      for (EdgeWPtr edge_wptr : edges_) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge->source().get() == this && edge->target() == tgt) {
            result = edge;
            break;
          }
          if (edge->target().get() == this && edge->source() == tgt) {
            result = edge;
            break;
          }
        }
      }
      return result;
    }

    /**
    * A vertex may be adjacent to the same facet more than once.
    * @deprecated
    */
    EdgeSPtr find_edge(const FacetSPtr& facet) const
    {
      EdgeSPtr result = EdgeSPtr();
      for (EdgeWPtr edge_wptr : edges_) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge->src(facet).get() == this) {
            result = edge;
            break;
          }
        }
      }
      return result;
    }

    void add_facet(const FacetSPtr& facet)
    {
      facets_.insert(facets_.end(), FacetWPtr(facet));
    }

    bool remove_facet(const FacetSPtr& facet)
    {
      bool result = false;
      typename std::list<FacetWPtr>::iterator it = facets_.begin();
      while (it != facets_.end()) {
        typename std::list<FacetWPtr>::iterator it_current = it;
        FacetWPtr facet_wptr = *it++;
        if (facet_wptr.lock() == facet) {
          facets_.erase(it_current);
          result = true;
          break;
        }
      }
      return result;
    }

    FacetSPtr first_facet() const
    {
      FacetSPtr result = FacetSPtr();
      for (FacetWPtr facet_wptr : facets_) {
        if (FacetSPtr facet = facet_wptr.lock()) {
          result = facet;
          break;
        }
      }
      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    FacetSPtr get_facet(unsigned int i)
    {
      FacetSPtr facet = first_facet();
      for (unsigned int j = 0; j < i; ++j) {
        facet = facet->next(this->shared_from_this());
      }
      return facet;
    }

    bool has_incident_edge(const EdgeSPtr& edge) const
    {
      EdgeWPtr edge_wptr = EdgeWPtr(edge);
      bool result = (edges_.end() !=
          STL_Extension::internal::weak_find(edges_.begin(), edges_.end(), edge_wptr));
      return result;
    }

    bool has_incident_facet(const FacetSPtr& facet) const
    {
      FacetWPtr facet_wptr = FacetWPtr(facet);
      bool result = (facets_.end() !=
          STL_Extension::internal::weak_find(facets_.begin(), facets_.end(), facet_wptr));
      return result;
    }

    void sort_edges()
    {
      typename std::list<EdgeWPtr>::iterator it_e = edges_.begin();
      while (it_e != edges_.end()) {
        typename std::list<EdgeWPtr>::iterator it_current = it_e;
        EdgeWPtr edge_wptr = *it_e++;
        if (edge_wptr.expired()) {
          edges_.erase(it_current);
        }
      }
      std::list<EdgeSPtr> tmp;
      typename std::list<EdgeSPtr>::iterator it_e_tmp = tmp.begin();
      if (edges_.size() > 0) {
        EdgeSPtr first = EdgeSPtr();
        EdgeSPtr edge = edges_.front().lock();
        CGAL_SS3_DEBUG_SPTR(edge);
        while (edge != first) {
          if (!first) {
            first = edge;
          }
          tmp.insert(tmp.end(), edge);
          edge = edge->next(this->shared_from_this());
        }
        while (it_e_tmp != tmp.end()) {
          EdgeSPtr edge = *it_e_tmp++;
          remove_edge(edge);
        }
      }
      edges_.clear();
      it_e_tmp = tmp.begin();
      while (it_e_tmp != tmp.end()) {
        EdgeSPtr edge = *it_e_tmp++;
        add_edge(edge);
      }
    }

    void sort_facets()
    {
      typename std::list<FacetWPtr>::iterator it_f = facets_.begin();
      while (it_f != facets_.end()) {
        typename std::list<FacetWPtr>::iterator it_current = it_f;
        FacetWPtr facet_wptr = *it_f++;
        if (facet_wptr.expired()) {
          facets_.erase(it_current);
        }
      }
      std::list<FacetSPtr> tmp;
      typename std::list<FacetSPtr>::iterator it_f_tmp = tmp.begin();
      if (facets_.size() > 0) {
        FacetSPtr facet = FacetSPtr(facets_.front());
        EdgeSPtr edge_first = find_edge(facet);
        EdgeSPtr edge;
        while (edge != edge_first) {
          if (!edge) {
            edge = edge_first;
          }
          tmp.insert(tmp.end(), facet);
          edge = edge->next(this->shared_from_this());
          facet = edge->other(facet);
        }
        while (it_f_tmp != tmp.end()) {
          FacetSPtr facet = *it_f_tmp++;
          remove_facet(facet);
        }
      }
      facets_.clear();
      it_f_tmp = tmp.begin();
      while (it_f_tmp != tmp.end()) {
        FacetSPtr facet = *it_f_tmp++;
        add_facet(facet);
      }
    }

    void sort()
    {
      typename std::list<EdgeWPtr>::iterator it_e = edges_.begin();
      while (it_e != edges_.end()) {
        typename std::list<EdgeWPtr>::iterator it_current = it_e;
        EdgeWPtr edge_wptr = *it_e++;
        if (edge_wptr.expired()) {
          edges_.erase(it_current);
        }
      }
      typename std::list<FacetWPtr>::iterator it_f = facets_.begin();
      while (it_f != facets_.end()) {
        typename std::list<FacetWPtr>::iterator it_current = it_f;
        FacetWPtr facet_wptr = *it_f++;
        if (facet_wptr.expired()) {
          facets_.erase(it_current);
        }
      }
      std::list<EdgeSPtr> edges_tmp;
      typename std::list<EdgeSPtr>::iterator it_e_tmp = edges_tmp.begin();
      std::list<FacetSPtr> facets_tmp;
      typename std::list<FacetSPtr>::iterator it_f_tmp = facets_tmp.begin();
      if (edges_.size() > 0) {
        EdgeSPtr edge_first = edges_.front().lock();
        CGAL_SS3_DEBUG_SPTR(edge_first);
        EdgeSPtr edge;
        FacetSPtr facet = edge_first->get_facet_L();
        if (edge_first->target() == this->shared_from_this()) {
          facet = edge_first->get_facet_R();
        }
        while (edge != edge_first) {
          if (!edge) {
            edge = edge_first;
          }
          edges_tmp.insert(edges_tmp.end(), edge);
          facets_tmp.insert(facets_tmp.end(), facet);
          edge = edge->next(this->shared_from_this());
          facet = edge->other(facet);
        }
        while (it_e_tmp != edges_tmp.end()) {
          EdgeSPtr edge = *it_e_tmp++;
          remove_edge(edge);
        }
        while (it_f_tmp != facets_tmp.end()) {
          FacetSPtr facet = *it_f_tmp++;
          remove_facet(facet);
        }
      }
      edges_.clear();
      it_e_tmp = edges_tmp.begin();
      while (it_e_tmp != edges_tmp.end()) {
        EdgeSPtr edge = *it_e_tmp++;
        add_edge(edge);
      }
      facets_.clear();
      it_f_tmp = facets_tmp.begin();
      while (it_f_tmp != facets_tmp.end()) {
        FacetSPtr facet = *it_f_tmp++;
        add_facet(facet);
      }
    }

    VertexSPtr next(const FacetSPtr& facet) const
    {
      VertexSPtr result = VertexSPtr();
      CGAL_assertion(facet->vertices().size() > 1);

      for (EdgeWPtr edge_wptr : edges_) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge->src(facet) == this->shared_from_this()) {
            result = edge->tgt(facet);
            break;
          }
        }
      }

      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    VertexSPtr prev(const FacetSPtr& facet) const
    {
      VertexSPtr result = VertexSPtr();
      CGAL_assertion(facet->vertices().size() > 1);

      for (EdgeWPtr edge_wptr : edges_) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge->tgt(facet) == this->shared_from_this()) {
            result = edge->src(facet);
            break;
          }
        }
      }

      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    /**
    * An edge will be created.
    * The target vertex will be returned.
    */
    VertexSPtr split(const FacetSPtr& facet_left,
                     FacetSPtr facet_right)
    {
      VertexSPtr result = VertexSPtr();
      if (!has_incident_facet(facet_left) || !has_incident_facet(facet_right)) {
        return result;
      }
    //  if (facet_left->find_edge(facet_right)) {
    //      return result;
    //  }
      result = Vertex::create(point());

      // 1. select edges and facets for result
      std::list<FacetSPtr> facets;
      std::list<EdgeSPtr> edges;
      FacetSPtr poly_curr = facet_right;
      EdgeSPtr edge_curr = find_edge(poly_curr);
      while (poly_curr != facet_left) {
        facets.insert(facets.end(), poly_curr);
        edge_curr = edge_curr->next(this->shared_from_this());
        edges.insert(edges.end(), edge_curr);
        poly_curr = edge_curr->other(poly_curr);
      }
      facets.insert(facets.end(), facet_left);

      // 2. split vertex
      typename std::list<EdgeSPtr>::iterator it_e = edges.begin();
      while (it_e != edges.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->source() == this->shared_from_this()) {
          edge->replace_vertex_src(result);
        } else if (edge->target() == this->shared_from_this()) {
          edge->replace_vertex_tgt(result);
        }
      }
      typename std::list<FacetSPtr>::iterator it_f = facets.begin();
      while (it_f != facets.end()) {
        FacetSPtr facet = *it_f++;
        if (facet != facet_left && facet != facet_right) {
          facet->remove_vertex(this->shared_from_this());
        }
        facet->add_vertex(result);
      }

      // 3. insert connecting edge
      EdgeSPtr edge = Edge::create(this->shared_from_this(), result);
      edge->set_facet_L(facet_left);
      edge->set_facet_R(facet_right);
      facet_left->add_edge(edge);
      facet_right->add_edge(edge);

      if (PolyhedronSPtr polyhedron = polyhedron_.lock()) {
        polyhedron->add_vertex(result);
        polyhedron->add_edge(edge);
      }
      return result;
    }

    std::string to_string() const
    {
      std::string result("Vertex(");
      result += "id=" + IO::String_factory::fromInteger(id_) + ", ";
      // result += "addr=" + IO::String_factory::fromPointer(this) +", ";
      result += "<" + IO::String_factory::fromDouble(CGAL::to_double(point().x())) + " ";
      result += IO::String_factory::fromDouble(CGAL::to_double(point().y())) + " ";
      result += IO::String_factory::fromDouble(CGAL::to_double(point().z())) + ">";
      result += ", Edges:" + IO::String_factory::fromInteger(edges_.size()) + " {";
      for (EdgeWPtr edge_wptr : edges_) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          result += " " + std::to_string(edge->id());
        }
      }
      result += " }";
      result += ", Facets:" + IO::String_factory::fromInteger(facets_.size()) + " {";
      for (FacetWPtr facet_wptr : facets_) {
        if (FacetSPtr facet = facet_wptr.lock()) {
          result += " " + std::to_string(facet->id());
        }
      }
      result += " }";
      result += ")";
      return result;
    }

  protected:
    Point_3 point_;

    std::list<EdgeWPtr> edges_;
    std::list<FacetWPtr> facets_;
    PolyhedronWPtr polyhedron_;
    typename std::list<VertexSPtr>::iterator polyhedron_list_it_;

    int id_;
    VertexDataSPtr data_;
  };

public:
  class Edge
    : public std::enable_shared_from_this<Edge>
  {
  public:
    class Edge_data
    {
      using EdgeWPtr = std::weak_ptr<Edge>;
      using EdgeSPtr = std::shared_ptr<Edge>;

      using EdgeDataSPtr = std::shared_ptr<Edge_data>;

    public:
      Edge_data() { /*intentionally does nothing*/ }
      virtual ~Edge_data() { /*intentionally does nothing*/ }

      static EdgeDataSPtr create(const EdgeSPtr& edge)
      {
        EdgeDataSPtr result = std::make_shared<Edge_data>();
        edge->set_data(result);
        return result;
      }

      virtual EdgeDataSPtr clone() const {
        return std::make_shared<Edge_data>(*this);
      }
    };

    class Skeleton_edge_data
      : public Edge_data
    {
      using EdgeWPtr = std::weak_ptr<Edge>;
      using EdgeSPtr = std::shared_ptr<Edge>;
      using FacetWPtr = std::weak_ptr<Facet>;
      using FacetSPtr = std::shared_ptr<Facet>;

      using SheetWPtr = typename CGAL::Straight_skeleton_3<Geom_traits>::SheetWPtr;
      using SheetSPtr = typename CGAL::Straight_skeleton_3<Geom_traits>::SheetSPtr;

      using EdgeDataSPtr = std::shared_ptr<Edge_data>;
      using SkelEdgeDataSPtr = std::shared_ptr<Skeleton_edge_data>;

    public:
      Skeleton_edge_data() : is_vanish_time_known_(false) { }
      virtual ~Skeleton_edge_data() { /*intentionally does nothing*/ }

      static SkelEdgeDataSPtr create(const EdgeSPtr& edge)
      {
        SkelEdgeDataSPtr result = std::make_shared<Skeleton_edge_data>();
        edge->set_data(result);
        return result;
      }

      virtual EdgeDataSPtr clone() const override {
        return std::make_shared<Skeleton_edge_data>(*this);
      }

      const std::optional<FT>& get_vanish_time()
      {
        CGAL_assertion(is_vanish_time_known_);
        return vanish_time_;
      }

      void set_vanish_time(const std::optional<FT>& vanish_time)
      {
        vanish_time_ = vanish_time;
        is_vanish_time_known_ = true;
      }

      SheetSPtr get_sheet() const
      {
        CGAL_SS3_DEBUG_WPTR(sheet_);
        return this->sheet_.lock();
      }

      void set_sheet(const SheetSPtr& sheet)
      {
        this->sheet_ = sheet;
      }

    protected:
      // The 'optional' + Boolean is a way to handle "unknown", "known and computed",
      // and "known and computed, but outside of bounds"
      std::optional<FT> vanish_time_;
      bool is_vanish_time_known_;
      SheetWPtr sheet_;
    };

  private:
    // using VertexWPtr = std::weak_ptr<Vertex>;
    using VertexSPtr = std::shared_ptr<Vertex>;
    using EdgeWPtr = std::weak_ptr<Edge>;
    using EdgeSPtr = std::shared_ptr<Edge>;
    using FacetWPtr = std::weak_ptr<Facet>;
    using FacetSPtr = std::shared_ptr<Facet>;

    using EdgeDataSPtr = std::shared_ptr<Edge_data>;
    using SkelEdgeDataSPtr = std::shared_ptr<Skeleton_edge_data>;

  public:
    Edge(const VertexSPtr& src, const VertexSPtr& tgt)
      : vertex_src_(src), vertex_tgt_(tgt), id_(-1)
    { }

    Edge(const Edge& edge)
      : vertex_src_(edge.vertex_src_->clone()),
        vertex_tgt_(edge.vertex_tgt_->clone()),
        id_(-1)
    { }

    virtual ~Edge()
    {
      vertex_src_.reset();
      vertex_tgt_.reset();
    }

    static EdgeSPtr create(const VertexSPtr& src, const VertexSPtr& tgt)
    {
      EdgeSPtr result = std::make_shared<Edge>(src, tgt);
      src->add_edge(result);
      tgt->add_edge(result);
      return result;
    }

    EdgeSPtr clone() const
    {
      EdgeSPtr result = std::make_shared<Edge>();
      CGAL_SS3_DEBUG_SPTR(result);
      result->vertex_src_->add_edge(result);
      result->vertex_tgt_->add_edge(result);
      result->set_id(id());
      if (has_data()) {
        result->set_data(get_data()->clone());
      }
      return result;
    }

    VertexSPtr source() const
    {
      CGAL_SS3_DEBUG_SPTR(this->vertex_src_);
      return this->vertex_src_;
    }

    void set_vertex_src(const VertexSPtr& src)
    {
      this->vertex_src_ = src;
    }

    typename std::list<EdgeWPtr>::iterator getVertexSrcListIt() const
    {
      return this->vertex_src_list_it_;
    }

    void setVertexSrcListIt(typename std::list<EdgeWPtr>::iterator list_it)
    {
      this->vertex_src_list_it_ = list_it;
    }

    VertexSPtr target() const
    {
      CGAL_SS3_DEBUG_SPTR(this->vertex_tgt_);
      return this->vertex_tgt_;
    }

    void set_vertex_tgt(const VertexSPtr& tgt)
    {
      this->vertex_tgt_ = tgt;
    }

    typename std::list<EdgeWPtr>::iterator getVertexTgtListIt() const
    {
      return this->vertex_tgt_list_it_;
    }

    void setVertexTgtListIt(typename std::list<EdgeWPtr>::iterator list_it)
    {
      this->vertex_tgt_list_it_ = list_it;
    }

    bool has_vertex(const VertexSPtr& vertex) const
    {
      return (this->vertex_src_ == vertex || this->vertex_tgt_ == vertex);
    }

    FacetSPtr get_facet_L() const
    {
      return this->facet_l_.lock();
    }

    void set_facet_L(const FacetSPtr& facet)
    {
      this->facet_l_ = facet;
      this->cachedReflexStatus_ = std::nullopt;
    }

    typename std::list<EdgeSPtr>::iterator getFacetLListIt() const
    {
      return this->facet_l_list_it_;
    }

    void setFacetLListIt(typename std::list<EdgeSPtr>::iterator list_it)
    {
      this->facet_l_list_it_ = list_it;
    }

    FacetSPtr get_facet_R() const
    {
      return this->facet_r_.lock();
    }

    void set_facet_R(const FacetSPtr& facet)
    {
      this->facet_r_ = facet;
      this->cachedReflexStatus_ = std::nullopt;
    }

    typename std::list<EdgeSPtr>::iterator getFacetRListIt() const
    {
      return this->facet_r_list_it_;
    }

    void setFacetRListIt(typename std::list<EdgeSPtr>::iterator list_it)
    {
      this->facet_r_list_it_ = list_it;
    }

    FacetSPtr get_facet_src() const
    {
      FacetSPtr result = FacetSPtr();
      VertexSPtr vertex_src = source();
      if (vertex_src->degree() == 3) {
        if (get_facet_L()) {
          result = get_facet_L()->next(vertex_src);
        }
      }
      return result;
    }

    FacetSPtr get_facet_tgt() const
    {
      FacetSPtr result = FacetSPtr();
      VertexSPtr vertex_tgt = target();
      if (vertex_tgt->degree() == 3) {
        if (get_facet_R()) {
          result = get_facet_R()->next(vertex_tgt);
        }
      }
      return result;
    }

    PolyhedronSPtr get_polyhedron() const
    {
      return this->polyhedron_.lock();
    }

    void set_polyhedron(const PolyhedronSPtr& polyhedron)
    {
      this->polyhedron_ = polyhedron;
    }

    typename std::list<EdgeSPtr>::iterator getPolyhedronListIt() const
    {
      return this->polyhedron_list_it_;
    }

    void setPolyhedronListIt(typename std::list<EdgeSPtr>::iterator list_it)

    {
      this->polyhedron_list_it_ = list_it;
    }

    int id() const
    {
      return this->id_;
    }

    void set_id(const int id)
    {
      this->id_ = id;
    }

    EdgeDataSPtr get_data() const
    {
      CGAL_SS3_DEBUG_SPTR(this->data_);
      return this->data_;
    }

    void set_data(const EdgeDataSPtr& data)
    {
      this->data_ = data;
    }

    bool has_data() const
    {
      bool result = false;
      if (data_) {
        result = true;
      }
      return result;
    }

    std::optional<bool> reflex_status() const
    {
      return cachedReflexStatus_;
    }

    Segment_3 segment() const
    {
      return { vertex_src_->point(), vertex_tgt_->point() };
    }

    Line_3 line() const
    {
      return { vertex_src_->point(), vertex_tgt_->point() };
    }

    VertexSPtr other(const VertexSPtr& vertex) const
    {
      CGAL_precondition(has_vertex(vertex));
      VertexSPtr result = VertexSPtr();
      if (vertex == vertex_src_) {
        result = vertex_tgt_;
      } else if (vertex == vertex_tgt_) {
        result = vertex_src_;
      }
      return result;
    }

    FacetSPtr other(const FacetSPtr& facet) const
    {
      FacetSPtr result = FacetSPtr();
      if (facet == facet_l_.lock()) {
        result = facet_r_.lock();
      } else if (facet == facet_r_.lock()) {
        result = facet_l_.lock();
      }
      return result;
    }

    VertexSPtr src(const FacetSPtr& facet_l) const
    {
      VertexSPtr result = VertexSPtr();
      if (facet_l == facet_l_.lock()) {
        result = vertex_src_;
      } else if (facet_l == facet_r_.lock()) {
        result = vertex_tgt_;
      }
      return result;
    }

    VertexSPtr tgt(const FacetSPtr& facet_l) const
    {
      VertexSPtr result = VertexSPtr();
      if (facet_l == facet_l_.lock()) {
        result = vertex_tgt_;
      } else if (facet_l == facet_r_.lock()) {
        result = vertex_src_;
      }
      return result;
    }

    FacetSPtr left(const VertexSPtr& vertex_src) const
    {
      FacetSPtr result = FacetSPtr();
      if (vertex_src == vertex_src_) {
        result = facet_l_.lock();
      } else if (vertex_src == vertex_tgt_) {
        result = facet_r_.lock();
      }
      return result;
    }

    FacetSPtr right(const VertexSPtr& vertex_src) const
    {
      FacetSPtr result = FacetSPtr();
      if (vertex_src == vertex_src_) {
        result = facet_r_.lock();
      } else if (vertex_src == vertex_tgt_) {
        result = facet_l_.lock();
      }
      return result;
    }

    EdgeSPtr next(const FacetSPtr& facet) const
    {
      EdgeSPtr result = EdgeSPtr();
      VertexSPtr vertex_tgt = tgt(facet);
      for (EdgeWPtr edge_wptr : vertex_tgt->edges()) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge->src(facet) == vertex_tgt) {
            result = edge;
            break;
          }
        }
      }
      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    EdgeSPtr prev(const FacetSPtr& facet) const
    {
      EdgeSPtr result = EdgeSPtr();
      VertexSPtr vertex_src = src(facet);
      for (EdgeWPtr edge_wptr : vertex_src->edges()) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge->tgt(facet) == vertex_src) {
            result = edge;
            break;
          }
        }
      }
      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    /**
      * counter clockwise from outside
      */
    EdgeSPtr next(const VertexSPtr& vertex) const
    {
      EdgeSPtr result = EdgeSPtr();
      FacetSPtr facet;
      if (vertex == vertex_src_) {
        facet = FacetSPtr(facet_l_);
      } else if (vertex == vertex_tgt_) {
        facet = FacetSPtr(facet_r_);
      }

      CGAL_SS3_DEBUG_SPTR(facet);

      for (EdgeWPtr edge_wptr : vertex->edges()) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge.get() == this) {
            continue;
          }
          if (edge->tgt(facet) == vertex) {
            result = edge;
            break;
          }
        }
      }

      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    EdgeSPtr prev(const VertexSPtr& vertex) const
    {
      EdgeSPtr result = EdgeSPtr();
      FacetSPtr facet;
      if (vertex == vertex_src_) {
        facet = FacetSPtr(facet_r_);
      } else if (vertex == vertex_tgt_) {
        facet = FacetSPtr(facet_l_);
      }

      CGAL_SS3_DEBUG_SPTR(facet);

      for (EdgeWPtr edge_wptr : vertex->edges()) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge.get() == this) {
            continue;
          }
          if (edge->src(facet) == vertex) {
            result = edge;
            break;
          }
        }
      }

      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    /**
    * Swaps src and tgt vertex and left and right facet
    */
    void invert()
    {
      VertexSPtr vertex_tmp_ = vertex_src_;
      typename std::list<EdgeWPtr>::iterator vertex_tmp_list_it_ = vertex_src_list_it_;
      FacetWPtr facet_tmp_ = facet_l_;
      typename std::list<EdgeSPtr>::iterator facet_tmp_list_it_ = facet_l_list_it_;
      vertex_src_ = vertex_tgt_;
      vertex_src_list_it_ = vertex_tgt_list_it_;
      facet_l_ = facet_r_;
      facet_l_list_it_ = facet_r_list_it_;
      vertex_tgt_ = vertex_tmp_;
      vertex_tgt_list_it_ = vertex_tmp_list_it_;
      facet_r_ = facet_tmp_;
      facet_r_list_it_ = facet_tmp_list_it_;
    }

    /**
      * Splits this edge into 2 edges.
      * The target vertex of this edge is set to the given middle vertex.
      * It returns newly created edge.
      */
    EdgeSPtr split(const VertexSPtr& middle)
    {
      EdgeSPtr result = Edge::create(middle, vertex_tgt_);
      if (FacetSPtr facet_l = facet_l_.lock()) {
        result->set_facet_L(facet_l);
        facet_l->add_edge(result);
      }
      if (FacetSPtr facet_r = facet_r_.lock()) {
        result->set_facet_R(facet_r);
        facet_r->add_edge(result);
      }
      vertex_tgt_->remove_edge(this->shared_from_this());
      vertex_tgt_ = middle;
      middle->add_edge(this->shared_from_this());
      if (PolyhedronSPtr polyhedron = polyhedron_.lock()) {
        polyhedron->add_edge(result);
      }
      return result;
    }

    /**
      * More than just a simple set method.
      */
    void replace_vertex_src(const VertexSPtr& vertex_src)
    {
      vertex_src_->remove_edge(this->shared_from_this());
      vertex_src_ = vertex_src;
      vertex_src->add_edge(this->shared_from_this());
    }

    void replace_vertex_tgt(const VertexSPtr& vertex_tgt)
    {
      vertex_tgt_->remove_edge(this->shared_from_this());
      vertex_tgt_ = vertex_tgt;
      vertex_tgt->add_edge(this->shared_from_this());
    }

    void replace_facet_L(const FacetSPtr& facet_l)
    {
      if (FacetSPtr facet = facet_l_.lock()) {
        facet->remove_edge(this->shared_from_this());
      }
      set_facet_L(facet_l);
      facet_l->add_edge(this->shared_from_this());
    }

    void replace_facet_R(const FacetSPtr& facet_r)
    {
      if (FacetSPtr facet = facet_r_.lock()) {
        facet->remove_edge(this->shared_from_this());
      }
      set_facet_R(facet_r);
      facet_r->add_edge(this->shared_from_this());
    }

    bool has_same_facets(const EdgeSPtr& edge) const
    {
      bool result = (facet_l_.lock() == edge->get_facet_L() &&
                     facet_r_.lock() == edge->get_facet_R()) ||
                    (facet_r_.lock() == edge->get_facet_L() &&
                     facet_l_.lock() == edge->get_facet_R());
      return result;
    }

    bool is_reflex() const
    {
      CGAL_precondition(vertex_src_->point() != vertex_tgt_->point());
      if (cachedReflexStatus_) {
        return *cachedReflexStatus_;
      }
      bool result = false;
      FacetSPtr facet_l = this->get_facet_L();
      FacetSPtr facet_r = this->get_facet_R();
      CGAL_SS3_DEBUG_SPTR(facet_l);
      CGAL_SS3_DEBUG_SPTR(facet_r);
      const Plane_3& plane_l = facet_l->get_plane();
      const Vector_3 normal_l = plane_l.orthogonal_vector();
      CGAL_assertion(normal_l != CGAL::NULL_VECTOR);
      const Vector_3 dir = line().to_vector();
      CGAL_assertion(dir != CGAL::NULL_VECTOR);
      const Point_3& p_src = vertex_src_->point();
      Point_3 p = p_src + CGAL::cross_product(normal_l, dir);
      const Plane_3& plane_r = facet_r->get_plane();
      if (plane_r.oriented_side(p) == CGAL::ON_POSITIVE_SIDE) {
        result = true;
      }
      cachedReflexStatus_ = result;
      return result;
    }

    std::string to_string() const
    {
      std::string result("Edge(");
      result += "id=" + IO::String_factory::fromInteger(id_);
      // result += ", addr=" + IO::String_factory::fromPointer(this);
      result += ", l=";
      if (FacetSPtr facet_l = get_facet_L()) {
          if (facet_l->id() != -1) {
            result += IO::String_factory::fromInteger(get_facet_L()->id());
          } else {
            // result += IO::String_factory::fromPointer(get_facet_L().get());
          }
      } else {
        result += "expired";
      }
      result += ", r=";
      if (FacetSPtr facet_r = get_facet_R()) {
          if (facet_r->id() != -1) {
            result += IO::String_factory::fromInteger(get_facet_R()->id());
          } else {
            // result += IO::String_factory::fromPointer(get_facet_R().get());
          }
      } else {
        result += "expired";
      }
      // src and tgt faces
      result += ", s=";
      if (get_facet_src()) {
          if (get_facet_src()->id() != -1) {
            result += IO::String_factory::fromInteger(get_facet_src()->id());
          } else {
            // result += IO::String_factory::fromPointer(get_facet_src().get());
          }
      } else {
        result += "nullptr";
      }
      result += ", d=";
      if (get_facet_tgt()) {
          if (get_facet_tgt()->id() != -1) {
            result += IO::String_factory::fromInteger(get_facet_tgt()->id());
          } else {
            // result += IO::String_factory::fromPointer(get_facet_tgt().get());
          }
      } else {
        result += "nullptr";
      }
      if (vertex_src_) {
        result += ",\n     src=" + vertex_src_->to_string();
      }
      if (vertex_tgt_) {
        result += ",\n     tgt=" + vertex_tgt_->to_string();
      }
      result += ")";
      return result;
    }

  protected:
    VertexSPtr vertex_src_;
    typename std::list<EdgeWPtr>::iterator vertex_src_list_it_;
    VertexSPtr vertex_tgt_;
    typename std::list<EdgeWPtr>::iterator vertex_tgt_list_it_;
    FacetWPtr facet_l_;
    typename std::list<EdgeSPtr>::iterator facet_l_list_it_;
    FacetWPtr facet_r_;
    typename std::list<EdgeSPtr>::iterator facet_r_list_it_;
    PolyhedronWPtr polyhedron_;
    typename std::list<EdgeSPtr>::iterator polyhedron_list_it_;

    int id_;
    mutable std::optional<bool> cachedReflexStatus_;
    EdgeDataSPtr data_;
  };

public:
  class Facet
    : public std::enable_shared_from_this<Facet>
  {
    friend class Polyhedron;

  public:
    class Facet_data
    {
      using FacetWPtr = std::weak_ptr<Facet>;
      using FacetSPtr = std::shared_ptr<Facet>;

      using FacetDataSPtr = std::shared_ptr<Facet_data>;

    public:
      Facet_data() { /*intentionally does nothing*/ }
      virtual ~Facet_data() { /*intentionally does nothing*/ }

      static FacetDataSPtr create(const FacetSPtr& facet)
      {
        FacetDataSPtr result = std::make_shared<Facet_data>();
        facet->set_data(result);
        return result;
      }

      virtual FacetDataSPtr clone(const FacetSPtr&) const {
        return std::make_shared<Facet_data>(*this);
      }
    };

    class Skeleton_facet_data
      : public Facet_data
    {
      using FacetWPtr = std::weak_ptr<Facet>;
      using FacetSPtr = std::shared_ptr<Facet>;

      using FacetDataSPtr = std::shared_ptr<Facet_data>;
      using SkelFacetDataSPtr = std::shared_ptr<Skeleton_facet_data>;

    public:
      Skeleton_facet_data() { /*intentionally does nothing*/ }
      virtual ~Skeleton_facet_data() { /*intentionally does nothing*/ }

      static SkelFacetDataSPtr create(const FacetSPtr& facet)
      {
        SkelFacetDataSPtr result = std::make_shared<Skeleton_facet_data>();
        result->set_facet_origin(facet);
        facet->set_data(result);
        return result;
      }

      virtual FacetDataSPtr clone(const FacetSPtr& origin) const override {
        SkelFacetDataSPtr result = std::make_shared<Skeleton_facet_data>(*this);
        result->set_facet_origin(origin);
        return result;
      }

      FacetSPtr get_facet_origin() const
      {
        CGAL_SS3_DEBUG_WPTR(facet_origin_);
        return this->facet_origin_.lock();
      }

      void set_facet_origin(const FacetSPtr& facet_origin)
      {
        this->facet_origin_ = facet_origin;
      }

      std::size_t get_input_face_id() const
      {
        return this->input_face_id_;
      }

      void set_input_face_id(const std::size_t input_face_id)
      {
        this->input_face_id_ = input_face_id;
      }

      const FT& get_speed() const
      {
        CGAL_assertion(speed_ >= 0);
        return speed_;
      }

      void set_speed(const FT& speed)
      {
        speed_ = speed;
      }

      const Plane_3& get_base_plane() const
      {
        return this->base_plane_;
      }

      void set_base_plane(const Plane_3& plane)
      {
        this->base_plane_ = plane;
      }

      bool has_final_plane() const
      {
        return this->final_plane_.has_value();
      }

      const Plane_3& get_final_plane() const
      {
        return *(this->final_plane_);
      }

      void set_final_plane(const std::optional<Plane_3>& plane)
      {
        this->final_plane_ = plane;
      }

    protected:
      FacetWPtr facet_origin_;
      std::size_t input_face_id_; // descriptor at std::next(faces(input).begin(), input_face_id_)

      FT speed_;

      Plane_3 base_plane_;
      std::optional<Plane_3> final_plane_;
    };

private:
    // using VertexWPtr = std::weak_ptr<Vertex>;
    using VertexSPtr = std::shared_ptr<Vertex>;
    // using EdgeWPtr = std::weak_ptr<Edge>;
     using EdgeSPtr = std::shared_ptr<Edge>;
    // using FacetWPtr = std::weak_ptr<Facet>;
    using FacetSPtr = std::shared_ptr<Facet>;

    using FacetDataSPtr = std::shared_ptr<Facet_data>;
    using SkelFacetDataSPtr = std::shared_ptr<Skeleton_facet_data>;

  public:
    Facet()
      : id_(-1)
    { }

    virtual ~Facet()
    {
      edges_.clear();
      vertices_.clear();
    }

    static FacetSPtr create()
    {
      return std::make_shared<Facet>();
    }

    static FacetSPtr create(const std::vector<VertexSPtr>& vertices)
    {
      FacetSPtr result = std::make_shared<Facet>();

      const std::size_t nv = vertices.size();
      for (std::size_t i=0; i<nv; ++i) {
        result->add_vertex(vertices[i]);
      }
      for (std::size_t i=0; i<nv; ++i) {
        EdgeSPtr edge = vertices[i]->find_edge(vertices[(i+1)%nv]);
        if (!edge) {
          edge = Edge::create(vertices[i], vertices[(i+1)%nv]);
        }
        result->add_edge(edge);
      }
      return result;
    }

    /**
      * Each edge may have a facet to the left and a facet to the right.
      * Left and right is defined by a view outside the polyhedron.
      * In case an edge is used the first time, it is assumed that the left
      * side of the edge points to the interior of the facet.
      * In case the edge is already part of a facet,
      * it is assumed that the right side of the edge points to
      * the interior of the facet that is created.
      */
    static FacetSPtr create(const std::vector<EdgeSPtr>& edges)
    {
      FacetSPtr result = std::make_shared<Facet>();
      for (std::size_t i=0; i<edges.size(); ++i) {
        result->add_edge(edges[i]);
      }
      return result;
    }

    FacetSPtr clone() const
    {
      FacetSPtr result = Facet::create();
      std::map<VertexSPtr, VertexSPtr> old_to_new;
      typename std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
      while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        VertexSPtr vertex_c = vertex->clone();
        old_to_new[vertex] = vertex_c;
        result->add_vertex(vertex_c);
      }
      typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
      while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr src = old_to_new.at(edge->source());
        VertexSPtr tgt = old_to_new.at(edge->target());
        EdgeSPtr edge_c = Edge::create(src, tgt);
        CGAL_assertion(edge_c->source() == src && edge_c->target() == tgt);
        if (edge->get_facet_L() == this->shared_from_this()) {
          edge_c->set_facet_L(result);
        }
        if (edge->get_facet_R() == this->shared_from_this()) {
          edge_c->set_facet_R(result);
        }
        result->add_edge(edge_c);
      }
      result->set_id(id());
      if (has_data()) {
        result->set_data(get_data()->clone(std::const_pointer_cast<Facet>(this->shared_from_this())));
      }
      return result;
    }

    void add_vertex(const VertexSPtr& vertex)
    {
      vertices_.insert(vertices_.end(), vertex);
      vertex->add_facet(this->shared_from_this());
    }

    bool remove_vertex(const VertexSPtr& vertex)
    {
      bool result = false;
      typename std::list<VertexSPtr>::iterator it_v =
        std::find(vertices_.begin(), vertices_.end(), vertex);
      if (it_v != vertices_.end()) {
        vertices_.erase(it_v);
        vertex->remove_facet(this->shared_from_this());
        result = true;
      }
      return result;
    }

    bool has_vertex(const VertexSPtr& vertex)
    {
      typename std::list<VertexSPtr>::iterator it_v =
        std::find(vertices_.begin(), vertices_.end(), vertex);
      return (it_v != vertices_.end());
    }

    bool is_triangle() const
    {
      return (vertices_.size() == 3);
    }

    void add_edge(const EdgeSPtr& edge)
    {
      typename std::list<EdgeSPtr>::iterator it = edges_.insert(edges_.end(), edge);
      FacetSPtr facet_l = edge->get_facet_L();
      FacetSPtr facet_r = edge->get_facet_R();
      if (facet_l == this->shared_from_this() && facet_r == this->shared_from_this()) {
        typename std::list<EdgeSPtr>::iterator it_e =
          std::find(edges_.begin(), edges_.end(), edge);
        if (it_e == edge->getFacetLListIt()) {
          edge->setFacetRListIt(it);
        } else {
          edge->setFacetLListIt(it);
        }
      } else if (facet_l == this->shared_from_this()) {
        edge->setFacetLListIt(it);
      } else if (facet_r == this->shared_from_this()) {
        edge->setFacetRListIt(it);
      } else if (!facet_l) {
        edge->set_facet_L(this->shared_from_this());
        edge->setFacetLListIt(it);
      } else if (!facet_r) {
        edge->set_facet_R(this->shared_from_this());
        edge->setFacetRListIt(it);
      } else {
        CGAL_SS3_HDS_TRACE(edge->to_string());
        throw std::runtime_error("The given edge already has a left and a right facet.");
      }
      VertexSPtr vertex_src = edge->src(this->shared_from_this());
      if (!contains_vertex(vertex_src)) {
        add_vertex(vertex_src);
      }
      VertexSPtr vertex_tgt = edge->tgt(this->shared_from_this());
      if (!contains_vertex(vertex_tgt)) {
        add_vertex(vertex_tgt);
      }
    }

    bool remove_edge(const EdgeSPtr& edge)
    {
      bool result = false;
      if (edge->get_facet_L() == this->shared_from_this()) {
        edges_.erase(edge->getFacetLListIt());
        edge->set_facet_L(FacetSPtr());
        edge->setFacetLListIt(typename std::list<EdgeSPtr>::iterator());
        result = true;
      } else if (edge->get_facet_R() == this->shared_from_this()) {
        edges_.erase(edge->getFacetRListIt());
        edge->set_facet_R(FacetSPtr());
        edge->setFacetRListIt(typename std::list<EdgeSPtr>::iterator());
        result = true;
      }
      return result;
    }

    /**
    * Searches for an edge to the given facet.
    * Does not work when there is more than one edge adjacent to both facets.
    * @deprecated
    */
    EdgeSPtr find_edge(const FacetSPtr& facet) const
    {
      EdgeSPtr result = EdgeSPtr();
      typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
      while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->get_facet_L() == facet || edge->get_facet_R() == facet) {
          result = edge;
          break;
        }
      }
      return result;
    }

    std::list<EdgeSPtr> find_edges(const FacetSPtr& facet) const
    {
      std::list<EdgeSPtr> result;
      typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
      while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->get_facet_L() == facet || edge->get_facet_R() == facet) {
          result.push_back(edge);
        }
      }
      return result;
    }

    bool contains_vertex(const VertexSPtr& vertex) const
    {
      bool result = (vertices_.end() != std::find(vertices_.begin(), vertices_.end(), vertex));
      return result;
    }

    bool has_incident_edge(const EdgeSPtr& edge) const
    {
      bool result = (edges_.end() != std::find(edges_.begin(), edges_.end(), edge));
      return result;
    }

    void sort_vertices()
    {
      FacetSPtr self(this->shared_from_this());
      std::list<VertexSPtr> tmp;
      typename std::list<VertexSPtr>::iterator it_v_tmp;
      while (vertices_.size() > 0) {
        VertexSPtr first = VertexSPtr();
        VertexSPtr vertex = vertices_.front();
        EdgeSPtr edge = EdgeSPtr();
        typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
        while (it_e != edges_.end()) {
          EdgeSPtr my_edge = *it_e++;
          if (vertex == my_edge->src(self)) {
            edge = my_edge;
            break;
          }
        }
        while (vertex != first) {
          if (!first) {
            first = vertex;
            it_v_tmp = tmp.insert(tmp.end(), vertex);
          } else {
            tmp.insert(tmp.end(), vertex);
          }
          vertex = edge->tgt(self);
          edge = edge->next(self);
        }
        while (it_v_tmp != tmp.end()) {
          VertexSPtr vertex = *it_v_tmp++;
          typename std::list<VertexSPtr>::iterator it_v =
            std::find(vertices_.begin(), vertices_.end(), vertex);
          if (it_v != vertices_.end()) {
            vertices_.erase(it_v);
          }
        }
      }
      vertices_.clear();
      it_v_tmp = tmp.begin();
      while (it_v_tmp != tmp.end()) {
        VertexSPtr vertex = *it_v_tmp++;
        vertices_.insert(vertices_.end(), vertex);
      }
    }

    void sort_edges()
    {
      std::list<EdgeSPtr> tmp;
      typename std::list<EdgeSPtr>::iterator it_e_tmp;
      while (edges_.size() > 0) {
        EdgeSPtr first = EdgeSPtr();
        EdgeSPtr edge = edges_.front();
        while (edge != first) {
          if (!first) {
            first = edge;
            it_e_tmp = tmp.insert(tmp.end(), edge);
          } else {
            tmp.insert(tmp.end(), edge);
          }
          edge = edge->next(this->shared_from_this());
        }
        while (it_e_tmp != tmp.end()) {
          EdgeSPtr edge = *it_e_tmp++;
          remove_edge(edge);
        }
      }
      edges_.clear();
      it_e_tmp = tmp.begin();
      while (it_e_tmp != tmp.end()) {
        EdgeSPtr edge = *it_e_tmp++;
        add_edge(edge);
      }
    }

    PolyhedronSPtr get_polyhedron() const
    {
      return this->polyhedron_.lock();
    }

    void set_polyhedron(const PolyhedronSPtr& polyhedron)
    {
      this->polyhedron_ = polyhedron;
    }

    typename std::list<FacetSPtr>::iterator getPolyhedronListIt() const
    {
      return this->polyhedron_list_it_;
    }

    void setPolyhedronListIt(typename std::list<FacetSPtr>::iterator list_it)
    {
      this->polyhedron_list_it_ = list_it;
    }

    FacetDataSPtr get_data() const
    {
      CGAL_SS3_DEBUG_SPTR(this->data_);
      return this->data_;
    }

    void set_data(const FacetDataSPtr& data)
    {
      this->data_ = data;
    }

    bool has_data() const
    {
      bool result = false;
      if (data_) {
        result = true;
      }
      return result;
    }

    std::list<VertexSPtr>& vertices()
    {
      return this->vertices_;
    }

    std::list<EdgeSPtr>& edges()
    {
      return this->edges_;
    }

    FacetSPtr next(const VertexSPtr& vertex) const
    {
      FacetSPtr result = FacetSPtr();
      typename std::list<FacetWPtr>::const_iterator it_f = vertex->facets().begin();
      while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f;
        if (FacetSPtr facet = facet_wptr.lock()) {
          if (facet == this->shared_from_this()) {
            if (vertex->degree() == 1) {
              result = facet;
            }
            break;
          }
        }
        it_f++;
      }
      if (it_f != vertex->facets().end()) {
        typename std::list<FacetWPtr>::const_iterator it_f_begin = it_f++;
        if (it_f == vertex->facets().end()) {
          it_f = vertex->facets().begin();
        }
        while (it_f != it_f_begin) {
          FacetWPtr facet_wptr = *it_f++;
          if (it_f == vertex->facets().end()) {
            it_f = vertex->facets().begin();
          }
          if (FacetSPtr facet = facet_wptr.lock()) {
            typename std::list<EdgeWPtr>::const_iterator it_e = vertex->edges().begin();
            while (it_e != vertex->edges().end()) {
              EdgeWPtr edge_wptr = *it_e++;
              if (EdgeSPtr edge = edge_wptr.lock()) {
                FacetSPtr facet_l = edge->get_facet_L();
                FacetSPtr facet_r = edge->get_facet_R();
                if ((facet_l == this->shared_from_this() &&
                     facet_r == facet &&
                     edge->target() == vertex) ||
                    (facet_r == this->shared_from_this() &&
                     facet_l == facet &&
                     edge->source() == vertex)) {
                  result = facet;
                  break;
                }
              }
            }
          }
        }
      }
      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    FacetSPtr prev(const VertexSPtr& vertex) const
    {
      FacetSPtr result = FacetSPtr();
      typename std::list<FacetWPtr>::const_reverse_iterator it_f = vertex->facets().rbegin();
      while (it_f != vertex->facets().rend()) {
        FacetWPtr facet_wptr = *it_f;
        if (FacetSPtr facet = facet_wptr.lock()) {
          if (facet == this->shared_from_this()) {
            if (vertex->degree() == 1) {
              result = facet;
            }
            break;
          }
        }
        it_f++;
      }
      if (it_f != vertex->facets().rend()) {
        typename std::list<FacetWPtr>::const_reverse_iterator it_f_begin = it_f++;
        if (it_f == vertex->facets().rend()) {
          it_f = vertex->facets().rbegin();
        }
        while (it_f != it_f_begin) {
          FacetWPtr facet_wptr = *it_f++;
          if (it_f == vertex->facets().rend()) {
            it_f = vertex->facets().rbegin();
          }
          if (FacetSPtr facet = facet_wptr.lock()) {
            typename std::list<EdgeWPtr>::const_iterator it_e = vertex->edges().begin();
            while (it_e != vertex->edges().end()) {
              EdgeWPtr edge_wptr = *it_e++;
              if (EdgeSPtr edge = edge_wptr.lock()) {
                FacetSPtr facet_l = edge->get_facet_L();
                FacetSPtr facet_r = edge->get_facet_R();
                if ((facet_l == this->shared_from_this() &&
                     facet_r == facet &&
                     edge->source() == vertex) ||
                    (facet_r == this->shared_from_this() &&
                     facet_l == facet &&
                     edge->target() == vertex)) {
                  result = facet;
                  break;
                }
              }
            }
          }
        }
      }
      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    /**
    * merge 'facet' into this facet.
    */
    void merge(const FacetSPtr& facet)
    {
      typename std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
      while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        facet->remove_vertex(vertex);
        if (!contains_vertex(vertex)) {
          add_vertex(vertex);
        }
      }
      typename std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
      while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        FacetSPtr facet_l = edge->get_facet_L();
        FacetSPtr facet_r = edge->get_facet_R();
        if ((facet_l == facet && facet_r == this->shared_from_this()) ||
            (facet_r == facet && facet_l == this->shared_from_this())) {
          facet->remove_edge(edge);
          remove_edge(edge);
          edge->get_polyhedron()->remove_edge(edge);
        } else {
          if (facet_l == facet) {
            facet->remove_edge(edge);
            edge->set_facet_L(this->shared_from_this());
            add_edge(edge);
          }
          if (facet_r == facet) {
            facet->remove_edge(edge);
            edge->set_facet_R(this->shared_from_this());
            add_edge(edge);
          }
        }
      }

      CGAL_assertion(get_polyhedron()->is_consistent());
    }

    int id() const
    {
      return this->id_;
    }

    void set_id(const int id)
    {
      this->id_ = id;
    }

    /**
    * The direction of the normal points to the outside.
    */
    const Plane_3& get_plane() const
    {
      return this->plane_;
    }

    void set_plane(const Plane_3& plane)
    {
      this->plane_ = plane;
    }

    /**
    * First vertices have to form a triangle that is inside.
    */
    bool init_plane()
    {
      bool result = false;

      CGAL_SS3_HDS_TRACE_V(32, "initialize plane of F" << id());

      const Point_3* point_prev = nullptr;
      std::vector<const Point_3*> points;
      typename std::list<VertexSPtr>::iterator it_v = vertices_.begin();
      while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        const Point_3& point = vertex->point();
        if (!point_prev || point != *point_prev) {
          points.push_back(&point);
        }
        point_prev = &point;
      }

      CGAL_SS3_HDS_TRACE_V(64, "computing normals from " << points.size() << " points");
      CGAL_SS3_HDS_TRACE_CODE(for (std::size_t i=0; i<points.size(); ++i))
      CGAL_SS3_HDS_TRACE_V(64, "point " << i << ": " << *points[i]);

      if (points.size() >= 3) {
        const Point_3& p0 = *(points[0]);
        Vector_3 normal = CGAL::NULL_VECTOR;
        std::size_t last_i = points.size() - 1;
        for (std::size_t i=1; i<last_i; ++i) {
          const Point_3& p1 = *(points[i]);
          const Point_3& p2 = *(points[i+1]);
          normal += 0.5 * CGAL::cross_product(p2 - p1, p0 - p1);
        }

        if (normal != CGAL::NULL_VECTOR) {
          plane_ = Plane_3 { p0, normal };
          result = true;
        }
      }

      CGAL_postcondition(result);

      return result;
    }

    bool make_first_convex()
    {
      bool result = false;
      EdgeSPtr edge_begin;
      FacetSPtr self(this->shared_from_this());
      const Vector_3 normal = plane_.orthogonal_vector();
      EdgeSPtr edge = edges_.front();
      EdgeSPtr first = EdgeSPtr();
      while (edge != first) {
        if (!first) {
          first = edge;
        }
        EdgeSPtr edge_next = edge->next(self);
        const Point_3* points[3];
        points[0] = &(edge->src(self)->point());
        points[1] = &(edge->tgt(self)->point());
        points[2] = &(edge_next->tgt(self)->point());

        // @fixme is this correct? the CGAL_PI/4.0 below is confusing...
        // Was it supposed to be CGAL_PI/2.0?
        if (!CGAL::collinear(*(points[0]), *(points[1]), *(points[2]))) {
          if (CGAL::angle(*(points[0]), *(points[1]), *(points[2]), *normal) == CGAL::ACUTE) {
            edge_begin = edge;
            result = true;
            break;
          }
        }
        edge = edge_next;
      }

      if (edge_begin) {
        typename std::list<EdgeSPtr>::iterator it_e = edges_.insert(edges_.begin(), edge_begin);
        if (edge_begin->get_facet_L() == this->shared_from_this()) {
          edges_.erase(edge->getFacetLListIt());
          edge_begin->setFacetLListIt(it_e);
        } else if (edge_begin->get_facet_R() == this->shared_from_this()) {
          edges_.erase(edge->getFacetRListIt());
          edge_begin->setFacetRListIt(it_e);
        }
        sort_edges();
        VertexSPtr vertex_begin = edge_begin->src(this->shared_from_this());
        typename std::list<VertexSPtr>::iterator it_v =
          std::find(vertices_.begin(), vertices_.end(), vertex_begin);
        vertices_.erase(it_v);
        vertices_.insert(vertices_.begin(), vertex_begin);
        sort_vertices();
      }
      if (!result) {
        CGAL_SS3_HDS_TRACE("Warning: Unable to make first 3 vertices convex.");
        CGAL_SS3_HDS_TRACE(to_string());
      }

      return result;
    }

    /**
    * returns the number of vertices whose degree is higher than 3
    */
    int num_high_degree_vertices() const
    {
      int result = 0;
      // get min degree
      for (const VertexSPtr& v : vertices_) {
        if (v->degree() > 3) {
          ++result;
        }
      }
      return result;
    }

    std::string to_string() const
    {
      std::stringstream sstr;
      sstr << "Facet(";
      sstr << "id=" << IO::String_factory::fromInteger(id_) << ", ";
      // sstr << "addr=" << IO::String_factory::fromPointer(this) << ", ";
      sstr << "Plane: <" << IO::String_factory::fromDouble(CGAL::to_double(plane_.a())) << ", "
            << IO::String_factory::fromDouble(CGAL::to_double(plane_.b())) << ", "
            << IO::String_factory::fromDouble(CGAL::to_double(plane_.c())) << ", "
            << IO::String_factory::fromDouble(CGAL::to_double(plane_.d())) << ">, ";

      if (has_data()) {
        sstr << "Speed: " << std::dynamic_pointer_cast<Skeleton_facet_data>(get_data())->get_speed() << ", ";
      }

      sstr << "Vertices:" + IO::String_factory::fromInteger(vertices_.size()) + ", ";
      sstr << "Edges:" + IO::String_factory::fromInteger(edges_.size()) + ",";
      if (vertices_.size() > 0) {
        sstr << "\n";
        typename std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
        while (it_v != vertices_.end()) {
          VertexSPtr vertex = *it_v++;
          sstr << vertex->to_string() << "\n";
        }
      }
      if (edges_.size() > 0) {
        sstr << "\n";
        typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
          while (it_e != edges_.end()) {
          EdgeSPtr edge = *it_e++;
          sstr << edge->to_string() << "\n";
        }
      }
      sstr << ") END FACET to_string()\n";

      return sstr.str();
    }

  protected:
    Plane_3 plane_;
    int id_;
    FacetDataSPtr data_;

    std::list<VertexSPtr> vertices_;
    std::list<EdgeSPtr> edges_;
    PolyhedronWPtr polyhedron_;
    typename std::list<FacetSPtr>::iterator polyhedron_list_it_;
  };

public:
  using VertexWPtr = std::weak_ptr<Vertex>;
  using VertexSPtr = std::shared_ptr<Vertex>;
  using EdgeWPtr = std::weak_ptr<Edge>;
  using EdgeSPtr = std::shared_ptr<Edge>;
  using FacetWPtr = std::weak_ptr<Facet>;
  using FacetSPtr = std::shared_ptr<Facet>;

  using Vertex_data = typename Vertex::Vertex_data;
  using VertexDataSPtr = std::shared_ptr<Vertex_data>;
  using Skeleton_vertex_data = typename Vertex::Skeleton_vertex_data;
  using SkelVertexDataSPtr = std::shared_ptr<Skeleton_vertex_data>;
  using Edge_data = typename Edge::Edge_data;
  using EdgeDataSPtr = std::shared_ptr<Edge_data>;
  using Skeleton_edge_data = typename Edge::Skeleton_edge_data;
  using SkelEdgeDataSPtr = std::shared_ptr<Skeleton_edge_data>;
  using Facet_data = typename Facet::Facet_data;
  using FacetDataSPtr = std::shared_ptr<Facet_data>;
  using Skeleton_facet_data = typename Facet::Skeleton_facet_data;
  using SkelFacetDataSPtr = std::shared_ptr<Skeleton_facet_data>;

public:
  Polyhedron()
    : id_(-1)
  { }

  virtual ~Polyhedron()
  {
    facets_.clear();
    edges_.clear();
    vertices_.clear();
  }

  static PolyhedronSPtr create()
  {
    return std::make_shared<Polyhedron>();
  }

  static PolyhedronSPtr create(const std::vector<FacetSPtr>& facets)
  {
    PolyhedronSPtr result = std::make_shared<Polyhedron>();
    for (std::size_t i=0; i<facets.size(); ++i) {
      result->add_facet(facets[i]);
    }
    return result;
  }

  PolyhedronSPtr clone() const
  {
    std::map<VertexSPtr, VertexSPtr> vertices_c;
    std::map<EdgeSPtr, EdgeSPtr> edges_c;
    PolyhedronSPtr result = Polyhedron::create();
    typename std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
      VertexSPtr vertex = *it_v++;
      const Point_3& point = vertex->point();
      VertexSPtr vertex_c = Vertex::create(point);
      vertex_c->set_id(vertex->id());
      if (vertex->has_data()) {
        vertex_c->set_data(vertex->get_data()->clone());
      }
      result->add_vertex(vertex_c);
      vertices_c[vertex] = vertex_c;
    }
    typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
      EdgeSPtr edge = *it_e++;
      VertexSPtr src = vertices_c[edge->source()];
      VertexSPtr tgt = vertices_c[edge->target()];
      EdgeSPtr edge_c = Edge::create(src, tgt);
      edge_c->set_id(edge->id());
      if (edge->has_data()) {
        edge_c->set_data(edge->get_data()->clone());
      }
      result->add_edge(edge_c);
      edges_c[edge] = edge_c;
    }
    typename std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
      FacetSPtr facet = *it_f++;
      FacetSPtr facet_c = Facet::create();
      facet_c->set_plane(facet->get_plane());
      facet_c->set_id(facet->id());
      typename std::list<VertexSPtr>::const_iterator it_v = facet->vertices().begin();
      while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        facet_c->add_vertex(vertices_c[vertex]);
      }
      typename std::list<EdgeSPtr>::const_iterator it_e = facet->edges().begin();
      while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        EdgeSPtr edge_c = edges_c[edge];
        if (edge->get_facet_L() == facet) {
          edge_c->set_facet_L(facet_c);
        }
        if (edge->get_facet_R() == facet) {
          edge_c->set_facet_R(facet_c);
        }
        facet_c->add_edge(edge_c);
      }
      if (facet->has_data()) {
        facet_c->set_data(facet->get_data()->clone(facet));
      }
      result->add_facet(facet_c);
    }
    return result;
  }

  bool empty()
  {
    return vertices_.empty() && edges_.empty() && facets_.empty();
  }

  int id() const
  {
    return this->id_;
  }

  void set_id(const int id)
  {
    this->id_ = id;
  }

  void initialize_all_IDs()
  {
    reset_all_IDs();

    for (const VertexSPtr& vertex : vertices_) {
      vertex->set_id(next_vertex_id_++);
    }
    for (const EdgeSPtr& edge : edges_) {
      edge->set_id(next_edge_id_++);
    }
    for (const FacetSPtr& facet : facets_) {
      facet->set_id(next_facet_id_++);
    }
    set_id(-1);
  }

  void reset_all_IDs()
  {
    for (const VertexSPtr& vertex : vertices_) {
      vertex->set_id(-1);
    }
    for (const EdgeSPtr& edge : edges_) {
      edge->set_id(-1);
    }
    for (const FacetSPtr& facet : facets_) {
      facet->set_id(-1);
    }
    set_id(-1);

    next_vertex_id_ = 0;
    next_edge_id_ = 0;
    next_facet_id_ = 0;
  }

  std::list<VertexSPtr>& vertices()
  {
    return this->vertices_;
  }

  std::list<EdgeSPtr>& edges()
  {
    return this->edges_;
  }

  std::list<FacetSPtr>& facets()
  {
    return this->facets_;
  }

  bool has_vertex(const VertexSPtr& vertex)
  {
    typename std::list<VertexSPtr>::iterator it_v =
      std::find(vertices_.begin(), vertices_.end(), vertex);
    return (it_v != vertices_.end());
  }

  bool has_edge(const EdgeSPtr& edge)
  {
    typename std::list<EdgeSPtr>::iterator it_e =
      std::find(edges_.begin(), edges_.end(), edge);
    return (it_e != edges_.end());
  }

  bool has_facet(const FacetSPtr& facet)
  {
    typename std::list<FacetSPtr>::iterator it_f =
      std::find(facets_.begin(), facets_.end(), facet);
    return (it_f != facets_.end());
  }

  int maxVertexDegree() const
  {
    unsigned int max_degree = 0;
    typename std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
      VertexSPtr vertex = *it_v++;
      if (vertex->degree() > max_degree) {
        max_degree = vertex->degree();
      }
    }
    return max_degree;
  }

  void add_vertex(const VertexSPtr& vertex)
  {
    vertex->set_id(next_vertex_id_++);
    typename std::list<VertexSPtr>::iterator it = vertices_.insert(vertices_.end(), vertex);
    vertex->set_polyhedron(this->shared_from_this());
    vertex->setPolyhedronListIt(it);
  }

  bool remove_vertex(const VertexSPtr& vertex)
  {
    bool result = false;
    if (vertex->get_polyhedron() == this->shared_from_this()) {
      vertices_.erase(vertex->getPolyhedronListIt());
      vertex->setPolyhedronListIt(typename std::list<VertexSPtr>::iterator());
      vertex->set_polyhedron(PolyhedronSPtr());
      typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
      while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          this->remove_facet(facet);
        }
      }
      typename std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
      while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (EdgeSPtr edge = edge_wptr.lock()) {
          this->remove_edge(edge);
        }
      }
      result = true;
    }
    return result;
  }

  /**
   * Searches for a vertex with the same coordinates as the given vertex.
   */
  VertexSPtr find_vertex(const VertexSPtr& needle)
  {
    VertexSPtr result = VertexSPtr();
    typename std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
      VertexSPtr vertex = *it_v++;
      if (vertex->point() == needle->point()) {
        result = vertex;
        break;
      }
    }
    return result;
  }

  void add_edge(const EdgeSPtr& edge)
  {
    edge->set_id(next_edge_id_++);
    typename std::list<EdgeSPtr>::iterator it = edges_.insert(edges_.end(), edge);
    edge->set_polyhedron(this->shared_from_this());
    edge->setPolyhedronListIt(it);
    VertexSPtr vertex = edge->source();
    if (vertex->get_polyhedron() != this->shared_from_this()) {
      this->add_vertex(vertex);
    }
    vertex = edge->target();
    if (vertex->get_polyhedron() != this->shared_from_this()) {
      this->add_vertex(vertex);
    }
  }

  bool remove_edge(const EdgeSPtr& edge)
  {
    bool result = false;
    if (edge->get_polyhedron() == this->shared_from_this()) {
      edges_.erase(edge->getPolyhedronListIt());
      edge->setPolyhedronListIt(typename std::list<EdgeSPtr>::iterator());
      edge->set_polyhedron(PolyhedronSPtr());
      FacetSPtr facet = edge->get_facet_L();
      if (facet) {
        this->remove_facet(facet);
      }
      facet = edge->get_facet_R();
      if (facet) {
        this->remove_facet(facet);
      }
      edge->source()->remove_edge(edge);
      edge->target()->remove_edge(edge);
      result = true;
    }
    return result;
  }

  /**
   * Searches for an edge with the same coordinates as the given edge.
   * The orientation of the edge is ignored.
   */
  EdgeSPtr find_edge(const EdgeSPtr& needle)
  {
    EdgeSPtr result = EdgeSPtr();
    typename std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
      EdgeSPtr edge = *it_e++;
      if (edge->source()->point() == needle->source()->point() &&
          edge->target()->point() == needle->target()->point()) {
        result = edge;
        break;
      }
      if (edge->source()->point() == needle->target()->point() &&
          edge->target()->point() == needle->source()->point()) {
        result = edge;
        break;
      }
    }
    return result;
  }

  void add_facet(const FacetSPtr& facet)
  {
    facet->set_id(next_facet_id_++);
    typename std::list<FacetSPtr>::iterator it_f = facets_.insert(facets_.end(), facet);
    facet->setPolyhedronListIt(it_f);
    facet->set_polyhedron(this->shared_from_this());
    // add content of facet
    typename std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
    while (it_v != facet->vertices().end()) {
      VertexSPtr vertex = *it_v++;
      if (vertex->get_polyhedron() != this->shared_from_this()) {
        this->add_vertex(vertex);
      }
    }
    typename std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
      EdgeSPtr edge = *it_e++;
      if (edge->get_polyhedron() != this->shared_from_this()) {
        this->add_edge(edge);
      }
    }
  }

  bool remove_facet(const FacetSPtr& facet)
  {
    bool result = false;
    if (facet->get_polyhedron() == this->shared_from_this()) {
      facets_.erase(facet->getPolyhedronListIt());
      facet->set_polyhedron(PolyhedronSPtr());
      facet->setPolyhedronListIt(typename std::list<FacetSPtr>::iterator());
      typename std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
      while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        facet->remove_edge(edge);  // clears list iterators
      }
      typename std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
      while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        facet->remove_vertex(vertex);
      }
      result = true;
    }
    return result;
  }

  void init_planes()
  {
    for (const FacetSPtr& facet : facets_) {
      facet->init_plane();
    }
  }

  void clear_data()
  {
    for (const VertexSPtr& vertex : vertices_) {
      vertex->set_data(VertexDataSPtr());
    }
    for (const EdgeSPtr& edge : edges_) {
      edge->set_data(EdgeDataSPtr());
    }
    for (const FacetSPtr& facet : facets_) {
      facet->set_data(FacetDataSPtr());
    }
  }

  bool is_consistent() const
  {
    bool result = true;

    typename std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
      VertexSPtr vertex = *it_v++;
      if (vertex->get_polyhedron() != this->shared_from_this()) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->to_string());
        result = false;
        break;
      }
      typename std::list<EdgeWPtr>::const_iterator it_e = vertex->edges().begin();
      while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (vertex != edge->source() && vertex != edge->target()) {
            CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->to_string());
            CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
            result = false;
            break;
          }

          if (edge->source() == edge->target())
          {
            CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->source());
            result = false;
            break;
          }
        } else {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->to_string());
        }
      }

      typename std::list<FacetWPtr>::const_iterator it_f = vertex->facets().begin();
      while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          if (!facet->contains_vertex(vertex)) {
            CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->to_string());
            CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << facet->to_string());
            result = false;
            break;
          }
        } else {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->to_string());
        }
      }
    }

    typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
      EdgeSPtr edge = *it_e++;
      if (edge->get_polyhedron() != this->shared_from_this()) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
        result = false;
        break;
      }
      if (!edge->source()->has_incident_edge(edge)) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->source()->to_string());
        result = false;
        break;
      }
      if (!edge->target()->has_incident_edge(edge)) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->target()->to_string());
        result = false;
        break;
      }
      EdgeWPtr edge_wptr;
      edge_wptr = *(edge->getVertexSrcListIt());
      if (edge_wptr.lock() != edge) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
      }
      edge_wptr = *(edge->getVertexTgtListIt());
      if (edge_wptr.lock() != edge) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
      }
      if (edge->get_facet_L()) {
        if (!edge->get_facet_L()->has_incident_edge(edge)) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->get_facet_L()->to_string());
          result = false;
          break;
        }
        edge_wptr = *(edge->getFacetLListIt());
        if (edge_wptr.lock() != edge) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
        }
      } else {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
        result = false;
      }
      if (edge->get_facet_R()) {
        if (!edge->get_facet_R()->has_incident_edge(edge)) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->get_facet_R()->to_string());
          result = false;
          break;
        }
        edge_wptr = *(edge->getFacetRListIt());
        if (edge_wptr.lock() != edge) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
        }
      } else {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
        result = false;
      }
    }

    typename std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
      FacetSPtr facet = *it_f++;
      if (facet->get_polyhedron() != this->shared_from_this()) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << facet->to_string());
        result = false;
        break;
      }
      typename std::list<VertexSPtr>::const_iterator it_v = facet->vertices().begin();
      while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (!vertex->has_incident_facet(facet)) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << facet->to_string());
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->to_string());
          result = false;
          break;
        }
      }
      typename std::list<EdgeSPtr>::const_iterator it_e = facet->edges().begin();
      while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->get_facet_L() != facet && edge->get_facet_R() != facet) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << facet->to_string());
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
          result = false;
          break;
        }
        if (!facet->contains_vertex(edge->source())) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << facet->to_string());
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
          result = false;
          break;
        }
        if (!facet->contains_vertex(edge->target())) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << facet->to_string());
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->to_string());
          result = false;
          break;
        }
      }
    }

    return result;
  }

  void clear()
  {
    typename std::list<FacetSPtr>::iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
      FacetSPtr facet = *it_f++;
      this->remove_facet(facet);
    }
    typename std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
      EdgeSPtr edge = *it_e++;
      this->remove_edge(edge);
    }
    typename std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
      VertexSPtr vertex = *it_v++;
      this->remove_vertex(vertex);
    }
  }

  std::string to_string() const
  {
    std::stringstream sstr;
    sstr << "Polyhedron(";
    if (id_ != -1) {
      sstr << "id=" << IO::String_factory::fromInteger(id_) << ", ";
    } else {
      // sstr << IO::String_factory::fromPointer(this) << ", ";
    }
    sstr << "Vertices:" + IO::String_factory::fromInteger(vertices_.size()) + ", ";
    sstr << "Edges:" + IO::String_factory::fromInteger(edges_.size()) + ", ";
    sstr << "Facets:" + IO::String_factory::fromInteger(facets_.size()) + ",\n";

    typename std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
      FacetSPtr facet = *it_f++;
      sstr << facet->to_string() << "\n";
    }
    sstr << ")";

  /*typename std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
      EdgeSPtr edge = *it_e++;
      sstr << edge->to_string() << "\n";
    }
    typename std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
      VertexSPtr vertex = *it_v++;
      sstr << vertex->to_string() << "\n";
    }*/
    return sstr.str();
  }

  void dump_edges(const std::string& base_filename) const
  {
    unsigned int fi = 0;
    typename std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
      FacetSPtr facet = *it_f++;

      std::ofstream edge_out(base_filename + "_" + std::to_string(fi++) + ".polylines.txt");
      edge_out.precision(17);

      typename std::list<EdgeSPtr>::const_iterator it_e = facet->edges().begin();
      while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        edge_out << "2 " << edge->source()->point() << " "
                         << edge->target()->point() << std::endl;
      }
    }
  }

private:
  std::list<VertexSPtr> vertices_;
  std::list<EdgeSPtr> edges_;
  std::list<FacetSPtr> facets_;

  int id_;

  int next_vertex_id_ = 0;
  int next_edge_id_ = 0;
  int next_facet_id_ = 0;
};

} // namespace HDS
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_INTERNAL_HDS_POLYHEDRON_H */
