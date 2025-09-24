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

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/internal/weak_find.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>
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
namespace Straight_skeletons_3 {
namespace internal {
namespace SDS {

template <typename Traits>
class Node;

template <typename Traits>
class Arc;

template <typename Traits>
class Sheet;

} // namespace SDS

namespace HDS {

template <typename Traits_>
class Polyhedron
  : public std::enable_shared_from_this<Polyhedron<Traits_> >
{
public:
  using Traits = Traits_;

  using PolyhedronSPtr = std::shared_ptr<Polyhedron>;
  using PolyhedronWPtr = std::weak_ptr<Polyhedron>;

private:
  using KernelFactory = kernel::KernelFactory<Traits>;

  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;
  using Segment_3 = typename Traits::Segment_3;
  using Vector_3 = typename Traits::Vector_3;
  using Line_3 = typename Traits::Line_3;
  using Plane_3 = typename Traits::Plane_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Segment3SPtr = std::shared_ptr<Segment_3>;
  using Vector3SPtr = std::shared_ptr<Vector_3>;
  using Line3SPtr = std::shared_ptr<Line_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;

public:
  template <typename GT>
  class Vertex;

  template <typename GT>
  class Edge;

  template <typename GT>
  class Facet;

public:
  template <typename GT>
  class Vertex
    : public std::enable_shared_from_this<Vertex<GT> >
  {
  public:
    class VertexData
    {
      using VertexWPtr = std::weak_ptr<Vertex<GT> >;
      using VertexSPtr = std::shared_ptr<Vertex<GT> >;

      using VertexDataSPtr = std::shared_ptr<VertexData>;

    public:
      VertexData() { /*intentionally does nothing*/ }
      virtual ~VertexData() { /*intentionally does nothing*/ }

      static VertexDataSPtr create(const VertexSPtr& vertex)
      {
        CGAL_SS3_DEBUG_SPTR(vertex);
        VertexDataSPtr result = std::make_shared<VertexData>();
        vertex->setData(result);
        return result;
      }

      virtual VertexDataSPtr clone() const {
        return std::make_shared<VertexData>(*this);
      }
    };

    class SkelVertexData
      : public VertexData
    {
      using VertexWPtr = std::weak_ptr<Vertex<GT> >;
      using VertexSPtr = std::shared_ptr<Vertex<GT> >;

      using NodeWPtr = std::weak_ptr<SDS::Node<GT> >;
      using NodeSPtr = std::shared_ptr<SDS::Node<GT> >;
      using ArcWPtr = std::weak_ptr<SDS::Arc<GT> >;
      using ArcSPtr = std::shared_ptr<SDS::Arc<GT> >;

      using VertexDataSPtr = std::shared_ptr<VertexData>;
      using SkelVertexDataSPtr = std::shared_ptr<SkelVertexData>;

    public:
      SkelVertexData() { /*intentionally does nothing*/ }
      virtual ~SkelVertexData() { /*intentionally does nothing*/ }

      static SkelVertexDataSPtr create(const VertexSPtr& vertex)
      {
        CGAL_SS3_DEBUG_SPTR(vertex);
        SkelVertexDataSPtr result = std::make_shared<SkelVertexData>();
        vertex->setData(result);
        return result;
      }

      virtual VertexDataSPtr clone() const override {
        return std::make_shared<SkelVertexData>(*this);
      }

      ArcSPtr getArc() const
      {
        CGAL_SS3_DEBUG_WPTR(arc_);
        return this->arc_.lock();
      }

      void setArc(ArcSPtr arc)
      {
        CGAL_SS3_DEBUG_SPTR(arc);
        this->arc_ = arc;
      }

      NodeSPtr getNode() const
      {
        CGAL_SS3_DEBUG_WPTR(node_);
        return this->node_.lock();
      }

      void setNode(NodeSPtr node)
      {
        CGAL_SS3_DEBUG_SPTR(node);
        this->node_ = node;
      }

      bool hasFinalPoint() const
      {
        return bool(final_point_);
      }

      Point3SPtr getFinalPoint() const
      {
        CGAL_SS3_DEBUG_SPTR(final_point_);
        return final_point_;
      }

      void setFinalPoint(Point3SPtr point)
      {
        final_point_ = point;
      }

    protected:
      ArcWPtr arc_;
      NodeWPtr node_;
      Point3SPtr final_point_;
    };

  private:
    using FT = typename GT::FT;
    using Point_3 = typename GT::Point_3;

    using VertexSPtr = std::shared_ptr<Vertex<GT> >;
    using EdgeWPtr = std::weak_ptr<Edge<GT> >;
    using EdgeSPtr = std::shared_ptr<Edge<GT> >;
    using FacetWPtr = std::weak_ptr<Facet<GT> >;
    using FacetSPtr = std::shared_ptr<Facet<GT> >;

    using VertexDataSPtr = std::shared_ptr<VertexData>;
    using SkelVertexDataSPtr = std::shared_ptr<SkelVertexData>;

  public:
    Vertex(Point3SPtr point)
      : point_(point), id_(-1)
    { }

    virtual ~Vertex()
    {
      facets_.clear();
      edges_.clear();
    }

    static VertexSPtr create(Point3SPtr point)
    {
      CGAL_SS3_DEBUG_SPTR(point);
      return std::make_shared<Vertex>(point);
    }

    VertexSPtr clone() const
    {
      VertexSPtr result = std::make_shared<Vertex>(getPoint());
      CGAL_SS3_DEBUG_SPTR(result);
      result->setID(getID());
      if (hasData()) {
        result->setData(getData()->clone());
      }
      return result;
    }

    Point3SPtr getPoint() const
    {
      CGAL_SS3_DEBUG_SPTR(point_);
      return point_;
    }

    void setPoint(Point3SPtr point)
    {
      CGAL_SS3_DEBUG_SPTR(point);
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

    PolyhedronSPtr getPolyhedron() const
    {
      return polyhedron_.lock();
    }

    void setPolyhedron(PolyhedronSPtr polyhedron)
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

    int getID() const
    {
      return id_;
    }

    void setID(int id)
    {
      id_ = id;
    }

    VertexDataSPtr getData() const
    {
      CGAL_SS3_DEBUG_SPTR(this->data_);
      return this->data_;
    }

    void setData(VertexDataSPtr data)
    {
      this->data_ = data;
    }

    bool hasData() const
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

    void addEdge(const EdgeSPtr& edge)
    {
      EdgeWPtr edge_wptr(edge);
      typename std::list<EdgeWPtr>::iterator it = edges_.insert(edges_.end(), edge_wptr);
      VertexSPtr vertex_src = edge->getVertexSrc();
      VertexSPtr vertex_dst = edge->getVertexDst();
      if (vertex_src == this->shared_from_this() && vertex_dst == this->shared_from_this()) {
        typename std::list<EdgeWPtr>::iterator it_e =
          STL_Extension::internal::weak_find(edges_.begin(), edges_.end(), edge_wptr);
        if (it_e == edge->getVertexSrcListIt()) {
          edge->setVertexDstListIt(it);
        } else {
          edge->setVertexSrcListIt(it);
        }
      } else if (vertex_src == this->shared_from_this()) {
        edge->setVertexSrcListIt(it);
      } else if (vertex_dst == this->shared_from_this()) {
        edge->setVertexDstListIt(it);
      }
    }

    bool removeEdge(const EdgeSPtr& edge)
    {
      bool result = false;
      if (edge->getVertexSrc() == this->shared_from_this()) {
        edges_.erase(edge->getVertexSrcListIt());
        edge->setVertexSrcListIt(typename std::list<EdgeWPtr>::iterator());
        result = true;
      } else if (edge->getVertexDst() == this->shared_from_this()) {
        edges_.erase(edge->getVertexDstListIt());
        edge->setVertexDstListIt(typename std::list<EdgeWPtr>::iterator());
        result = true;
      }
      return result;
    }

    EdgeSPtr firstEdge() const
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

    EdgeSPtr getEdge(unsigned int i)
    {
      EdgeSPtr edge = firstEdge();
      for (unsigned int j = 0; j < i; ++j) {
        edge = edge->next(this->shared_from_this());
      }
      return edge;
    }

    /**
      * Searches for an edge to the given destination.
      * The orientation of the edge is ignored.
      */
    EdgeSPtr findEdge(const VertexSPtr& dst) const
    {
      EdgeSPtr result = EdgeSPtr();
      for (EdgeWPtr edge_wptr : edges_) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge->getVertexSrc().get() == this && edge->getVertexDst() == dst) {
            result = edge;
            break;
          }
          if (edge->getVertexDst().get() == this && edge->getVertexSrc() == dst) {
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
    EdgeSPtr findEdge(const FacetSPtr& facet) const
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

    void addFacet(const FacetSPtr& facet)
    {
      facets_.insert(facets_.end(), FacetWPtr(facet));
    }

    bool removeFacet(const FacetSPtr& facet)
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

    FacetSPtr firstFacet() const
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

    FacetSPtr getFacet(unsigned int i)
    {
      FacetSPtr facet = firstFacet();
      for (unsigned int j = 0; j < i; ++j) {
        facet = facet->next(this->shared_from_this());
      }
      return facet;
    }

    bool containsEdge(const EdgeSPtr& edge) const
    {
      EdgeWPtr edge_wptr = EdgeWPtr(edge);
      bool result = (edges_.end() !=
          STL_Extension::internal::weak_find(edges_.begin(), edges_.end(), edge_wptr));
      return result;
    }

    bool containsFacet(const FacetSPtr& facet) const
    {
      FacetWPtr facet_wptr = FacetWPtr(facet);
      bool result = (facets_.end() !=
          STL_Extension::internal::weak_find(facets_.begin(), facets_.end(), facet_wptr));
      return result;
    }

    void sortEdges()
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
        CGAL_assertion(bool(edge));
        while (edge != first) {
          if (!first) {
            first = edge;
          }
          tmp.insert(tmp.end(), edge);
          edge = edge->next(this->shared_from_this());
        }
        while (it_e_tmp != tmp.end()) {
          EdgeSPtr edge = *it_e_tmp++;
          removeEdge(edge);
        }
      }
      edges_.clear();
      it_e_tmp = tmp.begin();
      while (it_e_tmp != tmp.end()) {
        EdgeSPtr edge = *it_e_tmp++;
        addEdge(edge);
      }
    }

    void sortFacets()
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
        EdgeSPtr edge_first = findEdge(facet);
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
          removeFacet(facet);
        }
      }
      facets_.clear();
      it_f_tmp = tmp.begin();
      while (it_f_tmp != tmp.end()) {
        FacetSPtr facet = *it_f_tmp++;
        addFacet(facet);
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
        CGAL_assertion(bool(edge_first));
        EdgeSPtr edge;
        FacetSPtr facet = edge_first->getFacetL();
        if (edge_first->getVertexDst() == this->shared_from_this()) {
          facet = edge_first->getFacetR();
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
          removeEdge(edge);
        }
        while (it_f_tmp != facets_tmp.end()) {
          FacetSPtr facet = *it_f_tmp++;
          removeFacet(facet);
        }
      }
      edges_.clear();
      it_e_tmp = edges_tmp.begin();
      while (it_e_tmp != edges_tmp.end()) {
        EdgeSPtr edge = *it_e_tmp++;
        addEdge(edge);
      }
      facets_.clear();
      it_f_tmp = facets_tmp.begin();
      while (it_f_tmp != facets_tmp.end()) {
        FacetSPtr facet = *it_f_tmp++;
        addFacet(facet);
      }
    }

    VertexSPtr next(const FacetSPtr& facet) const
    {
      VertexSPtr result = VertexSPtr();
      CGAL_assertion(facet->vertices().size() > 1);

      for (EdgeWPtr edge_wptr : edges_) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge->src(facet) == this->shared_from_this()) {
            result = edge->dst(facet);
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
          if (edge->dst(facet) == this->shared_from_this()) {
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
    * The destination vertex will be returned.
    */
    VertexSPtr split(const FacetSPtr& facet_left,
                     FacetSPtr facet_right)
    {
      VertexSPtr result = VertexSPtr();
      if (!containsFacet(facet_left) || !containsFacet(facet_right)) {
        return result;
      }
    //  if (facet_left->findEdge(facet_right)) {
    //      return result;
    //  }
      result = Vertex::create(getPoint());

      // 1. select edges and facets for result
      std::list<FacetSPtr> facets;
      std::list<EdgeSPtr> edges;
      FacetSPtr poly_curr = facet_right;
      EdgeSPtr edge_curr = findEdge(poly_curr);
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
        if (edge->getVertexSrc() == this->shared_from_this()) {
          edge->replaceVertexSrc(result);
        } else if (edge->getVertexDst() == this->shared_from_this()) {
          edge->replaceVertexDst(result);
        }
      }
      typename std::list<FacetSPtr>::iterator it_f = facets.begin();
      while (it_f != facets.end()) {
        FacetSPtr facet = *it_f++;
        if (facet != facet_left && facet != facet_right) {
          facet->removeVertex(this->shared_from_this());
        }
        facet->addVertex(result);
      }

      // 3. insert connecting edge
      EdgeSPtr edge = Edge<GT>::create(this->shared_from_this(), result);
      edge->setFacetL(facet_left);
      edge->setFacetR(facet_right);
      facet_left->addEdge(edge);
      facet_right->addEdge(edge);

      if (PolyhedronSPtr polyhedron = polyhedron_.lock()) {
        polyhedron->addVertex(result);
        polyhedron->addEdge(edge);
      }
      return result;
    }

    std::string toString() const
    {
      std::string result("Vertex(");
      result += "id=" + IO::StringFactory::fromInteger(id_) + ", ";
      // result += "addr=" + IO::StringFactory::fromPointer(this) +", ";
      result += "<" + IO::StringFactory::fromDouble(CGAL::to_double(getPoint()->x())) + " ";
      result += IO::StringFactory::fromDouble(CGAL::to_double(getPoint()->y())) + " ";
      result += IO::StringFactory::fromDouble(CGAL::to_double(getPoint()->z())) + ">";
      result += ", Edges:" + IO::StringFactory::fromInteger(edges_.size());
      result += ", Facets:" + IO::StringFactory::fromInteger(facets_.size()) + " {";
      for (FacetWPtr facet_wptr : facets_) {
        if (FacetSPtr facet = facet_wptr.lock()) {
          result += " " + std::to_string(facet->getID());
        }
      }
      result += " }";
      result += ")";
      return result;
    }

  protected:
    Point3SPtr point_;

    std::list<EdgeWPtr> edges_;
    std::list<FacetWPtr> facets_;
    PolyhedronWPtr polyhedron_;
    typename std::list<VertexSPtr>::iterator polyhedron_list_it_;

    int id_;
    VertexDataSPtr data_;
  };

public:
  template <typename GT>
  class Edge
    : public std::enable_shared_from_this<Edge<GT> >
  {
  public:
    class EdgeData
    {
      using EdgeWPtr = std::weak_ptr<Edge<GT> >;
      using EdgeSPtr = std::shared_ptr<Edge<GT> >;

      using EdgeDataSPtr = std::shared_ptr<EdgeData>;

    public:
      EdgeData() { /*intentionally does nothing*/ }
      virtual ~EdgeData() { /*intentionally does nothing*/ }

      static EdgeDataSPtr create(const EdgeSPtr& edge)
      {
        EdgeDataSPtr result = std::make_shared<EdgeData>();
        edge->setData(result);
        return result;
      }

      virtual EdgeDataSPtr clone() const {
        return std::make_shared<EdgeData>(*this);
      }
    };

    class SkelEdgeData
      : public EdgeData
    {
      using EdgeWPtr = std::weak_ptr<Edge<GT> >;
      using EdgeSPtr = std::shared_ptr<Edge<GT> >;
      using FacetWPtr = std::weak_ptr<Facet<GT> >;
      using FacetSPtr = std::shared_ptr<Facet<GT> >;

      using SheetWPtr = std::weak_ptr<SDS::Sheet<GT> >;
      using SheetSPtr = std::shared_ptr<SDS::Sheet<GT> >;

      using EdgeDataSPtr = std::shared_ptr<EdgeData>;
      using SkelEdgeDataSPtr = std::shared_ptr<SkelEdgeData>;

    public:
      SkelEdgeData() { /*intentionally does nothing*/ }
      virtual ~SkelEdgeData() { /*intentionally does nothing*/ }

      static SkelEdgeDataSPtr create(const EdgeSPtr& edge)
      {
        SkelEdgeDataSPtr result = std::make_shared<SkelEdgeData>();
        edge->setData(result);
        return result;
      }

      virtual EdgeDataSPtr clone() const override {
        return std::make_shared<SkelEdgeData>(*this);
      }

      SheetSPtr getSheet() const
      {
        CGAL_SS3_DEBUG_WPTR(sheet_);
        return this->sheet_.lock();
      }

      void setSheet(SheetSPtr sheet)
      {
        this->sheet_ = sheet;
      }

    protected:
      SheetWPtr sheet_;
    };

  private:
    using FT = typename GT::FT;
    using Segment_3 = typename GT::Segment_3;
    using Vector_3 = typename GT::Vector_3;
    using Line_3 = typename GT::Line_3;

    using Segment3SPtr = std::shared_ptr<Segment_3>;
    using Vector3SPtr = std::shared_ptr<Vector_3>;
    using Line3SPtr = std::shared_ptr<Line_3>;

    // using VertexWPtr = std::weak_ptr<Vertex<GT> >;
    using VertexSPtr = std::shared_ptr<Vertex<GT> >;
    using EdgeWPtr = std::weak_ptr<Edge<GT> >;
    using EdgeSPtr = std::shared_ptr<Edge<GT> >;
    using FacetWPtr = std::weak_ptr<Facet<GT> >;
    using FacetSPtr = std::shared_ptr<Facet<GT> >;

    using EdgeDataSPtr = std::shared_ptr<EdgeData>;
    using SkelEdgeDataSPtr = std::shared_ptr<SkelEdgeData>;

  public:
    Edge(const VertexSPtr& src, const VertexSPtr& dst)
      : vertex_src_(src), vertex_dst_(dst), id_(-1)
    { }

    Edge(const Edge& edge)
      : vertex_src_(edge.vertex_src_->clone()),
        vertex_dst_(edge.vertex_dst_->clone()),
        id_(-1)
    { }

    virtual ~Edge()
    {
      vertex_src_.reset();
      vertex_dst_.reset();
    }

    static EdgeSPtr create(const VertexSPtr& src, const VertexSPtr& dst)
    {
      EdgeSPtr result = std::make_shared<Edge>(src, dst);
      src->addEdge(result);
      dst->addEdge(result);
      return result;
    }

    EdgeSPtr clone() const
    {
      EdgeSPtr result = std::make_shared<Edge>();
      CGAL_SS3_DEBUG_SPTR(result);
      result->vertex_src_->addEdge(result);
      result->vertex_dst_->addEdge(result);
      result->setID(getID());
      if (hasData()) {
        result->setData(getData()->clone());
      }
      return result;
    }

    VertexSPtr getVertexSrc() const
    {
      CGAL_SS3_DEBUG_SPTR(this->vertex_src_);
      return this->vertex_src_;
    }

    void setVertexSrc(const VertexSPtr& src)
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

    VertexSPtr getVertexDst() const
    {
      CGAL_SS3_DEBUG_SPTR(this->vertex_dst_);
      return this->vertex_dst_;
    }

    void setVertexDst(const VertexSPtr& dst)
    {
      this->vertex_dst_ = dst;
    }

    typename std::list<EdgeWPtr>::iterator getVertexDstListIt() const
    {
      return this->vertex_dst_list_it_;
    }

    void setVertexDstListIt(typename std::list<EdgeWPtr>::iterator list_it)
    {
      this->vertex_dst_list_it_ = list_it;
    }

    FacetSPtr getFacetL() const
    {
      return this->facet_l_.lock();
    }

    void setFacetL(const FacetSPtr& facet)
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

    FacetSPtr getFacetR() const
    {
      return this->facet_r_.lock();
    }

    void setFacetR(const FacetSPtr& facet)
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

    FacetSPtr getFacetSrc() const
    {
      FacetSPtr result = FacetSPtr();
      VertexSPtr vertex_src = getVertexSrc();
      if (vertex_src->degree() == 3) {
        if (getFacetL()) {
          result = getFacetL()->next(vertex_src);
        }
      }
      return result;
    }

    FacetSPtr getFacetDst() const
    {
      FacetSPtr result = FacetSPtr();
      VertexSPtr vertex_dst = getVertexDst();
      if (vertex_dst->degree() == 3) {
        if (getFacetR()) {
          result = getFacetR()->next(vertex_dst);
        }
      }
      return result;
    }

    PolyhedronSPtr getPolyhedron() const
    {
      return this->polyhedron_.lock();
    }

    void setPolyhedron(PolyhedronSPtr polyhedron)
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

    int getID() const
    {
      return this->id_;
    }

    void setID(int id)
    {
      this->id_ = id;
    }

    EdgeDataSPtr getData() const
    {
      CGAL_SS3_DEBUG_SPTR(this->data_);
      return this->data_;
    }

    void setData(EdgeDataSPtr data)
    {
      this->data_ = data;
    }

    bool hasData() const
    {
      bool result = false;
      if (data_) {
        result = true;
      }
      return result;
    }

    std::optional<bool> getReflexStatus() const
    {
      return cachedReflexStatus_;
    }

    Segment3SPtr segment() const
    {
      return KernelFactory::createSegment3(vertex_src_->getPoint(), vertex_dst_->getPoint());
    }

    Line3SPtr line() const
    {
      return KernelFactory::createLine3(vertex_src_->getPoint(), vertex_dst_->getPoint());
    }

    VertexSPtr other(const VertexSPtr& vertex) const
    {
      VertexSPtr result = VertexSPtr();
      if (vertex == vertex_src_) {
        result = vertex_dst_;
      } else if (vertex == vertex_dst_) {
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
        result = vertex_dst_;
      }
      return result;
    }

    VertexSPtr dst(const FacetSPtr& facet_l) const
    {
      VertexSPtr result = VertexSPtr();
      if (facet_l == facet_l_.lock()) {
        result = vertex_dst_;
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
      } else if (vertex_src == vertex_dst_) {
        result = facet_r_.lock();
      }
      return result;
    }

    FacetSPtr right(const VertexSPtr& vertex_src) const
    {
      FacetSPtr result = FacetSPtr();
      if (vertex_src == vertex_src_) {
        result = facet_r_.lock();
      } else if (vertex_src == vertex_dst_) {
        result = facet_l_.lock();
      }
      return result;
    }

    EdgeSPtr next(const FacetSPtr& facet) const
    {
      EdgeSPtr result = EdgeSPtr();
      VertexSPtr vertex_dst = dst(facet);
      for (EdgeWPtr edge_wptr : vertex_dst->edges()) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge->src(facet) == vertex_dst) {
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
          if (edge->dst(facet) == vertex_src) {
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
      } else if (vertex == vertex_dst_) {
        facet = FacetSPtr(facet_r_);
      }

      CGAL_SS3_DEBUG_SPTR(facet);

      for (EdgeWPtr edge_wptr : vertex->edges()) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge.get() == this) {
            continue;
          }
          if (edge->dst(facet) == vertex) {
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
      } else if (vertex == vertex_dst_) {
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
    * Swaps src and dst vertex and left and right facet
    */
    void invert()
    {
      VertexSPtr vertex_tmp_ = vertex_src_;
      typename std::list<EdgeWPtr>::iterator vertex_tmp_list_it_ = vertex_src_list_it_;
      FacetWPtr facet_tmp_ = facet_l_;
      typename std::list<EdgeSPtr>::iterator facet_tmp_list_it_ = facet_l_list_it_;
      vertex_src_ = vertex_dst_;
      vertex_src_list_it_ = vertex_dst_list_it_;
      facet_l_ = facet_r_;
      facet_l_list_it_ = facet_r_list_it_;
      vertex_dst_ = vertex_tmp_;
      vertex_dst_list_it_ = vertex_tmp_list_it_;
      facet_r_ = facet_tmp_;
      facet_r_list_it_ = facet_tmp_list_it_;
    }

    /**
      * Splits this edge into 2 edges.
      * The destination vertex of this edge is set to the given middle vertex.
      * It returns newly created edge.
      */
    EdgeSPtr split(const VertexSPtr& middle)
    {
      EdgeSPtr result = Edge::create(middle, vertex_dst_);
      if (FacetSPtr facet_l = facet_l_.lock()) {
        result->setFacetL(facet_l);
        facet_l->addEdge(result);
      }
      if (FacetSPtr facet_r = facet_r_.lock()) {
        result->setFacetR(facet_r);
        facet_r->addEdge(result);
      }
      vertex_dst_->removeEdge(this->shared_from_this());
      vertex_dst_ = middle;
      middle->addEdge(this->shared_from_this());
      if (PolyhedronSPtr polyhedron = polyhedron_.lock()) {
        polyhedron->addEdge(result);
      }
      return result;
    }

    /**
      * More than just a simple set method.
      */
    void replaceVertexSrc(const VertexSPtr& vertex_src)
    {
      vertex_src_->removeEdge(this->shared_from_this());
      vertex_src_ = vertex_src;
      vertex_src->addEdge(this->shared_from_this());
    }

    void replaceVertexDst(const VertexSPtr& vertex_dst)
    {
      vertex_dst_->removeEdge(this->shared_from_this());
      vertex_dst_ = vertex_dst;
      vertex_dst->addEdge(this->shared_from_this());
    }

    void replaceFacetL(const FacetSPtr& facet_l)
    {
      if (FacetSPtr facet = facet_l_.lock()) {
        facet->removeEdge(this->shared_from_this());
      }
      setFacetL(facet_l);
      facet_l->addEdge(this->shared_from_this());
    }

    void replaceFacetR(const FacetSPtr& facet_r)
    {
      if (FacetSPtr facet = facet_r_.lock()) {
        facet->removeEdge(this->shared_from_this());
      }
      setFacetR(facet_r);
      facet_r->addEdge(this->shared_from_this());
    }

    bool hasSameFacets(const EdgeSPtr& edge) const
    {
      bool result = (facet_l_.lock() == edge->getFacetL() &&
                     facet_r_.lock() == edge->getFacetR()) ||
                    (facet_r_.lock() == edge->getFacetL() &&
                     facet_l_.lock() == edge->getFacetR());
      return result;
    }

    bool isReflex() const
    {
      CGAL_precondition(*(vertex_src_->getPoint()) != *(vertex_dst_->getPoint()));
      if (cachedReflexStatus_) {
        return *cachedReflexStatus_;
      }
      bool result = false;
      FacetSPtr facet_l = this->getFacetL();
      FacetSPtr facet_r = this->getFacetR();
      CGAL_precondition(bool(facet_l));
      CGAL_precondition(bool(facet_r));
      Plane3SPtr plane_l = facet_l->plane();
      Vector3SPtr normal_l = KernelFactory::createVector3(plane_l);
      Vector3SPtr dir = KernelFactory::createVector3(line());
      Point3SPtr p_src = vertex_src_->getPoint();
      CGAL_assertion(*normal_l != CGAL::NULL_VECTOR);
      CGAL_assertion(*dir != CGAL::NULL_VECTOR);
      Point_3 p = (*p_src) + CGAL::cross_product(*normal_l, *dir);
      Plane3SPtr plane_r = facet_r->plane();
      if (plane_r->oriented_side(p) == CGAL::ON_POSITIVE_SIDE) {
        result = true;
      }
      cachedReflexStatus_ = result;
      return result;
    }

    std::string toString() const
    {
      std::string result("Edge(");
      result += "id=" + IO::StringFactory::fromInteger(id_);
      // result += ", addr=" + IO::StringFactory::fromPointer(this);
      result += ", l=";
      if (FacetSPtr facet_l = getFacetL()) {
          if (facet_l->getID() != -1) {
            result += IO::StringFactory::fromInteger(getFacetL()->getID());
          } else {
            // result += IO::StringFactory::fromPointer(getFacetL().get());
          }
      } else {
        result += "expired";
      }
      result += ", r=";
      if (FacetSPtr facet_r = getFacetR()) {
          if (facet_r->getID() != -1) {
            result += IO::StringFactory::fromInteger(getFacetR()->getID());
          } else {
            // result += IO::StringFactory::fromPointer(getFacetR().get());
          }
      } else {
        result += "expired";
      }
      // src and dst faces
      result += ", s=";
      if (getFacetSrc()) {
          if (getFacetSrc()->getID() != -1) {
            result += IO::StringFactory::fromInteger(getFacetSrc()->getID());
          } else {
            // result += IO::StringFactory::fromPointer(getFacetSrc().get());
          }
      } else {
        result += "nullptr";
      }
      result += ", d=";
      if (getFacetDst()) {
          if (getFacetDst()->getID() != -1) {
            result += IO::StringFactory::fromInteger(getFacetDst()->getID());
          } else {
            // result += IO::StringFactory::fromPointer(getFacetDst().get());
          }
      } else {
        result += "nullptr";
      }
      if (vertex_src_) {
        result += ",\n     src=" + vertex_src_->toString();
      }
      if (vertex_dst_) {
        result += ",\n     dst=" + vertex_dst_->toString();
      }
      result += ")";
      return result;
    }

  protected:
    VertexSPtr vertex_src_;
    typename std::list<EdgeWPtr>::iterator vertex_src_list_it_;
    VertexSPtr vertex_dst_;
    typename std::list<EdgeWPtr>::iterator vertex_dst_list_it_;
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
  template <typename GT>
  class Facet
    : public std::enable_shared_from_this<Facet<GT> >
  {
    friend class Polyhedron;

  public:
    class FacetData
    {
      using FacetWPtr = std::weak_ptr<Facet<GT> >;
      using FacetSPtr = std::shared_ptr<Facet<GT> >;

      using FacetDataSPtr = std::shared_ptr<FacetData>;

    public:
      FacetData() { /*intentionally does nothing*/ }
      virtual ~FacetData() { /*intentionally does nothing*/ }

      static FacetDataSPtr create(const FacetSPtr& facet)
      {
        FacetDataSPtr result = std::make_shared<FacetData>();
        facet->setData(result);
        return result;
      }

      virtual FacetDataSPtr clone() const {
        return std::make_shared<FacetData>(*this);
      }
    };

    class SkelFacetData
      : public FacetData
    {
      using FT = typename GT::FT;

      using FacetWPtr = std::weak_ptr<Facet<GT> >;
      using FacetSPtr = std::shared_ptr<Facet<GT> >;

      using FacetDataSPtr = std::shared_ptr<FacetData>;
      using SkelFacetDataSPtr = std::shared_ptr<SkelFacetData>;

    public:
      SkelFacetData() { /*intentionally does nothing*/ }
      virtual ~SkelFacetData() { /*intentionally does nothing*/ }

      static SkelFacetDataSPtr create(const FacetSPtr& facet)
      {
        SkelFacetDataSPtr result = std::make_shared<SkelFacetData>();
        result->setFacetOrigin(facet);
        facet->setData(result);
        return result;
      }

      virtual FacetDataSPtr clone() const override {
        return std::make_shared<SkelFacetData>(*this);
      }

      FacetSPtr getFacetOrigin() const
      {
        CGAL_SS3_DEBUG_WPTR(facet_origin_);
        return this->facet_origin_.lock();
      }

      void setFacetOrigin(const FacetSPtr& facet_origin)
      {
        this->facet_origin_ = facet_origin;
      }

      const FT& getSpeed() const
      {
        CGAL_assertion(speed_ > 0);
        return speed_;
      }

      void setSpeed(const FT& speed)
      {
        speed_ = speed;
      }

      Plane3SPtr getBasePlane() const
      {
        CGAL_precondition(bool(this->base_plane_));
        return this->base_plane_;
      }

      void setBasePlane(Plane3SPtr plane)
      {
        this->base_plane_ = plane;
      }

      bool hasFinalPlane() const
      {
        return bool(this->final_plane_);
      }

      Plane3SPtr getFinalPlane() const
      {
        CGAL_precondition(bool(this->final_plane_));
        return this->final_plane_;
      }

      void setFinalPlane(Plane3SPtr plane)
      {
        this->final_plane_ = plane;
      }

    protected:
      FacetWPtr facet_origin_;
      FT speed_;

      Plane3SPtr base_plane_;
      Plane3SPtr final_plane_;
    };

private:
    using FT = typename GT::FT;
    using Point_3 = typename GT::Point_3;
    using Plane_3 = typename GT::Plane_3;

    // using VertexWPtr = std::weak_ptr<Vertex<GT> >;
    using VertexSPtr = std::shared_ptr<Vertex<GT> >;
    // using EdgeWPtr = std::weak_ptr<Edge<GT> >;
     using EdgeSPtr = std::shared_ptr<Edge<GT> >;
    // using FacetWPtr = std::weak_ptr<Facet<GT> >;
    using FacetSPtr = std::shared_ptr<Facet<GT> >;

    using FacetDataSPtr = std::shared_ptr<FacetData>;
    using SkelFacetDataSPtr = std::shared_ptr<SkelFacetData>;

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
        result->addVertex(vertices[i]);
      }
      for (std::size_t i=0; i<nv; ++i) {
        EdgeSPtr edge = vertices[i]->findEdge(vertices[(i+1)%nv]);
        if (!edge) {
          edge = Edge<GT>::create(vertices[i], vertices[(i+1)%nv]);
        }
        result->addEdge(edge);
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
        result->addEdge(edges[i]);
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
        result->addVertex(vertex_c);
      }
      typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
      while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr src = old_to_new.at(edge->getVertexSrc());
        VertexSPtr dst = old_to_new.at(edge->getVertexDst());
        EdgeSPtr edge_c = Edge<GT>::create(src, dst);
        CGAL_assertion(edge_c->getVertexSrc() == src && edge_c->getVertexDst() == dst);
        if (edge->getFacetL() == this->shared_from_this()) {
          edge_c->setFacetL(result);
        }
        if (edge->getFacetR() == this->shared_from_this()) {
          edge_c->setFacetR(result);
        }
        result->addEdge(edge_c);
      }
      result->setID(getID());
      if (hasData()) {
        result->setData(getData()->clone());
      }
      return result;
    }

    void addVertex(const VertexSPtr& vertex)
    {
      vertices_.insert(vertices_.end(), vertex);
      vertex->addFacet(this->shared_from_this());
    }

    bool removeVertex(const VertexSPtr& vertex)
    {
      bool result = false;
      typename std::list<VertexSPtr>::iterator it_v =
        std::find(vertices_.begin(), vertices_.end(), vertex);
      if (it_v != vertices_.end()) {
        vertices_.erase(it_v);
        vertex->removeFacet(this->shared_from_this());
        result = true;
      }
      return result;
    }

    bool hasVertex(const VertexSPtr& vertex)
    {
      typename std::list<VertexSPtr>::iterator it_v =
        std::find(vertices_.begin(), vertices_.end(), vertex);
      return (it_v != vertices_.end());
    }

    bool isTriangle() const
    {
      return (vertices_.size() == 3);
    }

    void addEdge(const EdgeSPtr& edge)
    {
      typename std::list<EdgeSPtr>::iterator it = edges_.insert(edges_.end(), edge);
      FacetSPtr facet_l = edge->getFacetL();
      FacetSPtr facet_r = edge->getFacetR();
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
        edge->setFacetL(this->shared_from_this());
        edge->setFacetLListIt(it);
      } else if (!facet_r) {
        edge->setFacetR(this->shared_from_this());
        edge->setFacetRListIt(it);
      } else {
        CGAL_SS3_HDS_TRACE(edge->toString());
        throw std::runtime_error("The given edge already has a left and a right facet.");
      }
      VertexSPtr vertex_src = edge->src(this->shared_from_this());
      if (!containsVertex(vertex_src)) {
        addVertex(vertex_src);
      }
      VertexSPtr vertex_dst = edge->dst(this->shared_from_this());
      if (!containsVertex(vertex_dst)) {
        addVertex(vertex_dst);
      }
    }

    bool removeEdge(const EdgeSPtr& edge)
    {
      bool result = false;
      if (edge->getFacetL() == this->shared_from_this()) {
        edges_.erase(edge->getFacetLListIt());
        edge->setFacetL(FacetSPtr());
        edge->setFacetLListIt(typename std::list<EdgeSPtr>::iterator());
        result = true;
      } else if (edge->getFacetR() == this->shared_from_this()) {
        edges_.erase(edge->getFacetRListIt());
        edge->setFacetR(FacetSPtr());
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
    EdgeSPtr findEdge(const FacetSPtr& facet) const
    {
      EdgeSPtr result = EdgeSPtr();
      typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
      while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->getFacetL() == facet || edge->getFacetR() == facet) {
          result = edge;
          break;
        }
      }
      return result;
    }

    std::list<EdgeSPtr> findEdges(const FacetSPtr& facet) const
    {
      std::list<EdgeSPtr> result;
      typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
      while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->getFacetL() == facet || edge->getFacetR() == facet) {
          result.push_back(edge);
        }
      }
      return result;
    }

    bool containsVertex(const VertexSPtr& vertex) const
    {
      bool result = (vertices_.end() != std::find(vertices_.begin(), vertices_.end(), vertex));
      return result;
    }

    bool containsEdge(const EdgeSPtr& edge) const
    {
      bool result = (edges_.end() != std::find(edges_.begin(), edges_.end(), edge));
      return result;
    }

    void sortVertices()
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
          vertex = edge->dst(self);
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

    void sortEdges()
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
          removeEdge(edge);
        }
      }
      edges_.clear();
      it_e_tmp = tmp.begin();
      while (it_e_tmp != tmp.end()) {
        EdgeSPtr edge = *it_e_tmp++;
        addEdge(edge);
      }
    }

    PolyhedronSPtr getPolyhedron() const
    {
      return this->polyhedron_.lock();
    }

    void setPolyhedron(PolyhedronSPtr polyhedron)
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

    FacetDataSPtr getData() const
    {
      CGAL_SS3_DEBUG_SPTR(this->data_);
      return this->data_;
    }

    void setData(FacetDataSPtr data)
    {
      this->data_ = data;
    }

    bool hasData() const
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
                FacetSPtr facet_l = edge->getFacetL();
                FacetSPtr facet_r = edge->getFacetR();
                if ((facet_l == this->shared_from_this() &&
                     facet_r == facet &&
                     edge->getVertexDst() == vertex) ||
                    (facet_r == this->shared_from_this() &&
                     facet_l == facet &&
                     edge->getVertexSrc() == vertex)) {
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
                FacetSPtr facet_l = edge->getFacetL();
                FacetSPtr facet_r = edge->getFacetR();
                if ((facet_l == this->shared_from_this() &&
                     facet_r == facet &&
                     edge->getVertexSrc() == vertex) ||
                    (facet_r == this->shared_from_this() &&
                     facet_l == facet &&
                     edge->getVertexDst() == vertex)) {
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
        facet->removeVertex(vertex);
        if (!containsVertex(vertex)) {
          addVertex(vertex);
        }
      }
      typename std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
      while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if ((facet_l == facet && facet_r == this->shared_from_this()) ||
            (facet_r == facet && facet_l == this->shared_from_this())) {
          facet->removeEdge(edge);
          removeEdge(edge);
          edge->getPolyhedron()->removeEdge(edge);
        } else {
          if (facet_l == facet) {
            facet->removeEdge(edge);
            edge->setFacetL(this->shared_from_this());
            addEdge(edge);
          }
          if (facet_r == facet) {
            facet->removeEdge(edge);
            edge->setFacetR(this->shared_from_this());
            addEdge(edge);
          }
        }
      }

      CGAL_assertion(getPolyhedron()->isConsistent());
    }

    int getID() const
    {
      return this->id_;
    }

    void setID(int id)
    {
      this->id_ = id;
    }

    /**
    * The direction of the normal points to the outside.
    */
    Plane3SPtr getPlane() const
    {
      CGAL_precondition(bool(this->plane_));
      return this->plane_;
    }

    void setPlane(Plane3SPtr plane)
    {
      this->plane_ = plane;
    }

    /**
    * First vertices have to form a triangle that is inside.
    */
    bool initPlane()
    {
      bool result = false;

      CGAL_SS3_HDS_TRACE("initPlane() of " << getID());

      Point3SPtr point_prev;
      std::vector<Point3SPtr> points;
      typename std::list<VertexSPtr>::iterator it_v = vertices_.begin();
      while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr point = vertex->getPoint();
        if (point_prev != point) {
          points.push_back(point);
        }
        point_prev = point;
      }

      CGAL_SS3_HDS_TRACE("computing normals from " << points.size() << " points");
      CGAL_SS3_HDS_TRACE_CODE(for (std::size_t i=0; i<points.size(); ++i))
      CGAL_SS3_HDS_TRACE("point " << i << ": " << *points[i]);

      if (points.size() >= 3) {
        Point3SPtr p0 = points[0];
        Vector3SPtr normal = KernelFactory::createVector3(0,0,0);
        std::size_t last_i = points.size() - 1;
        for (std::size_t i=1; i<last_i; ++i) {
          Point3SPtr p1 = points[i];
          Point3SPtr p2 = points[i+1];
          CGAL_assertion(p1 && p2);
          *normal += 0.5 * CGAL::cross_product(*p2 - *p1, *p0 - *p1);
        }

        if (*normal != CGAL::NULL_VECTOR) {
          plane_ = KernelFactory::createPlane3(p0, normal);
          result = true;
        }
      }

      CGAL_postcondition(result);

      return result;
    }

    Plane3SPtr plane()
    {
      if (!this->plane_) {
        this->initPlane();
      }
      CGAL_SS3_DEBUG_SPTR(this->plane_);
      return this->plane_;
    }

    bool makeFirstConvex()
    {
      bool result = false;
      if (!plane_) {
        return false;
      }
      EdgeSPtr edge_begin;
      FacetSPtr self(this->shared_from_this());
      Vector3SPtr normal = KernelFactory::createVector3(plane_);
      EdgeSPtr edge = edges_.front();
      EdgeSPtr first = EdgeSPtr();
      while (edge != first) {
        if (!first) {
          first = edge;
        }
        EdgeSPtr edge_next = edge->next(self);
        Point3SPtr points[3];
        points[0] = edge->src(self)->getPoint();
        points[1] = edge->dst(self)->getPoint();
        points[2] = edge_next->dst(self)->getPoint();

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
        if (edge_begin->getFacetL() == this->shared_from_this()) {
          edges_.erase(edge->getFacetLListIt());
          edge_begin->setFacetLListIt(it_e);
        } else if (edge_begin->getFacetR() == this->shared_from_this()) {
          edges_.erase(edge->getFacetRListIt());
          edge_begin->setFacetRListIt(it_e);
        }
        sortEdges();
        VertexSPtr vertex_begin = edge_begin->src(this->shared_from_this());
        typename std::list<VertexSPtr>::iterator it_v =
          std::find(vertices_.begin(), vertices_.end(), vertex_begin);
        vertices_.erase(it_v);
        vertices_.insert(vertices_.begin(), vertex_begin);
        sortVertices();
      }
      if (!result) {
        CGAL_SS3_HDS_TRACE("Warning: Unable to make first 3 vertices convex.");
        CGAL_SS3_HDS_TRACE(toString());
      }

      return result;
    }

    /**
    * returns the number of vertices whose degree is higher than 3
    */
    int numHighDegreeVertices() const
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

    std::string toString() const
    {
      std::stringstream sstr;
      sstr << "Facet(";
      sstr << "id=" << IO::StringFactory::fromInteger(id_) << ", ";
      // sstr << "addr=" << IO::StringFactory::fromPointer(this) << ", ";
      if (plane_) {
        sstr << "Plane: <" << IO::StringFactory::fromDouble(CGAL::to_double(plane_->a())) << ", "
             << IO::StringFactory::fromDouble(CGAL::to_double(plane_->b())) << ", "
             << IO::StringFactory::fromDouble(CGAL::to_double(plane_->c())) << ", "
             << IO::StringFactory::fromDouble(CGAL::to_double(plane_->d())) << ">, ";
      }

      if (hasData()) {
        sstr << "Speed: " << std::dynamic_pointer_cast<SkelFacetData>(getData())->getSpeed() << ", ";
      }

      sstr << "Vertices:" + IO::StringFactory::fromInteger(vertices_.size()) + ", ";
      sstr << "Edges:" + IO::StringFactory::fromInteger(edges_.size()) + ",";
      if (vertices_.size() > 0) {
        sstr << "\n";
        typename std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
        while (it_v != vertices_.end()) {
          VertexSPtr vertex = *it_v++;
          sstr << vertex->toString() << "\n";
        }
      }
      if (edges_.size() > 0) {
        sstr << "\n";
        typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
          while (it_e != edges_.end()) {
          EdgeSPtr edge = *it_e++;
          sstr << edge->toString() << "\n";
        }
      }
      sstr << ") END FACET TOSTRING()\n";

      return sstr.str();
    }

  protected:
    Plane3SPtr plane_;
    int id_;
    FacetDataSPtr data_;

    std::list<VertexSPtr> vertices_;
    std::list<EdgeSPtr> edges_;
    PolyhedronWPtr polyhedron_;
    typename std::list<FacetSPtr>::iterator polyhedron_list_it_;
  };

public:
  using VertexWPtr = std::weak_ptr<Vertex<Traits> >;
  using VertexSPtr = std::shared_ptr<Vertex<Traits> >;
  using EdgeWPtr = std::weak_ptr<Edge<Traits> >;
  using EdgeSPtr = std::shared_ptr<Edge<Traits> >;
  using FacetWPtr = std::weak_ptr<Facet<Traits> >;
  using FacetSPtr = std::shared_ptr<Facet<Traits> >;

  using VertexData = typename Vertex<Traits>::VertexData;
  using VertexDataSPtr = std::shared_ptr<VertexData>;
  using SkelVertexData = typename Vertex<Traits>::SkelVertexData;
  using SkelVertexDataSPtr = std::shared_ptr<SkelVertexData>;
  using EdgeData = typename Edge<Traits>::EdgeData;
  using EdgeDataSPtr = std::shared_ptr<EdgeData>;
  using SkelEdgeData = typename Edge<Traits>::SkelEdgeData;
  using SkelEdgeDataSPtr = std::shared_ptr<SkelEdgeData>;
  using FacetData = typename Facet<Traits>::FacetData;
  using FacetDataSPtr = std::shared_ptr<FacetData>;
  using SkelFacetData = typename Facet<Traits>::SkelFacetData;
  using SkelFacetDataSPtr = std::shared_ptr<SkelFacetData>;

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
      result->addFacet(facets[i]);
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
      Point3SPtr point = vertex->getPoint();
      VertexSPtr vertex_c = Vertex<Traits>::create(point);
      vertex_c->setID(vertex->getID());
      if (vertex->hasData()) {
        vertex_c->setData(vertex->getData()->clone());
      }
      result->addVertex(vertex_c);
      vertices_c[vertex] = vertex_c;
    }
    typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
      EdgeSPtr edge = *it_e++;
      VertexSPtr src = vertices_c[edge->getVertexSrc()];
      VertexSPtr dst = vertices_c[edge->getVertexDst()];
      EdgeSPtr edge_c = Edge<Traits>::create(src, dst);
      edge_c->setID(edge->getID());
      if (edge->hasData()) {
        edge_c->setData(edge->getData()->clone());
      }
      result->addEdge(edge_c);
      edges_c[edge] = edge_c;
    }
    typename std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
      FacetSPtr facet = *it_f++;
      FacetSPtr facet_c = Facet<Traits>::create();
      facet_c->setPlane(facet->getPlane());
      facet_c->setID(facet->getID());
      typename std::list<VertexSPtr>::const_iterator it_v = facet->vertices().begin();
      while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        facet_c->addVertex(vertices_c[vertex]);
      }
      typename std::list<EdgeSPtr>::const_iterator it_e = facet->edges().begin();
      while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        EdgeSPtr edge_c = edges_c[edge];
        if (edge->getFacetL() == facet) {
          edge_c->setFacetL(facet_c);
        }
        if (edge->getFacetR() == facet) {
          edge_c->setFacetR(facet_c);
        }
        facet_c->addEdge(edge_c);
      }
      if (facet->hasData()) {
        facet_c->setData(facet->getData()->clone());
      }
      result->addFacet(facet_c);
    }
    return result;
  }

  void addVertex(const VertexSPtr& vertex)
  {
    vertex->setID(next_vertex_id_++);
    typename std::list<VertexSPtr>::iterator it = vertices_.insert(vertices_.end(), vertex);
    vertex->setPolyhedron(this->shared_from_this());
    vertex->setPolyhedronListIt(it);
  }

  bool removeVertex(const VertexSPtr& vertex)
  {
    bool result = false;
    if (vertex->getPolyhedron() == this->shared_from_this()) {
      vertices_.erase(vertex->getPolyhedronListIt());
      vertex->setPolyhedronListIt(typename std::list<VertexSPtr>::iterator());
      vertex->setPolyhedron(PolyhedronSPtr());
      typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
      while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          this->removeFacet(facet);
        }
      }
      typename std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
      while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (EdgeSPtr edge = edge_wptr.lock()) {
          this->removeEdge(edge);
        }
      }
      result = true;
    }
    return result;
  }

  /**
   * Searches for a vertex with the same coordinates as the given vertex.
   */
  VertexSPtr findVertex(const VertexSPtr& needle)
  {
    VertexSPtr result = VertexSPtr();
    typename std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
      VertexSPtr vertex = *it_v++;
      if (vertex->getPoint() == needle->getPoint()) {
        result = vertex;
        break;
      }
    }
    return result;
  }

  void addEdge(const EdgeSPtr& edge)
  {
    edge->setID(next_edge_id_++);
    typename std::list<EdgeSPtr>::iterator it = edges_.insert(edges_.end(), edge);
    edge->setPolyhedron(this->shared_from_this());
    edge->setPolyhedronListIt(it);
    VertexSPtr vertex = edge->getVertexSrc();
    if (vertex->getPolyhedron() != this->shared_from_this()) {
      this->addVertex(vertex);
    }
    vertex = edge->getVertexDst();
    if (vertex->getPolyhedron() != this->shared_from_this()) {
      this->addVertex(vertex);
    }
  }

  bool removeEdge(const EdgeSPtr& edge)
  {
    bool result = false;
    if (edge->getPolyhedron() == this->shared_from_this()) {
      edges_.erase(edge->getPolyhedronListIt());
      edge->setPolyhedronListIt(typename std::list<EdgeSPtr>::iterator());
      edge->setPolyhedron(PolyhedronSPtr());
      FacetSPtr facet = edge->getFacetL();
      if (facet) {
        this->removeFacet(facet);
      }
      facet = edge->getFacetR();
      if (facet) {
        this->removeFacet(facet);
      }
      edge->getVertexSrc()->removeEdge(edge);
      edge->getVertexDst()->removeEdge(edge);
      result = true;
    }
    return result;
  }

  /**
   * Searches for an edge with the same coordinates as the given edge.
   * The orientation of the edge is ignored.
   */
  EdgeSPtr findEdge(const EdgeSPtr& needle)
  {
    EdgeSPtr result = EdgeSPtr();
    typename std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
      EdgeSPtr edge = *it_e++;
      if (edge->getVertexSrc()->getPoint() == needle->getVertexSrc()->getPoint() &&
          edge->getVertexDst()->getPoint() == needle->getVertexDst()->getPoint()) {
        result = edge;
        break;
      }
      if (edge->getVertexSrc()->getPoint() == needle->getVertexDst()->getPoint() &&
          edge->getVertexDst()->getPoint() == needle->getVertexSrc()->getPoint()) {
        result = edge;
        break;
      }
    }
    return result;
  }

  void addFacet(const FacetSPtr& facet)
  {
    facet->setID(next_facet_id_++);
    typename std::list<FacetSPtr>::iterator it_f = facets_.insert(facets_.end(), facet);
    facet->setPolyhedronListIt(it_f);
    facet->setPolyhedron(this->shared_from_this());
    // add content of facet
    typename std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
    while (it_v != facet->vertices().end()) {
      VertexSPtr vertex = *it_v++;
      if (vertex->getPolyhedron() != this->shared_from_this()) {
        this->addVertex(vertex);
      }
    }
    typename std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
      EdgeSPtr edge = *it_e++;
      if (edge->getPolyhedron() != this->shared_from_this()) {
        this->addEdge(edge);
      }
    }
  }

  bool removeFacet(const FacetSPtr& facet)
  {
    bool result = false;
    if (facet->getPolyhedron() == this->shared_from_this()) {
      facets_.erase(facet->getPolyhedronListIt());
      facet->setPolyhedron(PolyhedronSPtr());
      facet->setPolyhedronListIt(typename std::list<FacetSPtr>::iterator());
      typename std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
      while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        facet->removeEdge(edge);  // clears list iterators
      }
      typename std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
      while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        facet->removeVertex(vertex);
      }
      result = true;
    }
    return result;
  }

  void initPlanes()
  {
    for (const FacetSPtr& facet : facets_) {
      facet->initPlane();
    }
  }

  void clearData()
  {
    for (const VertexSPtr& vertex : vertices_) {
      vertex->setData(VertexDataSPtr());
    }
    for (const EdgeSPtr& edge : edges_) {
      edge->setData(EdgeDataSPtr());
    }
    for (const FacetSPtr& facet : facets_) {
      facet->setData(FacetDataSPtr());
    }
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

  bool isConsistent() const
  {
    bool result = true;

    typename std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
      VertexSPtr vertex = *it_v++;
      if (vertex->getPolyhedron() != this->shared_from_this()) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->toString());
        result = false;
        break;
      }
      typename std::list<EdgeWPtr>::const_iterator it_e = vertex->edges().begin();
      while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (vertex != edge->getVertexSrc() && vertex != edge->getVertexDst()) {
            CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->toString());
            CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
            result = false;
            break;
          }

          if (edge->getVertexSrc() == edge->getVertexDst())
          {
            CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->getVertexSrc());
            result = false;
            break;
          }
        } else {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->toString());
        }
      }

      typename std::list<FacetWPtr>::const_iterator it_f = vertex->facets().begin();
      while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          if (!facet->containsVertex(vertex)) {
            CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->toString());
            CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << facet->toString());
            result = false;
            break;
          }
        } else {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->toString());
        }
      }
    }

    typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
      EdgeSPtr edge = *it_e++;
      if (edge->getPolyhedron() != this->shared_from_this()) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
        result = false;
        break;
      }
      if (!edge->getVertexSrc()->containsEdge(edge)) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->getVertexSrc()->toString());
        result = false;
        break;
      }
      if (!edge->getVertexDst()->containsEdge(edge)) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->getVertexDst()->toString());
        result = false;
        break;
      }
      EdgeWPtr edge_wptr;
      edge_wptr = *(edge->getVertexSrcListIt());
      if (edge_wptr.lock() != edge) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
      }
      edge_wptr = *(edge->getVertexDstListIt());
      if (edge_wptr.lock() != edge) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
      }
      if (edge->getFacetL()) {
        if (!edge->getFacetL()->containsEdge(edge)) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->getFacetL()->toString());
          result = false;
          break;
        }
        edge_wptr = *(edge->getFacetLListIt());
        if (edge_wptr.lock() != edge) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
        }
      } else {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
        result = false;
      }
      if (edge->getFacetR()) {
        if (!edge->getFacetR()->containsEdge(edge)) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->getFacetR()->toString());
          result = false;
          break;
        }
        edge_wptr = *(edge->getFacetRListIt());
        if (edge_wptr.lock() != edge) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
        }
      } else {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
        result = false;
      }
    }

    typename std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
      FacetSPtr facet = *it_f++;
      if (facet->getPolyhedron() != this->shared_from_this()) {
        CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << facet->toString());
        result = false;
        break;
      }
      typename std::list<VertexSPtr>::const_iterator it_v = facet->vertices().begin();
      while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (!vertex->containsFacet(facet)) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << facet->toString());
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << vertex->toString());
          result = false;
          break;
        }
      }
      typename std::list<EdgeSPtr>::const_iterator it_e = facet->edges().begin();
      while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->getFacetL() != facet && edge->getFacetR() != facet) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << facet->toString());
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
          result = false;
          break;
        }
        if (!facet->containsVertex(edge->getVertexSrc()) ||
            !facet->containsVertex(edge->getVertexDst())) {
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << facet->toString());
          CGAL_SS3_HDS_TRACE("Inconsistency @ L" << __LINE__ << "\n" << edge->toString());
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
      this->removeFacet(facet);
    }
    typename std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
      EdgeSPtr edge = *it_e++;
      this->removeEdge(edge);
    }
    typename std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
      VertexSPtr vertex = *it_v++;
      this->removeVertex(vertex);
    }
  }

  bool empty()
  {
    return vertices_.empty() && edges_.empty() && facets_.empty();
  }

  int getID() const
  {
    return this->id_;
  }

  void setID(int id)
  {
    this->id_ = id;
  }

  void initializeAllIDs()
  {
    resetAllIDs();

    for (const VertexSPtr& vertex : vertices_) {
      vertex->setID(next_vertex_id_++);
    }
    for (const EdgeSPtr& edge : edges_) {
      edge->setID(next_edge_id_++);
    }
    for (const FacetSPtr& facet : facets_) {
      facet->setID(next_facet_id_++);
    }
    setID(-1);
  }

  void resetAllIDs()
  {
    for (const VertexSPtr& vertex : vertices_) {
      vertex->setID(-1);
    }
    for (const EdgeSPtr& edge : edges_) {
      edge->setID(-1);
    }
    for (const FacetSPtr& facet : facets_) {
      facet->setID(-1);
    }
    setID(-1);

    next_vertex_id_ = 0;
    next_edge_id_ = 0;
    next_facet_id_ = 0;
  }

  std::string toString() const
  {
    std::stringstream sstr;
    sstr << "Polyhedron(";
    if (id_ != -1) {
      sstr << "id=" << IO::StringFactory::fromInteger(id_) << ", ";
    } else {
      // sstr << IO::StringFactory::fromPointer(this) << ", ";
    }
    sstr << "Vertices:" + IO::StringFactory::fromInteger(vertices_.size()) + ", ";
    sstr << "Edges:" + IO::StringFactory::fromInteger(edges_.size()) + ", ";
    sstr << "Facets:" + IO::StringFactory::fromInteger(facets_.size()) + ",\n";

    typename std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
      FacetSPtr facet = *it_f++;
      sstr << facet->toString() << "\n";
    }
    sstr << ")";

  /*typename std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
      EdgeSPtr edge = *it_e++;
      sstr << edge->toString() << "\n";
    }
    typename std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
      VertexSPtr vertex = *it_v++;
      sstr << vertex->toString() << "\n";
    }*/
    return sstr.str();
  }

  void dumpEdges(const std::string& base_filename) const
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
        edge_out << "2 " << *(edge->getVertexSrc()->getPoint()) << " "
                         << *(edge->getVertexDst()->getPoint()) << std::endl;
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
