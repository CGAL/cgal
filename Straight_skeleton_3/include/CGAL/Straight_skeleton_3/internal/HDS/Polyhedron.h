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
 * @file   data/3d/Polyhedron.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef CGAL_STRAIGHT_SKELETON_INTERNAL_HDS_POLYHEDRON_H
#define CGAL_STRAIGHT_SKELETON_INTERNAL_HDS_POLYHEDRON_H

#include <CGAL/Straight_skeleton_3/internal/debug.h>

#include <CGAL/Straight_skeleton_3/internal/weak_find.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>
#include <CGAL/Straight_skeleton_3/Configuration.h>

#ifdef CGAL_SS3_DUMP_FILES
# include <CGAL/Straight_skeleton_3/IO/OBJ.h>
#endif

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
  using Sphere_3 = typename Traits::Sphere_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Segment3SPtr = std::shared_ptr<Segment_3>;
  using Vector3SPtr = std::shared_ptr<Vector_3>;
  using Line3SPtr = std::shared_ptr<Line_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;
  using Sphere3SPtr = std::shared_ptr<Sphere_3>;

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

      static VertexDataSPtr create(VertexSPtr vertex)
      {
        VertexDataSPtr result = std::make_shared<VertexData>();
        result->setVertex(vertex);
        vertex->setData(result);
        return result;
      }

      VertexSPtr getVertex() const
      {
        return this->vertex_.lock();
      }

      void setVertex(VertexSPtr vertex)
      {
        this->vertex_ = vertex;
      }

    protected:
      VertexWPtr vertex_;
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

      using SkelVertexDataSPtr = std::shared_ptr<SkelVertexData>;

    public:
      SkelVertexData() { /*intentionally does nothing*/ }
      virtual ~SkelVertexData() { /*intentionally does nothing*/ }

      static SkelVertexDataSPtr create(VertexSPtr vertex)
      {
        SkelVertexDataSPtr result = std::make_shared<SkelVertexData>();
        result->setVertex(vertex);
        vertex->setData(result);
        return result;
      }

      ArcSPtr getArc() const
      {
        CGAL_SS3_DEBUG_WPTR(arc_);
        return this->arc_.lock();
      }

      void setArc(ArcSPtr arc)
      {
        this->arc_ = arc;
      }

      NodeSPtr getNode() const
      {
        return this->node_.lock();
      }

      void setNode(NodeSPtr node)
      {
        this->node_ = node;
      }

      VertexSPtr getOffsetVertex() const
      {
        CGAL_SS3_DEBUG_WPTR(offset_vertex_);
        return this->offset_vertex_.lock();
      }

      void setOffsetVertex(VertexSPtr offset_vertex)
      {
        this->offset_vertex_ = offset_vertex;
      }

    protected:
      ArcWPtr arc_;
      NodeWPtr node_;
      VertexWPtr offset_vertex_;
    };

  private:
    using FT = typename GT::FT;
    using Point_3 = typename GT::Point_3;

    using Point3SPtr = std::shared_ptr<Point_3>;

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
      return std::make_shared<Vertex>(point);
    }

    VertexSPtr clone() const
    {
      return std::make_shared<Vertex>(getPoint());
    }

    Point3SPtr getPoint() const
    {
      CGAL_SS3_DEBUG_SPTR(point_);
      return point_;
    }

    void setPoint(Point3SPtr point)
    {
      point_ = point;
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

    int getID() const
    {
      return id_;
    }

    void setID(int id)
    {
      id_ = id;
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

    bool isReflex() const
    {
      if (degree() == 0) {
          return false;
      }
      bool result = true;
      for (EdgeWPtr edge_wptr : edges_) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (!edge->isReflex()) {
            result = false;
          }
        }
      }
      return result;
    }

    bool isConvex() const
    {
      if (degree() == 0) {
        return false;
      }
      bool result = true;
      for (EdgeWPtr edge_wptr : edges_) {
        if (EdgeSPtr edge = edge_wptr.lock()) {
          if (edge->isReflex()) {
            result = false;
          }
        }
      }
      return result;
    }

    void addEdge(EdgeSPtr edge)
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

    bool removeEdge(EdgeSPtr edge)
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
      EdgeSPtr result;
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
    EdgeSPtr findEdge(VertexSPtr dst) const
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
    EdgeSPtr findEdge(FacetSPtr facet) const
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

    void addFacet(FacetSPtr facet)
    {
      facets_.insert(facets_.end(), FacetWPtr(facet));
    }

    bool removeFacet(FacetSPtr facet)
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
      FacetSPtr result;
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

    bool containsEdge(EdgeSPtr edge) const
    {
      EdgeWPtr edge_wptr = EdgeWPtr(edge);
      bool result = (edges_.end() !=
          STL_Extension::internal::weak_find(edges_.begin(), edges_.end(), edge_wptr));
      return result;
    }

    bool containsFacet(FacetSPtr facet) const
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

    std::list<EdgeWPtr>& edges()
    {
      return edges_;
    }

    std::list<FacetWPtr>& facets()
    {
      return facets_;
    }

    VertexSPtr next(FacetSPtr facet) const
    {
      VertexSPtr result = VertexSPtr(); // @fixme nullptr & whatnot

      if (facet->vertices().size() == 1) { // @fixme what's the point of this loop...
        VertexSPtr vertex = *facet->vertices().begin();
        CGAL_assertion(vertex == this->shared_from_this());
        return vertex;
      }

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

# if 0
      VertexSPtr result = VertexSPtr();
      std::list<VertexSPtr>::const_iterator it_v = facet->vertices().begin();
      while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v;
        if (vertex == this->shared_from_this()) {
          if (facet->vertices().size() == 1) { // @fixme what's the point of this loop...
            result = vertex;
          }
          break;
        }
        it_v++;
      }
      if (it_v != facet->vertices().end()) {
        std::list<VertexSPtr>::const_iterator it_v_begin = it_v++;
        if (it_v == facet->vertices().end()) {
          it_v = facet->vertices().begin();
        }
        while (it_v != it_v_begin) {
          VertexSPtr vertex = *it_v++;
          if (it_v == facet->vertices().end()) {
            it_v = facet->vertices().begin();
          }
          EdgeSPtr edge = findEdge(vertex);
          if (edge) {
            if (edge->dst(facet) == vertex) {
              result = vertex;
              break;
            }
          }
        }
      }
      CGAL_SS3_DEBUG_SPTR(result);
      return result;
#endif
  }

    VertexSPtr prev(FacetSPtr facet) const
    {
      VertexSPtr result = VertexSPtr();
      typename std::list<VertexSPtr>::const_reverse_iterator it_v = facet->vertices().rbegin();
      while (it_v != facet->vertices().rend()) {
        VertexSPtr vertex = *it_v;
        if (vertex == this->shared_from_this()) {
          if (facet->vertices().size() == 1) {
            result = vertex;
          }
          break;
        }
        ++it_v;
      }
      if (it_v != facet->vertices().rend()) {
        typename std::list<VertexSPtr>::const_reverse_iterator it_v_begin = it_v++;
        if (it_v == facet->vertices().rend()) {
          it_v = facet->vertices().rbegin();
        }
        while (it_v != it_v_begin) {
          VertexSPtr vertex = *it_v++;
          if (it_v == facet->vertices().rend()) {
            it_v = facet->vertices().rbegin();
          }
          EdgeSPtr edge = findEdge(vertex);
          if (edge) {
            if (edge->src(facet) == vertex) {
              result = vertex;
              break;
            }
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
    VertexSPtr split(FacetSPtr facet_left,
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
    Point3SPtr final_point_;

    std::list<EdgeWPtr> edges_;
    std::list<FacetWPtr> facets_;
    PolyhedronWPtr polyhedron_;
    typename std::list<VertexSPtr>::iterator polyhedron_list_it_;

    VertexDataSPtr data_;
    int id_;
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

      static EdgeDataSPtr create(EdgeSPtr edge)
      {
        EdgeDataSPtr result = std::make_shared<EdgeData>();
        result->setEdge(edge);
        edge->setData(result);
        return result;
      }

      EdgeSPtr getEdge() const
      {
        return this->edge_.lock();
      }
      void setEdge(EdgeSPtr edge)
      {
        this->edge_ = edge;
      }

    protected:
      EdgeWPtr edge_;
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

      using SkelEdgeDataSPtr = std::shared_ptr<SkelEdgeData>;

    public:
      SkelEdgeData() { /*intentionally does nothing*/ }
      virtual ~SkelEdgeData() { /*intentionally does nothing*/ }

      static SkelEdgeDataSPtr create(EdgeSPtr edge)
      {
        SkelEdgeDataSPtr result = std::make_shared<SkelEdgeData>();
        result->setEdge(edge);
        edge->setData(result);
        return result;
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

      EdgeSPtr getOffsetEdge() const
      {
        CGAL_SS3_DEBUG_WPTR(offset_edge_);
        return this->offset_edge_.lock();
      }

      void setOffsetEdge(EdgeSPtr offset_edge)
      {
        this->offset_edge_ = offset_edge;
      }

      /**
      * Used by WeightVertexSplitter (intersection with a plane)
      */
      FacetSPtr getFacetOrigin() const
      {
        CGAL_SS3_DEBUG_WPTR(facet_origin_);
        return this->facet_origin_.lock();
      }

      void setFacetOrigin(FacetSPtr facet_origin)
      {
        this->facet_origin_ = facet_origin;
      }

    protected:
      SheetWPtr sheet_;
      EdgeWPtr offset_edge_;
      FacetWPtr facet_origin_;
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
    Edge(VertexSPtr src, VertexSPtr dst)
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

    static EdgeSPtr create(VertexSPtr src, VertexSPtr dst)
    {
      EdgeSPtr result = std::make_shared<Edge>(src, dst);
      src->addEdge(result);
      dst->addEdge(result);
      return result;
    }

    EdgeSPtr clone() const
    {
      EdgeSPtr result = EdgeSPtr(new Edge(*this));
      result->vertex_src_->addEdge(result);
      result->vertex_dst_->addEdge(result);
      return result;
    }

    VertexSPtr getVertexSrc() const
    {
      CGAL_SS3_DEBUG_SPTR(this->vertex_src_);
      return this->vertex_src_;
    }

    void setVertexSrc(VertexSPtr src)
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

    void setVertexDst(VertexSPtr dst)
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

    void setFacetL(FacetSPtr facet)
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

    void setFacetR(FacetSPtr facet)
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

    void setPolyhedronListIt(typename std::list<EdgeSPtr>::iterator list_it) {
      this->polyhedron_list_it_ = list_it;
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

    Segment3SPtr segment() const
    {
      return KernelFactory::createSegment3(vertex_src_->getPoint(), vertex_dst_->getPoint());
    }

    Line3SPtr line() const
    {
      return KernelFactory::createLine3(vertex_src_->getPoint(), vertex_dst_->getPoint());
    }

    VertexSPtr other(VertexSPtr vertex) const
    {
      VertexSPtr result = VertexSPtr();
      if (vertex == vertex_src_) {
        result = vertex_dst_;
      } else if (vertex == vertex_dst_) {
        result = vertex_src_;
      }
      return result;
    }

    FacetSPtr other(FacetSPtr facet) const
    {
      FacetSPtr result = FacetSPtr();
      if (facet == facet_l_.lock()) {
        result = facet_r_.lock();
      } else if (facet == facet_r_.lock()) {
        result = facet_l_.lock();
      }
      return result;
    }

    VertexSPtr src(FacetSPtr facet_l) const
    {
      VertexSPtr result = VertexSPtr();
      if (facet_l == facet_l_.lock()) {
        result = vertex_src_;
      } else if (facet_l == facet_r_.lock()) {
        result = vertex_dst_;
      }
      return result;
    }

    VertexSPtr dst(FacetSPtr facet_l) const
    {
      VertexSPtr result = VertexSPtr();
      if (facet_l == facet_l_.lock()) {
        result = vertex_dst_;
      } else if (facet_l == facet_r_.lock()) {
        result = vertex_src_;
      }
      return result;
    }

    FacetSPtr left(VertexSPtr vertex_src) const
    {
      FacetSPtr result = FacetSPtr();
      if (vertex_src == vertex_src_) {
        result = facet_l_.lock();
      } else if (vertex_src == vertex_dst_) {
        result = facet_r_.lock();
      }
      return result;
    }

    FacetSPtr right(VertexSPtr vertex_src) const
    {
      FacetSPtr result = FacetSPtr();
      if (vertex_src == vertex_src_) {
        result = facet_r_.lock();
      } else if (vertex_src == vertex_dst_) {
        result = facet_l_.lock();
      }
      return result;
    }

    EdgeSPtr next(FacetSPtr facet) const
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

#if 0
      EdgeSPtr result = EdgeSPtr();
      std::list<EdgeSPtr>::const_iterator it_e = facet->edges().begin();
      while (it_e != facet->edges().end()) {
          EdgeSPtr edge = *it_e;
          if (edge == this->shared_from_this()) {
            result = edge;
            break;
          }
          ++it_e;
      }
      if (it_e != facet->edges().end()) {
        VertexSPtr vertex_dst = this->dst(facet);
        std::list<EdgeSPtr>::const_iterator it_e_begin = it_e++;
        if (it_e == facet->edges().end()) {
          it_e = facet->edges().begin();
        }
        while (it_e != it_e_begin) {
          EdgeSPtr edge = *it_e++;
          if (it_e == facet->edges().end()) {
            it_e = facet->edges().begin();
          }
          if (vertex_dst == edge->src(facet)) {
            result = edge;
            break;
          }
        }
      }
      CGAL_SS3_DEBUG_SPTR(result);
      return result;
#endif
    }

    EdgeSPtr prev(FacetSPtr facet) const
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

#if 0
      EdgeSPtr result = EdgeSPtr();
      std::list<EdgeSPtr>::const_reverse_iterator it_e = facet->edges().rbegin();
      while (it_e != facet->edges().rend()) {
        EdgeSPtr edge = *it_e;
        if (edge == this->shared_from_this()) {
          result = edge;
          break;
        }
        ++it_e;
      }
      if (it_e != facet->edges().rend()) {
        VertexSPtr vertex_src = this->src(facet);
        std::list<EdgeSPtr>::const_reverse_iterator it_e_begin = it_e++;
        if (it_e == facet->edges().rend()) {
            it_e = facet->edges().rbegin();
        }
        while (it_e != it_e_begin) {
          EdgeSPtr edge = *it_e++;
          if (it_e == facet->edges().rend()) {
            it_e = facet->edges().rbegin();
          }
          if (vertex_src == edge->dst(facet)) {
            result = edge;
            break;
          }
        }
      }
      CGAL_SS3_DEBUG_SPTR(result);
      return result;
#endif
    }

    /**
      * counter clockwise from outside
      */
    EdgeSPtr next(VertexSPtr vertex) const
    {
      EdgeSPtr result = EdgeSPtr();
      FacetSPtr facet;
      if (vertex == vertex_src_) {
        facet = FacetSPtr(facet_l_);
      } else if (vertex == vertex_dst_) {
        facet = FacetSPtr(facet_r_);
      }
      if (facet) {
        std::list<EdgeSPtr> edges_possible;
        for (EdgeWPtr edge_wptr : vertex->edges()) {
          if (EdgeSPtr edge = edge_wptr.lock()) {
            if (edge.get() == this) {
              continue;
            }
            if (edge->dst(facet) == vertex) {
              edges_possible.push_back(edge);
            }
          }
        }
        if (edges_possible.size() == 1) {
          result = edges_possible.front();
        } else {
          double angle_min = 2*CGAL_PI;
          typename std::list<EdgeSPtr>::iterator it_e = edges_possible.begin();
          while (it_e != edges_possible.end()) {
            EdgeSPtr edge = *it_e++;
            double angle = angleTo(edge);
            if (angle == angle_min) {
              CGAL_SS3_HDS_TRACE("Warning: Not able to distinguish possible next edges.");
            }
            if (angle <= angle_min) {
              result = edge;
              angle_min = angle;
            }
          }
        }
      }
      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    EdgeSPtr prev(VertexSPtr vertex) const
    {
      EdgeSPtr result = EdgeSPtr();
      FacetSPtr facet;
      if (vertex == vertex_src_) {
        facet = FacetSPtr(facet_r_);
      } else if (vertex == vertex_dst_) {
        facet = FacetSPtr(facet_l_);
      }
      if (facet) {
        std::list<EdgeSPtr> edges_possible;
        for (EdgeWPtr edge_wptr : vertex->edges()) {
          if (EdgeSPtr edge = edge_wptr.lock()) {
            if (edge.get() == this) {
              continue;
            }
            if (edge->src(facet) == vertex) {
              edges_possible.push_back(edge);
            }
          }
        }
        if (edges_possible.size() == 1) {
          result = edges_possible.front();
        } else {
          double angle_max = 0.0;
          typename std::list<EdgeSPtr>::iterator it_e = edges_possible.begin();
          while (it_e != edges_possible.end()) {
            EdgeSPtr edge = *it_e++;
            double angle = angleTo(edge);
            if (angle == angle_max) {
              CGAL_SS3_HDS_TRACE("Warning: Not able to distinguish possible next edges.");
            }
            if (angle >= angle_max) {
              result = edge;
              angle_max = angle;
            }
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
    EdgeSPtr split(VertexSPtr middle)
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
    void replaceVertexSrc(VertexSPtr vertex_src)
    {
      vertex_src_->removeEdge(this->shared_from_this());
      vertex_src_ = vertex_src;
      vertex_src->addEdge(this->shared_from_this());
    }

    void replaceVertexDst(VertexSPtr vertex_dst)
    {
      vertex_dst_->removeEdge(this->shared_from_this());
      vertex_dst_ = vertex_dst;
      vertex_dst->addEdge(this->shared_from_this());
    }

    void replaceFacetL(FacetSPtr facet_l)
    {
      if (FacetSPtr facet = facet_l_.lock()) {
        facet->removeEdge(this->shared_from_this());
      }
      setFacetL(facet_l);
      facet_l->addEdge(this->shared_from_this());
    }

    void replaceFacetR(FacetSPtr facet_r)
    {
      if (FacetSPtr facet = facet_r_.lock()) {
        facet->removeEdge(this->shared_from_this());
      }
      setFacetR(facet_r);
      facet_r->addEdge(this->shared_from_this());
    }

    bool hasSameFacets(EdgeSPtr edge) const
    {
      bool result = (facet_l_.lock() == edge->getFacetL() &&
                     facet_r_.lock() == edge->getFacetR()) ||
                    (facet_r_.lock() == edge->getFacetL() &&
                     facet_l_.lock() == edge->getFacetR());
      return result;
    }

    int getID() const
    {
      return this->id_;
    }

    void setID(int id)
    {
      this->id_ = id;
    }

    double angle() const
    {
      double result = 0.0;
      FacetSPtr facet_l, facet_r;
      if ((facet_l = getFacetL()) && (facet_r = getFacetR())) {
        Vector3SPtr v1 = KernelFactory::createVector3(facet_l->plane());
        Vector3SPtr v2 = KernelFactory::createVector3(facet_r->plane());
  # ifdef USE_CGAL
        result = acos(CGAL::to_double(((*v1) * (*v2)) /
                  CGAL::disallowed_sqrt(v1->squared_length() * v2->squared_length())));
  # else
        result = acos(((*v1) * (*v2)) /
                  sqrt(v1->squared_length() * v2->squared_length()));
  # endif
        result = CGAL_PI - result;
        if (isReflex()) {
          result = 2.0*CGAL_PI - result;
        }
      } else {
        CGAL_SS3_HDS_TRACE("Warning: Not able to determine angle.");
        CGAL_SS3_HDS_TRACE(toString());
      }
      return result;
    }

    std::optional<bool> getReflexStatus() const
    {
      return cachedReflexStatus_;
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

    double angleTo(EdgeSPtr edge) const
    {
      double result = 0.0;
      VertexSPtr vertex;
      if (vertex_src_ == edge->getVertexSrc() ||
          vertex_src_ == edge->getVertexDst()) {
        vertex = vertex_src_;
      } else if (vertex_dst_ == edge->getVertexSrc() ||
                 vertex_dst_ == edge->getVertexDst()) {
        vertex = vertex_dst_;
      }
      FacetSPtr facet;
      FacetSPtr facet_l = getFacetL();
      FacetSPtr facet_r = getFacetR();
      if (facet_l == edge->getFacetL() ||
          facet_l == edge->getFacetR()) {
        facet = facet_l;
      } else if (facet_r == edge->getFacetL() ||
                 facet_r == edge->getFacetR()) {
        facet = facet_r;
      }
      if (vertex && facet) {
        Vector3SPtr normal = KernelFactory::createVector3(facet->plane());
        Vector3SPtr dir_self;
        if (vertex == vertex_src_) {
          dir_self = KernelFactory::createVector3(*(vertex_dst_->getPoint()) - *(vertex_src_->getPoint()));
        } else {
          dir_self = KernelFactory::createVector3(*(vertex_src_->getPoint()) - *(vertex_dst_->getPoint()));
        }
        Vector3SPtr dir_other;
        if (vertex == edge->getVertexSrc()) {
          dir_other = KernelFactory::createVector3(*(edge->getVertexDst()->getPoint()) - *(edge->getVertexSrc()->getPoint()));
        } else {
          dir_other = KernelFactory::createVector3(*(edge->getVertexSrc()->getPoint()) - *(edge->getVertexDst()->getPoint()));
        }
        double angle = acos(CGAL::to_double(((*dir_self) * (*dir_other)) /
          CGAL::disallowed_sqrt(dir_self->squared_length() * dir_other->squared_length())));

        if ((facet_l == edge->getFacetR() && facet_r == edge->getFacetL()) ||
            (facet_l == edge->getFacetL() && facet_r == edge->getFacetR())) {
          if (angle < CGAL_PI/2.0) {
            result = 0.0;
          } else {
            result = CGAL_PI;
          }
        } else {
          result = angle;
          double angle_normal = 0.0;

          Vector_3 crossprod = CGAL::cross_product(*dir_self, *dir_other);
          angle_normal = acos(CGAL::to_double(((*normal) * crossprod) /
            CGAL::disallowed_sqrt(normal->squared_length() * crossprod.squared_length())));

          if (angle_normal > CGAL_PI/2.0) {
            result += CGAL_PI;
          }
        }
      } else {
        CGAL_SS3_HDS_TRACE("Warning: Edges do not have a shared Vertex and a shared Facet.");
      }
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

    EdgeDataSPtr data_;

    mutable std::optional<bool> cachedReflexStatus_;

    int id_;
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

      static FacetDataSPtr create(FacetSPtr facet)
      {
        FacetDataSPtr result = std::make_shared<FacetData>();
        result->setFacet(facet);
        facet->setData(result);
        return result;
      }

      FacetSPtr getFacet() const
      {
        return this->facet_.lock();
      }

      void setFacet(FacetSPtr facet)
      {
        this->facet_ = facet;
      }

    protected:
      FacetWPtr facet_;
    };

    class SkelFacetData
      : public FacetData
    {
      using FT = typename GT::FT;

      using FacetWPtr = std::weak_ptr<Facet<GT> >;
      using FacetSPtr = std::shared_ptr<Facet<GT> >;

      using SkelFacetDataSPtr = std::shared_ptr<SkelFacetData>;

    public:
      SkelFacetData() { /*intentionally does nothing*/ }
      virtual ~SkelFacetData() { /*intentionally does nothing*/ }

      static SkelFacetDataSPtr create(FacetSPtr facet)
      {
        SkelFacetDataSPtr result = std::make_shared<SkelFacetData>();
        result->setFacet(facet);
        result->setFacetOrigin(facet);
        facet->setData(result);
        return result;
      }

      FacetSPtr getOffsetFacet() const
      {
        CGAL_SS3_DEBUG_WPTR(offset_facet_);
        return this->offset_facet_.lock();
      }

      void setOffsetFacet(FacetSPtr offset_facet)
      {
        this->offset_facet_ = offset_facet;
      }

      FacetSPtr getFacetOrigin() const
      {
        CGAL_SS3_DEBUG_WPTR(facet_origin_);
        return this->facet_origin_.lock();
      }

      void setFacetOrigin(FacetSPtr facet_origin)
      {
        this->facet_origin_ = facet_origin;
      }

      const FT& getSpeed() const
      {
        CGAL_assertion(speed_ != 0);
        return speed_;
      }

      void setSpeed(const FT& speed)
      {
        speed_ = speed;
      }

    protected:
      FacetWPtr offset_facet_;
      FacetWPtr facet_origin_;
      FT speed_;
    };

private:
    using FT = typename GT::FT;
    using Point_3 = typename GT::Point_3;
    using Plane_3 = typename GT::Plane_3;

    using Point3SPtr = std::shared_ptr<Point_3>;
    using Plane3SPtr = std::shared_ptr<Plane_3>;

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

      result->plane_ = this->plane_;
      result->base_plane_ = this->base_plane_;
      result->final_plane_ = this->final_plane_;

      result->cachedSpeed_ = this->cachedSpeed_;
      result->cachedPlane_ = this->cachedPlane_;

      return result;
    }

    void addVertex(VertexSPtr vertex)
    {
      vertices_.insert(vertices_.end(), vertex);
      vertex->addFacet(this->shared_from_this());
    }

    bool removeVertex(VertexSPtr vertex)
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

    bool hasVertex(VertexSPtr vertex)
    {
      typename std::list<VertexSPtr>::iterator it_v =
        std::find(vertices_.begin(), vertices_.end(), vertex);
      return (it_v != vertices_.end());
    }

    bool isTriangle() const
    {
      return (vertices_.size() == 3);
    }

    void addEdge(EdgeSPtr edge)
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

    bool removeEdge(EdgeSPtr edge)
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
    EdgeSPtr findEdge(FacetSPtr facet) const
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

    std::list<EdgeSPtr> findEdges(FacetSPtr facet) const
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

    bool containsVertex(VertexSPtr vertex) const
    {
      bool result = (vertices_.end() != std::find(vertices_.begin(), vertices_.end(), vertex));
      return result;
    }

    bool containsEdge(EdgeSPtr edge) const
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

    FacetSPtr next(VertexSPtr vertex) const
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

    FacetSPtr prev(VertexSPtr vertex) const
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
    void merge(FacetSPtr facet)
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

    /**
    * Normalize the plane coefficients to obtain a canonical plane representation
    */
    void normalizePlaneCoefficients()
    {
      CGAL_precondition(bool(this->plane_));

      const FT& a = plane_->a();
      const FT& b = plane_->b();
      const FT& c = plane_->c();
      const FT& d = plane_->d();
      // this should be the only place with unavoidable SQRTs
      const FT n = CGAL::approximate_sqrt(CGAL::square(a) + CGAL::square(b) + CGAL::square(c));

      if (!is_zero(n)) {
        // @todo to_double() it here too?
        plane_ = KernelFactory::createPlane3(a/n, b/n, c/n, d/n);
      }
    }

    /**
    * Check if the plane is normalized
    */
    bool isNormalizedPlane()
    {
      CGAL_precondition(bool(this->plane_));
      const FT& a = this->plane_->a();
      const FT& b = this->plane_->b();
      const FT& c = this->plane_->c();
      return (a*a + b*b + c*c - 1) <= 1e-5;
    }

    enum class PerturbationType
    {
      NUDGE,
      STEPS,
      EXACT,
      HIGH_DEGREES
    };

    /**
    * Store the current plane ahead of perturbation
    */
    void storePlaneCoefficients()
    {
      CGAL_precondition(bool(this->plane_));

      // need a different shared ptr here because plane_ will change with the perturbation
      cachedPlane_ = KernelFactory::createPlane3(*(plane_));

      CGAL_SS3_TRANSF_TRACE("caching plane of Facet " << this->id_ << " : [" << *cachedPlane_ << "]");
    }

    /**
    * Nudge the plane coefficients by a random value in the range [low, high].
    */
    void perturbPlaneCoefficientsNudge(const double range)
    {
      CGAL_precondition(isNormalizedPlane());

      // @todo storing coefficients is only useful if we plan on untilting at the end.
      // storePlaneCoefficients();

      CGAL_SS3_TRANSF_TRACE_V(16, "Nudging (Nudge) Face " << this->getID());
      CGAL_SS3_TRANSF_TRACE_V(16, "  From coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                          << plane_->c() << " " << plane_->d() << "]");

      auto nudge = [&](const FT& v) {
        static std::random_device rd;
        unsigned int s = 0; // rd()
        // CGAL_SS3_TRANSF_TRACE("seed = " << s);
        static std::mt19937 gen(s);
        static std::uniform_real_distribution<> rdist(-range, range);

        // Since we are perturbing, we might as well collapse the DAG of 'v'.
        // the point is also that once 'nv' is a double, its interval will be a singleton,
        // and we will have access to static filters
        double step = rdist(gen);
        double nv = CGAL::to_double(v) + step;
        return nv;
      };

      double na = nudge(plane_->a());
      double nb = nudge(plane_->b());
      double nc = nudge(plane_->c());
      double nd = nudge(plane_->d()); // @todo do not nudge 'd'? (mind the 'to_double()')

      double n = CGAL::approximate_sqrt(CGAL::square(na) + CGAL::square(nb) + CGAL::square(nc));
      CGAL_assertion(n != 0); // should not happen since we have normalized and the shift is tiny

      // below doesn't seem to matter? Probably need specific static filters...
#if 0
      plane_ = KernelFactory::createPlane3(na/n, nb/n, nc/n, nd/n);
#else
      // cast to_double() *after* the normalization to have double coordinates in the planes
      // the downside is that we won't have a^2 + b^2 + c^2 == 1,
      // but then again, who does...
      const double a = CGAL::to_double(na/n);
      const double b = CGAL::to_double(nb/n);
      const double c = CGAL::to_double(nc/n);
      const double d = CGAL::to_double(nd/n);
      plane_ = KernelFactory::createPlane3(a, b, c, d);
#endif

      CGAL_SS3_TRANSF_TRACE_V(16, "  To coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                        << plane_->c() << " " << plane_->d() << "]");

      CGAL_postcondition(isNormalizedPlane());
    }

    /**
    * Nudge the plane coefficients by a number of nextafter steps
    * in a random direction.
    */
    // I can go lower
    void perturbPlaneCoefficientsSteps(int steps)
    {
      CGAL_precondition(isNormalizedPlane());

      // @todo storing coefficients is only useful if we plan on untilting at the end.
      // storePlaneCoefficients();

      CGAL_SS3_TRANSF_TRACE_V(16, "Nudging (Steps) Face " << this->getID());
      CGAL_SS3_TRANSF_TRACE_V(16, "  From coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                          << plane_->c() << " " << plane_->d() << "]");

      auto nudge = [&](const FT& v) {
        std::random_device rd;
        unsigned int s = 0; // rd()
        // CGAL_SS3_TRANSF_TRACE("seed = " << s);
        std::mt19937 gen(s);
        std::uniform_int_distribution<> step_dis(1, steps);
        std::uniform_int_distribution<> dir_dis(0, 1); // 0: negative, 1: positive
        int n = step_dis(gen);
        int dir = dir_dis(gen);
        double direction = dir ? INFINITY : -INFINITY;
        double nv = CGAL::to_double(v);
        for (int i = 0; i < n; ++i) {
          nv = std::nextafter(nv, direction);
        }
        return nv;
      };

      double na = nudge(plane_->a());
      double nb = nudge(plane_->b());
      double nc = nudge(plane_->c());
      double nd = nudge(plane_->d()); // @todo do not nudge 'd'? (mind the 'to_double()')

      double n = CGAL::approximate_sqrt(CGAL::square(na) + CGAL::square(nb) + CGAL::square(nc));
      CGAL_assertion(n != 0); // should not happen since we have normalized and the shift is tiny

      // below doesn't seem to matter? Probably need specific static filters...
  #if 0
      plane_ = KernelFactory::createPlane3(na/n, nb/n, nc/n, nd/n);
  #else
      // cast to_double() *after* the normalization to have double coordinates in the planes
      // the downside is that we won't have a^2 + b^2 + c^2 == 1,
      // but then again, who does...
      plane_ = KernelFactory::createPlane3(CGAL::to_double(na/n),
                                           CGAL::to_double(nb/n),
                                           CGAL::to_double(nc/n),
                                           CGAL::to_double(nd/n));
  #endif

      CGAL_SS3_TRANSF_TRACE_V(16, "  To coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                        << plane_->c() << " " << plane_->d() << "]");

      CGAL_postcondition(isNormalizedPlane());
    }

    /**
    * Nudge the plane coefficients by a random value in the range [0, 1] / den.
    */
    // I can go lower
    void perturbPlaneCoefficientsExact(const FT& den)
    {
      CGAL_precondition(isNormalizedPlane());

      // @todo storing coefficients is only useful if we plan on untilting at the end.
      // storePlaneCoefficients();

      CGAL_SS3_TRANSF_TRACE_V(16, "Nudging (Exact) Face " << this->getID());
      CGAL_SS3_TRANSF_TRACE_V(16, "  From coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                          << plane_->c() << " " << plane_->d() << "]");

      auto nudge = [&](const FT& v) {
        static std::random_device rd;
        unsigned int s = 0; // rd()
        // CGAL_SS3_TRANSF_TRACE("seed = " << s);
        static std::mt19937 gen(s);
        static std::uniform_real_distribution<> rdist(0, 1);
        FT step = rdist(gen);
        return v + step / CGAL::square(den);
      };

      FT na = nudge(plane_->a());
      FT nb = nudge(plane_->b());
      FT nc = nudge(plane_->c());
      FT nd = nudge(plane_->d()); // @todo do not nudge 'd'?

      // so small, needless to normalize (@todo should we even normalize others?)
      plane_ = KernelFactory::createPlane3(na, nb, nc, nd);

      CGAL_SS3_TRANSF_TRACE_V(16, "  To coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                        << plane_->c() << " " << plane_->d() << "]");

      CGAL_postcondition(isNormalizedPlane());
    }

    /**
    * Nudge the plane coefficients but ensure that the perturbed plane goes through 0, 1, or 2 fixed points.
    * If 0 points: nudge all coefficients independently.
    * If 1 point: nudge (a, b, c), recompute d so the plane passes through the point.
    * If 2 points: nudge (a, b, c) with the constraint that the new plane passes through both points.
    */
    void perturbPlaneCoefficientsFixedPoints(const double range,
                                            const std::vector<Point3SPtr>& fixed_points)
    {
      CGAL_precondition(isNormalizedPlane());
      CGAL_precondition(fixed_points.size() <= 2);

      // @todo storing coefficients is only useful if we plan on untilting at the end.
      // storePlaneCoefficients();

      CGAL_SS3_TRANSF_TRACE_V(16, "Nudging Face " << this->getID());
      CGAL_SS3_TRANSF_TRACE_V(16, "  From coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                          << plane_->c() << " " << plane_->d() << "]");
      CGAL_SS3_TRANSF_TRACE_V(16, "  with " << fixed_points.size() << " fixed points");
      CGAL_SS3_TRANSF_TRACE_CODE(for (Point3SPtr fp : fixed_points))
      CGAL_SS3_TRANSF_TRACE_V(16, "    " << *fp);

      static std::random_device rd;
      unsigned int s = 0; // rd()
      // CGAL_SS3_TRANSF_TRACE("seed = " << s);
      static std::mt19937 gen(s);
      static std::uniform_real_distribution<> rdist(-range, range);

      auto nudge = [&](const FT& v) {
        // Since we are perturbing, we might as well collapse the DAG of 'v'.
        // the point is also that once 'nv' is a double, its interval will be a singleton,
        // and we will have access to static filters
        double step = rdist(gen);
        double nv = CGAL::to_double(v) + step;
        return nv;
      };

  #ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
      auto nudge_to_simplest_rational_in_interval = [&](const FT& v) {
        double d1 = nudge(v);
        double d2 = nudge(v);
        if (d2 < d1) {
          std::swap(d1, d2);
        }
        FT nv = CGAL::simplest_rational_in_interval<CGAL::K::Exact_kernel::FT>(d1, d2);
        return nv;
      };
  #endif

      // 0 fixed points: nudge all coefficients independently
      if (fixed_points.size() == 0) {
  #ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        FT na = nudge_to_simplest_rational_in_interval(plane_->a());
        FT nb = nudge_to_simplest_rational_in_interval(plane_->b());
        FT nc = nudge_to_simplest_rational_in_interval(plane_->c());
        FT nd = nudge_to_simplest_rational_in_interval(plane_->d());
  #else
        double na = nudge(plane_->a());
        double nb = nudge(plane_->b());
        double nc = nudge(plane_->c());
        double nd = nudge(plane_->d());
  #endif
        plane_ = KernelFactory::createPlane3(na, nb, nc, nd);
      } else if (fixed_points.size() == 1) {
        // 1 fixed point: nudge (a, b, c), recompute d so the plane passes through the point
        Point3SPtr p0 = fixed_points[0];

  #ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        FT na = nudge_to_simplest_rational_in_interval(plane_->a());
        FT nb = nudge_to_simplest_rational_in_interval(plane_->b());
        FT nc = nudge_to_simplest_rational_in_interval(plane_->c());
  #else
        double na = nudge(plane_->a());
        double nb = nudge(plane_->b());
        double nc = nudge(plane_->c());
  #endif
        const FT& x0 = p0->x();
        const FT& y0 = p0->y();
        const FT& z0 = p0->z();
        FT d = - (na * x0 + nb * y0 + nc * z0);
        plane_ = KernelFactory::createPlane3(na, nb, nc, d);
        CGAL_postcondition(plane_->has_on(*p0));
      } else if (fixed_points.size() == 2) {
        // 2 fixed points: construct a plane through both points, nudge the normal within the allowed family
        Point3SPtr p0 = fixed_points[0];
        Point3SPtr p1 = fixed_points[1];
        CGAL_assertion(*p0 != *p1);

        const FT& p0x = p0->x();
        const FT& p0y = p0->y();
        const FT& p0z = p0->z();
        const FT& p1x = p1->x();
        const FT& p1y = p1->y();
        const FT& p1z = p1->z();

        // Step 1: Direction vector between points
        FT ux = p1x - p0x;
        FT uy = p1y - p0y;
        FT uz = p1z - p0z;
        FT uu = ux*ux + uy*uy + uz*uz;

        // Step 2: Original normal
        const FT& a0 = plane_->a();
        const FT& b0 = plane_->b();
        const FT& c0 = plane_->c();

        // Step 3: Project original normal onto plane orthogonal to u
        FT dot = a0*ux + b0*uy + c0*uz;
        FT ab = a0 - dot * ux / uu;
        FT bb = b0 - dot * uy / uu;
        FT cb = c0 - dot * uz / uu;

        // Step 4: Find a direction to nudge (cross product)
        FT vx = uy * cb - uz * bb;
        FT vy = uz * ab - ux * cb;
        FT vz = ux * bb - uy * ab;

        // Step 5: Nudge the normal
  #ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        FT epsilon = nudge_to_simplest_rational_in_interval(rdist(gen));
  #else
        double epsilon = rdist(gen);
  #endif

        FT a1 = ab + epsilon * vx;
        FT b1 = bb + epsilon * vy;
        FT c1 = cb + epsilon * vz;

        // Step 6: Compute d so plane passes through p0
        FT d1 = - (a1 * p0x + b1 * p0y + c1 * p0z);
        plane_ = KernelFactory::createPlane3(a1, b1, c1, d1);

        CGAL_postcondition(plane_->has_on(*p0));
        CGAL_postcondition(plane_->has_on(*p1));
      } else {
        CGAL_SS3_TRANSF_TRACE("Error: called fixed point facet perturbation with > 2 fixed points");
      }

      CGAL_SS3_TRANSF_TRACE_V(16, "  To coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                        << plane_->c() << " " << plane_->d() << "]");

      CGAL_postcondition(isNormalizedPlane());
    }

    void perturbPlaneCoefficientsHighDegrees(const double range)
    {
      std::vector<Point3SPtr> high_degree_points;
      for (VertexSPtr v : vertices_) {
        if (v->degree() > 3) {
          high_degree_points.push_back(v->getPoint());
        }
      }
      return perturbPlaneCoefficientsFixedPoints(range, high_degree_points);
    }

    /**
    * Nudge the plane coefficients to get rid of simultaneous events
    */
    void perturbPlaneCoefficients(PerturbationType type = PerturbationType::HIGH_DEGREES)
    {
      double range = 1e-10;
      ConfigurationSPtr config = Configuration::getInstance();
      if (config->isLoaded()) {
        range = config->getDouble("main", "rand_move_points_range");
      }

      if (type == PerturbationType::NUDGE) {
        perturbPlaneCoefficientsNudge(range);
      } else if (type == PerturbationType::STEPS) {
        perturbPlaneCoefficientsSteps(10);
      } else if (type == PerturbationType::EXACT) {
        perturbPlaneCoefficientsExact(FT(1e100));
      } else if (type == PerturbationType::HIGH_DEGREES) {
        perturbPlaneCoefficientsHighDegrees(range);
      } else {
        CGAL_error_msg("Unknown perturbation type");
      }
    }

    /**
    * Restore the plane coefficients to the previous value, updating 'd' so that the plane
    * matches the desired offset.
    */
    void restorePlaneCoefficients(const FT& perturbationOffset,
                                  const FT& perturbationEndOffset)
    {
      CGAL_SS3_TRANSF_TRACE("plane of Facet " << this->id_ << " is [" << *plane_ << "]");

      if (!cachedPlane_) {
        std::cerr << "Warning: no plane coefficients to restore" << std::endl;
        return;
      }

      FT speed = 1.0;
      if (hasData()) {
        speed = std::dynamic_pointer_cast<SkelFacetData>(getData())->getSpeed();
      }

      CGAL_SS3_TRANSF_TRACE("OLD d = " << cachedPlane_->d());
      CGAL_SS3_TRANSF_TRACE("perturbationOffset = " << perturbationOffset);
      CGAL_SS3_TRANSF_TRACE("perturbationEndOffset = " << perturbationEndOffset);

      // The minus sign "d - ..." is because we shrink, so the plane needs to be offset
      // by the difference of offsets, but in the direction opposite of its normal.
      //
      // This is similar to when we call, e.g.:
      //   Plane3SPtr offset_plane_l = KernelWrapper::offsetPlane(plane_l, - speed_l);
      //                                                                  ^^^
      FT d = cachedPlane_->d() - speed * (perturbationEndOffset - perturbationOffset);

      plane_ = KernelFactory::createPlane3(cachedPlane_->a(), cachedPlane_->b(), cachedPlane_->c(), d);
      CGAL_assertion_code(FT sq_n = CGAL::square(plane_->a()) + CGAL::square(plane_->b()) + CGAL::square(plane_->c()));
      CGAL_assertion((sq_n - 1) < 1e-5);

      CGAL_SS3_TRANSF_TRACE("plane of Facet " << this->id_ << " restored to [" << *plane_ << "]");
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

#if 1 // @fixme is this correct? the CGAL_PI/4.0 below is confusing... Was it supposed to be CGAL_PI/2.0?
        if (!CGAL::collinear(*(points[0]), *(points[1]), *(points[2]))) {
          if (CGAL::angle(*(points[0]), *(points[1]), *(points[2]), *normal) == CGAL::ACUTE) {
            edge_begin = edge;
            result = true;
            break;
          }
        }
#else // old code; has issues with collinear points
        if (points[0] != points[1] &&
            points[1] != points[2] &&
            points[2] != points[0]) {
          Plane3SPtr plane_current = KernelFactory::createPlane3(
                  points[0], points[1], points[2]);
          Vector3SPtr normal_current = KernelFactory::createVector3(plane_current);
          double angle = 0.0;
          double arg = 0.0;
          arg = CGAL::to_double(((*normal)*(*normal_current)) /
                  CGAL::disallowed_sqrt(normal->squared_length() * normal_current->squared_length()));

          // fixes issues with floating point precision
          if (arg <= -1.0) {
            angle = CGAL_PI;
          } else if (arg >= 1.0) {
            angle = 0.0;
          } else {
            angle = acos(arg);
          }
          if (angle < CGAL_PI/4.0) {
            edge_begin = edge;
            result = true;
            break;
          }
        }
#endif
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
      for (VertexSPtr v : vertices_) {
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

  public: // @tmp should be 'protected', but tmp public for shiftfacets() and handleEdgeEvent (which should shift in place degree 1 vertices)
    std::list<VertexSPtr> vertices_;
    std::list<EdgeSPtr> edges_;
    PolyhedronWPtr polyhedron_;
    typename std::list<FacetSPtr>::iterator polyhedron_list_it_;
    FacetDataSPtr data_;

    Plane3SPtr plane_;
    Plane3SPtr base_plane_;
    Plane3SPtr final_plane_;

    int id_;

    Plane3SPtr cachedPlane_;
    FT cachedSpeed_;
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
      VertexSPtr vertex_c = vertex->clone();
      result->addVertex(vertex_c);
      vertices_c[vertex] = vertex_c;
    }
    typename std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
      EdgeSPtr edge = *it_e++;
      VertexSPtr src = vertices_c[edge->getVertexSrc()];
      VertexSPtr dst = vertices_c[edge->getVertexDst()];
      EdgeSPtr edge_c = Edge<Traits>::create(src, dst);
      result->addEdge(edge_c);
      edges_c[edge] = edge_c;
    }
    typename std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
      FacetSPtr facet = *it_f++;
      FacetSPtr facet_c = Facet<Traits>::create();
      facet_c->plane_ = facet->plane_;
      facet_c->base_plane_ = facet->base_plane_;
      facet_c->final_plane_ = facet->final_plane_;

      if (facet->hasData()) {
        FacetDataSPtr data = facet->getData();
        SkelFacetDataSPtr data_c = SkelFacetData::create(facet_c);
        data_c->setSpeed(std::dynamic_pointer_cast<SkelFacetData>(data)->getSpeed());
      }
      if (facet->getID() != -1) { // @todo remove this? Anyway the other IDs are not copied...
        facet_c->setID(facet->getID());
      }
      facet_c->cachedPlane_ = facet->cachedPlane_;
      facet_c->cachedSpeed_ = facet->cachedSpeed_;

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
      result->addFacet(facet_c);
    }
    return result;
  }

  void addVertex(VertexSPtr vertex)
  {
    vertex->setID(next_vertex_id_++);
    typename std::list<VertexSPtr>::iterator it = vertices_.insert(vertices_.end(), vertex);
    vertex->setPolyhedron(this->shared_from_this());
    vertex->setPolyhedronListIt(it);
  }

  bool removeVertex(VertexSPtr vertex)
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
  VertexSPtr findVertex(VertexSPtr needle)
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

  void addEdge(EdgeSPtr edge)
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

  // @todo should this remove incident degree 1 vertices?
  bool removeEdge(EdgeSPtr edge)
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
  EdgeSPtr findEdge(EdgeSPtr needle)
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

  void addFacet(FacetSPtr facet)
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

  bool removeFacet(FacetSPtr facet)
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
    for (FacetSPtr facet : facets_) {
      facet->initPlane();
    }
  }

  void clearData()
  {
    for (VertexSPtr vertex : vertices_) {
      vertex->setData(VertexDataSPtr());
    }
    for (EdgeSPtr edge : edges_) {
      edge->setData(EdgeDataSPtr());
    }
    for (FacetSPtr facet : facets_) {
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
        std::cerr << "no left facet?" << std::endl;
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
        std::cerr << "no right facet?" << std::endl;
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

    for (VertexSPtr vertex : vertices_) {
      vertex->setID(next_vertex_id_++);
    }
    for (EdgeSPtr edge : edges_) {
      edge->setID(next_edge_id_++);
    }
    for (FacetSPtr facet : facets_) {
      facet->setID(next_facet_id_++);
    }
    setID(-1);
  }

  void resetAllIDs()
  {
    for (VertexSPtr vertex : vertices_) {
      vertex->setID(-1);
    }
    for (EdgeSPtr edge : edges_) {
      edge->setID(-1);
    }
    for (FacetSPtr facet : facets_) {
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
