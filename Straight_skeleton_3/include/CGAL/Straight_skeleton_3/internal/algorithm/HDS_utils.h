// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_HDS_UTILS_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_HDS_UTILS_H

#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>

#include <CGAL/assertions.h>

#include <memory>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class PolyhedronTransformation;

template <typename Traits>
class HdsUtils
{
  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;
  using Vector_3 = typename Traits::Vector_3;
  using Line_3 = typename Traits::Line_3;
  using Plane_3 = typename Traits::Plane_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Vector3SPtr = std::shared_ptr<Vector_3>;
  using Line3SPtr = std::shared_ptr<Line_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Vertex = typename Polyhedron::template Vertex<Traits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using SkelVertexData = typename Polyhedron::SkelVertexData;
  using SkelVertexDataSPtr = typename Polyhedron::SkelVertexDataSPtr;
  using Edge = typename Polyhedron::template Edge<Traits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using SkelEdgeData = typename Polyhedron::SkelEdgeData;
  using SkelEdgeDataSPtr = typename Polyhedron::SkelEdgeDataSPtr;
  using Facet = typename Polyhedron::template Facet<Traits>;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;
  using SkelFacetData = typename Polyhedron::SkelFacetData;
  using SkelFacetDataSPtr = typename Polyhedron::SkelFacetDataSPtr;

private:
  using KernelWrapper = kernel::KernelWrapper<Traits>;
  using KernelFactory = kernel::KernelFactory<Traits>;
  using Transformation = algorithm::PolyhedronTransformation<Traits>;

public:
  static Line3SPtr line(const EdgeSPtr& edge)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    Line3SPtr result = Line3SPtr();
    VertexSPtr vertex_src = edge->getVertexSrc();
    VertexSPtr vertex_dst = edge->getVertexDst();
    if (*(vertex_src->getPoint()) == *(vertex_dst->getPoint())) {
      FacetSPtr facet_l = edge->getFacetL();
      FacetSPtr facet_r = edge->getFacetR();
      FacetSPtr facet_src = edge->getFacetSrc();
      FacetSPtr facet_dst = edge->getFacetDst();

      Plane3SPtr offset_plane_l = Transformation::shiftPlane(facet_l, -1);
      Plane3SPtr offset_plane_r = Transformation::shiftPlane(facet_r, -1);
      Plane3SPtr offset_plane_src = Transformation::shiftPlane(facet_src, -1);
      Plane3SPtr offset_plane_dst = Transformation::shiftPlane(facet_dst, -1);

      Point3SPtr p_src = KernelWrapper::intersection(offset_plane_src,
              offset_plane_l, offset_plane_r);
      Point3SPtr p_dst = KernelWrapper::intersection(offset_plane_dst,
              offset_plane_l, offset_plane_r);
      Vector3SPtr v_dir = KernelFactory::createVector3((*p_dst) - (*p_src));
      result = KernelFactory::createLine3(vertex_src->getPoint(), v_dir);
    } else {
      result = edge->line();
    }
    return result;
  }

public:
  // just for convenience to avoid all the verbosity of extracting the data
  static const FT& getSpeed(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->hasData());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->getSpeed();
  }

  static void setSpeed(const FacetSPtr& facet, const FT& speed)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->hasData());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
    CGAL_SS3_DEBUG_SPTR(data);
    data->setSpeed(speed);
  }

  static Plane3SPtr getBasePlane(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->hasData());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->getBasePlane();
  }

  static void setBasePlane(const FacetSPtr& facet, const Plane3SPtr& plane)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_SS3_DEBUG_SPTR(plane);
    CGAL_precondition(facet->hasData());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
    CGAL_SS3_DEBUG_SPTR(data);
    data->setBasePlane(plane);
    CGAL_SS3_DEBUG_SPTR(getBasePlane(facet));
  }

  static bool hasFinalPlane(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->hasData());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->hasFinalPlane();
  }

  static Plane3SPtr getFinalPlane(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->hasData());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->getFinalPlane();
  }

  static void setFinalPlane(const FacetSPtr& facet, const Plane3SPtr& plane)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_SS3_DEBUG_SPTR(plane);
    CGAL_precondition(facet->hasData());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->setFinalPlane(plane);
  }

    /**
    * updates the final coefficients of the plane of a facet
    */
  // resets using the first 3 final planes, even if there are more
  static void resetFinalPlane(const FacetSPtr& facet,
                              const FT& offset_future_bound)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    Plane3SPtr plane = Transformation::shiftPlane(facet, offset_future_bound);

    CGAL_precondition(facet->hasData());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
    CGAL_SS3_DEBUG_SPTR(data);
    data->setFinalPlane(plane);
  }

  static Plane3SPtr getFinalPlane(const FacetSPtr& facet,
                                  const FT& offset_future_bound)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->hasData());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
    CGAL_SS3_DEBUG_SPTR(data);

    // @todo technically we should check if it's the correct plane
    if (!data->hasFinalPlane()) {
      resetFinalPlane(facet, offset_future_bound);
    }
    return data->getFinalPlane();
  }

public:
  static bool hasFinalPoint(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->hasData());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->hasFinalPoint();
  }

  static Point3SPtr getFinalPoint(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->hasData());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->getFinalPoint();
  }

  static void setFinalPoint(const VertexSPtr& vertex, const Point3SPtr& point)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_SS3_DEBUG_SPTR(point);
    CGAL_precondition(vertex->hasData());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->setFinalPoint(point);
  }


  /**
    * updates the final position of the vertex of a polyhedron,
    * computed from the planes of its incident faces.
    */
  static bool resetFinalPoint(const VertexSPtr& vertex,
                              const FT& offset_future_bound)
  {
    CGAL_SS3_TRANSF_TRACE_V(16, "resetFinalPoint() of " << vertex->toString());
    CGAL_SS3_DEBUG_SPTR(vertex);

    std::array<Plane3SPtr, 3> final_planes;
    unsigned int i = 0;
    typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (i < 3 && it_f != vertex->facets().end()) {
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        if (!hasFinalPlane(facet)) {
          resetFinalPlane(facet, offset_future_bound);
        }

        final_planes[i++] = getFinalPlane(facet);

        CGAL_SS3_TRANSF_TRACE_V(32,"  Facet " << facet->getID() << " [" << *(getFinalPlane(facet)) << "]");
      }
    }
    CGAL_postcondition(i == 3);

    Point3SPtr point = KernelWrapper::intersection(final_planes[0], final_planes[1], final_planes[2]);
    if (!point) {
      CGAL_SS3_TRANSF_TRACE_V(1,"Warning: triplet of planes does not define a point!");
      Point3SPtr result = Point3SPtr();
      CGAL_SS3_DEBUG_SPTR(result);
      return false;
    }

    CGAL_precondition(vertex->hasData());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
    CGAL_SS3_DEBUG_SPTR(data);

    data->setFinalPoint(point);
    CGAL_SS3_TRANSF_TRACE_V(32,"  New final point = " << *point);

    return true;
  }

  static Point3SPtr getFinalPoint(const VertexSPtr& vertex,
                                  const FT& offset_future_bound)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->hasData());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
    CGAL_SS3_DEBUG_SPTR(data);

    // @todo technically we should check if it's the correct point
    if (!data->hasFinalPoint()) {
      resetFinalPoint(vertex, offset_future_bound);
    }
    return data->getFinalPoint();
  }

public:
  static bool isReflex(const EdgeSPtr& edge,
                       const bool future_facing = true)
  {
    CGAL_SS3_DEBUG_SPTR(edge);

    if (edge->getReflexStatus()) {
      return *(edge->getReflexStatus());
    }

    bool result = false;

    VertexSPtr vertex_src = edge->getVertexSrc();
    VertexSPtr vertex_dst = edge->getVertexDst();

    // @todo is degenerate edge? pointer comparison or actual position comparison?
    if (*(vertex_src->getPoint()) == *(vertex_dst->getPoint())) {
      FacetSPtr facet_l = edge->getFacetL();
      FacetSPtr facet_r = edge->getFacetR();
      FacetSPtr facet_src = edge->getFacetSrc();
      FacetSPtr facet_dst = edge->getFacetDst();

      FT od = future_facing ? -1 : 1;
      Plane3SPtr offset_plane_l = Transformation::shiftPlane(facet_l, od);
      Plane3SPtr offset_plane_r = Transformation::shiftPlane(facet_r, od);
      Plane3SPtr offset_plane_src = Transformation::shiftPlane(facet_src, od);
      Plane3SPtr offset_plane_dst = Transformation::shiftPlane(facet_dst, od);

      Point3SPtr p_src = KernelWrapper::intersection(offset_plane_src, offset_plane_l, offset_plane_r);
      Point3SPtr p_dst = KernelWrapper::intersection(offset_plane_dst, offset_plane_l, offset_plane_r);
      CGAL_assertion(p_src && p_dst);

      Vector3SPtr v_dir = KernelFactory::createVector3((*p_dst) - (*p_src));
      CGAL_assertion(*v_dir != CGAL::NULL_VECTOR);

      Vector3SPtr normal_l = KernelFactory::createVector3(offset_plane_l);
      CGAL_assertion(*normal_l != CGAL::NULL_VECTOR);

      Vector3SPtr v_cross = KernelWrapper::cross(normal_l, v_dir);
      Point3SPtr p = KernelFactory::createPoint3((*p_src) + (*v_cross));
      if (KernelWrapper::side(offset_plane_r, p) > 0) {
        result = true;
      }
    } else {
      result = edge->isReflex();
    }
    return result;
  }

  static bool isReflex(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    if (vertex->degree() == 0) {
      return false;
    }
    bool result = true;
    for (EdgeWPtr edge_wptr : vertex->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if (!isReflex(edge)) {
          result = false;
          break;
        }
      }
    }
    return result;
  }

  static bool isConvex(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    if (vertex->degree() == 0) {
      return false;
    }
    bool result = true;
    for (EdgeWPtr edge_wptr : vertex->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if (isReflex(edge)) {
          result = false;
          break;
        }
      }
    }
    return result;
  }

  /**
    * Checks if the given edge is part of a triangle
    * on the surface of the polyhedron.
    * Used by: nextEdgeEvent, nextTriangleEvent
    */
  static bool isTriangle(const FacetSPtr& facet,
                         EdgeSPtr edge_begin)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_SS3_DEBUG_SPTR(edge_begin);
    bool result = false;
    if (!facet->containsEdge(edge_begin)) {
      return false;
    }
    unsigned int i = 0;
    // VertexSPtr vertices[3];
    EdgeSPtr edge = edge_begin;
    for (i = 0; i < 3; ++i) {
      // vertices[i] = edge->src(facet);
      edge = edge->next(facet);
      if (!edge) {
        break;
      }
    }
    if (i == 3 && edge_begin == edge) {
      result = true;
      // careful about ptr equality if below is uncommented
      // if (vertices[0]->getPoint() != vertices[1]->getPoint() &&
      //         vertices[1]->getPoint() != vertices[2]->getPoint() &&
      //         vertices[2]->getPoint() != vertices[0]->getPoint()) {
      //     // check if triangle is a hole
      //     Plane3SPtr plane = KernelFactory::createPlane3(
      //             vertices[0]->getPoint(),
      //             vertices[1]->getPoint(),
      //             vertices[2]->getPoint());
      //     Vector3SPtr v1 = KernelFactory::createVector3(plane);
      //     Vector3SPtr v2 = KernelFactory::createVector3(facet->getPlane());
      //     if (((*v1) * (*v2)) > 0.0) {
      //         // angle between orthogonal vectors of planes < CGAL_PI/2.0
      //         // not a hole
      //         result = true;
      //     }
      // }
    }
    return result;
  }

  /**
    * Checks if the given edge is part of a tetrahedron.
    * Used by: nextTriangleEvent, nextTetrahedronEvent
    */
  static bool isTetrahedron(const EdgeSPtr& edge_begin)
  {
    CGAL_SS3_DEBUG_SPTR(edge_begin);
    bool result = true;
    VertexSPtr vertices[4];
    for (unsigned int i = 0; i < 4; ++i) {
      vertices[i] = VertexSPtr();
    }
    vertices[0] = edge_begin->getVertexSrc();
    vertices[1] = edge_begin->getVertexDst();
    EdgeSPtr edge_l = edge_begin->next(edge_begin->getFacetL());
    vertices[2] = edge_l->dst(edge_begin->getFacetL());
    EdgeSPtr edge_r = edge_begin->next(edge_begin->getFacetR());
    vertices[3] = edge_r->dst(edge_begin->getFacetR());
    for (unsigned int i = 0; i < 4; ++i) {
      if (!vertices[i]) {
        result = false;
        return result;
      }
    }
    for (unsigned int i = 0; i < 4; ++i) {
      for (unsigned int j = 1; j < 4; ++j) {
        EdgeSPtr edge = vertices[i]->findEdge(vertices[(i+j)%4]);
        if (!edge) {
          result = false;
        }
      }
    }
    return result;
  }

};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_HDS_UTILS_H */
