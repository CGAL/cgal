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
 * @file   algo/3d/AbstractVertexSplitter.h
 * @author Gernot Walzl
 * @date   2012-10-17
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_ABSTRACT_VERTEX_SPLITTER_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_ABSTRACT_VERTEX_SPLITTER_H

#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_self_intersection.h>

#include <list>
#include <memory>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class AbstractVertexSplitter
{
  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;
  using Plane_3 = typename Traits::Plane_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  using SkelVertexData = typename Polyhedron::SkelVertexData;
  using SkelVertexDataSPtr = typename Polyhedron::SkelVertexDataSPtr;
  using SkelFacetData = typename Polyhedron::SkelFacetData;

private:
  using KernelWrapper = kernel::KernelWrapper<Traits>;
  using PolyhedronTransformation = algorithm::PolyhedronTransformation<Traits>;
  using SelfIntersection = algorithm::SelfIntersection<Traits>;

public:
  // static const int ANGLE_VERTEX_SPLITTER = -1;  // does not work
  static const int COMBI_VERTEX_SPLITTER = 1;
  static const int CONVEX_VERTEX_SPLITTER = 2;
  // static const int VOLUME_VERTEX_SPLITTER = 3;
  // static const int WEIGHT_VERTEX_SPLITTER = 4;
  // static const int SPHERE_VERTEX_SPLITTER = 5;

public:
  AbstractVertexSplitter(int type = 0)
    : type_(type)
  { }

  virtual ~AbstractVertexSplitter() { /*intentionally does nothing*/ }

  virtual PolyhedronSPtr splitVertex(VertexSPtr vertex) = 0;  // abstract

  virtual int getType() const
  {
    return type_;
  }

  static PolyhedronSPtr splitConvexVertex(VertexSPtr vertex)
  {
    assert(vertex->isConvex());
    PolyhedronSPtr result = vertex->getPolyhedron();
    while (vertex->degree() > 3) {
      FT sq_speed_max = 0.0;
      FacetSPtr facet_1;
      FacetSPtr facet_2;
      typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
      while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          FacetSPtr facet_prev = facet->prev(vertex);
          FacetSPtr facet_next = facet->next(vertex);
          FT speed = 1.0;
          if (facet->hasData()) {
            speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
          }
          FT speed_prev = 1.0;
          if (facet_prev->hasData()) {
            speed_prev = std::dynamic_pointer_cast<SkelFacetData>(facet_prev->getData())->getSpeed();
          }
          FT speed_next = 1.0;
          if (facet_next->hasData()) {
            speed_next = std::dynamic_pointer_cast<SkelFacetData>(facet_next->getData())->getSpeed();
          }

          // @todo I don't think the quick split of this function works for weighted cases
          CGAL_assertion(speed == 1.0 && speed_prev == 1.0 && speed_next == 1.0);

          Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->plane(), - speed);
          Plane3SPtr offset_plane_prev = KernelWrapper::offsetPlane(facet_prev->plane(), -speed_prev);
          Plane3SPtr offset_plane_next = KernelWrapper::offsetPlane(facet_next->plane(), -speed_next);

          Point3SPtr offset_point = KernelWrapper::intersection(offset_plane_prev, offset_plane, offset_plane_next);

          FT sq_speed = KernelWrapper::squared_distance(vertex->getPoint(), offset_point);
          if (sq_speed > sq_speed_max) {
            facet_1 = facet_next;
            facet_2 = facet_prev;
            sq_speed_max = sq_speed;
          }
        }
      }
      VertexSPtr vertex_splitted = vertex->split(facet_1, facet_2);
      if (vertex->hasData()) {
        SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
        SkelVertexDataSPtr data_splitted = SkelVertexData::create(vertex_splitted);
        data_splitted->setNode(data->getNode());
      }
    }
    return result;
  }

  static PolyhedronSPtr splitReflexVertex(VertexSPtr vertex)
  {
    CGAL_precondition(vertex->isReflex());

    PolyhedronSPtr result = vertex->getPolyhedron();
    while (vertex->degree() > 3) {
      FT speed_min = std::numeric_limits<double>::max(); // do not put FT
      FacetSPtr facet_1;
      FacetSPtr facet_2;
      typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
      while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          FacetSPtr facet_prev = facet->prev(vertex);
          FacetSPtr facet_next = facet->next(vertex);
          FT speed = 1.0;
          if (facet->hasData()) {
            speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
          }
          FT speed_prev = 1.0;
          if (facet_prev->hasData()) {
            speed_prev = std::dynamic_pointer_cast<SkelFacetData>(facet_prev->getData())->getSpeed();
          }
          FT speed_next = 1.0;
          if (facet_next->hasData()) {
            speed_next = std::dynamic_pointer_cast<SkelFacetData>(facet_next->getData())->getSpeed();
          }

          // @todo I don't think the quick split of this function works for weighted cases
          CGAL_assertion(speed == 1.0 && speed_prev == 1.0 && speed_next == 1.0);

          Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->plane(), -speed);
          Plane3SPtr offset_plane_prev = KernelWrapper::offsetPlane(facet_prev->plane(), -speed_prev);
          Plane3SPtr offset_plane_next = KernelWrapper::offsetPlane(facet_next->plane(), -speed_next);

          Point3SPtr offset_point = KernelWrapper::intersection(
                  offset_plane_prev, offset_plane, offset_plane_next);

          FT v_speed = KernelWrapper::distance(vertex->getPoint(), offset_point);
          if (v_speed < speed_min) {
            facet_1 = facet_next;
            facet_2 = facet_prev;
            speed_min = v_speed;
          }
        }
      }
      VertexSPtr vertex_splitted = vertex->split(facet_1, facet_2);
      if (vertex->hasData()) {
        SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
        SkelVertexDataSPtr data_splitted = SkelVertexData::create(vertex_splitted);
        data_splitted->setNode(data->getNode());
      }
    }
    return result;
  }

  static bool checkSplitted(PolyhedronSPtr polyhedron)
  {
    bool result = false;
    PolyhedronSPtr polyhedron_offset = PolyhedronTransformation::shiftFacets(polyhedron, -1.0);
    if (polyhedron_offset) {
      result = !SelfIntersection::hasSelfIntersectingSurface(polyhedron_offset);
    }
    return result;
  }

  virtual std::string toString() const
  {
    std::string result;
    switch (getType()) {
      case COMBI_VERTEX_SPLITTER:
        result = "CombiVertexSplitter";
        break;
      case CONVEX_VERTEX_SPLITTER:
        result = "ConvexVertexSplitter";
        break;
      default:
        result = "AbstractVertexSplitter";
    }
    return result;
  }

protected:
  int type_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_ABSTRACT_VERTEX_SPLITTER_H */
