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
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_self_intersection.h>

#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class AbstractVertexSplitter
{
  using Polyhedron = HDS::Polyhedron<Traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using VertexSPtr = typename Polyhedron::VertexSPtr;

private:
  using Transformation = algorithm::PolyhedronTransformation<Traits>;
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

  virtual PolyhedronSPtr splitVertex(const VertexSPtr& vertex) = 0;  // abstract

  virtual int getType() const
  {
    return type_;
  }

  static bool checkSplitted(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    bool result = false;
    PolyhedronSPtr polyhedron_cpy = polyhedron->clone();
    result = Transformation::shiftFacetsDegree1(polyhedron_cpy, -1.0);
    result = result && !SelfIntersection::hasSelfIntersectingSurface(polyhedron_cpy);
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
