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
 * file   algo/3d/AbstractVertexSplitter.h
 * author Gernot Walzl
 * date   2012-10-17
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

template <typename GeomTraits>
class Abstract_vertex_splitter
{
  using Polyhedron = HDS::Polyhedron<GeomTraits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using VertexSPtr = typename Polyhedron::VertexSPtr;

private:
  using Transformation = algorithm::Polyhedron_transformation<GeomTraits>;
  using Self_intersection = algorithm::Self_intersection<GeomTraits>;

public:
  // static const int ANGLE_VERTEX_SPLITTER = -1;  // does not work
  static const int COMBI_VERTEX_SPLITTER = 1;
  static const int CONVEX_VERTEX_SPLITTER = 2;
  // static const int VOLUME_VERTEX_SPLITTER = 3;
  // static const int WEIGHT_VERTEX_SPLITTER = 4;
  // static const int SPHERE_VERTEX_SPLITTER = 5;

public:
  Abstract_vertex_splitter(int type = 0)
    : type_(type)
  { }

  virtual ~Abstract_vertex_splitter() { /*intentionally does nothing*/ }

  virtual PolyhedronSPtr split_vertex(const VertexSPtr& vertex) = 0; // abstract

  virtual int getType() const
  {
    return type_;
  }

  static bool check_splitted(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    bool result = false;
    PolyhedronSPtr polyhedron_cpy = polyhedron->clone();
    result = Transformation::shift_facets_deg1(polyhedron_cpy, -1.0);
    result = result && !Self_intersection::has_self_intersecting_surface(polyhedron_cpy);
    return result;
  }

  virtual std::string to_string() const
  {
    std::string result;
    switch (getType()) {
      case COMBI_VERTEX_SPLITTER:
        result = "Combi_vertex_splitter";
        break;
      case CONVEX_VERTEX_SPLITTER:
        result = "Convex_vertex_splitter";
        break;
      default:
        result = "Abstract_vertex_splitter";
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
