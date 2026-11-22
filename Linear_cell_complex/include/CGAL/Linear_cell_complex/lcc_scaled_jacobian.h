// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef CGAL_LCC_SCALED_JACOBIAN_H
#define CGAL_LCC_SCALED_JACOBIAN_H

#include <CGAL/assertions.h>
#include <CGAL/Kernel/global_functions.h>
#include <algorithm>
#include <array>
#include <cmath>

namespace CGAL {
namespace internal {
////////////////////////////////////////////////////////////////////////////////
// Scaled Jacobian at one corner: determinant of the three unit edge vectors.
// 0 for a degenerate corner (an edge of length 0).
template<typename Vector>
double scaled_jacobian_of_corner(const Vector& v1, const Vector& v2,
                                 const Vector& v3)
{
  const double l=CGAL::to_double(v1.squared_length())*
                 CGAL::to_double(v2.squared_length())*
                 CGAL::to_double(v3.squared_length());
  if(l==0.) { return 0.; }
  return CGAL::to_double(CGAL::determinant(v1, v2, v3))/std::sqrt(l);
}
////////////////////////////////////////////////////////////////////////////////
/* Scaled Jacobian of a hexahedron given by its 8 vertices in this order,
 * i.e. the minimum over its 8 corners.
 *       4----7
 *      /|   /|
 *     5----6 |
 *     | 3--|-2
 *     |/   |/
 *     0----1
 */
template<typename LCC>
double scaled_jacobian_of_hexahedron_vertices
(const LCC& lcc,
 const std::array<typename LCC::Vertex_attribute_const_descriptor, 8>& vertices)
{
  // Each corner, followed by its three neighbors.
  static const int corners[8][4]={{0,1,3,5}, {3,0,2,4}, {2,3,1,7}, {1,2,0,6},
                                  {4,7,5,3}, {5,4,6,0}, {6,5,7,1}, {7,6,4,2}};
  using Vector=typename LCC::Vector;
  double res=1.;
  for(const auto& c: corners)
  {
    const auto& p=lcc.point_of_vertex_attribute(vertices[c[0]]);
    res=(std::min)(res, scaled_jacobian_of_corner
                   (Vector(p, lcc.point_of_vertex_attribute(vertices[c[1]])),
                    Vector(p, lcc.point_of_vertex_attribute(vertices[c[2]])),
                    Vector(p, lcc.point_of_vertex_attribute(vertices[c[3]]))));
  }
  return res;
}
////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
double jacobian_tri_hex_for_dart_v1(const LCC& lcc,
                                    typename LCC::Dart_const_descriptor d)
{
  // Get the vertices of virtual hexahedra, considering virtually merging
  // face(d) and face(beta2(d)), and face(beta1,2(d)) and face(beta1,2,1,2(d)).
  std::array<typename LCC::Vertex_attribute_const_descriptor, 8> vertices;
  vertices[0]=lcc.vertex_attribute(lcc.beta(d, 0));
  vertices[1]=lcc.vertex_attribute(d);
  vertices[2]=lcc.vertex_attribute(lcc.beta(d, 2, 0));
  vertices[3]=lcc.vertex_attribute(lcc.beta(d, 1));

  vertices[4]=lcc.vertex_attribute(lcc.beta(d, 1, 2, 1, 2, 0));
  vertices[5]=lcc.vertex_attribute(lcc.beta(d, 1, 2, 0));
  vertices[6]=lcc.vertex_attribute(lcc.beta(d, 0, 2, 0));
  if(vertices[6]==vertices[5])
  { vertices[6]=lcc.vertex_attribute(lcc.beta(d, 0, 2, 0, 2, 0)); }
  vertices[7]=lcc.vertex_attribute(lcc.beta(d, 2, 1, 2, 0));
  if(vertices[7]==vertices[6])
  { vertices[7]=lcc.vertex_attribute(lcc.beta(d, 2, 1, 2, 0, 2, 0)); }

  return scaled_jacobian_of_hexahedron_vertices(lcc, vertices);
}
////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
double jacobian_tri_hex_for_dart_v2(const LCC& lcc,
                                    typename LCC::Dart_const_descriptor d)
{
  // Get the vertices of virtual hexahedra, considering virtually merging
  // face(d) and face(beta2(d)), and face(beta1,2(d)) and face(beta1,2,0,2(d)).
  std::array<typename LCC::Vertex_attribute_const_descriptor, 8> vertices;
  vertices[0]=lcc.vertex_attribute(lcc.beta(d, 0));
  vertices[1]=lcc.vertex_attribute(d);
  vertices[2]=lcc.vertex_attribute(lcc.beta(d, 2, 0));
  vertices[3]=lcc.vertex_attribute(lcc.beta(d, 1));

  vertices[4]=lcc.vertex_attribute(lcc.beta(d, 1, 2, 0));
  vertices[5]=lcc.vertex_attribute(lcc.beta(d, 1, 2, 0, 2, 0));
  vertices[6]=lcc.vertex_attribute(lcc.beta(d, 0, 2, 0));
  if(vertices[6]==vertices[5])
  { vertices[6]=lcc.vertex_attribute(lcc.beta(d, 0, 2, 0, 2, 0)); }
  vertices[7]=lcc.vertex_attribute(lcc.beta(d, 2, 1, 2, 0));
  if(vertices[7]==vertices[6])
  { vertices[7]=lcc.vertex_attribute(lcc.beta(d, 2, 1, 2, 0, 2, 0)); }

  return scaled_jacobian_of_hexahedron_vertices(lcc, vertices);
}
////////////////////////////////////////////////////////////////////////////////
// Scaled Jacobian of a hexahedron whose faces are triangulated.
template<typename LCC>
double scaled_jacobian_of_triangulated_hexahedron
(const LCC& lcc, typename LCC::Dart_const_descriptor d)
{
  return (std::max)({jacobian_tri_hex_for_dart_v1(lcc, d),
                     jacobian_tri_hex_for_dart_v1(lcc, lcc.beta(d, 1)),
                     jacobian_tri_hex_for_dart_v1(lcc, lcc.beta(d, 0)),
                     jacobian_tri_hex_for_dart_v2(lcc, d),
                     jacobian_tri_hex_for_dart_v2(lcc, lcc.beta(d, 1)),
                     jacobian_tri_hex_for_dart_v2(lcc, lcc.beta(d, 0))});
}
////////////////////////////////////////////////////////////////////////////////
} // namespace internal

////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
double scaled_jacobian_of_hexahedron(const LCC& lcc,
                                     typename LCC::Dart_const_descriptor d)
{
  CGAL_precondition(lcc.is_volume_combinatorial_hexahedron(d));

  std::array<typename LCC::Vertex_attribute_const_descriptor, 8> vertices;
  vertices[0]=lcc.vertex_attribute(d);
  vertices[1]=lcc.vertex_attribute(lcc.beta(d, 0));
  vertices[2]=lcc.vertex_attribute(lcc.beta(d, 0, 2, 1, 1));
  vertices[3]=lcc.vertex_attribute(lcc.beta(d, 2, 1, 1));

  vertices[4]=lcc.vertex_attribute(lcc.beta(d, 2, 0));
  vertices[5]=lcc.vertex_attribute(lcc.beta(d, 2));
  vertices[6]=lcc.vertex_attribute(lcc.beta(d, 1, 1));
  vertices[7]=lcc.vertex_attribute(lcc.beta(d, 1, 2, 0));

  return internal::scaled_jacobian_of_hexahedron_vertices(lcc, vertices);
}
////////////////////////////////////////////////////////////////////////////////
} // namespace CGAL

#endif // CGAL_LCC_SCALED_JACOBIAN_H
