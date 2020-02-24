// Copyright (c) 2014 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_BARYCENTRIC_H
#define CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_BARYCENTRIC_H

#include <CGAL/license/Surface_mesh_shortest_path.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/number_utils.h>

#include <utility>
#include <cstddef>

namespace CGAL {

namespace Surface_mesh_shortest_paths_3 {

template <class K, class B, class Construct_barycentric_coordinates>
class Construct_barycentric_coordinates_in_triangle_2
{
public:
  typedef typename K::FT FT;
  typedef typename K::Triangle_2 Triangle_2;
  typedef typename K::Point_2 Point_2;
  typedef typename K::Vector_2 Vector_2;

  typedef typename K::Construct_vector_2 Construct_vector_2;
  typedef typename K::Compute_scalar_product_2 Compute_scalar_product_2;

  typedef B result_type;

private:
  Construct_barycentric_coordinates m_construct_barycentric_coordinates;
  Construct_vector_2 m_construct_vector_2;
  Compute_scalar_product_2 m_compute_scalar_product_2;

public:
  Construct_barycentric_coordinates_in_triangle_2()
  {
  }

  Construct_barycentric_coordinates_in_triangle_2(const Construct_barycentric_coordinates& cbc,
                                                  const Construct_vector_2& cv2,
                                                  const Compute_scalar_product_2& csp2)
    : m_construct_barycentric_coordinates(cbc)
    , m_construct_vector_2(cv2)
    , m_compute_scalar_product_2(csp2)
  {
  }

  result_type operator()(const Triangle_2& t, const Point_2& p) const
  {
    Vector_2 v0 = m_construct_vector_2(t[0], t[1]);
    Vector_2 v1 = m_construct_vector_2(t[0], t[2]);
    Vector_2 v2 = m_construct_vector_2(t[0], p);

    FT d00 = m_compute_scalar_product_2(v0, v0);
    FT d01 = m_compute_scalar_product_2(v0, v1);
    FT d11 = m_compute_scalar_product_2(v1, v1);
    FT d20 = m_compute_scalar_product_2(v2, v0);
    FT d21 = m_compute_scalar_product_2(v2, v1);

    FT denom = d00 * d11 - d01 * d01;

    FT v = (d11 * d20 - d01 * d21) / denom;
    FT w = (d00 * d21 - d01 * d20) / denom;
    return m_construct_barycentric_coordinates(FT(1) - v - w, v, w);
  }
};

template <class K, class B, class Construct_barycentric_coordinates>
class Construct_barycentric_coordinates_in_triangle_3
{
public:
  typedef typename K::FT FT;
  typedef typename K::Triangle_3 Triangle_3;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;

  typedef typename K::Construct_vector_3 Construct_vector_3;
  typedef typename K::Compute_scalar_product_3 Compute_scalar_product_3;

  typedef B result_type;

private:
  Construct_barycentric_coordinates m_construct_barycentric_coordinates;
  Construct_vector_3 m_construct_vector_3;
  Compute_scalar_product_3 m_compute_scalar_product_3;

public:
  Construct_barycentric_coordinates_in_triangle_3()
  {
  }

  Construct_barycentric_coordinates_in_triangle_3(const Construct_barycentric_coordinates& cbc,
                                                  const Construct_vector_3& cv3,
                                                  const Compute_scalar_product_3& csp3)
    : m_construct_barycentric_coordinates(cbc)
    , m_construct_vector_3(cv3)
    , m_compute_scalar_product_3(csp3)
  {
  }

  result_type operator()(const Triangle_3& t, const Point_3& p) const
  {
    Vector_3 v0 = m_construct_vector_3(t[0], t[1]);
    Vector_3 v1 = m_construct_vector_3(t[0], t[2]);
    Vector_3 v2 = m_construct_vector_3(t[0], p);

    FT d00 = m_compute_scalar_product_3(v0, v0);
    FT d01 = m_compute_scalar_product_3(v0, v1);
    FT d11 = m_compute_scalar_product_3(v1, v1);
    FT d20 = m_compute_scalar_product_3(v2, v0);
    FT d21 = m_compute_scalar_product_3(v2, v1);

    FT denom = d00 * d11 - d01 * d01;

    FT v = (d11 * d20 - d01 * d21) / denom;
    FT w = (d00 * d21 - d01 * d20) / denom;
    return m_construct_barycentric_coordinates(FT(1) - v - w, v, w);
  }
};

/*!
  \brief Classification of the location of a 3-tuple barycentric coordinates in a triangle
*/
enum Barycentric_coordinates_type
{
  /// If the coordinates is invalid
  BARYCENTRIC_COORDINATES_INVALID = 0,
  /// if the coordinates has exactly one non-zero weight equal to 1, and the rest are zero
  BARYCENTRIC_COORDINATES_ON_VERTEX,
  ///if the coordinates has exactly one zero weight, and the rest sum to 1
  BARYCENTRIC_COORDINATES_ON_BOUNDARY,
  /// if the coordinates has no non-zero weight, and they all sum to 1
  BARYCENTRIC_COORDINATES_ON_BOUNDED_SIDE,
  /// if the weights of the coordinates do not sum to 1
  BARYCENTRIC_COORDINATES_ON_UNBOUNDED_SIDE
  #ifndef CGAL_NO_DEPRECATED_CODE
  // deprecated in CGAL 4.10
  , BARYCENTRIC_COORDINATE_INVALID = 0
  , BARYCENTRIC_COORDINATE_ON_VERTEX
  , BARYCENTRIC_COORDINATE_ON_BOUNDARY
  , BARYCENTRIC_COORDINATE_ON_BOUNDED_SIDE
  , BARYCENTRIC_COORDINATE_ON_UNBOUNDED_SIDE
  #endif
};

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated in CGAL 4.10
typedef Barycentric_coordinates_type Barycentric_coordinate_type;
#endif

template <class B, class Construct_barycentric_coordinates_weight>
class Classify_barycentric_coordinates
{
public:
  typedef B Barycentric_coordinates;
  typedef std::pair<Barycentric_coordinates_type, std::size_t> result_type;

private:
  Construct_barycentric_coordinates_weight m_construct_barycentric_coordinates_weight;

public:

  Classify_barycentric_coordinates()
  {
  }

  Classify_barycentric_coordinates(const Construct_barycentric_coordinates_weight& construct_barycentric_coordinates_weight)
    : m_construct_barycentric_coordinates_weight(construct_barycentric_coordinates_weight)
  {
  }

  result_type operator()(const Barycentric_coordinates& baryCoords)
  {
    Construct_barycentric_coordinates_weight cbcw;

    bool nonZero[3];
    std::size_t numNonZero = 0;

    if (cbcw(baryCoords, 0) + cbcw(baryCoords, 1) + cbcw(baryCoords, 2) > 1.00001 ||
        cbcw(baryCoords, 0) + cbcw(baryCoords, 1) + cbcw(baryCoords, 2) < 0.99999)
    {
      return std::make_pair(BARYCENTRIC_COORDINATES_ON_UNBOUNDED_SIDE, 0);
    }

    for (std::size_t i = 0; i < 3; ++i)
    {
      nonZero[i] = !CGAL::is_zero(cbcw(baryCoords, i));

      if (nonZero[i])
      {
        ++numNonZero;
      }
    }

    if (numNonZero == 3)
    {
      return std::make_pair(BARYCENTRIC_COORDINATES_ON_BOUNDED_SIDE, 0);
    }
    else if (numNonZero == 2)
    {
      std::size_t associatedEdge = 3;

      if (!nonZero[0])
      {
        associatedEdge = 1;
      }
      else if (!nonZero[1])
      {
        associatedEdge = 2;
      }
      else
      {
        associatedEdge = 0;
      }

      return std::make_pair(BARYCENTRIC_COORDINATES_ON_BOUNDARY, associatedEdge);
    }
    else if (numNonZero == 1)
    {
      std::size_t associatedEdge = 3;

      if (nonZero[0])
      {
        associatedEdge = 0;
      }
      else if (nonZero[1])
      {
        associatedEdge = 1;
      }
      else
      {
        associatedEdge = 2;
      }

      return std::make_pair(BARYCENTRIC_COORDINATES_ON_VERTEX, associatedEdge);
    }
    else
    {
      return std::make_pair(BARYCENTRIC_COORDINATES_INVALID, 0);
    }
  }

};

} // namespace Surface_mesh_shortest_paths_3

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_SHORTEST_PATHS_3_BARYCENTRIC_H
