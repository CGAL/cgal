// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_BARYCENTRIC_H
#define CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_BARYCENTRIC_H

#include <CGAL/number_utils.h>

namespace CGAL {

namespace internal {

template <class FT, class B, class P, class V, class T>
inline B Construct_barycentric_coordinate_any(const T& t, const P& p)
{
  V v0 = t[1] - t[0];
  V v1 = t[2] - t[0];
  V v2 = p - t[0];

  FT d00 = v0 * v0;
  FT d01 = v0 * v1;
  FT d11 = v1 * v1;
  FT d20 = v2 * v0;
  FT d21 = v2 * v1;
  FT denom = d00 * d11 - d01 * d01;

  FT v = (d11 * d20 - d01 * d21) / denom;
  FT w = (d00 * d21 - d01 * d20) / denom;
  return B(FT(1.0) - v - w, v, w);
}

template <class FT, class B, class P, class T>
inline P Construct_triangle_location_any(const T& t, const B& a)
{
  return P(ORIGIN) + (((t[0] - P(ORIGIN)) * a[0]) + ((t[1] - P(ORIGIN)) * a[1]) + ((t[2] - P(ORIGIN)) * a[2]));
}

enum Barycentric_coordinate_type
{
  BARYCENTRIC_COORDINATE_INVALID = 0,
  BARYCENTRIC_COORDINATE_VERTEX,
  BARYCENTRIC_COORDINATE_EDGE,
  BARYCENTRIC_COORDINATE_INTERNAL,
  BARYCENTRIC_COORDINATE_EXTERNAL,
};

template <class B>
Barycentric_coordinate_type classify_barycentric_coordinate(const B& baryCoords, size_t& associatedEdge)
{
  bool nonZero[3];
  size_t numNonZero = 0;

  if (CGAL::is_positive(baryCoords * baryCoords - 1.0))
  {
    return BARYCENTRIC_COORDINATE_EXTERNAL;
  }
  
  for (size_t i = 0; i < 3; ++i)
  {
    nonZero[i] = !CGAL::is_zero(baryCoords[i]);
    
    if (nonZero[i])
    {
      ++numNonZero;
    }
  }

  if (numNonZero == 3)
  {
    associatedEdge = 3;
    return BARYCENTRIC_COORDINATE_INTERNAL;
  }
  else if (numNonZero == 2)
  {
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

    return BARYCENTRIC_COORDINATE_EDGE;
  }
  else if (numNonZero == 1)
  {
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

    return BARYCENTRIC_COORDINATE_VERTEX;
  }
  else
  {
    return BARYCENTRIC_COORDINATE_INVALID;
  }
}

template <class B>
Barycentric_coordinate_type classify_barycentric_coordinate(const B& baryCoords)
{
  size_t dummy;
  return classify_barycentric_coordinate(baryCoords, dummy);
}

} // namespace Internal

} // namespace CGAL

#endif // CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_BARYCENTRIC_H
