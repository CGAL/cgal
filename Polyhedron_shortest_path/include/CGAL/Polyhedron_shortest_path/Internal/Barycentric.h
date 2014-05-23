// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_BARYCENTRIC_H
#define CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_BARYCENTRIC_H

namespace CGAL {

namespace Internal {

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

} // namespace Internal

} // namespace CGAL

#endif // CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_BARYCENTRIC_H
