// Copyright (c) 2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Georg Osang
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_2_TRIANGULATION_2_LATTICE_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_2_LATTICE_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_offset_2.h>

#include <algorithm>
#include <array>
#include <iostream>

namespace CGAL {

template <typename K>
class Lattice_2
{
public:
  typedef typename K::FT                                 FT;
  typedef typename K::Vector_2                           Vector;
  typedef typename K::Point_2                            Point;

  typedef typename CGAL::Periodic_2_offset_2             Offset;

  typedef std::array<Vector, 2>                          Basis;
  typedef std::array<Vector, 3>                          Voronoi_face_normals;

  // Constructors
  Lattice_2() { }

  Lattice_2(const Vector& v0, const Vector& v1, const K& k = K())
    : _origin(CGAL::ORIGIN), _basis(CGAL::make_array(v0, v1)), _k(k)
  {
    initialize();
  }

  Lattice_2(const Basis& basis, const K& k = K())
    : _origin(CGAL::ORIGIN), _basis(basis), _k(k)
  {
    initialize();
  }

  Lattice_2(const Point& origin, const Vector& v0, const Vector& v1, const K& k = K())
    : _origin(origin), _basis(CGAL::make_array(v0, v1)), _k(k)
  {
    initialize();
  }

  Lattice_2(const Point& origin, const Basis& basis, const K& k = K())
    : _origin(origin), _basis(basis), _k(k)
  {
    initialize();
  }

  Lattice_2(const Lattice_2& other) = default;

  // Access
  Basis& basis() { return _basis; }
  const Basis& basis() const { return _basis; }

  const std::array<std::array<int, 2>, 12>& overlapping_offsets() const { return _overlapping_offsets; }

  const Voronoi_face_normals& Voronoi_vectors() const { return _Vfn; }
  FT systole_sq_length() const { return _systole_sq_length; }

  // Initialization
  void initialize()
  {
    reduce_basis();
    construct_Voronoi_face_normals();
  }

  void reduce_basis()
  {
    bool reduced = false;
    while(!reduced)
    {
      FT c01 = _k.compute_scalar_product_2_object()(_basis[0], _basis[1]);
      FT c00 = _k.compute_scalar_product_2_object()(_basis[0], _basis[0]);
      FT c11 = _k.compute_scalar_product_2_object()(_basis[1], _basis[1]);
      if(c11 < c00)
      {
        std::swap(_basis[0], _basis[1]);
        std::swap(c00, c11);
      }

      if(4*c01*c01 <= c00*c00)
      {
        // Basis is Lagrange-reduced.
        if(c01 > FT(0))
        {
          // Negate b1 if necessary to ensure obtuse angle between b0 and b1.
          _basis[1] = _k.construct_opposite_vector_2_object()(_basis[1]);
        }
        reduced = true;
      }
      else
      {
        // Basis is not Lagrange-reduced.
        if(c01 > 0)
        {
          // b1 -= b0
          _basis[1] = _k.construct_sum_of_vectors_2_object()(_basis[1],
                        _k.construct_opposite_vector_2_object()(_basis[0]));
        }
        else
        {
          // b1 += b0
          _basis[1] = _k.construct_sum_of_vectors_2_object()(_basis[1], _basis[0]);
        }
      }
    }

    CGAL_postcondition(_basisis_reduced());
  }

  // Only used in assertion check.
  bool _basisis_reduced()
  {
    Vector ext = _k.construct_opposite_vector_2_object()(
                   _k.construct_sum_of_vectors_2_object()(_basis[0], _basis[1]));
    return _k.compute_scalar_product_2_object()(_basis[0], _basis[1]) <= 0 &&
           _k.compute_scalar_product_2_object()(_basis[0], ext) <= 0 &&
           _k.compute_scalar_product_2_object()(_basis[1], ext) <= 0;
  }

  void construct_Voronoi_face_normals()
  {
    Vector third = _k.construct_opposite_vector_2_object()(
                     _k.construct_sum_of_vectors_2_object()(_basis[0], _basis[1]));

    _Vfn = CGAL::make_array(_basis[0], _basis[1], third);

    // @todo check the reduction algorithm, might not be needed to also check the third's length
    _systole_sq_length = (std::min)((std::min)(_basis[0].squared_length(),
                                               _basis[1].squared_length()),
                                    third.squared_length());
  }

  // Canonicalization
  // @exact
  bool is_in_scaled_domain(const Point& p,
                           const FT scaling_factor = 1) const
  {
    for(int i=0; i<3; ++i)
    {
      const Vector& vfn = _Vfn[i];
      const Vector ptv(_origin, p);

      const FT sp = FT(2) * _k.compute_scalar_product_2_object()(ptv, vfn) /
                              _k.compute_scalar_product_2_object()(vfn, vfn);

      if(scaling_factor <= sp || sp < -scaling_factor)
        return false;
    }

    return true;
  }

  // @exact
  Point construct_canonical_point(const Point& p) const
  {
    Point cp = p;

    // @todo factorize it with 'is_in_domain()' (somehow)
    int vfn_pos = 0;
    while(vfn_pos < 3)
    {
      const Vector& vfn = _Vfn[vfn_pos];
      const Vector ptv(_origin, cp);

      const FT sp = _k.compute_scalar_product_2_object()(ptv, vfn) /
                      _k.compute_scalar_product_2_object()(vfn, vfn);

      if(-0.5 <= sp && sp < 0.5)
      {
        ++vfn_pos;
      }
      else
      {
        Vector tv = vfn;
        tv = _k.construct_scaled_vector_2_object()(tv, - std::floor(CGAL::to_double(sp + 0.5)));
        cp = _k.construct_translated_point_2_object()(cp, tv);
        vfn_pos = 0;
      }
    }

    return cp;
  }

  Point translate_by_offset(const Point& p, const Offset o) const
  {
//    std::cout << "translate_by_offset(" << p << " Off: " << o << ") = ";

//    std::cout << "Reduced basis b[0] = " << _basis[0] << std::endl;
//    std::cout << "Reduced basis b[1] = " << _basis[1] << std::endl;

    Vector translation = _k.construct_sum_of_vectors_2_object()(
                           _k.construct_scaled_vector_2_object()(_basis[0], o.x()),
                           _k.construct_scaled_vector_2_object()(_basis[1], o.y()));

    return _k.construct_translated_point_2_object()(p, translation);
  }

public:
  // A list of those offsets such that the domain translated along the offset
  // overlaps the scaled domain.
  std::array<std::array<int, 2>, 12> _overlapping_offsets =
  {{
    // entirely contained in scaled domains
    {-1, -1}, {0, 1}, {1, 0}, {-1, 0}, {0, -1}, {1, 1},
    // intersecting the scaled domain
    {-1, -2}, {1, 2}, {-2, -1}, {2, 1}, {-1, 1}, {1, -1}
  }};

private:
  Point _origin;
  std::array<Vector, 2> _basis;
  std::array<Vector, 3> _Vfn;
  FT _systole_sq_length;

  K _k;
};

} // namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_2_LATTICE_2_H
