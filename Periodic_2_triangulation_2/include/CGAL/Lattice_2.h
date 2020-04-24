// Copyright (c) 2019-2020 XXXXXXXX
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     :

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
    : basis_(CGAL::make_array(v0, v1)), k(k)
  {
    initialize();
  }

  Lattice_2(const Basis& basis, const K& k = K())
    : basis_(basis), k(k)
  {
    initialize();
  }

  Lattice_2(const Lattice_2& other) = default;

  // Access
  Basis& basis() { return basis_; }
  const Basis& basis() const { return basis_; }

  const Voronoi_face_normals& Voronoi_vectors() const { return Vfn_; }
  FT systole_sq_length() const { return systole_sq_length_; }

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
      FT c01 = k.compute_scalar_product_2_object()(basis_[0], basis_[1]);
      FT c00 = k.compute_scalar_product_2_object()(basis_[0], basis_[0]);
      FT c11 = k.compute_scalar_product_2_object()(basis_[1], basis_[1]);
      if(c11 < c00)
      {
        std::swap(basis_[0], basis_[1]);
        std::swap(c00, c11);
      }

      if(4*c01*c01 <= c00*c00)
      {
        // Basis is Lagrange-reduced.
        if(c01 > 0)
        {
          // Negate b1 if necessary to ensure obtuse angle between b0 and b1.
          basis_[1] = k.construct_opposite_vector_2_object()(basis_[1]);
        }
        reduced = true;
      }
      else
      {
        // Basis is not Lagrange-reduced.
        if(c01 > 0)
        {
          // b1 -= b0
          basis_[1] = k.construct_sum_of_vectors_2_object()(basis_[1],
                        k.construct_opposite_vector_2_object()(basis_[0]));
        }
        else
        {
          // b1 += b0
          basis_[1] = k.construct_sum_of_vectors_2_object()(basis_[1], basis_[0]);
        }
      }
    }

    CGAL_assertion(basis_is_reduced());
  }

  // Only used in assertion check.
  bool basis_is_reduced()
  {
    Vector ext = k.construct_opposite_vector_2_object()(
                   k.construct_sum_of_vectors_2_object()(basis_[0], basis_[1]));
    return k.compute_scalar_product_2_object()(basis_[0], basis_[1]) <= 0 &&
             k.compute_scalar_product_2_object()(basis_[0], ext) <= 0 &&
             k.compute_scalar_product_2_object()(basis_[1], ext) <= 0;
  }

  void construct_Voronoi_face_normals()
  {
    // @tmp is this really needed or can things be done with predicates?
    Vector third = k.construct_opposite_vector_2_object()(
                     k.construct_sum_of_vectors_2_object()(basis_[0], basis_[1]));

    Vfn_ = CGAL::make_array(basis_[0], basis_[1], third);

    // @todo check the reduction algorithm, might not be needed to also check the third's length
    systole_sq_length_ = (std::min)((std::min)(basis_[0].squared_length(),
                                               basis_[1].squared_length()),
                                    third.squared_length());

    std::cout << "Vectors:" << basis_[0] << "|" << basis_[1] << "|" << third << "|" << systole_sq_length_ << std::endl;
  }

  // Canonicalization
  // @fixme, this shouldn't take the offsetted point, but the canonical point and
  // the offset (to obtain an exact predicate)
  bool is_in_scaled_domain(const Point& p,
                           const FT scaling_factor = 1) const
  {
    for(int i=0; i<3; ++i)
    {
      const Vector& vfn = Vfn_[i];
      const Vector ptv(CGAL::ORIGIN, p);

      const FT sp = k.compute_scalar_product_2_object()(ptv, vfn) /
                      k.compute_scalar_product_2_object()(vfn, vfn);

      if(!(-0.5 * scaling_factor <= sp && sp < 0.5 * scaling_factor))
        return false;
    }

    return true;
  }

  Point construct_canonical_point(const Point& p) const
  {
    // @check It is fine to do constructions here because an approximation
    // of the exact canonical position of 'p' is fine: we only care about
    // consistency between that approximate position and its offsetted positions

    Point cp = p;

    // @todo factorize it with 'is_in_domain()' (somehow)
    int vfn_pos = 0;
    while(vfn_pos < 3)
    {
      const Vector& vfn = Vfn_[vfn_pos]; // @todo operator(int)
      const Vector ptv(CGAL::ORIGIN, cp);

      const FT sp = k.compute_scalar_product_2_object()(ptv, vfn) /
                      k.compute_scalar_product_2_object()(vfn, vfn);

      if(-0.5 <= sp && sp < 0.5)
      {
        ++vfn_pos;
      }
      else
      {
        Vector tv = vfn;
        tv = k.construct_scaled_vector_2_object()(tv, - std::floor(CGAL::to_double(sp + 0.5) ));
        cp = k.construct_translated_point_2_object()(cp, tv);
        vfn_pos = 0;
      }
    }

    return cp;
  }

  Point translate_by_offset(const Point& p, const Offset o) const
  {
//    std::cout << "Reduced basis b[0] = " << basis_[0] << std::endl;
//    std::cout << "Reduced basis b[1] = " << basis_[1] << std::endl;

//    std::cout << "translate_by_offset(" << p << " Off: " << o << ") = ";

    Vector translation = k.construct_sum_of_vectors_2_object()(
                           k.construct_scaled_vector_2_object()(basis_[0], o.x()),
                           k.construct_scaled_vector_2_object()(basis_[1], o.y()));

//    std::cout << k.construct_translated_point_2_object()(p, translation) << std::endl;

    return k.construct_translated_point_2_object()(p, translation);
  }

private:
  FT systole_sq_length_;
  std::array<Vector, 2> basis_;
  std::array<Vector, 3> Vfn_;

  K k;
};

namespace Periodic_2_triangulations_2 {
namespace internal {

template < typename K_,
           typename Construct_point_2_base_ = typename K_::Construct_point_2>
class Lattice_construct_point_2
  : public Construct_point_2_base_
{
  typedef Construct_point_2_base_            Base;
  typedef K_                                 Kernel;

  typedef typename Kernel::Point_2           Point;
  typedef CGAL::Periodic_2_offset_2          Offset;

  typedef Lattice_2<K_>                      Lattice;

public:
  Lattice_construct_point_2(const Lattice* lattice, const Base& cp)
    : Base(cp), lattice_(lattice)
  { }

  using Base::operator();

  Point operator()(const Point& p, const Offset& o) const
  {
    CGAL_assertion(lattice_ != nullptr);
    return lattice_->translate_by_offset(p, o);
  }

private:
  const Lattice* lattice_;
};

} // namespace internal
} // namespace Periodic_2_triangulations_2
} // namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_2_LATTICE_2_H
