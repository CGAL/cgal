// Copyright (c) 2019-2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//                 Georg Osang

#ifndef CGAL_PERIODIC_3_TRIANGULATION_3_LATTICE_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_3_LATTICE_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/Periodic_3_offset_3.h>

#include <algorithm>
#include <array>
#include <iostream>

namespace CGAL {

template <typename K>
class Lattice_3
{
public:
  typedef typename K::FT                                 FT;
  typedef typename K::Vector_3                           Vector;
  typedef typename K::Point_3                            Point;

  typedef typename CGAL::Periodic_3_offset_3             Offset;

  typedef std::array<Vector, 3>                          Basis;
  typedef std::array<Vector, 7>                          Voronoi_face_normals;

  // Constructors
  Lattice_3() { }

  Lattice_3(const Vector& v0, const Vector& v1, const Vector& v2, const K& k = K())
    : basis_(CGAL::make_array(v0, v1, v2)), k(k)
  {
    initialize();
  }

  Lattice_3(const Basis& basis, const K& k = K())
    : basis_(basis), k(k)
  {
    initialize();
  }

  Lattice_3(const Lattice_3& other) = default;

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
    const std::vector<std::array<int, 4> > pair_partitions {{0,1,2,3},{0,2,1,3},{0,3,1,2},{1,2,0,3},{1,3,0,2},{2,3,0,1}};

    // List of basis vectors with their negative sum appended at the end
    std::vector<Vector> ext_basis = std::vector<Vector>(basis_.begin(), basis_.end());
    ext_basis.push_back(construct_extended_basis_vector());

    // Selling's algorithm to reduce a basis in 3D
    bool reduced = false;
    while (!reduced) {
      // lattice is reduced if all pair-wise scalar products are non-positive
      reduced = true;
      for (const std::array<int, 4>& pp : pair_partitions) {
        const FT sp = gt_.compute_scalar_product_3_object()(ext_basis[pp[0]], ext_basis[pp[1]]);
        if (sp > 0) {
          reduced = false;
          ext_basis[pp[2]] = k.construct_sum_of_vectors_3_object()(ext_basis[pp[2]], ext_basis[pp[0]]);
          ext_basis[pp[3]] = k.construct_sum_of_vectors_3_object()(ext_basis[pp[3]], ext_basis[pp[0]]);
          ext_basis[pp[0]] = k.construct_opposite_vector_3_object()(ext_basis[pp[0]]);
          break;
        }
      }
    }

    basis_ = CGAL::make_array(ext_basis[0], ext_basis[1], ext_basis[2]);
    CGAL_assertion(basis_is_reduced());
  }

  // Only used in assertion check.
  bool basis_is_reduced()
  {
    Vector ext = construct_extended_basis_vector();
    return k.compute_scalar_product_3_object()(basis_[0], basis_[1]) <= 0 &&
           k.compute_scalar_product_3_object()(basis_[0], basis_[2]) <= 0 &&
           k.compute_scalar_product_3_object()(basis_[1], basis_[2]) <= 0 &&
           k.compute_scalar_product_3_object()(basis_[0], ext) <= 0 &&
           k.compute_scalar_product_3_object()(basis_[1], ext) <= 0 &&
           k.compute_scalar_product_3_object()(basis_[2], ext) <= 0;
  }

  void construct_Voronoi_face_normals()
  {
    Vector ext = construct_extended_basis_vector();
    Vector fn1 = k.construct_sum_of_vectors_3_object()(basis_[0] + basis_[1]);
    Vector fn2 = k.construct_sum_of_vectors_3_object()(basis_[0] + basis_[2]);
    Vector fn3 = k.construct_sum_of_vectors_3_object()(basis_[1] + basis_[2]);

    Vfn_ = CGAL::make_array(basis_[0], basis_[1], basis_[2], extb, fn1, fn2, fn3);

    systole_sq_length_ = basis_[0].squared_length();
    for (int i = 1; i < 7; ++i) {
      systole_sq_length_ = (std::min)(systole_sq_length_, Vfn_[i].squared_length());
    }
  }

  // Canonicalization
  // @fixme, this shouldn't take the offsetted point, but the canonical point and
  // the offset (to obtain an exact predicate)
  bool is_in_scaled_domain(const Point& p,
                           const FT scaling_factor = 1) const
  {
    for(int i=0; i<7; ++i)
    {
      const Vector& vfn = Vfn_[i];
      const Vector ptv(CGAL::ORIGIN, p);

      const FT sp = k.compute_scalar_product_3_object()(ptv, vfn) /
                      k.compute_scalar_product_3_object()(vfn, vfn);

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
    while(vfn_pos < 7)
    {
      const Vector& vfn = Vfn_[vfn_pos]; // @todo operator(int)
      const Vector ptv(CGAL::ORIGIN, cp);

      const FT sp = k.compute_scalar_product_3_object()(ptv, vfn) /
                      k.compute_scalar_product_3_object()(vfn, vfn);

      if(-0.5 <= sp && sp < 0.5)
      {
        ++vfn_pos;
      }
      else
      {
        Vector tv = vfn;
        tv = k.construct_scaled_vector_3_object()(tv, - std::floor(CGAL::to_double(sp + 0.5) ));
        cp = k.construct_translated_point_3_object()(cp, tv);
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

    Vector translation = k.construct_sum_of_vectors_3_object()(
                           k.construct_sum_of_vectors_3_object()(
                             k.construct_scaled_vector_3_object()(basis_[0], o.x()),
                             k.construct_scaled_vector_3_object()(basis_[1], o.y())),
                           k.construct_scaled_vector_3_object()(basis_[2], o.z()));

//    std::cout << k.construct_translated_point_3_object()(p, translation) << std::endl;

    return k.construct_translated_point_3_object()(p, translation);
  }
private:
  Vector construct_extended_basis_vector() {
    return k.construct_opposite_vector_3_object()(
             k.construct_sum_of_vectors_3_object()(
               basis_[0],
               k.construct_sum_of_vectors_3_object()(basis_[1], basis_[2])
             ));
  }

private:
  FT systole_sq_length_;
  std::array<Vector, 3> basis_;
  std::array<Vector, 7> Vfn_;

  K k;
};

namespace Periodic_3_triangulations_3 {
namespace internal {

template < typename K_,
           typename Construct_point_3_base_ = typename K_::Construct_point_3>
class Lattice_construct_point_3
  : public Construct_point_3_base_
{
  typedef Construct_point_3_base_            Base;
  typedef K_                                 Kernel;

  typedef typename Kernel::Point_3           Point;
  typedef CGAL::Periodic_3_offset_3          Offset;

  typedef Lattice_3<K_>                      Lattice;

public:
  Lattice_construct_point_3(const Lattice* lattice, const Base& cp)
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
} // namespace Periodic_3_triangulations_3
} // namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_3_LATTICE_3_H
