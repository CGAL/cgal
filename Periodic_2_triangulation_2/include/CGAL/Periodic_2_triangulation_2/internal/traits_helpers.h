// Copyright (c) 2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_2_TRIANGULATION_2_TRAITS_HELPERS_H
#define CGAL_PERIODIC_2_TRIANGULATION_2_TRAITS_HELPERS_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_triangulation_2/Square_flat_torus_2/Construct_point_on_square_flat_torus_2.h>
#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Construct_point_on_lattice_2.h>
#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Lattice_2.h>

#include <CGAL/Has_conversion.h>
#include <CGAL/triangulation_assertions.h>

#include <type_traits>

namespace CGAL {
namespace P2T2 {
namespace internal {

template <typename K, typename O, typename D>
struct Construct_point_getter
{
  typedef void type;
};

// Iso_rectangle_2
template <typename K, typename O>
struct Construct_point_getter<K, O, typename K::Iso_rectangle_2>
{
  typedef Construct_point_on_square_flat_torus_2<K, O> type;
};

// Lattice_2
template <typename K, typename O>
struct Construct_point_getter<K, O, Lattice_2<K> >
{
  typedef Construct_point_on_lattice_2<K, O> type;
};

/// Default domains for traits constructors

// Iso_rectangle_2
template <typename K, typename D>
D get_default_domain(typename std::enable_if<std::is_same<D, typename K::Iso_rectangle_2>::value>::type* = nullptr)
{
  return D(0, 0, 1, 1);
}

// Lattice_2
template <typename K, typename D>
D get_default_domain(typename std::enable_if<std::is_same<D, Lattice_2<K> >::value>::type* = nullptr)
{
  typename K::Vector_2 v0(1, 0), v1(0, 1);
  return Lattice_2<K>(v0, v1); // centered at the origin
}

template <typename K, typename D>
D get_default_domain()
{
  CGAL_triangulation_assertion(false);
  return D();
}

/// Domain conversions (to exact/approx) for filtered traits

template <typename K1, typename D, typename K2>
struct Exact_domain_getter
{
  typedef void type;
};

// Iso_rectangle_2
template <typename K1, typename K2>
struct Exact_domain_getter<K1, typename K1::Iso_rectangle_2, K2>
{
  typedef typename K2::Iso_rectangle_2 type;
};

// Lattice_2
template <typename K1, typename K2>
struct Exact_domain_getter<K1, Lattice_2<K1>, K2>
{
  typedef Lattice_2<K2> type;
};

template <typename K1, typename K2>
Lattice_2<K2> convert_domain(const Lattice_2<K1>& /*lattice*/)
{
  // @exact
  // - what to do? Simply convert the members? What if the values (e.g. the basis reduction) obtained
  // with the input kernel are not correct (inexactness)
  // - Do the first computation with an exact kernel back into the double? There might be a loss
  // of precision on the way back to EK (AK->EK->AK-here->EK)
  // --> just do lattice computations in exact, and store the approximate and exact stuff in the lattice_2?
  CGAL_assertion(false);
  return Lattice_2<K2>();
}

template <typename K1, typename K2>
typename K2::Iso_rectangle_2 convert_domain(const typename K1::Iso_rectangle_2& domain)
{
  typename CGAL::internal::Converter_selector<K1, K2>::type converter;
  return converter(domain);
}

} // namespace internal
} // namespace P2T2

template <class Kernel_,
          class Offset_ = CGAL::Periodic_2_offset_2,
          class Domain_ = typename Kernel_::Iso_rectangle_2>
class Periodic_2_triangulation_traits_base_2;

namespace P2T2 {
namespace internal {

/// Offset conversions for filtered traits

// The Offset_converter is parametrized by a usual kernel converter,
// and adds the conversions for Offsets.
template <typename Converter_>
struct Offset_converter_2
  : public Converter_
{
  typedef Converter_                        Base;

  typedef typename Base::Source_kernel      Source_kernel;
  typedef typename Base::Target_kernel      Target_kernel;

  typedef typename Periodic_2_triangulation_traits_base_2<Source_kernel>::Offset   Source_off;
  typedef typename Periodic_2_triangulation_traits_base_2<Source_kernel>::Point_2  Source_pt;

  typedef typename Periodic_2_triangulation_traits_base_2<Target_kernel>::Offset   Target_off;
  typedef typename Periodic_2_triangulation_traits_base_2<Target_kernel>::Point_2  Target_pt;

  using Base::operator();

  Target_off operator()(const Source_off &off) const { return off; }
};

} // namespace internal
} // namespace P2T2
} // namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_2_TRAITS_HELPERS_H
