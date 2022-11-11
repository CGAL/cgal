// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef CGAL_ARRANGEMENTS_ARR_TRAITS_ADAPTOR_H
#define CGAL_ARRANGEMENTS_ARR_TRAITS_ADAPTOR_H

#include "ForwardDeclarations.h"

/**
 * Support for new ArrTraits should specify types:
 *
 * Kernel - a not-necessarily-exact kernel to represent the arrangement
 * graphically. We'll use the Point_2 type provided by this kernel for
 * computing distances
 * Point_2 - the point type used in the particular arrangement
 * CoordinateType - the coordinate type used by the point type
 */
template <typename ArrTraits>
class ArrTraitsAdaptor
{
};

template <typename Kernel_>
class ArrTraitsAdaptor<CGAL::Arr_segment_traits_2<Kernel_>>
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_segment_traits_2<Kernel> ArrTraits;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename Kernel_>
class ArrTraitsAdaptor<CGAL::Arr_linear_traits_2<Kernel_>>
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_linear_traits_2<Kernel> ArrTraits;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename SegmentTraits>
class ArrTraitsAdaptor<CGAL::Arr_polyline_traits_2<SegmentTraits>>
{
public:
  typedef CGAL::Arr_polyline_traits_2<SegmentTraits> ArrTraits;
  typedef typename SegmentTraits::Kernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class ArrTraitsAdaptor<CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>
{
public:
  typedef CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits> ArrTraits;
  typedef AlgKernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
class ArrTraitsAdaptor<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>>
{
public:
  typedef CGAL::Arr_Bezier_curve_traits_2<
    RatKernel, AlgKernel, NtTraits, BoundingTraits>
    ArrTraits;
  typedef RatKernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename Coefficient_>
class ArrTraitsAdaptor<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>
{
public:
  typedef Coefficient_ Coefficient;
  typedef typename CGAL::Arr_algebraic_segment_traits_2<Coefficient> ArrTraits;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename ArrTraits::Algebraic_real_1 CoordinateType;
  typedef CGAL::Cartesian<typename ArrTraits::Bound> Kernel;
};

template <typename AlgebraicKernel_d_1>
class ArrTraitsAdaptor<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>
{
public:
  typedef typename CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>
    ArrTraits;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename ArrTraits::Algebraic_real_1 CoordinateType;
  typedef CGAL::Cartesian<typename ArrTraits::Bound> Kernel;
};

#endif
