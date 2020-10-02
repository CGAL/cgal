// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#ifndef CGAL_ARRANGEMENT_DEMO_FORWARD_DECLARATION_H
#define CGAL_ARRANGEMENT_DEMO_FORWARD_DECLARATION_H

namespace CGAL
{

// Arrangement Traits
template <typename Kernel_>
class Arr_segment_traits_2;

template <typename Kernel_>
class Arr_linear_traits_2;

template <typename SegmentTraits_2>
class Arr_polyline_traits_2;

template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class Arr_conic_traits_2;

template <
  typename RatKernel_, typename AlgKernel_, typename NtTraits_,
  typename BoundingTraits_>
class Arr_Bezier_curve_traits_2;

template <class Coefficient_>
class Arr_algebraic_segment_traits_2;

template <typename AlgebraicKernel_d_1>
class Arr_rational_function_traits_2;

template <
  class RatKernel_, class AlgKernel_, class NtTraits_, class BoundingTraits_>
class _Bezier_point_2;

// Tags
struct Arr_boundary_side_tag;
struct Arr_oblivious_side_tag;
struct Arr_open_side_tag;
struct Arr_closed_side_tag;
struct Arr_contracted_side_tag;
struct Arr_identified_side_tag;

template <typename T>
class Rational_traits;

template <typename FT_>
struct Cartesian;

class Object;

namespace Qt
{
class CurveInputMethod;
}

} // namespace CGAL

namespace demo_types
{
enum class TraitsType : int;
struct DemoTypes;
}

#endif
