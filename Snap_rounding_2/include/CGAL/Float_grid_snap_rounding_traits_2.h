// Copyright (c) 2025 Geometry Factory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// author(s)     : Léo Valque

#ifndef CGAL_FLOAT_GRID_SNAP_ROUNDING_TRAITS_2_H
#define CGAL_FLOAT_GRID_SNAP_ROUNDING_TRAITS_2_H

#include <CGAL/license/Snap_rounding_2.h>

#include <CGAL/Arr_segment_traits_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <type_traits>

#include <CGAL/Named_function_parameters.h>

namespace CGAL {

#if DOXYGEN_RUNNING
/*!
\ingroup Snap_rounding_vertical_slab_grp

The class `Float_grid_snap_rounding_traits_2<InputKernel, ExactKernel>`
is a model of the `VerticalSlabSnapRoundingTraits_2` concept. It is identical to `Double_grid_snap_rounding_traits_2<InputKernel, ExactKernel>`,
except that points are rounded to single-precision floating-point coordinates.

\tparam InputKernel specifies the kernel used for the input
geometric objects. It must model the \cgal `Kernel` concept. These objects must be convertible to and from the types defined by
`ExactKernel` using `CGAL::Cartesian_converter`.

\tparam ExactKernel specifies the exact kernel used internally
to compute intersections and subdivision points before to round them. It must model the \cgal `Kernel` concept.

\cgalModels{VerticalSlabSnapRoundingTraits_2}

\sa{CGAL::snap_rounding_2()}
\sa{CGAL::vertical_slab_snap_rounding_2()}
\sa{CGAL::Double_grid_snap_rounding_traits_2}
\sa{CGAL::Integer_grid_snap_rounding_traits_2}
*/
template<typename InputKernel, typename ExactKernel = Exact_predicates_exact_constructions_kernel>
struct Float_grid_snap_rounding_traits_2
  :  Arr_segment_traits_2<ExactKernel>
{
  using Base = Arr_segment_traits_2<ExactKernel>;
#else

template<typename InputKernel, typename ExactKernel = Exact_predicates_exact_constructions_kernel, typename BaseTraits = Arr_segment_traits_2<ExactKernel> >
struct Float_grid_snap_rounding_traits_2
  : BaseTraits
{
  using Base = BaseTraits;
#endif

  using FT = typename Base::FT;
  using Target_FT = float;
  using Point_2   = typename Base::Point_2;
  using Segment_2 = typename Base::Segment_2;

  using Less_xy_2 = typename Base::Less_xy_2;
  using Less_y_2  = typename Base::Less_y_2;
  using Equal_2   = typename Base::Equal_2;

  using Construct_point_2   = typename Base::Construct_point_2;
  using Construct_source_2  = typename Base::Construct_source_2;
  using Construct_target_2  = typename Base::Construct_target_2;
  using Construct_segment_2 = typename Base::Construct_segment_2;

  using Evaluation_tag = Tag_true;
#if DOXYGEN_RUNNING
  using Evaluate = unspecified_type;
#else
  using Evaluate = internal::Evaluate<FT>;
#endif

  typedef Cartesian_converter<InputKernel, ExactKernel> Converter_to_exact;
  typedef Cartesian_converter<ExactKernel, InputKernel> Converter_from_exact;

  struct Evaluation{
    void operator()(const Point_2 &p) const{
      internal::Evaluate<FT>()(p);
    }
  };

  struct Construct_point_at_x_on_segment_2{
    Point_2 operator()(const Segment_2 &seg, const FT &x) const{
      FT y= (seg.supporting_line().y_at_x(x));
      return Base().construct_point_2_object()(x, y);
    }
  };

  struct Compute_squared_round_bound_2{
    double operator()(const FT &x) const{
      double inf = std::numeric_limits<double>::infinity();
      double b=std::nextafter(std::nextafterf((float) to_interval(x).second, -inf) - std::nextafterf((float) to_interval(x).first, inf), inf);
      return b*b;
    }
    double operator()(const Point_2 &p) const{
      return (*this)(p.x())+(*this)(p.y());
    }
    double operator()(const Segment_2 &seg) const{
      return (std::max)( (*this)(Base().construct_source_2_object()(seg)),
                         (*this)(Base().construct_target_2_object()(seg)));
    }
  };

  struct Construct_rounded_point_2{
    Target_FT operator()(const FT &x) const{
      return (float) to_double(x);
    }
    Point_2 operator()(const Point_2 &p) const{
      return Point_2((*this)(p.x()),(*this)(p.y()));
    }
  };

  Evaluate evaluate_object() const{ return Evaluate(); }
  Converter_to_exact converter_to_exact_object() const{ return Converter_to_exact(); }
  Converter_from_exact converter_from_exact_object() const{ return Converter_from_exact(); }
  Construct_point_at_x_on_segment_2 construct_point_at_x_on_segment_2_object() const{ return Construct_point_at_x_on_segment_2(); }
  Compute_squared_round_bound_2 compute_squared_round_bound_2_object() const{ return Compute_squared_round_bound_2(); }
  Construct_rounded_point_2 construct_rounded_point_2_object() const{ return Construct_rounded_point_2(); }
};

} //namespace CGAL


#endif
