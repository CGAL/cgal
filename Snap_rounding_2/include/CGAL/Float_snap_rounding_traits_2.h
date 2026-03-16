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

#ifndef CGAL_FLOAT_SNAP_ROUNDING_TRAITS_2_H
#define CGAL_FLOAT_SNAP_ROUNDING_TRAITS_2_H

#include <CGAL/license/Snap_rounding_2.h>

#include <CGAL/Arr_segment_traits_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <type_traits>

#include <CGAL/Named_function_parameters.h>

namespace CGAL {

namespace internal {
template<typename Input_Kernel, typename Exact_Kernel = Exact_predicates_exact_constructions_kernel, typename BaseTraits = Arr_segment_traits_2<Exact_Kernel> >
struct Float_snap_rounding_traits_base_2: BaseTraits{
  using Base = BaseTraits;

  using FT = typename Base::FT;
  using Target_FT = double;
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
  using Evaluate = internal::Evaluate<FT>;

  typedef Cartesian_converter<Input_Kernel, Exact_Kernel> Converter_to_exact;
  typedef Cartesian_converter<Exact_Kernel, Input_Kernel> Converter_from_exact;

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

  Evaluate evaluate_object() const{ return Evaluate(); }
  Converter_to_exact converter_to_exact_object() const{ return Converter_to_exact(); }
  Converter_from_exact converter_from_exact_object() const{ return Converter_from_exact(); }

  Construct_point_at_x_on_segment_2 construct_point_at_x_on_segment_2_object() const{ return Construct_point_at_x_on_segment_2(); }
};
}

/*!
\ingroup PkgFloatSnapRounding2Ref

The class `Double_snap_rounding_traits_2<Kernel>` is a model of the
`FloatSnapRoundingTraits_2` concept. The points are rounded on double precision floating point coordinates.

This class should be instantiated with a kernel with points and segments that are convertible to
`Exact_predicates_exact_constructions_kernel`, the user must specified to used another exact geometric kernel that conforms
to the \cgal kernel-concept, such as `Cartesian<Gmpq>`.

\cgalModels{SnapRoundingTraits_2}
*/
template<typename Input_Kernel, typename Exact_Kernel = Exact_predicates_exact_constructions_kernel, typename BaseTraits = Arr_segment_traits_2<Exact_Kernel> >
struct Double_snap_rounding_traits_2: internal::Float_snap_rounding_traits_base_2<Input_Kernel, Exact_Kernel, BaseTraits>{
  using Target_FT = double;
  using Base = internal::Float_snap_rounding_traits_base_2<Input_Kernel, Exact_Kernel, BaseTraits>;
  using FT = typename Base::FT;
  using Point_2 = typename Base::Point_2;
  using Segment_2 = typename Base::Segment_2;

  struct Compute_squared_round_bound_2{
    double operator()(const FT &x) const{
      double b=std::nextafter(to_interval(x).second - to_interval(x).first, std::numeric_limits<double>::infinity());
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
      return to_double(x);
    }
    Point_2 operator()(const Point_2 &p) const{
      return Point_2((*this)(p.x()),(*this)(p.y()));
    }
  };
  Compute_squared_round_bound_2 compute_squared_round_bound_2_object() const{ return Compute_squared_round_bound_2(); }
  Construct_rounded_point_2 construct_rounded_point_2_object() const{ return Construct_rounded_point_2(); }
};

/*!
\ingroup PkgFloatSnapRounding2Ref

The class `Float_snap_rounding_traits_2<Kernel>` is a model of the `FloatSnapRoundingTraits_2` concept. It is identical to `Double_snap_rounding_traits_2<Kernel>`,
except that points are rounded to single-precision floating-point coordinates.

\cgalModels{SnapRoundingTraits_2}
*/
template<typename Input_Kernel, typename Exact_Kernel = Exact_predicates_exact_constructions_kernel, typename BaseTraits = Arr_segment_traits_2<Exact_Kernel> >
struct Float_snap_rounding_traits_2: internal::Float_snap_rounding_traits_base_2<Input_Kernel, Exact_Kernel, BaseTraits>{
  using Target_FT = float;
  using Base = internal::Float_snap_rounding_traits_base_2<Input_Kernel, Exact_Kernel, BaseTraits>;
  using FT = typename Base::FT;
  using Point_2 = typename Base::Point_2;
  using Segment_2 = typename Base::Segment_2;

  struct Compute_squared_round_bound_2{
    Target_FT operator()(const FT &x) const{
      Target_FT inf = std::numeric_limits<double>::infinity();
      Target_FT b=std::nextafter(std::nextafterf(to_interval(x).second, -inf) - std::nextafterf(to_interval(x).first, inf), inf);
      return b*b;
    }
    Target_FT operator()(const Point_2 &p) const{
      return (*this)(p.x())+(*this)(p.y());
    }
    Target_FT operator()(const Segment_2 &seg) const{
      return (std::max)( (*this)(Base().construct_source_2_object()(seg)),
                         (*this)(Base().construct_target_2_object()(seg)));
    }
  };

  struct Construct_rounded_point_2{
    float operator()(const FT &x) const{
      return (float) to_double(x);
    }
    Point_2 operator()(const Point_2 &p) const{
      return Point_2((*this)(p.x()),(*this)(p.y()));
    }
  };
  Compute_squared_round_bound_2 compute_squared_round_bound_2_object() const{ return Compute_squared_round_bound_2(); }
  Construct_rounded_point_2 construct_rounded_point_2_object() const{ return Construct_rounded_point_2(); }
};

} //namespace CGAL

#endif
