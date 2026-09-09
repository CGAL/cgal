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

#ifndef CGAL_INTEGER_GRID_SNAP_ROUNDING_TRAITS_2_H
#define CGAL_INTEGER_GRID_SNAP_ROUNDING_TRAITS_2_H

#include <CGAL/license/Snap_rounding_2.h>

#include <CGAL/Arr_segment_traits_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <type_traits>

#include <CGAL/Named_function_parameters.h>

namespace CGAL {

namespace internal::vertical_slab_snap_rounding_impl {

// Duplicate from PMP::triangle_soup_snap_rounding to avoid dependencies. TODO Factorize

// Certified ceil function for exact number types
template <class NT> double double_ceil(const Lazy_exact_nt< NT > &x);
template <class NT> double double_ceil(const NT &x);

template <class NT>
double double_ceil(const Lazy_exact_nt< NT > &x){
  // If both sides are in the same ceil, return this ceil
  double ceil_left=std::ceil(to_interval(x).first);
  if(ceil_left==std::ceil(to_interval(x).second))
    return ceil_left;
  // If not refine the interval by contracting the DAG and try again
  x.exact();
  ceil_left=std::ceil(to_interval(x).first);
  if(ceil_left==std::ceil(to_interval(x).second))
    return ceil_left;
  // If not return the ceil of the exact value
  return double_ceil( x.exact());
};

template <class NT>
double double_ceil(const NT &x){
  using FT = Fraction_traits<NT>;
  if constexpr(FT::Is_fraction::value){
    // If NT is a fraction, the ceil value is the result of the Euclidean division of the numerator and the denominator.
    typename FT::Numerator_type num, r, e;
    typename FT::Denominator_type denom;
    typename FT::Decompose()(x,num,denom);
    div_mod(num, denom, r, e);
    if((r>=0) && e!=0) //If the result is positive, the ceil value is one above
      return to_double(r+1);
    return to_double(r);
  } else {
    // Return the ceil of the approximation
    return std::ceil(to_double(x));
  }
};


}

#if DOXYGEN_RUNNING
/*!
\ingroup Snap_rounding_vertical_slab_grp

The class `Integer_grid_snap_rounding_traits_2<InputKernel, ExactKernel>`
is a model of the `VerticalSlabSnapRoundingTraits_2` concept. It is identical to `Double_grid_snap_rounding_traits_2<InputKernel, ExactKernel>`,
except that points are rounded to the closest point with integer coordinates.

\tparam InputKernel specifies the kernel used for the input
geometric objects. It must model the \cgal `Kernel` concept. These objects must be convertible to and from the types defined by
`ExactKernel` using `CGAL::Cartesian_converter`.

\tparam ExactKernel specifies the exact kernel used internally
to compute intersections and subdivision points before to round them. It must model the \cgal `Kernel` concept.

\cgalModels{VerticalSlabSnapRoundingTraits_2}

\sa{CGAL::snap_rounding_2()}
\sa{CGAL::vertical_slab_snap_rounding_2()}
\sa{CGAL::Double_grid_snap_rounding_traits_2}
\sa{CGAL::Float_grid_snap_rounding_traits_2}
*/
template<typename InputKernel, typename ExactKernel = Exact_predicates_exact_constructions_kernel>
struct Integer_grid_snap_rounding_traits_2
  :  Arr_segment_traits_2<ExactKernel>
{
  using Base = Arr_segment_traits_2<ExactKernel>;
#else

template<typename InputKernel, typename ExactKernel = Exact_predicates_exact_constructions_kernel, typename BaseTraits = Arr_segment_traits_2<ExactKernel> >
struct Integer_grid_snap_rounding_traits_2
  : BaseTraits
{
  using Base = BaseTraits;
#endif

  using FT = typename Base::FT;
  using Target_FT = double; // For precision
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

  Integer_grid_snap_rounding_traits_2(): m_pixel_size(1.0){}
  Integer_grid_snap_rounding_traits_2(double pixel_size): m_pixel_size(pixel_size){}

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
    Compute_squared_round_bound_2(): m_round_bound(1./4.){}
    Compute_squared_round_bound_2(double pixel_size): m_round_bound((pixel_size*pixel_size)/4.){}

    double operator()(const FT& /*x*/) const{
      return m_round_bound;
    }
    double operator()(const Point_2& /*p*/) const{
      return 2*m_round_bound;
    }
    double operator()(const Segment_2& /*seg*/) const{
      return 2*m_round_bound;
    }

  private:
    const double m_round_bound;
  };

  struct Construct_rounded_point_2{
    Construct_rounded_point_2(): m_pixel_size(1.0){}
    Construct_rounded_point_2(double pixel_size): m_pixel_size(pixel_size){}

    Target_FT operator()(const FT &x) const{
      return internal::vertical_slab_snap_rounding_impl::double_ceil((x / m_pixel_size) - 0.5) * m_pixel_size;
    }
    Point_2 operator()(const Point_2 &p) const{
      return Point_2((*this)(p.x()),(*this)(p.y()));
    }
  private:
    const double m_pixel_size;
  };

  Evaluate evaluate_object() const{ return Evaluate(); }
  Converter_to_exact converter_to_exact_object() const{ return Converter_to_exact(); }
  Converter_from_exact converter_from_exact_object() const{ return Converter_from_exact(); }
  Construct_point_at_x_on_segment_2 construct_point_at_x_on_segment_2_object() const{ return Construct_point_at_x_on_segment_2(); }
  Compute_squared_round_bound_2 compute_squared_round_bound_2_object() const{ return Compute_squared_round_bound_2(m_pixel_size); }
  Construct_rounded_point_2 construct_rounded_point_2_object() const{ return Construct_rounded_point_2(m_pixel_size); }

private:
  const double m_pixel_size;
};

} //namespace CGAL

#endif
