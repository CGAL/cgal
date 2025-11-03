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

#include <CGAL/Arr_segment_traits_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Named_function_parameters.h>

#include  <CGAL/mutex.h>

namespace CGAL {

/*!
\ingroup PkgFloatSnapRounding2Ref

The class `Float_snap_rounding_traits_2<Kernel>` is a model of the
`FloatSnapRoundingTraits_2` concept, and is the only traits class supplied
with the package.
This class should be instantiated with a kernel with points and segments that are convertible to
`Exact_predicates_exact_constructions_kernel`, the user must specified to used another exact geometric kernel that conforms
to the \cgal kernel-concept, such as `Cartesian<Gmpq>`.

\cgalModels{SnapRoundingTraits_2}

*/
template<typename Input_Kernel, typename Exact_Kernel = Exact_predicates_exact_constructions_kernel>
struct Float_snap_rounding_traits_2: Arr_segment_traits_2<Exact_Kernel>{
  typedef Arr_segment_traits_2<Exact_Kernel> Base;

  typedef typename Base::FT        FT;
  typedef typename Base::Point_2   Point_2;
  typedef typename Base::Segment_2 Segment_2;

  typedef typename Base::Less_xy_2 Less_xy_2;
  typedef typename Base::Less_y_2  Less_y_2;
  typedef typename Base::Equal_2   Equal_2;

  typedef typename Base::Construct_point_2   Construct_point_2;
  typedef typename Base::Construct_source_2  Construct_source_2;
  typedef typename Base::Construct_target_2  Construct_target_2;
  typedef typename Base::Construct_segment_2 Construct_segment_2;

  typedef Cartesian_converter<Input_Kernel, Exact_Kernel> Converter_to_exact;
  typedef Cartesian_converter<Exact_Kernel, Input_Kernel> Converter_from_exact;

  // Return an upperbound of the squared distance between a point and its rounded value
  struct Squared_round_bound_2{
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

  struct Construct_round_point_2{
    double operator()(const FT &x) const{
      return to_double(x);
    }
    Point_2 operator()(const Point_2 &p) const{
      return Point_2((*this)(p.x()),(*this)(p.y()));
    }
  };

  struct Construct_point_at_x_on_segment_2{
    Point_2 operator()(const Segment_2 &seg, const FT &x) const{
      FT y= (seg.supporting_line().y_at_x(x));
      return Base().construct_point_2_object()(x, y);
    }
  };

  Converter_to_exact converter_to_exact_object() const{
    return Converter_to_exact();
  }

  Converter_from_exact converter_from_exact_object() const{
    return Converter_from_exact();
  }

  Squared_round_bound_2 squared_round_bound_2_object() const{
    return Squared_round_bound_2();
  }

  Construct_round_point_2 construct_round_point_2_object() const{
    return Construct_round_point_2();
  }

  Construct_point_at_x_on_segment_2 construct_point_at_x_on_segment_2_object() const{
    return Construct_point_at_x_on_segment_2();
  }
};

} //namespace CGAL

#endif
