// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_EPICK_D_H
#define CGAL_EPICK_D_H

#include <CGAL/disable_warnings.h>

#include <CGAL/NewKernel_d/Cartesian_base.h>
#include <CGAL/NewKernel_d/Cartesian_static_filters.h>
#include <CGAL/NewKernel_d/Cartesian_filter_K.h>
#include <CGAL/NewKernel_d/Wrapper/Cartesian_wrap.h>
#include <CGAL/NewKernel_d/Kernel_d_interface.h>
#include <CGAL/internal/Exact_type_selector.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/NewKernel_d/Types/Weighted_point.h>


namespace CGAL {
#define CGAL_BASE \
  Cartesian_filter_K< \
    Cartesian_base_d<double, Dim>, \
    Cartesian_base_d<Interval_nt_advanced, Dim>, \
    Cartesian_base_d<internal::Exact_field_selector<double>::Type, Dim> \
  >
template<class Dim>
struct Epick_d_help1
: CGAL_BASE
{
  constexpr Epick_d_help1(){}
  constexpr Epick_d_help1(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
#define CGAL_BASE \
  Cartesian_filter_K< \
    Epick_d_help1<Dim>, \
    Cartesian_base_d<Interval_nt_advanced, Dim>, \
    Cartesian_base_d<internal::Exact_ring_selector<double>::Type, Dim>, \
    typename Functors_without_division<Dim>::type \
  >
template<class Dim>
struct Epick_d_help2
: CGAL_BASE
{
  constexpr Epick_d_help2(){}
  constexpr Epick_d_help2(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
#define CGAL_BASE \
  Cartesian_static_filters<Dim,Epick_d_help2<Dim>,Epick_d_help3<Dim> >

template<class Dim>
struct Epick_d_help3
: CGAL_BASE
{
  constexpr Epick_d_help3(){}
  constexpr Epick_d_help3(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
#define CGAL_BASE \
  Kernel_d_interface< \
    Cartesian_wrap< \
      Epick_d_help3<Dim>, \
      Epick_d<Dim> > >
template<class Dim>
struct Epick_d
: CGAL_BASE
{
  constexpr Epick_d(){}
  constexpr Epick_d(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
}

#include <CGAL/enable_warnings.h>

#endif
