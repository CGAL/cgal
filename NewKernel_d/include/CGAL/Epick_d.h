// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Marc Glisse

#ifndef CGAL_EPICK_D_H
#define CGAL_EPICK_D_H
#include <CGAL/NewKernel_d/Cartesian_base.h>
#include <CGAL/NewKernel_d/Cartesian_static_filters.h>
#include <CGAL/NewKernel_d/Cartesian_filter_K.h>
#include <CGAL/NewKernel_d/Wrapper/Cartesian_wrap.h>
#include <CGAL/NewKernel_d/Kernel_d_interface.h>
#include <CGAL/internal/Exact_type_selector.h>
#include <CGAL/Interval_nt.h>


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
  CGAL_CONSTEXPR Epick_d_help1(){}
  CGAL_CONSTEXPR Epick_d_help1(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
#define CGAL_BASE \
  Cartesian_static_filters<Dim,Epick_d_help1<Dim>,Epick_d_help2<Dim> >
template<class Dim>
struct Epick_d_help2
: CGAL_BASE
{
  CGAL_CONSTEXPR Epick_d_help2(){}
  CGAL_CONSTEXPR Epick_d_help2(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
#define CGAL_BASE \
  Kernel_d_interface< \
    Cartesian_wrap< \
      Epick_d_help2<Dim>, \
      Epick_d<Dim> > >
template<class Dim>
struct Epick_d
: CGAL_BASE
{
  CGAL_CONSTEXPR Epick_d(){}
  CGAL_CONSTEXPR Epick_d(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
}
#endif
