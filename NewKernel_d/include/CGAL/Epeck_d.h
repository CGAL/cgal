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

#ifndef CGAL_EPECK_D_H
#define CGAL_EPECK_D_H

#include <CGAL/disable_warnings.h>

#include <CGAL/NewKernel_d/Cartesian_base.h>
#include <CGAL/NewKernel_d/Wrapper/Cartesian_wrap.h>
#include <CGAL/NewKernel_d/Kernel_d_interface.h>
#include <CGAL/NewKernel_d/Lazy_cartesian.h>
#include <CGAL/NewKernel_d/KernelD_converter.h>
#include <CGAL/internal/Exact_type_selector.h>
#include <CGAL/NewKernel_d/Types/Weighted_point.h>

// TODO: add static filters somewhere
namespace CGAL {
#define CGAL_KA Cartesian_base_d<Interval_nt_advanced,Dim>
#define CGAL_KE Cartesian_base_d<internal::Exact_field_selector<double>::Type, Dim>
template<class Dim> using Epeck_d_help1 = Lazy_cartesian<CGAL_KE, CGAL_KA, KernelD_converter<CGAL_KE, CGAL_KA>>;
#undef CGAL_KE
#undef CGAL_KA
#undef CGAL_BASE
#define CGAL_BASE \
  Kernel_d_interface< \
    Cartesian_wrap< \
      Epeck_d_help1<Dim>, \
      Epeck_d<Dim> > >
template<class Dim>
struct Epeck_d
: CGAL_BASE
{
  constexpr Epeck_d(){}
  constexpr Epeck_d(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
}

#include <CGAL/enable_warnings.h>

#endif
