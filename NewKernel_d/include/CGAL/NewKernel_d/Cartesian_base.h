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

#ifndef CGAL_KERNEL_D_CARTESIAN_BASE_H
#define CGAL_KERNEL_D_CARTESIAN_BASE_H

#include <CGAL/basic.h>
#include <CGAL/NewKernel_d/Cartesian_complete.h>
#include <CGAL/NewKernel_d/Cartesian_LA_base.h>

namespace CGAL {
#define CGAL_BASE \
  Cartesian_LA_base_d< FT_, Dim_ >
template < typename FT_, typename Dim_, typename Derived_=Default>
struct Cartesian_base_d : public CGAL_BASE
{
    constexpr Cartesian_base_d(){}
    constexpr Cartesian_base_d(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_BASE_H
