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
// SPDX-License-Identifier: LGPL-3.0+
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
    CGAL_CONSTEXPR Cartesian_base_d(){}
    CGAL_CONSTEXPR Cartesian_base_d(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_BASE_H
