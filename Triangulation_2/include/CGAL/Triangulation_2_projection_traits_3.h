// Copyright (c) 2009  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H
#define CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/Triangulation_2/internal/Triangulation_2_filtered_projection_traits_3.h>

namespace CGAL{

// This declaration is needed to break the cyclic dependency.
template < class Filtered_kernel >
class Triangulation_2_filtered_projection_traits_3;

template <class Kernel, bool Has_filtered_predicates=Kernel::Has_filtered_predicates>
class Triangulation_2_projection_traits_3
  : public Triangulation_2_projection_traits_base_3<Kernel>
{
public:
  explicit
  Triangulation_2_projection_traits_3(const typename Kernel::Vector_3& n_)
    : Triangulation_2_projection_traits_base_3<Kernel>(n_)
  {}
};

template <class Kernel>
class Triangulation_2_projection_traits_3<Kernel, true>
  : public Triangulation_2_filtered_projection_traits_3<Kernel>
{
public:
  explicit
  Triangulation_2_projection_traits_3(const typename Kernel::Vector_3& n_)
    : Triangulation_2_filtered_projection_traits_3<Kernel>(n_)
  {}
};

} // end namespace CGAL

#endif // CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H
