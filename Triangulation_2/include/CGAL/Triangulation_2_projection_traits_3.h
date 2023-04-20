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

#define CGAL_DEPRECATED_HEADER "<CGAL/Triangulation_2_projection_traits_3.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Projection_traits_3.h>"
#include <CGAL/Installation/internal/deprecation_warning.h>

#ifndef CGAL_NO_DEPRECATED_CODE

#include <CGAL/Kernel_23/internal/Filtered_projection_traits_3.h>

namespace CGAL {

template < class Filtered_kernel >
class Filtered_projection_traits_3;

template <class Kernel, bool Has_filtered_predicates = Kernel::Has_filtered_predicates>
class Triangulation_2_projection_traits_3
  : public Projection_traits_base_3<Kernel>
{
public:
  explicit Triangulation_2_projection_traits_3(const typename Kernel::Vector_3& n_)
    : Projection_traits_base_3<Kernel>(n_)
  {}
};

template <class Kernel>
class Triangulation_2_projection_traits_3<Kernel, true>
  : public Filtered_projection_traits_3<Kernel>
{
public:
  explicit Triangulation_2_projection_traits_3(const typename Kernel::Vector_3& n_)
    : Filtered_projection_traits_3<Kernel>(n_)
  {}
};

} // end namespace CGAL

#endif //CGAL_NO_DEPRECATED_CODE

#endif // CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H
