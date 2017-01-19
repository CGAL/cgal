// Copyright (c) 2009  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H
#define CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/internal/Triangulation_2_filtered_projection_traits_3.h>

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
