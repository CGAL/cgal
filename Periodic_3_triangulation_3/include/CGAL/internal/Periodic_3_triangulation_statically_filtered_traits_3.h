// Copyright (c) 2001,2004,2008-2009   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_TRIANGULATION_STATICALLY_FILTERED_TRAITS_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_STATICALLY_FILTERED_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/internal/Static_filters/Periodic_3_orientation_3.h>
#include <CGAL/internal/Periodic_3_triangulation_filtered_traits_3.h>

namespace CGAL {

// The `Traits` argument is supposed to provide exact primitives.
template < typename Traits >
class Periodic_3_triangulation_statically_filtered_traits_3
    : public Traits
{
  typedef Periodic_3_triangulation_statically_filtered_traits_3<Traits> Self;

public:
  typedef internal::Static_filters_predicates::Periodic_3_orientation_3<Traits>
    Orientation_3;

  Orientation_3 orientation_3_object() const {
    return Orientation_3(&this->_domain,&this->_domain_e,&this->_domain_f);
  }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_STATICALLY_FILTERED_TRAITS_3_H
