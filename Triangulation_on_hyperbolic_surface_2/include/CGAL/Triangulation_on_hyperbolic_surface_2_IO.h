// Copyright (c) 2024
// INRIA Nancy (France), and Université Gustave Eiffel Marne-la-Vallee (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Loïc Dubois, Marc Pouget, Monique Teillaud

// This file contains the declaration and the implementation of the input/output
// functions for the package Triangulation_on_hyperbolic_surface_2

#ifndef CGAL_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2_IO_H
#define CGAL_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2_IO_H

#include <CGAL/license/Triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/basic.h>
#include <iostream>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////
template<class Traits>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_fundamental_domain_2<Traits>& domain){
  CGAL_precondition(domain.is_valid());
  return domain.to_stream(s);
}

template<class Traits>
std::istream& operator>>(std::istream& s, Hyperbolic_fundamental_domain_2<Traits>& domain){
  return domain.from_stream(s);
}

////////////////////////////////////////////////////////////////////////////////
template<class Traits>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_isometry_2<Traits>& isometry){
  for (int k=0; k<4; k++){
    s << isometry.get_coefficient(k);
  }
  return s;
}

////////////////////////////////////////////////////////////////////////////////
template<class Traits, class Attributes>
std::ostream& operator<<(std::ostream& s, const Triangulation_on_hyperbolic_surface_2<Traits, Attributes>& triangulation){
  triangulation.to_stream(s);
  return s;
}

template<class Traits, class Attributes>
void operator>>(std::istream& s, Triangulation_on_hyperbolic_surface_2<Traits, Attributes>& triangulation){
  triangulation.from_stream(s);
}

} // namespace CGAL

#endif // CGAL_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2_IO_H
