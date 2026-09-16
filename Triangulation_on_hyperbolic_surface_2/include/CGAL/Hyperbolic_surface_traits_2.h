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

#ifndef CGAL_HYPERBOLIC_SURFACE_TRAITS_2
#define CGAL_HYPERBOLIC_SURFACE_TRAITS_2

#include <CGAL/license/Triangulation_on_hyperbolic_surface_2.h>

#include <CGAL/Complex_number.h>

namespace CGAL {

template<class HyperbolicTraitsClass>
class Hyperbolic_surface_traits_2
  : public HyperbolicTraitsClass
{
public:
  typedef typename HyperbolicTraitsClass::FT                          FT;
  typedef typename HyperbolicTraitsClass::Hyperbolic_point_2          Hyperbolic_point_2;
  typedef Complex_number<FT>                                          Complex;
};

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_SURFACE_TRAITS_2
