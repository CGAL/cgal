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
#include <CGAL/Root_of_traits.h>

namespace CGAL {

  /*
  Marc: Convert cosh_hd to an object fct?
  remove "static" and construct the traits each time it is used?
  (may follow the kernel cgal style CGAL_DISTANCE_2_POINT_2_POINT_2_H that uses inline?? it may not be an object function but a global kernel function seems ok since there are no construction)
 it is used in CGAL_HYPERBOLIC_FUNDAMENTAL_DOMAIN_2_H
  */

  //TODO DOC add cosh of hyperbolic distance
template<class HyperbolicTraitsClass>
class Hyperbolic_surface_traits_2
  : public HyperbolicTraitsClass
{
public:
  typedef typename HyperbolicTraitsClass::FT                          FT;
  typedef typename HyperbolicTraitsClass::Hyperbolic_point_2          Hyperbolic_point_2;
  typedef Complex_number<FT>                                          Complex;

  /*
    Note: In the HyperbolicDelaunayTriangulationTraits_2 Concept, there is no
    computation of distances since all computations are about
    circumcircles, that are also Euclidean circles, and only incircle or
    orientation tests are used.
  */
template <typename Number = FT, typename Point = Hyperbolic_point_2>
static Number cosh_hd(Point const& u, Hyperbolic_point_2 const& v)
{
      Number num = (u.x() - v.x()) * (u.x() - v.x()) + (u.y() - v.y()) * (u.y() - v.y());
      Number den = (1 - (u.x() * u.x() + u.y() * u.y())) * (1 - (v.x() * v.x() + v.y() * v.y()));
      return 2 * num / den + 1;
}

};

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_SURFACE_TRAITS_2
