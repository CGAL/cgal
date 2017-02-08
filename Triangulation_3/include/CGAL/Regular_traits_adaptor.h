// Copyright (c) 1999-2004,2006-2009,2014-2015   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Aymeric Pellé <Aymeric.Pelle@sophia.inria.fr>
#ifndef CGAL_REGULAR_TRAITS_ADAPTOR_H
#define CGAL_REGULAR_TRAITS_ADAPTOR_H

#include <CGAL/license/Triangulation_3.h>

namespace CGAL {

  template < class RTT, class ConstructPoint, class Functor_>
class Regular_traits_adaptor
{
  const ConstructPoint& cp;
  const Functor_& f;

  typedef RTT                                              RTraits;
  typedef Functor_                                         Functor;

  typedef typename RTraits::FT                             FT;
  typedef typename RTraits::Point_3                        Point_3;
  typedef typename RTraits::Weighted_point_3               Weighted_point_3;

public:
  typedef typename Functor::result_type                     result_type;

public:
  Regular_traits_adaptor (const ConstructPoint& cp, const Functor& f)
    : cp(cp), f(f)
  { }
  
public:



public:
  // with offset ---------------------------------------------------------------
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1) const
  {
    return f(cp(p0), cp(p1));
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2) const
  {
    return f(cp(p0), cp(p1), cp(p2));
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3) const
  {
    return f(cp(p0), cp(p1), cp(p2), cp(p3));
  }
  
  result_type operator() (const Point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2) const
  {
    return f(p0, cp(p1), cp(p2));
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const FT w) const
  {
    return f(cp(p0), cp(p1), cp(p2), cp(p3), w);
  }
  
};
 
}  // namespace CGAL

#endif /* CGAL_REGULAR_TRAITS_ADAPTOR_H */
