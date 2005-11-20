// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef BOOLEAN_SET_OPERATIONS_H
#define BOOLEAN_SET_OPERATIONS_H

#include <CGAL/Bso_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>

CGAL_BEGIN_NAMESPACE

template <class Polygon>
struct Default_bso_traits
{};


template <class Kernel, class Container>
struct Default_bso_traits<CGAL::Polygon_2<Kernel, Container> >
{
  typedef Bso_segment_traits_2<Kernel,
                               Container,
                               Arr_segment_traits_2<Kernel> >    Traits;
};

template <class Polygon>
struct Default_bso_traits<CGAL::General_polygon_with_holes_2<Polygon> >
{
  typedef typename Default_bso_traits<Polygon>::Traits Traits;
};



template <class Polygon_A, class Polygon_B>
inline bool do_intersect(const Polygon_A& pgn1, const Polygon_B& pgn2)
{
  typedef typename Default_bso_traits<Polygon_A>::Traits    Traits;
  
  Traits tr;
  return (do_intersect(pgn1, pgn2, tr));
}

template <class Traits>
inline bool do_intersect(const typename Traits::Polygon_2& pgn1,
                         const typename Traits::Polygon_2& pgn2,
                  Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  return (gps.do_intersect(pgn2));
}

template <class Traits>
inline bool do_intersect(const typename Traits::Polygon_2& pgn1,
                         const typename Traits::Polygon_with_holes_2& pgn2,
                  Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  return (gps.do_intersect(pgn2));
}

template <class Traits>
inline bool do_intersect(const typename Traits::Polygon_with_holes_2& pgn1,
                         const typename Traits::Polygon_2& pgn2,
                         Traits& tr)
{
  return do_intersect(pgn2, pgn1, tr);
}

template <class Traits>
inline bool do_intersect(const typename Traits::Polygon_with_holes_2& pgn1,
                  const typename Traits::Polygon_with_holes_2& pgn2,
                  Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  return (gps.do_intersect(pgn2));
}

CGAL_END_NAMESPACE

#endif
