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

#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Gps_traits_adaptor_2.h>
#include <CGAL/iterator.h> 

CGAL_BEGIN_NAMESPACE

template <class Polygon>
struct Gps_default_traits
{};


template <class Kernel, class Container>
struct Gps_default_traits<CGAL::Polygon_2<Kernel, Container> >
{
  typedef Gps_segment_traits_2<Kernel,
                               Container,
                               Arr_segment_traits_2<Kernel> >    Traits;
};

template <class Polygon>
struct Gps_default_traits<CGAL::General_polygon_with_holes_2<Polygon> >
{
  typedef typename Gps_default_traits<Polygon>::Traits Traits;
};

template <class Arr_traits>
struct Gps_default_traits<CGAL::General_polygon_2<Arr_traits> >
{
  typedef Gps_traits_adaptor_2<Arr_traits>    Traits;
};


/////////////////////
//  do_intersect  //
///////////////////
template <class Polygon_A, class Polygon_B>
inline bool do_intersect(const Polygon_A& pgn1, const Polygon_B& pgn2)
{
  typedef typename Gps_default_traits<Polygon_A>::Traits    Traits;
  
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

/////////////////////
//  intersection  //
///////////////////

template <class Polygon_A, class Polygon_B, class OutputIterator>
inline OutputIterator intersection(const Polygon_A& pgn1,
                                   const Polygon_B& pgn2,
                                   OutputIterator oi)
{
  typedef typename Gps_default_traits<Polygon_A>::Traits    Traits;
  
  Traits tr;
  return (intersection(pgn1, pgn2, oi, tr));
}

template <class Traits, class OutputIterator>
inline OutputIterator intersection(const typename Traits::Polygon_2& pgn1,
                                   const typename Traits::Polygon_2& pgn2,
                                   OutputIterator oi,
                                   Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.intersection(pgn2);
  return (gps.polygons_with_holes(oi));
}

template <class Traits, class OutputIterator>
inline OutputIterator intersection(const typename Traits::Polygon_2& pgn1,
                         const typename Traits::Polygon_with_holes_2& pgn2,
                         OutputIterator oi,
                         Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.intersection(pgn2);
  return (gps.polygons_with_holes(oi));
}

template <class Traits, class OutputIterator>
inline OutputIterator intersection(const typename Traits::Polygon_with_holes_2& pgn1,
                         const typename Traits::Polygon_2& pgn2,
                         OutputIterator oi,
                         Traits& tr)
{
  return intersection(pgn2, pgn1, oi, tr);
}

template <class Traits, class OutputIterator>
inline OutputIterator intersection(const typename Traits::Polygon_with_holes_2& pgn1,
                         const typename Traits::Polygon_with_holes_2& pgn2,
                         OutputIterator oi,
                         Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.intersection(pgn2);
  return (gps.polygons_with_holes(oi));
}


/////////////////////
//      join      //
///////////////////

template <class Polygon_A, class Polygon_B, class OutputPolygon>
inline bool join(const Polygon_A& pgn1,
                 const Polygon_B& pgn2,
                 OutputPolygon&   res)
{
  typedef typename Gps_default_traits<Polygon_A>::Traits    Traits;
  
  Traits tr;
  return (join(pgn1, pgn2, res, tr));
}

template <class Traits>
inline bool join(const typename Traits::Polygon_2& pgn1,
                 const typename Traits::Polygon_2& pgn2,
                 typename Traits::Polygon_with_holes_2&   res,
                 Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.join(pgn2);
  if(gps.number_of_polygons_with_holes() == 1)
  {
    Oneset_iterator<typename Traits::Polygon_with_holes_2> oi (res);
    gps.polygons_with_holes(oi);
    return true;
  }

  // the polygon doesnt intersect, the original pgn1, pgn2 contain the union
  return false;
}

template <class Traits>
inline bool join(const typename Traits::Polygon_2& pgn1,
                 const typename Traits::Polygon_with_holes_2& pgn2,
                 typename Traits::Polygon_with_holes_2&   res,
                 Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.join(pgn2);
  if(gps.number_of_polygons_with_holes() == 1)
  {
    Oneset_iterator<typename Traits::Polygon_with_holes_2> oi (res);
    gps.polygons_with_holes(oi);
    return true;
  }

  // the polygon doesnt intersect, the original pgn1, pgn2 contain the union
  return false;
}

template <class Traits>
inline bool join(const typename Traits::Polygon_with_holes_2& pgn1,
                 const typename Traits::Polygon_2& pgn2,
                 typename Traits::Polygon_with_holes_2&   res,
                 Traits& tr)
{
  return join(pgn2, pgn1, res, tr);
}

template <class Traits>
inline bool join(const typename Traits::Polygon_with_holes_2& pgn1,
                 const typename Traits::Polygon_with_holes_2& pgn2,
                 typename Traits::Polygon_with_holes_2&   res,
                 Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.join(pgn2);
  if(gps.number_of_polygons_with_holes() == 1)
  {
    Oneset_iterator<typename Traits::Polygon_with_holes_2> oi (res);
    gps.polygons_with_holes(oi);
    return true;
  }

  // the polygon doesnt intersect, the original pgn1, pgn2 contain the union
  return false;
}


/////////////////////
//  difference    //
///////////////////

template <class Polygon_A, class Polygon_B, class OutputIterator>
inline OutputIterator difference(const Polygon_A& pgn1,
                                 const Polygon_B& pgn2,
                                 OutputIterator oi)
{
  typedef typename Gps_default_traits<Polygon_A>::Traits    Traits;
  
  Traits tr;
  return (difference(pgn1, pgn2, oi, tr));
}

template <class Traits, class OutputIterator>
inline OutputIterator difference(const typename Traits::Polygon_2& pgn1,
                                 const typename Traits::Polygon_2& pgn2,
                                 OutputIterator oi,
                                 Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.difference(pgn2);
  return (gps.polygons_with_holes(oi));
}

template <class Traits, class OutputIterator>
inline OutputIterator difference(const typename Traits::Polygon_2& pgn1,
                                 const typename Traits::Polygon_with_holes_2& pgn2,
                                 OutputIterator oi,
                                 Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.difference(pgn2);
  return (gps.polygons_with_holes(oi));
}

template <class Traits, class OutputIterator>
inline OutputIterator difference(const typename Traits::Polygon_with_holes_2& pgn1,
                                 const typename Traits::Polygon_2& pgn2,
                                 OutputIterator oi,
                                 Traits& tr)
{
  return difference(pgn2, pgn1, oi, tr);
}

template <class Traits, class OutputIterator>
inline OutputIterator difference(const typename Traits::Polygon_with_holes_2& pgn1,
                                 const typename Traits::Polygon_with_holes_2& pgn2,
                                 OutputIterator oi,
                                 Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.difference(pgn2);
  return (gps.polygons_with_holes(oi));
}


/////////////////////////////
//  symmetric_difference  //
///////////////////////////

template <class Polygon_A, class Polygon_B, class OutputIterator>
inline OutputIterator symmetric_difference(const Polygon_A& pgn1,
                                           const Polygon_B& pgn2,
                                           OutputIterator oi)
{
  typedef typename Gps_default_traits<Polygon_A>::Traits    Traits;
  
  Traits tr;
  return (symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Traits, class OutputIterator>
inline OutputIterator symmetric_difference(const typename Traits::Polygon_2& pgn1,
                                           const typename Traits::Polygon_2& pgn2,
                                           OutputIterator oi,
                                           Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.symmetric_difference(pgn2);
  return (gps.polygons_with_holes(oi));
}

template <class Traits, class OutputIterator>
inline OutputIterator symmetric_difference(const typename Traits::Polygon_2& pgn1,
                                           const typename Traits::Polygon_with_holes_2& pgn2,
                                           OutputIterator oi,
                                           Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.symmetric_difference(pgn2);
  return (gps.polygons_with_holes(oi));
}

template <class Traits, class OutputIterator>
inline OutputIterator symmetric_difference(const typename Traits::Polygon_with_holes_2& pgn1,
                                           const typename Traits::Polygon_2& pgn2,
                                           OutputIterator oi,
                                           Traits& tr)
{
  return symmetric_difference(pgn2, pgn1, oi, tr);
}

template <class Traits, class OutputIterator>
inline OutputIterator symmetric_difference(const typename Traits::Polygon_with_holes_2& pgn1,
                                           const typename Traits::Polygon_with_holes_2& pgn2,
                                           OutputIterator oi,
                                           Traits& tr)
{
  General_polygon_set_2<Traits> gps(pgn1);
  gps.symmetric_difference(pgn2);
  return (gps.polygons_with_holes(oi));
}

CGAL_END_NAMESPACE

#endif
