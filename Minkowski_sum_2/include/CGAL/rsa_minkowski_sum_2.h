// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// $Source: $
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_RSA_MINKOWSKI_SUM_H
#define CGAL_RSA_MINKOWSKI_SUM_H

CGAL_BEGIN_NAMESPACE

/*!
 */
template <class ConicTraits, class Container>
typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
rsa_minkowski_sum_2
    (const ConicTraits& ,
     const Polygon_2<typename ConicTraits::Rat_kernel, Container>& pgn1,
     const typename Gps_circle_segment_traits_2<typename 
                              ConicTraits::Rat_kernel>::Polygon_2& pgn2)
{
  typedef ConicTraits                                   Conic_traits_2;
  typedef typename Conic_traits_2::Rat_kernel           Rat_kernel;
  typedef typename Rat_kernel::FT                       Rational;
  typedef typename Rat_kernel::Point_2                  Point_2;
  typedef typename Rat_kernel::
  typedef Polygon_2<Rat_kernel, Container>              Polygon_2;
  
  typedef typename Conic_traits_2::Alg_kernel           Alg_kernel;
  typedef typename Alg_kernel::FT                       Algebraic;
  typedef typename Alg_kernel::Point_2                  Alg_point_2;
  
  typedef Exact_offset_base_2<ConicTraits, Container>        Base;
  typedef Offset_by_convolution_2<Base>                      Exact_offset_2;
  typedef typename Exact_offset_2::Offset_polygon_2          Offset_polygon_2;

  Base                                               base;
  Exact_offset_2                                     exact_offset (base);
  Offset_polygon_2                                   offset_bound;
  std::list<Offset_polygon_2>                        offset_holes;

  exact_offset (pgn, r, 
                offset_bound, std::back_inserter(offset_holes));

  return (typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
          (offset_bound, offset_holes.begin(), offset_holes.end()));
}

CGAL_END_NAMESPACE

#endif
