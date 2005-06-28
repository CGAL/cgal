// Copyright (c) 1997  Tel-Aviv University (Israel).
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

#ifndef ARR_OVERLAY_H
#define ARR_OVERLAY_H

#include <CGAL/Arrangement_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_impl.h>
#include <CGAL/Sweep_line_2/Arr_sweep_line_event.h>
#include <CGAL/Arr_overlay_2/Overlay_subcurve.h>
#include <CGAL/Arr_overlay_2/Overlay_visitor.h>
#include <CGAL/Arr_overlay_2/Overlay_meta_traits.h>

#include <vector>

CGAL_BEGIN_NAMESPACE


/*! */
template < class Traits_,
           class Dcel1,
           class Dcel2,
           class ResDcel,
           class OverlayTraits >

void overlay (const Arrangement_2<Traits_, Dcel1>   & arr1,
              const Arrangement_2<Traits_, Dcel2>   & arr2,
              Arrangement_2<Traits_, ResDcel>       & res,
              OverlayTraits & traits)
{
  typedef Traits_                                     Base_Traits;
  typedef typename Base_Traits::X_monotone_curve_2    Base_X_monotone_curve_2;

  typedef Arrangement_2<Traits_, Dcel1>               Arrangement1;
  typedef Arrangement_2<Traits_, Dcel2>               Arrangement2;
  typedef Arrangement_2<Traits_, ResDcel>             Res_Arrangement;

  typedef typename Arrangement1::Halfedge_const_handle       
                                                      Halfedge_const_handle_1;
  typedef typename Arrangement1::Edge_const_iterator         
                                                      Edge_const_iterator_1;
  typedef typename Arrangement1::Face_const_handle
                                                      Face_handle1;


  typedef typename Arrangement2::Halfedge_const_handle       
                                                      Halfedge_const_handle_2;
  typedef typename Arrangement2::Edge_const_iterator  
                                                      Edge_const_iterator_2;
  typedef typename Arrangement2::Face_const_handle    
                                                      Face_handle2;

  typedef typename Res_Arrangement::Halfedge_handle   Halfedge_handle_res;
  typedef typename Res_Arrangement::Face_handle       Res_Face_handle;

  typedef Overlay_meta_traits<Base_Traits,
                              Arrangement1,
                              Arrangement2>            Traits;

  typedef typename Traits::X_monotone_curve_2          X_monotone_curve_2;

  typedef Overlay_subcurve<Traits>  Subcurve;
  typedef Arr_sweep_line_event<Traits,
                               Subcurve,
                               Halfedge_handle_res >   Event;
  typedef Overlay_visitor<Traits,
                          Arrangement1,
                          Arrangement2,
                          Res_Arrangement,
                          Event,
                          Subcurve,
                          OverlayTraits>                 Visitor;

  typedef Sweep_line_2<Traits,
                            Event,
                            Subcurve,
                            Visitor>                     Sweep_line;


                            

  std::vector<X_monotone_curve_2>   arr_curves;
  arr_curves.resize(arr1.number_of_edges() + arr2.number_of_edges());
  
  typename Base_Traits::Compare_xy_2    comp_xy =
    arr1.get_traits()->compare_xy_2_object();
    
  //iterate over arr1's edges and create X_monotone_curve_2 from each edge
  unsigned int i = 0;
  for(Edge_const_iterator_1 itr1 = arr1.edges_begin();
      itr1 != arr1.edges_end();
      ++itr1, ++i)
  {
    Halfedge_const_handle_1 he = *itr1;

    // Associate each x-monotone curve with the halfedge that represent it
    // that is directed from right to left.
    if(comp_xy(he.source().point(),
	             he.target().point()) == SMALLER)
       he = he.twin();
    const Base_X_monotone_curve_2& base_cv = he.curve();

    arr_curves[i] =
      X_monotone_curve_2(base_cv, he, Halfedge_const_handle_2());
  }

  //iterate over arr2's edges and create X_monotone_curve_2 from each edge
  for(Edge_const_iterator_2 itr2 = arr2.edges_begin();
      itr2 != arr2.edges_end();
      ++itr2, ++i)
  {
    Halfedge_const_handle_2 he = *itr2;

    // Associate each x-monotone curve with the halfedge that represent it
    // that is directed from right to left.
    if(comp_xy(he.source().point(),
	             he.target().point()) == SMALLER)
       he = he.twin();

    const Base_X_monotone_curve_2& base_cv = he.curve();

    arr_curves[i] =
      X_monotone_curve_2(base_cv, Halfedge_const_handle_1(), he);
  }

  Visitor visitor(arr1, arr2, res, traits);
  Sweep_line sweep_object(&visitor);
  sweep_object.sweep(arr_curves.begin(), arr_curves.end());
}



CGAL_END_NAMESPACE

#endif
