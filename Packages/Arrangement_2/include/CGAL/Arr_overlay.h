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


#include <CGAL/Sweep_line_2/Sweep_line_2_impl.h>
#include <CGAL/Sweep_line_2/Arr_sweep_line_event.h>
#include <CGAL/Arr_overlay_2/Overlay_subcurve.h>
#include <CGAL/Arr_overlay_2/Overlay_visitor.h>
#include <CGAL/Arr_overlay_2/Overlay_meta_traits.h>

#include <vector>

CGAL_BEGIN_NAMESPACE


/*! */
template < class Arrangement1,
           class Arrangement2,
           class Res_Arrangement,
           class OverlayTraits >

void arr_overlay (const Arrangement1  & arr1,
                  const Arrangement2  & arr2,
                  Res_Arrangement     & res,
                  const OverlayTraits & traits)
{
  typedef typename Arrangement1::Traits_2                Base_Traits;
  typedef typename Base_Traits::X_monotone_curve_2     Base_X_monotone_curve_2;

  typedef typename Arrangement1::Halfedge_const_handle Halfedge_const_handle_1;
  typedef typename Arrangement1::Edge_const_iterator   Edge_const_iterator_1;

  typedef typename Arrangement2::Halfedge_const_handle Halfedge_const_handle_2;
  typedef typename Arrangement2::Edge_const_iterator   Edge_const_iterator_2;

  typedef typename Res_Arrangement::Halfedge_handle    Halfedge_handle_res;

  typedef Overlay_meta_traits<Base_Traits,
                              Halfedge_const_handle_1,
                              Halfedge_const_handle_2> Traits;

  typedef typename Traits::X_monotone_curve_2          X_monotone_curve_2;

  typedef Overlay_subcurve<Traits, Halfedge_handle_res>  Subcurve;
  typedef Arr_sweep_line_event<Traits, Subcurve>         Event;
  typedef Overlay_visitor<Traits,
                          Arrangement1,
                          Arrangement2,
                          Res_Arrangement,
                          Event,
                          Subcurve,
                          OverlayTraits>                 Visitor;

  typedef Sweep_line_2_impl<Traits,
                            Event,
                            Subcurve,
                            Visitor>                     Sweep_line;


                            

  std::vector<X_monotone_curve_2>   arr1_curves;
  std::vector<X_monotone_curve_2>   arr2_curves;

  arr1_curves.resize(arr1.number_of_edges());
  arr2_curves.resize(arr2.number_of_edges());

    
    
  //iterate over arr1's edges and create X_monotone_curve_2 from each edge
  unsigned int i = 0;
  for(Edge_const_iterator_1 itr1 = arr1.edges_begin();
      itr1 != arr1.edges_end();
      ++itr1, ++i)
  {
    Halfedge_const_handle_1 he = *itr1;
    const Base_X_monotone_curve_2& base_cv = he.curve();

    arr1_curves[i] = X_monotone_curve_2(base_cv, he, Halfedge_const_handle_2());
  }

  //iterate over arr2's edges and create X_monotone_curve_2 from each edge
  i = 0;
  for(Edge_const_iterator_2 itr2 = arr2.edges_begin();
      itr2 != arr2.edges_end();
      ++itr2, ++i)
  {
    Halfedge_const_handle_2 he = *itr2;
    const Base_X_monotone_curve_2& base_cv = he.curve();

    arr2_curves[i] = X_monotone_curve_2(base_cv, Halfedge_const_handle_1(), he);
  }

  Visitor visitor(arr1, arr2, res, traits);
  Sweep_line sweep_object(&visitor);

  sweep_object.init_x_curves(arr1_curves.begin(), arr1_curves.end());
  sweep_object.init_x_curves(arr2_curves.begin(), arr2_curves.end());
  sweep_object.sweep();


}




CGAL_END_NAMESPACE


#endif
