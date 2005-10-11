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

#ifndef ARR_OVERLAY_H
#define ARR_OVERLAY_H

#include <CGAL/Arrangement_2.h>
#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Arr_construction_event.h>
#include <CGAL/Arr_overlay_2/Overlay_subcurve.h>
#include <CGAL/Arr_overlay_2/Overlay_visitor.h>
#include <CGAL/Arr_overlay_2/Overlay_meta_traits.h>
#include <CGAL/Object.h>

#include <vector>

CGAL_BEGIN_NAMESPACE

/*!
 * Compute the overlay of two input arrangements.
 * \param arr1 The first arrangement.
 * \param arr2 The second arrangement.
 * \param res Output: The resulting arrangement.
 * \param traits An overlay-traits class. As arr1, arr2 and res are all
 *               templated with the same arrangement-traits class but with
 *               different DCELs, the overlay-traits class defines the
 *               various overlay operations of pairs of DCEL features from
 *               Dcel1 and Dcel2 to the resulting ResDcel.
 */
template <class Traits_,
          class Dcel1,
          class Dcel2,
          class ResDcel,
          class OverlayTraits >
void overlay (const Arrangement_2<Traits_, Dcel1>& arr1,
              const Arrangement_2<Traits_, Dcel2>& arr2,
              Arrangement_2<Traits_, ResDcel>& res,
              OverlayTraits& traits)
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

  // in order to iterate over isolated vertices
  typedef typename Arrangement1::Vertex_const_iterator
                                                      Vertex_const_iterator_1;

  typedef typename Arrangement2::Vertex_const_iterator
                                                      Vertex_const_iterator_2;

  typedef Overlay_meta_traits<Base_Traits,
                              Arrangement1,
                              Arrangement2>            Traits;

  typedef typename Traits::X_monotone_curve_2          X_monotone_curve_2;
  typedef typename Traits::Point_2                     Point_2;

  typedef Overlay_subcurve<Traits>  Subcurve;
  typedef Arr_construction_event<Traits,
                                 Subcurve,
                                 Halfedge_handle_res > Event;
  typedef Overlay_visitor<Traits,
                          Arrangement1,
                          Arrangement2,
                          Res_Arrangement,
                          Event,
                          Subcurve,
                          OverlayTraits>               Visitor;

  typedef Sweep_line_2<Traits,
                       Visitor,
                       Subcurve,
                       Event>                          Sweep_line;


    
  //the res arrangement can't be the same one as one of the input arrangement
  CGAL_precondition((&res != &arr1) && (&res != &arr2));

  std::vector<X_monotone_curve_2>   arr_curves;
  std::list<Point_2>                iso_points;
  arr_curves.resize(arr1.number_of_edges() + arr2.number_of_edges());
  
  //iterate over arr1's edges and create X_monotone_curve_2 from each edge
  unsigned int            i = 0;
  Edge_const_iterator_1   itr1;
  Halfedge_const_handle_1 he1;

  for (itr1 = arr1.edges_begin(); itr1 != arr1.edges_end(); ++itr1, ++i)
  {
    he1 = itr1;

    // Associate each x-monotone curve with the halfedge that represent it
    // that is directed from right to left.
    if (he1->direction() == SMALLER)
      he1 = he1->twin();

    const Base_X_monotone_curve_2& base_cv = he1->curve();

    arr_curves[i] =
      X_monotone_curve_2 (base_cv, he1, Halfedge_const_handle_2());
  }

  //iterate over arr2's edges and create X_monotone_curve_2 from each edge
  Edge_const_iterator_2   itr2;
  Halfedge_const_handle_2 he2;

  for (itr2 = arr2.edges_begin(); itr2 != arr2.edges_end(); ++itr2, ++i)
  {
    he2 = itr2;

    // Associate each x-monotone curve with the halfedge that represent it
    // that is directed from right to left.
    if (he2->direction() == SMALLER)
       he2 = he2->twin();

    const Base_X_monotone_curve_2& base_cv = he2->curve();

    arr_curves[i] =
      X_monotone_curve_2 (base_cv, Halfedge_const_handle_1(), he2);
  }

  // iterate over arr1's vertices and associate each isolated point with
  // its vertex
  Vertex_const_iterator_1 v_itr1;

  for(v_itr1 = arr1.vertices_begin();
      v_itr1 != arr1.vertices_end();
      ++v_itr1)
  {
    if(v_itr1->is_isolated())
    {
      iso_points.push_back(Point_2(v_itr1->point(),
                           CGAL::make_object(v_itr1),
                           Object()));
    }
  }
      
  // iterate over arr2's vertices and associate each isolated point with
  // its vertex
  Vertex_const_iterator_2 v_itr2;

  for (v_itr2 = arr2.vertices_begin();
       v_itr2 != arr2.vertices_end();
       ++v_itr2)
  {
    if (v_itr2->is_isolated())
    {
      iso_points.push_back(Point_2(v_itr2->point(),
                                   Object(),
                                   CGAL::make_object(v_itr2)));
    }
  }

  // clear the result arrangement
  res.clear();

  //perform the sweep whith overlay visitor
  Visitor    visitor(arr1, arr2, res, traits);
  Traits     meta_traits(res.get_traits());
  Sweep_line sweep_object(&meta_traits, &visitor);
  sweep_object.sweep(arr_curves.begin(),
                     arr_curves.end(),
                     iso_points.begin(),
                     iso_points.end());
}



CGAL_END_NAMESPACE

#endif
