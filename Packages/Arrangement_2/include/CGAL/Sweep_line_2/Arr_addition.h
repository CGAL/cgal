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

#ifndef ARR_ADDITION_H
#define ARR_ADDITION_H

#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_utils.h>
#include <CGAL/Sweep_line_2/Arr_construction_event.h>
#include <CGAL/Sweep_line_2/Arr_construction_curve.h>
#include <CGAL/Sweep_line_2/Arr_addition_visitor.h>
#include <CGAL/Sweep_line_2/Arr_addition_traits.h>

#include <CGAL/assertions.h>
#include <list>
#include <vector>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template <class Arrangement_>
class Arr_addition 
{
  typedef Arrangement_                                     Arrangement;
  typedef typename Arrangement::Halfedge_handle            Halfedge_handle;
  typedef typename Arrangement::Edge_iterator              Edge_iterator;
  typedef typename Arrangement::Vertex_iterator            Vertex_iterator;
  typedef typename Arrangement::Face_handle                Face_handle;
  typedef typename Arrangement::Traits_2                   Base_traits;

  typedef Arr_addition_traits<Base_traits, Arrangement>    Traits;
  typedef Arr_construction_curve<Traits>                   Subcurve; 
  typedef Arr_construction_event<Traits,
                                 Subcurve,
                                 Halfedge_handle>          Event;
  
  typedef typename Traits::Curve_2                         Curve_2;
  typedef typename Traits::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Traits::Point_2                         Point_2;
  typedef Arr_addition_visitor <Traits,
                                Arrangement,
                                Event,
                                Subcurve>                  Visitor;

  
 
  typedef Sweep_line_2<Traits,
                       Visitor,
                       Subcurve,
                       Event>                              Sweep_line;
 


public:

  Arr_addition(Arrangement &arr):
      m_traits(new Traits(*(arr.get_traits()))),
      m_arr(&arr),
      m_visitor(&arr),
      m_sweep_line(m_traits, &m_visitor)
      {}



  template<class CurveInputIterator>
  void insert_curves(CurveInputIterator begin, 
                     CurveInputIterator end)
  {
    std::vector<X_monotone_curve_2>      xcurves_vec;
    std::vector<Point_2>                 iso_points;

    make_x_monotone (begin,
		     end,
		     std::back_inserter(xcurves_vec),
		     std::back_inserter(iso_points),
		     m_traits);

    typename Traits::Compare_xy_2  comp_xy = m_traits->compare_xy_2_object();
    for (Edge_iterator eit = m_arr->edges_begin();
         eit != m_arr->edges_end();
         ++eit) 
    {
      Halfedge_handle he;
      if(comp_xy(eit->source()->point(),
	       eit->target()->point()) == SMALLER)
         he = eit->handle()->twin();
      else
        he = eit->handle();

      xcurves_vec.push_back(X_monotone_curve_2(he->curve(), he));
    }

    for(Vertex_iterator v_itr = m_arr->vertices_begin();
        v_itr != m_arr->vertices_end();
        ++v_itr)
    {
      if(v_itr->is_isolated())
        iso_points.push_back(Point_2(v_itr->point(), v_itr->handle()));
    }


    m_sweep_line.sweep(xcurves_vec.begin(),
                       xcurves_vec.end(),
                       iso_points.begin(),
                       iso_points.end());
  }



  template<class XCurveInputIterator>
  void insert_x_curves(XCurveInputIterator begin,
                       XCurveInputIterator end)
  {
    std::vector<X_monotone_curve_2>      xcurves_vec;
    std::vector<Point_2>                 iso_points;

    typename Traits::Compare_xy_2  comp_xy = m_traits->compare_xy_2_object();
    for (Edge_iterator eit = m_arr->edges_begin();
         eit != m_arr->edges_end();
         ++eit) 
    {
      Halfedge_handle he;
      if(comp_xy(eit->source()->point(),
	       eit->target()->point()) == SMALLER)
         he = eit->handle()->twin();
      else
        he = eit->handle();

      xcurves_vec.push_back(X_monotone_curve_2(he->curve(), he));
    }

    for(Vertex_iterator v_itr = m_arr->vertices_begin();
        v_itr != m_arr->vertices_end();
        ++v_itr)
    {
      if(v_itr->is_isolated())
        iso_points.push_back(Point_2(v_itr->point(), v_itr->handle()));
    }
    std::copy(begin, end, std::back_inserter(xcurves_vec));
    m_sweep_line.sweep(xcurves_vec.begin(),
                       xcurves_vec.end(),
                       iso_points.begin(),
                       iso_points.end());
  }

 

  ~Arr_addition()
  {
    delete m_traits;
  }


              
protected:

  Traits*              m_traits;
  Arrangement*         m_arr;
  Visitor              m_visitor;
  Sweep_line           m_sweep_line;
};

CGAL_END_NAMESPACE

#endif
