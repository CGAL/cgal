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

#ifndef ARR_ADDITION_INSERT_H
#define ARR_ADDITION_INSERT_H

#include <CGAL/Sweep_line_2/Sweep_line_2_impl.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_utils.h>
#include <CGAL/Sweep_line_2/Arr_sweep_line_event.h>
#include <CGAL/Sweep_line_2/Arr_sweep_line_curve.h>
#include <CGAL/Sweep_line_2/Arr_sweep_line_visitor.h>
#include <CGAL/Sweep_line_2/Arr_agg_addition_visitor.h>
#include <CGAL/Sweep_line_2/Arr_meta_traits.h>

#include <CGAL/assertions.h>
#include <list>
#include <vector>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template <class Arr>
class Arr_addition_insert 
{
  typedef typename Arr::Halfedge_handle                    Halfedge_handle;
  typedef typename Arr::Edge_iterator                      Edge_iterator;
  typedef typename Arr::Face_handle                        Face_handle;
  typedef typename Arr::Traits_2                           Base_traits;
  typedef Arr_meta_traits<Base_traits, Halfedge_handle>    Traits;
  typedef Arr_sweep_line_curve<Traits>                     Subcurve; 
  typedef Arr_sweep_line_event<Traits,
                               Subcurve,
                               Halfedge_handle>            Event;
  
  typedef typename Traits::Curve_2                         Curve_2;
  typedef typename Traits::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Traits::Point_2                         Point_2;
  typedef Arr_agg_addition_visitor <Traits,
                                    Arr,
                                    Event,
                                    Subcurve>              Visitor;

  
 
  typedef Sweep_line_2<Traits,
                       Visitor,
                       Subcurve,
                       Event>                              Sweep_line;
 


public:

  Arr_addition_insert(Arr *arr):
      m_traits(new Traits(*(arr->get_traits()))),
      m_arr(arr),
      m_visitor(arr),
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

    m_sweep_line.sweep(xcurves_vec.begin(),
                       xcurves_vec.end(),
                       iso_points.begin(),
                       iso_points.end());
  }

  /*template<class XCurveInputIterator>
  void insert_x_curves(XCurveInputIterator begin,
                       XCurveInputIterator end)
  {
    std::vector<X_monotone_curve_2>      xcurves_vec;
    for (Edge_iterator eit = m_arr->edges_begin(); eit != m_arr->edges_end(); ++eit) 
    {
      xcurves_vec.push_back(eit->curve());
    }
    m_arr->clear();
    std::copy(begin, end, std::back_inserter(xcurves_vec));
    m_sweep_line.init_x_curves(xcurves_vec.begin(),xcurves_vec.end());
    m_sweep_line.sweep();
  }*/

 

  ~Arr_addition_insert()
  {
    delete m_traits;
  }


              
protected:

  Traits*              m_traits;
  Arr*                 m_arr;
  Visitor              m_visitor;
  Sweep_line           m_sweep_line;
};

CGAL_END_NAMESPACE

#endif
