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

#ifndef CGAL_ARR_AGGREGATE_INSERT_H
#define CGAL_ARR_AGGREGATE_INSERT_H

#include <CGAL/Sweep_line_2/Sweep_line_2_impl.h>
#include <CGAL/Sweep_line_2/Arr_sweep_line_event.h>
#include <CGAL/Sweep_line_2/Arr_sweep_line_curve.h>
#include <CGAL/Sweep_line_2/Arr_sweep_line_visitor.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_utils.h>
#include <CGAL/assertions.h>
#include <list>
#include <vector>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template <class Arr>
class Arr_aggregate_insert 
{
  typedef typename Arr::Halfedge_handle                      Halfedge_handle;
  typedef typename Arr::Edge_iterator                        Edge_iterator;
  typedef typename Arr::Traits_2                             Traits;
  typedef Arr_sweep_line_curve<Traits>                       Subcurve; 
  typedef Arr_sweep_line_event<Traits,
                               Subcurve,
                               Halfedge_handle>              Event;
  typedef typename Traits::X_monotone_curve_2                X_monotone_curve_2;
  typedef typename Traits::Point_2                           Point_2;
  typedef Arr_sweep_line_visitor<Traits,
                                 Arr,
                                 Event,
                                 Subcurve>                   Visitor;

  
 
  typedef Sweep_line_2<Traits,
                            Event,
                            Subcurve,
                            Visitor,
                            CGAL_ALLOCATOR(int)>       Sweep_line;
 


public:

  Arr_aggregate_insert(Arr *arr):
      m_traits(new Traits()),
      m_traits_owner(true),
      m_arr(arr),
      m_visitor(arr),
      m_sweep_line(m_traits, &m_visitor)
      {}


  Arr_aggregate_insert(const Traits *traits, Arr *arr):
      m_traits(traits),
      m_traits_owner(false),
      m_arr(arr),
      m_visitor(arr),
      m_sweep_line(m_traits, &m_visitor)
      {}

  template<class CurveInputIterator>
  void insert_curves(CurveInputIterator begin, 
                     CurveInputIterator end)
  {
    std::vector<X_monotone_curve_2>      curves_vec;
    curves_vec.reserve(std::distance(begin, end) + m_arr->number_of_edges());

    std::vector<Point_2> points_vec;
    make_x_monotone(begin,
                    end,
                    std::back_inserter(curves_vec),
                    std::back_inserter(points_vec),
                    m_traits);
   
    for (Edge_iterator eit = m_arr->edges_begin();
         eit != m_arr->edges_end();
         ++eit) 
    {
      curves_vec.push_back(eit->curve());
    }

    m_arr->clear();

    //Perform the sweep
    m_sweep_line.sweep(curves_vec.begin(),
                       curves_vec.end(),
                       points_vec.begin(),
                       points_vec.end());
   
  }

  template<class XCurveInputIterator>
  void insert_x_curves(XCurveInputIterator begin,
                       XCurveInputIterator end)
  {
    std::vector<X_monotone_curve_2>      curves_vec;
    curves_vec.reserve(std::distance(begin, end) + m_arr->number_of_edges());
    std::copy(begin, end, std::back_inserter(curves_vec));

    for (Edge_iterator eit = m_arr->edges_begin();
         eit != m_arr->edges_end();
         ++eit) 
    {
      curves_vec.push_back(eit->curve());
    }

    m_arr->clear();
    m_sweep_line.init(curves_vec.begin(), curves_vec.end());
    m_sweep_line.sweep();
  }

 

  ~Arr_aggregate_insert()
  {
    if(m_traits_owner)
      delete m_traits;
  }


              
protected:

  const Traits*        m_traits;
  bool                 m_traits_owner;
  Arr*                 m_arr;
  Visitor              m_visitor;
  Sweep_line           m_sweep_line;
};

CGAL_END_NAMESPACE

#endif
