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

#ifndef ARR_NON_X_AGGREGATE_INSERT_H
#define ARR_NON_X_AGGREGATE_INSERT_H

#include <CGAL/Sweep_line_2/Sweep_line_2_impl.h>
#include <CGAL/Sweep_line_2/Arr_sweep_line_event.h>
#include <CGAL/Sweep_line_2/Arr_sweep_line_curve.h>
#include <CGAL/Sweep_line_2/Arr_sweep_line_visitor.h>
#include <CGAL/Sweep_line_2/Arr_non_x_traits.h>
#include <CGAL/assertions.h>
#include <list>
#include <vector>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template <class Arr>
class Arr_non_x_aggregate_insert 
{
  typedef typename Arr::Halfedge_handle                      Halfedge_handle;
  typedef typename Arr::Edge_iterator                        Edge_iterator;
  typedef typename Arr::Traits_2                             Arr_traits;
  typedef Arr_non_x_traits<Arr_traits>                       Traits;
  typedef Arr_sweep_line_curve<Traits, Halfedge_handle>      Subcurve; 
  typedef Arr_sweep_line_event<Traits, Subcurve>             Event;
  typedef typename Traits::X_monotone_curve_2                X_monotone_curve_2;
  typedef Arr_sweep_line_visitor<Traits,
                                 Arr,
                                 Event,
                                 Subcurve>                   Visitor;

  
 
  typedef Sweep_line_2_impl<Traits,
                            Event,
                            Subcurve,
                            Visitor,
                            CGAL_ALLOCATOR(int)>       Sweep_line;
 


public:

  Arr_non_x_aggregate_insert(Arr *arr):
      m_traits(new Traits()),
      m_traits_owner(true),
      m_arr(arr),
      m_visitor(arr),
      m_sweep_line(m_traits, &m_visitor)
      {}


  Arr_non_x_aggregate_insert(const Arr_traits *traits, Arr *arr):
      m_traits(static_cast<const Traits*>(traits)), 
      m_traits_owner(false),
      m_arr(arr),
      m_visitor(arr),
      m_sweep_line(m_traits, &m_visitor)
      {}

  
  template<class XCurveInputIterator>
  void insert_curves(XCurveInputIterator  begin,
                     XCurveInputIterator  end)
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
  }

  

  ~Arr_non_x_aggregate_insert()
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
