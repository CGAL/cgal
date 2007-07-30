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
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_ADDITION_H
#define CGAL_ARR_ADDITION_H

#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_utils.h>
#include <CGAL/Sweep_line_2/Arr_construction_event.h>
#include <CGAL/Sweep_line_2/Arr_construction_subcurve.h>
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
  typedef Arr_construction_subcurve<Traits>                   Subcurve; 
  typedef Arr_construction_event<Traits,
                                 Subcurve,
                                 Arrangement>              Event;
  
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

  Arr_addition(Arrangement &arr) :
    m_arr(&arr),
    m_traits(new Traits(*(arr.get_traits()))),
    m_visitor(&arr),
    m_sweep_line(m_traits, &m_visitor)
  {}

  template<class CurveInputIterator>
  void insert_curves(CurveInputIterator begin, 
                     CurveInputIterator end)
  {
    // Subdivide all input curves into basic x-monotone curves and isolated
    // points.
    std::list<typename Base_traits::X_monotone_curve_2>   base_xcurves;
    std::list<typename Base_traits::Point_2>              base_points;

    make_x_monotone (begin,
                     end,
                     std::back_inserter(base_xcurves),
                     std::back_inserter(base_points),
                     m_arr->get_traits());

    // Covert the basic objects to the types defined by the addition traits.
    typename std::list<typename Base_traits::X_monotone_curve_2>::iterator xit;
    typename std::list<typename Base_traits::Point_2>::iterator            pit;
    std::vector<X_monotone_curve_2>   xcurves_vec (base_xcurves.size() +
                                                   m_arr->number_of_edges());
    std::vector<Point_2>              iso_points (base_points.size() +
                                         m_arr->number_of_isolated_vertices());
    int                               i_cv = 0, i_pt = 0;

    for (xit = base_xcurves.begin();
         xit != base_xcurves.end(); ++xit, i_cv++)
    {
      xcurves_vec[i_cv] = X_monotone_curve_2 (*xit);
    }

    for (pit = base_points.begin();
         pit != base_points.end(); ++pit, i_pt++)
    {
      iso_points[i_pt] = Point_2 (*pit);
    }

    // Add the x-montone curves and the isolated points from the arrangement.
    Edge_iterator                  eit;
    Halfedge_handle                he;

    for (eit = m_arr->edges_begin(); eit != m_arr->edges_end(); ++eit, i_cv++) 
    {
      if (eit->direction() == SMALLER)
        he = eit->twin();
      else
        he = eit;

      xcurves_vec[i_cv] = X_monotone_curve_2 (he->curve(), he);
    }

    Vertex_iterator                vit;

    for(vit = m_arr->vertices_begin(); vit != m_arr->vertices_end(); ++vit)
    {
      if (vit->is_isolated())
      {
        iso_points[i_pt] = Point_2 (vit->point(), vit);
        i_pt++;
      }
    }

    // Perform the sweep.
    m_sweep_line.sweep(xcurves_vec.begin(),
                       xcurves_vec.end(),
                       iso_points.begin(),
                       iso_points.end());
    return;
  }

  template<class XCurveInputIterator>
  void insert_x_curves(XCurveInputIterator begin,
                       XCurveInputIterator end)
  {
    std::vector<X_monotone_curve_2>      xcurves_vec;
    std::vector<Point_2>                 iso_points;

    Edge_iterator    eit;
    Halfedge_handle  he;

    for (eit = m_arr->edges_begin(); eit != m_arr->edges_end(); ++eit) 
    {
      if (eit->direction() == SMALLER)
        he = eit->twin();
      else
        he = eit;

      xcurves_vec.push_back(X_monotone_curve_2(he->curve(), he));
    }

    Vertex_iterator v_itr;

    for (v_itr = m_arr->vertices_begin(); 
         v_itr != m_arr->vertices_end(); ++v_itr)
    {
      if(v_itr->is_isolated())
        iso_points.push_back(Point_2(v_itr->point(), v_itr));
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

  Arrangement*         m_arr;
  Traits*              m_traits;
  Visitor              m_visitor;
  Sweep_line           m_sweep_line;
};

CGAL_END_NAMESPACE

#endif
