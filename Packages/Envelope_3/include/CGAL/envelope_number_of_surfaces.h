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
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_NUMBER_OF_SURFACES_H
#define CGAL_ENVELOPE_NUMBER_OF_SURFACES_H

#include <CGAL/Timer.h>

#include <iostream>
#include <cassert>
#include <list>
#include <set>

CGAL_BEGIN_NAMESPACE

// return the number of different xy-monotone surfaces that appear in the
// minimization diagram "arr"
template <class MinimizationDiagram_2>
std::size_t envelope_find_number_of_surfaces(MinimizationDiagram_2& arr)
{
  typedef MinimizationDiagram_2                                     Minimization_diagram_2;
  typedef typename Minimization_diagram_2::Traits_2                 Traits_2;
  typedef typename Traits_2::Xy_monotone_surface_3                  Xy_monotone_surface_3;

  typedef std::list<Xy_monotone_surface_3>                          Surfaces_list;
  
  Surfaces_list slist;
  envelope_find_unique_surfaces(arr, std::back_inserter(slist));
  
  return slist.size();
}

template <class MinimizationDiagram_2, class OutputIterator>
OutputIterator envelope_find_unique_surfaces(MinimizationDiagram_2& arr, OutputIterator o)
{
  typedef MinimizationDiagram_2                                     Minimization_diagram_2;
  typedef typename Minimization_diagram_2::Traits_2                 Traits_2;
  typedef typename Traits_2::Xy_monotone_surface_3                  Xy_monotone_surface_3;

  typedef typename Minimization_diagram_2::Halfedge_iterator        Halfedge_iterator;
  typedef typename Minimization_diagram_2::Face_iterator            Face_iterator;
  typedef typename Minimization_diagram_2::Vertex_iterator          Vertex_iterator;
  typedef typename Minimization_diagram_2::Dcel::Face_data_iterator Data_iterator;

  typedef std::set<Xy_monotone_surface_3>                           Surfaces_set;
  typedef typename std::set<Xy_monotone_surface_3>::iterator        Surfaces_set_it; 

  Surfaces_set sset;
  envelope_find_unique_surfaces_set(arr, sset);
  
  for(Surfaces_set_it it = sset.begin(); it != sset.end(); ++it)
  {
    *o = *it;
    ++o;
  }
  return o;
}

template <class MinimizationDiagram_2, class Surfaces_set>
void envelope_find_unique_surfaces_set(MinimizationDiagram_2& arr, Surfaces_set& sset)
{
  typedef MinimizationDiagram_2                                     Minimization_diagram_2;
  typedef typename Minimization_diagram_2::Traits_2                 Traits_2;
  typedef typename Traits_2::Xy_monotone_surface_3                  Xy_monotone_surface_3;

  typedef typename Minimization_diagram_2::Halfedge_iterator        Halfedge_iterator;
  typedef typename Minimization_diagram_2::Face_iterator            Face_iterator;
  typedef typename Minimization_diagram_2::Vertex_iterator          Vertex_iterator;
  typedef typename Minimization_diagram_2::Dcel::Face_data_iterator Data_iterator;

  Data_iterator di;
  // vertices
  Vertex_iterator vi = arr.vertices_begin();
  for(; vi != arr.vertices_end(); ++vi)
  {
    di = vi->begin_data();
    for(; di != vi->end_data(); ++di)
      sset.insert(*di);
  }
  // edges
  Halfedge_iterator hi = arr.halfedges_begin();
  for(; hi != arr.halfedges_end(); ++hi)
  {
    di = hi->begin_data();
    for(; di != hi->end_data(); ++di)
      sset.insert(*di);
  }

  // faces
  Face_iterator fi = arr.faces_begin();
  for(; fi != arr.faces_end(); ++fi)
  {
    di = fi->begin_data();
    for(; di != fi->end_data(); ++di)
      sset.insert(*di);
  }
}

CGAL_END_NAMESPACE

#endif 
