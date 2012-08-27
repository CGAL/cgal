// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//                 Baruch Zukerman        <baruchzu@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_OVERLAY_2_H
#define CGAL_ENVELOPE_OVERLAY_2_H

#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Envelope_3/Envelope_overlay_functor.h>

#include <iostream>

namespace CGAL {

template <class MinimizationDiagram_2, 
          class OverlayFunctor = Envelope_overlay_functor<MinimizationDiagram_2> >
class Envelope_overlay_2
{
public:
  typedef MinimizationDiagram_2                                  Minimization_diagram_2;
  
  typedef typename Minimization_diagram_2::Face_handle           Face_handle;
  typedef typename Minimization_diagram_2::Face_iterator         Face_iterator;

  typedef typename Minimization_diagram_2::Vertex_handle         Vertex_handle;
  typedef typename Minimization_diagram_2::Vertex_iterator       Vertex_iterator;

  typedef typename Minimization_diagram_2::Halfedge_handle       Halfedge_handle;
  typedef typename Minimization_diagram_2::Halfedge_iterator     Halfedge_iterator;

  typedef OverlayFunctor                                         Overlay_functor;
protected:
  typedef typename Minimization_diagram_2::Geometry_traits_2     Traits;
  typedef typename Traits::Xy_monotone_surface_3                 Xy_monotone_surface_3;

public:
  
  void operator()(Minimization_diagram_2& md1,
                  Minimization_diagram_2& md2,
                  Minimization_diagram_2& result)
  {
    CGAL_assertion(md1.is_valid());
    CGAL_assertion(md2.is_valid());

    Overlay_functor overlay_func(md1, md2, result);
    overlay(md1, md2, result, overlay_func);
        
    CGAL_assertion_code(post_test_assertions(result));
  }


public:
  
  /*
  void print_face(Face_handle fh)
  {
    std::cout << (fh->is_unbounded() ? "unbounded" : "bounded");
    
    if (fh->get_is_set())
    {
      std::cout << " #data= " << fh->number_of_data_objects();
      if (fh->number_of_data_objects() > 0)
        std::cout << " data= " << fh->get_data();
    }

    if (fh->get_aux_is_set(0))
    {
      std::cout << " #data1= " << get_number_of_aux_data_objects(fh, 0);
      if (get_number_of_aux_data_objects(fh, 0)>0)
        std::cout << " data#1= " << get_aux_data(fh, 0);
    }
    if (fh->get_aux_is_set(1))
    {
      std::cout << " #data2= " << get_number_of_aux_data_objects(fh, 1);
      if (get_number_of_aux_data_objects(fh, 1)>0)
        std::cout << " data#2= " << get_aux_data(fh, 1);
    }
    std::cout << std::endl;
  }

  // print the aux data in the faces of md
  void print_faces(Minimization_diagram_2& md)
  {
    Face_iterator fit = md.faces_begin();
    for(; fit != md.faces_end(); ++fit)
    {
      Face_handle fh = fit;
      print_face(fh);
    }
    std::cout << std::endl;
  }

  void print_vertices(Minimization_diagram_2& md)
  {
    Vertex_iterator vit = md.vertices_begin();
    for(; vit != md.vertices_end(); ++vit)
    {
      Vertex_handle vh = vit;
      std::cout << vh->point();

      if (vh->get_is_set())
      {
        std::cout << " #data= " << vh->number_of_data_objects();
        if (vh->number_of_data_objects() > 0)
          std::cout << " data= " << vh->get_data();
      }

      if (vh->get_aux_is_set(0))
      {
        std::cout << " #data1= " << get_number_of_aux_data_objects(vh, 0);
        if (get_number_of_aux_data_objects(vh, 0)>0)
          std::cout << " data#1= " << get_aux_data(vh, 0);
      }
      if (vh->get_aux_is_set(1))
      {
        std::cout << " #data2= " << get_number_of_aux_data_objects(vh, 1);
        if (get_number_of_aux_data_objects(vh, 1)>0)
          std::cout << " data#2= " << get_aux_data(vh, 1);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;  
  }

  void print_edges(Minimization_diagram_2& md)
  {
    Halfedge_iterator hit = md.halfedges_begin();
    for(; hit != md.halfedges_end(); ++hit, ++hit)
    {
      Halfedge_handle hh = hit;
      std::cout << hh->curve();

      if (hh->get_is_set())
      {
        std::cout << " #data= " << hh->number_of_data_objects();
        if (hh->number_of_data_objects() > 0)
          std::cout << " data= " << hh->get_data();
      }


      if (hh->get_aux_is_set(0))
      {
        std::cout << " #data1= " << get_number_of_aux_data_objects(hh, 0);
        if (get_number_of_aux_data_objects(hh, 0)>0)
          std::cout << " data#1= " << get_aux_data(hh, 0);
      }
      if (hh->get_aux_is_set(1))
      {
        std::cout << " #data2= " << get_number_of_aux_data_objects(hh, 1);

        if (get_number_of_aux_data_objects(hh, 1)>0)
          std::cout << " data#2= " << get_aux_data(hh, 1);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  */

  void post_test_assertions(Minimization_diagram_2& md)
  {
    // check that all data is filled in result
    Face_iterator fi = md.faces_begin();
    for(; fi != md.faces_end(); ++fi)
    {
      Face_handle fh = fi;
      CGAL_assertion_msg(fh->get_aux_is_set(0), "data from md1 on face is not set");
      CGAL_assertion_msg(fh->get_aux_is_set(1), "data from md2 on face is not set");
    }

    Halfedge_iterator hi = md.halfedges_begin();
    for(; hi != md.halfedges_end(); ++hi)
    {
      Halfedge_handle hh = hi;
      CGAL_assertion_msg(hh->get_aux_is_set(0), "data from md1 on halfedge is not set");
      CGAL_assertion_msg(hh->get_aux_is_set(1), "data from md2 on halfedge is not set");
    }

    Vertex_iterator vi = md.vertices_begin();
    for(; vi != md.vertices_end(); ++vi)
    {
      Vertex_handle vh = vi;
      CGAL_assertion_msg(vh->get_aux_is_set(0), "data from md1 on vertex is not set");
      CGAL_assertion_msg(vh->get_aux_is_set(1), "data from md2 on vertex is not set");
    }
  }
protected:
  // helper methods
  template <class FeatureHandle>
  Xy_monotone_surface_3 get_aux_data(FeatureHandle fh, unsigned int id)
  {
    const Object& o = fh->get_aux_source(id);
    Xy_monotone_surface_3 data;

    Halfedge_handle h;
    Vertex_handle v;
  	Face_handle f;
  	if (assign(v, o))
  	  data = v->get_data();
  	else if (assign(h, o))
  	  data = h->get_data();
  	else
  	{
  	  CGAL_assertion(assign(f, o));
      assign(f, o);
  	  data = f->get_data();
  	}
    return data;
  }
  template <class FeatureHandle>
  int get_number_of_aux_data_objects(FeatureHandle fh, unsigned int id)
  {
	  const Object& o = fh->get_aux_source(id);
    int data;

    Halfedge_handle h;
    Vertex_handle v;
  	Face_handle f;
  	if (assign(v, o))
  	  data = v->number_of_data_objects();
  	else if (assign(h, o))
  	  data = h->number_of_data_objects();
  	else
  	{
  	  CGAL_assertion(assign(f, o));
      assign(f, o);
  	  data = f->number_of_data_objects();
  	}
    return data;
  }

};

} //namespace CGAL

#endif
