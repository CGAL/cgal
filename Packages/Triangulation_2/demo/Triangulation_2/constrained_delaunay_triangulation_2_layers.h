// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : triangulation_2_constrained_layers.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifndef CGAL_TRIANGULATION_2_CONSTRAINED_LAYERS_H
#define CGAL_TRIANGULATION_2_CONSTRAINED_LAYERS_H


#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Triangulation_2.h>

template <class T>
class Qt_layer_show_constraints : public CGAL::Qt_widget_layer {
public:
  typedef typename T::Edge            Edge;
  typedef typename T::Finite_edges_iterator
                                      Finite_edges_iterator;

  Qt_layer_show_constraints(T &t) : tr(t){};

  void draw()
  {  
    Finite_edges_iterator it = tr.finite_edges_begin();
    *widget << CGAL::RED << CGAL::LineWidth(2);
    while(it != tr.finite_edges_end()) {
      if(tr.is_constrained(*it))
        *widget << tr.segment(it);
      ++it;
    }
  };
private:
  T	&tr;
  
};//end class 

template <class T>
class Qt_layer_show_points : public CGAL::Qt_widget_layer {
public:
  typedef typename T::Point           Point;
  typedef typename T::Segment         Segment;
  typedef typename T::Vertex          Vertex;
  typedef typename T::Vertex_iterator	Vertex_iterator;

  Qt_layer_show_points(T &t) : tr(t){};

  void draw()
  {  
    Vertex_iterator it = tr.vertices_begin(), 
		beyond = tr.vertices_end();
    *widget << CGAL::GREEN << CGAL::PointSize (3) 
		<< CGAL::PointStyle (CGAL::DISC);    
    while(it != beyond) {      
      *widget << (*it).point();
      ++it;
    }
  };
private:
  T	&tr;
  
};//end class 

template <class T>
class Qt_layer_show_triangulation : public CGAL::Qt_widget_layer
{
public:
	
  Qt_layer_show_triangulation(T &t) : tr(t){};


  void draw()
  {
    *widget << CGAL::BLUE; 
    *widget << tr;
  };
	
private:
  T &tr;
};//end class 


#endif
