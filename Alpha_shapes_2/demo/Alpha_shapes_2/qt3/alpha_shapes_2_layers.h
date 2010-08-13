// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Radu Ursu

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Triangulation_2.h>
#include <CGAL/IO/Qt_widget_Alpha_shape_2.h>
#include <qimage.h>

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
    *widget << CGAL::LineWidth(2) ;
    *widget << CGAL::BLUE;
    *widget << tr;
  };

private:
  T &tr;
};//end class

template <class T>
class Qt_layer_show_voronoi : public CGAL::Qt_widget_layer
{
public:
  Qt_layer_show_voronoi(T &t1) : tr(t1){};

  void draw()
  {
    *widget << CGAL::LineWidth(2) ;
    *widget << CGAL::RED ;
    tr.draw_dual(*widget) ;
  };

private:
  T	&tr;
};//end class

class Qt_layer_show_alpha_shape : public CGAL::Qt_widget_layer
{
public:
  Qt_layer_show_alpha_shape(Alpha_shape *a1) : a(a1){};

  void draw(){
    a->set_mode(Alpha_shape::GENERAL);
    *widget << CGAL::LineWidth(2) << CGAL::GREEN;
    *widget << *a;
  }
private:
  Alpha_shape *a;
};

class Qt_layer_show_image : public CGAL::Qt_widget_layer{
public:
  Qt_layer_show_image(QImage* i) : image(i){}
  ~Qt_layer_show_image(){};
  void draw(){
    widget->get_painter().drawImage(0, 0, *image, 0, 0, widget->width(), widget->height(), 0);
  }
  QImage *image;
};
