// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : include/CGAL/IO/Qt_layer_show_nearest_vertex.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_LAYER_SHOW_NEAREST_VERTEX_H
#define CGAL_QT_LAYER_SHOW_NEAREST_VERTEX_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <qobject.h>
#include <qcolor.h>

namespace CGAL {

template <class DT, class Line>
class Qt_layer_show_nearest_vertex : public Qt_widget_layer
{
  enum State {DRAW_NOTHING, DRAW_POINT, DRAW_LINE};

public:
  typedef typename DT::Point			Point;
  typedef typename DT::Segment			Segment;
  typedef typename DT::Face_handle		Face_handle;
  typedef typename DT::Vertex_handle		Vertex_handle;
  typedef typename DT::Geom_traits::FT		FT;
  typedef Qt_layer_show_nearest_vertex<DT, Line> Self;

  Qt_layer_show_nearest_vertex(const DT &t,
			       Self* twin_layer,
			       const QColor &point_color = Qt::green,
			       const int point_size = 10,
			       const PointStyle point_style = CIRCLE,
			       const QColor &line_color = Qt::green,
			       const int line_width = 1)
    : tr(t),
      _point_color(point_color),
      _point_size(point_size),
      _point_style(point_style),
      _line_color(line_color),
      _line_width(line_width),
      state(DRAW_NOTHING),
      other_layer(twin_layer)
    {};

  void draw(){
    if(state == DRAW_NOTHING) return;

    // save properties
    const QColor save_color = widget->color();
    const int save_point_size = widget->pointSize();
    const PointStyle save_point_style = widget->pointStyle();
    const int save_line_width = widget->lineWidth();

    widget->lock();

    // set properties
    switch ( state )
      {
      case DRAW_POINT: 
	{
	  widget->setColor(_point_color);
	  widget->setPointSize(_point_size);
	  widget->setPointStyle(_point_style);

	  *widget << point;
	  break;
	}
      case DRAW_LINE:
	{
	  widget->setColor(_line_color);
	  widget->setLineWidth(_line_width);
	  *widget << line;
	  break;
	}
      case DRAW_NOTHING:
	;
      }

    widget->unlock();

    // restore properties
    widget->setColor(save_color);
    widget->setPointSize(save_point_size);
    widget->setPointStyle(save_point_style);
    widget->setLineWidth(save_line_width);
  }

  void mouseMoveEvent(QMouseEvent *e)
  {
    if (tr.dimension()<1) return;

    FT x, y;
    widget->x_real(e->x(), x),
    widget->y_real(e->y(), y);
    Point p(x, y);

    state = DRAW_POINT;

    Vertex_handle v = tr.nearest_vertex(p);

    if( point != v->point() )
      {
	point = v->point();
	widget->redraw();
	other_layer->set_line( dual(point) );
      }
  }

  void leaveEvent(QEvent *e)
  {
    state = DRAW_NOTHING;
    other_layer->draw_nothing();
  }

  void set_twin(Self* twin_layer) { other_layer = twin_layer; };

  void set_line(const Line& l)
  {
    line = l;
    state = DRAW_LINE;
    widget->redraw();
  }

  void draw_nothing()
  {
    state = DRAW_NOTHING;
    widget->redraw();
  }

private:
  const DT &tr;
  QColor _point_color;
  int _point_size;
  PointStyle _point_style;
  QColor _line_color;
  int _line_width;

  State state;
  Point point;
  Line line;

  Self* other_layer;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_NEAREST_VERTEX_H
