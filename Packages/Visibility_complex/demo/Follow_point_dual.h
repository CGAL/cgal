#ifndef FOLLOW_POINT_DUAL_H
#define FOLLOW_POINT_DUAL_H

#include "types.h"
#include "utils.h"

#include <CGAL/IO/Qt_widget_layer.h>

class Follow_point_dual : public CGAL::Qt_widget_layer
{
  enum State {DRAW_NOTHING, DRAW_POINT, DRAW_LINE};

public:
  Follow_point_dual(Follow_point_dual* twin_layer = 0,
		    const QColor &point_color = Qt::gray,
		    const int point_size = 1,
		    const CGAL::PointStyle point_style = CGAL::CIRCLE,
		    const QColor &line_color = Qt::gray,
		    const int line_width = 1)
    : _point_color(point_color),
      _point_size(point_size),
      _point_style(point_style),
      _line_color(line_color),
      _line_width(line_width),
      other_layer(twin_layer)
    {}

  void draw(){
    if(state == DRAW_NOTHING) return;

    // save properties
    const QColor save_color = widget->color();
    const int save_point_size = widget->pointSize();
    const CGAL::PointStyle save_point_style = widget->pointStyle();
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

	  *widget << p;
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
    FT x, y;
    widget->x_real(e->x(), x),
    widget->y_real(e->y(), y);
    Point_2 p2(x, y);

    state = DRAW_POINT;

    if( p2 != p )
      {
	p = p2;
	widget->redraw();
	other_layer->set_line( dual(p) );
      }
  }

  void leaveEvent(QEvent*)
  {
    state = DRAW_NOTHING;
    other_layer->draw_nothing();
    widget->redraw();
  }

  void set_twin(Follow_point_dual* twin_layer) 
  { 
    other_layer = twin_layer; 
  };

  void set_line(const Line_2& l)
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
  QColor _point_color;
  int _point_size;
  CGAL::PointStyle _point_style;
  QColor _line_color;
  int _line_width;

  Follow_point_dual* other_layer;

  State state;
  Point_2 p;
  Line_2 line;
};

#endif
