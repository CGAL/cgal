#ifndef SHOW_NEAREST_VERTEX_H
#define SHOW_NEAREST_VERTEX_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <qobject.h>
#include <qcolor.h>

namespace CGAL {

template <class DT, class OtherDT, class Line, 
	  class Draw_other_things,
	  class Other_draw_other_things>
class Show_nearest_vertex : public Qt_widget_layer
{
  enum State {DRAW_NOTHING, DRAW_POINT, DRAW_LINE};

public:
  typedef typename DT::Point			Point;
  typedef typename DT::Segment			Segment;
  typedef typename DT::Face_handle		Face_handle;
  typedef typename DT::Vertex                   Vertex;
  typedef typename DT::Vertex_handle		Vertex_handle;
  typedef typename DT::Geom_traits::FT		FT;
  typedef Show_nearest_vertex<DT, OtherDT,
				       Line,
				       Draw_other_things,
				       Other_draw_other_things> Self;
  typedef Show_nearest_vertex<OtherDT, DT,
				       Line,
				       Other_draw_other_things,
				       Draw_other_things> Other;

  Show_nearest_vertex(const DT &t,
		      Other* twin_layer,
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

	  const Point& p = vh->point();

	  *widget << p;
	  Draw_other_things()(widget, vh);

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

    const Vertex_handle v = tr.nearest_vertex(p);

    if( vh != v )
      {
	vh = v;
	widget->redraw();
	other_layer->set_line( dual(vh->point()) );
      }
  }

  void leaveEvent(QEvent*)
  {
    state = DRAW_NOTHING;
    other_layer->draw_nothing();
    widget->redraw();
  }

  void set_twin(Other* twin_layer) 
  { 
    other_layer = twin_layer; 
  };

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
  Vertex_handle vh;
  Line line;

  Other* other_layer;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_NEAREST_VERTEX_H
