// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Radu Ursu
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file          : include/CGAL/IO/Qt_widget_get_arc.h

#ifndef CGAL_QT_WIDGET_GET_ARC_H
#define CGAL_QT_WIDGET_GET_ARC_H

#include <cmath>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/squared_distance_2.h>

#ifndef CGAL_QT_WIDGET_GET_POINT_BUTTON
#define CGAL_QT_WIDGET_GET_POINT_BUTTON Qt::LeftButton
#endif

namespace CGAL {

namespace CGALi {

class Qt_widget_get_arc_helper
  : public CGAL::Qt_widget_layer
{
Q_OBJECT

public:

  Qt_widget_get_arc_helper(QObject* parent = 0, const char* name = 0)
    : Qt_widget_layer(parent, name) {}

  virtual void got_it() { emit new_object_time(); }

signals:
  void new_object_time();
};

} // namespace CGALi

template <class R>
class Qt_widget_get_arc
  : public CGALi::Qt_widget_get_arc_helper
{
  typedef typename R::FT                         FT;
  typedef typename R::Circle_2                   Circle;
  typedef typename R::Circular_arc_2             Circular_arc_2;
  typedef typename R::Segment_2                  Segment;
  typedef typename R::Line_2                     Line;
  typedef typename R::Point_2                    Point;

public:

  Qt_widget_get_arc(const QCursor c=QCursor(Qt::crossCursor),
		       QObject* parent = 0, const char* name = 0)
     : CGALi::Qt_widget_get_arc_helper(parent, name),
       cursor(c), firstpoint(false),
       firsttime(true), secondpoint(false), thirdpoint(false),
       clockwise(false)
  {}

  void draw()
  {
    if (secondpoint) {
      *widget << CGAL::GREEN << introduced_circle;
      if (thirdpoint){
	*widget << CGAL::BLUE << Segment(Point(x1, y1), first_point);
      }
    }
    firsttime = true;
  }

private:

  bool is_pure(Qt::ButtonState s) {
    return ! ((s & Qt::ControlButton) ||
              (s & Qt::ShiftButton) ||
              (s & Qt::AltButton));
  }

  void mousePressEvent(QMouseEvent *e)
  {
    if (e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON
        && !firstpoint
        && is_pure(e->state()))
    {
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      x1 = x;
      y1 = y;
      x2 = x;
      y2 = y;
      firstpoint = true;
      center = Point(x, y);
    } else if (e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON
	       && !secondpoint && is_pure(e->state()))
    {
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      introduced_circle = Circle(Point(x1,y1),
	      squared_distance(Point(x1, y1), Point(x,y)));
      secondpoint = true;
      firsttime = true;
      on_the_boundary = Point(x, y);
    } else if (e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON
	       && !thirdpoint && is_pure(e->state())){
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      first_point = Point(x, y);
      thirdpoint = true;
      firsttime = true;
    } else if (e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON){
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      second_point = Point(x, y);
      firstpoint = false; secondpoint = false; thirdpoint = false;
      firsttime = true;
      got_it();
    }
  }

  void keyPressEvent(QKeyEvent *e)
  {
    switch ( e->key() ) {
      case Key_Escape:			// key_escape
         if (firstpoint && !secondpoint)
         {
	   firstpoint = false;
           RasterOp old_raster = widget->rasterOp();
           QColor old_color = widget->color();
           widget->lock();
           widget->setRasterOp(XorROP);
           *widget << CGAL::GREEN;
           *widget << Circle(Point(x1,y1),
                      squared_distance(Point(x1, y1), Point(x2,y2)));
           widget->setRasterOp(old_raster);
           widget->setColor(old_color);
           widget->unlock();
	   firsttime = true;
         }
         break;
      case Key_C:
	if(thirdpoint){
	  clockwise = !clockwise;
	  if(!firsttime){
	    //draw the arc
            QColor old_color = widget->color();
            RasterOp old_raster = widget->rasterOp();
	    //save the initial raster mode
            widget->setRasterOp(XorROP);
            widget->lock();
	    *widget << CGAL::RED;
	    int x_screen = widget->x_pixel(CGAL::to_double(center.x()));
	    int y_screen = widget->y_pixel(CGAL::to_double(center.y()));
	    int x_screen_b = widget->x_pixel(CGAL::to_double(on_the_boundary.x()));
	    int y_screen_b = widget->y_pixel(CGAL::to_double(on_the_boundary.y()));
	    int radius = (int)sqrt(CGAL::to_double((x_screen_b - x_screen) *
			  (x_screen_b - x_screen) +
			  (y_screen_b - y_screen) *
			  (y_screen_b - y_screen)));
	    double diff = 180/3.1415926;
	
	    double a = atan2( to_double(first_point.y() - center.y()),
			  to_double(first_point.x() - center.x()));
	    double a2p = atan2( to_double(y2 - center.y()),
			      to_double(x2 - center.x()));
            if(!clockwise){double tempa = a; a = a2p; a2p = tempa;}

	    double alen2;
	    if(a2p > a)
	      alen2 = a2p - a;
	    else
	      alen2 = 2 * 3.1415926 + a2p - a;

            widget->get_painter().drawArc(x_screen - radius,
				        y_screen - radius,
				        2 * radius, 2 * radius,
				        (int)(a * diff * 16),
					(int)(alen2 * diff * 16));
	    a = atan2( to_double(first_point.y() - center.y()),
			  to_double(first_point.x() - center.x()));
	    double a2 = atan2( to_double(y2 - center.y()),
			   to_double(x2 - center.x()));
	    if(clockwise){double tempa = a; a = a2; a2 = tempa;}
	    double alen;
	    if(a2 > a)
	      alen = a2 - a;
	    else
	      alen = 2 * 3.1415926 + a2 - a;

            widget->get_painter().drawArc(x_screen - radius,
				        y_screen - radius,
				        2 * radius, 2 * radius,
				        (int)(a * diff * 16),
					(int)(alen * diff * 16));
            widget->unlock();
            widget->setRasterOp(old_raster);
            widget->setColor(old_color);
	  }
	}
	break;
    }//endswitch
  }

  void leaveEvent(QEvent *)
  {
    if (firstpoint && !secondpoint)
    {
      QColor old_color = widget->color();
      RasterOp old_raster = widget->rasterOp(); //save the initial raster mode

      widget->lock();
      widget->setRasterOp(XorROP);
      *widget << CGAL::GREEN << Circle(Point(x1,y1),
                                       squared_distance(Point(x1, y1),
                                                        Point(x2,y2)));
      widget->unlock();
      widget->setRasterOp(old_raster);
      widget->setColor(old_color);
      firsttime = true;
    }
  }

  void mouseMoveEvent(QMouseEvent *e)
  {
    if (firstpoint==true && !secondpoint)
    {		
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      QColor old_color = widget->color();
      RasterOp old_raster = widget->rasterOp();//save the initial raster mode
      widget->setRasterOp(XorROP);
      widget->lock();
      *widget << CGAL::GREEN;
      if(!firsttime)
	*widget << Circle(Point(x1,y1),
		  squared_distance(Point(x1, y1), Point(x2,y2)));
      *widget << Circle(Point(x1,y1),
		  squared_distance(Point(x1, y1), Point(x,y)));
      widget->unlock();
      widget->setRasterOp(old_raster);
      widget->setColor(old_color);

      //save the last coordinates to redraw the screen
      x2 = x;
      y2 = y;	
      firsttime = false;
    } else if(secondpoint || thirdpoint){
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      QColor old_color = widget->color();
      RasterOp old_raster = widget->rasterOp();//save the initial raster mode
      widget->setRasterOp(XorROP);
      widget->lock();
      *widget << CGAL::BLUE;
      if(!firsttime)
        *widget << Segment(Point(x1,y1),Point(x2,y2));

      *widget << Segment(Point(x1,y1),Point(x,y));
      widget->unlock();
      widget->setRasterOp(old_raster);
      widget->setColor(old_color);
      //save the last coordinates to redraw the screen

      if(thirdpoint){
	//draw the arc
        QColor old_color = widget->color();
        RasterOp old_raster = widget->rasterOp();//save the initial raster mode
        widget->setRasterOp(XorROP);
        widget->lock();
	*widget << CGAL::RED;

	int x_screen = widget->x_pixel(CGAL::to_double(center.x()));
	int y_screen = widget->y_pixel(CGAL::to_double(center.y()));
	int x_screen_b = widget->x_pixel(CGAL::to_double(on_the_boundary.x()));
	int y_screen_b = widget->y_pixel(CGAL::to_double(on_the_boundary.y()));

	int radius = (int)sqrt(CGAL::to_double((x_screen_b - x_screen) *
			  (x_screen_b - x_screen) +
			  (y_screen_b - y_screen) *
			  (y_screen_b - y_screen)));
	
	double diff = 180/3.1415926;

	double a = atan2( to_double(first_point.y() - center.y()),
			  to_double(first_point.x() - center.x()));

	double a2 = atan2( to_double(y - center.y()),
			   to_double(x - center.x()));
	if(clockwise){double tempa = a; a = a2; a2 = tempa;}
	double alen;
	if(a2 > a)
	  alen = a2 - a;
	else
	  alen = 2 * 3.1415926 + a2 - a;

	if(!firsttime){
	  double a = atan2( to_double(first_point.y() - center.y()),
			  to_double(first_point.x() - center.x()));
	  double a2p = atan2( to_double(y2 - center.y()),
			      to_double(x2 - center.x()));
          if(clockwise){double tempa = a; a = a2p; a2p = tempa;}

	  double alen2;
	  if(a2p > a)
	    alen2 = a2p - a;
	  else
	    alen2 = 2 * 3.1415926 + a2p - a;

          widget->get_painter().drawArc(x_screen - radius,
				        y_screen - radius,
				        2 * radius, 2 * radius,
				        (int)(a * diff * 16),
					(int)(alen2 * diff * 16));
	}

	widget->get_painter().drawArc(x_screen - radius,
				      y_screen - radius,
				      2 * radius, 2 * radius,
				      (int)(a * diff * 16),
				      (int)(alen * diff * 16));
        widget->unlock();
        widget->setRasterOp(old_raster);
        widget->setColor(old_color);
      }
      x2 = x;
      y2 = y;
      firsttime = false;
    }
  }

public:

  Circular_arc_2 get_circular_arc() {
    //Circle circ (introduced_circle.center(),
                 //introduced_circle.squared_radius());
    Circle circ = introduced_circle;

    if (clockwise)
      std::swap(first_point, second_point);

    Line l1 (first_point, center);
    Line l2 (second_point, center);

    return Circular_arc_2 (circ,
     // FIXME which one is newer?	
     //l1, compare_lexicographically_xy(first_point,  center) == CGAL::SMALLER,
     //l2, compare_lexicographically_xy(second_point, center) == CGAL::SMALLER);
	       l1, CGAL::compare_xy(first_point,  center) == CGAL::SMALLER,
	       l2, CGAL::compare_xy(second_point, center) == CGAL::SMALLER);
  }

private:

  void activating()
  {
    oldpolicy = widget->focusPolicy();
    widget->setFocusPolicy(QWidget::StrongFocus);
    oldcursor = widget->cursor();
    widget->setCursor(cursor);
  }

  void deactivating()
  {
    widget->setFocusPolicy(oldpolicy);
    widget->setCursor(oldcursor);
    firstpoint = false;
  }

  QCursor cursor;
  QCursor oldcursor;

  FT    x1, //the X of the first point
        y1; //the Y of the first point
  FT	x2, //the old second point's X
        y2; //the old second point's Y
  bool	firstpoint, //true if the user left clicked once
	firsttime;  //true if the line is not drawn
  bool  secondpoint;//true after the secons point
  bool  thirdpoint; //true after the third point

  Circle introduced_circle;
  Point  center;
  Point  on_the_boundary;

  Point  first_point;
  Point  second_point;

  bool  clockwise;

  QWidget::FocusPolicy	oldpolicy;
};

} // namespace CGAL

#endif // CGAL_QT_WIDGET_GET_ARC_H
