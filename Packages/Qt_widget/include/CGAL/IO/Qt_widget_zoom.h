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
// file          : include/CGAL/IO/Qt_widget_zoom.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_ZOOM_H
#define CGAL_QT_WIDGET_ZOOM_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_tool.h>
#include <qcolor.h>

namespace CGAL {

class Qt_widget_zoom : public Qt_widget_standard_tool
{
private:
  int	x2, y2;
  bool	circle_already_drawn;
  bool	old_mouse_tracking;

public:
  Qt_widget_zoom() : circle_already_drawn(false) {};

  void draw_circle(int x,int y)
  {
    const int
      d=100, // diameter
      r=50;  // radius (should be d/2 :-)

    RasterOp oldRasterOp = widget->rasterOp();	//save the initial raster mode
    widget->setRasterOp(XorROP);
    QColor oldColor=widget->color(); // save the initial color
    widget->setColor(Qt::green);
    widget->painter().drawEllipse(x-r, y-r, d, d);
    widget->setColor(oldColor);
    widget->setRasterOp(oldRasterOp);
    widget->do_paint();
  };

  void leaveEvent(QEvent*)
  {
    if(circle_already_drawn)
      draw_circle(x2,y2);
    circle_already_drawn=false;
  };

  void widget_repainted(){
    circle_already_drawn=false;
  };

  void mousePressEvent(QMouseEvent *e)
  {
    const double ratios=2.0;
    double
      x=widget->x_real(e->x()),
      y=widget->y_real(e->y());
		
    if(e->button() == Qt::LeftButton)
			widget->zoom_in(ratios, x, y);
    if(e->button() == Qt::RightButton)
			widget->zoom_out(ratios, x, y);
		
    widget->redraw();
    emit(redraw()); 
  };

  void mouseMoveEvent(QMouseEvent *e)
  {
    int
      x = e->x(),
      y = e->y();
    
    widget->lock();
    draw_circle(x,y); // draw the new circle
    if(circle_already_drawn) // erase the old one if needed
      draw_circle(x2,y2);
    widget->unlock();
		
    //save the last coordinates to redraw the screen
    x2 = x;
    y2 = y;
    circle_already_drawn=true;
  };

  void attaching()
  {
    old_mouse_tracking=widget->hasMouseTracking();
    widget->setMouseTracking(TRUE);
    oldcursor = widget->cursor();
    widget->setCursor(crossCursor);
    circle_already_drawn=false;
  };

  void detaching()
  {
    if(circle_already_drawn)
      draw_circle(x2,y2); // erase the circle if needed
    widget->setCursor(oldcursor);
    widget->setMouseTracking(old_mouse_tracking);
  };
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WIDGET_ZOOM_H
