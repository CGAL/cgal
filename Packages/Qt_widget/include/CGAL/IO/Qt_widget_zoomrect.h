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
// file          : include/CGAL/IO/Qt_widget_zoomrect.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_ZOOMRECT_H
#define CGAL_QT_WIDGET_ZOOMRECT_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <qrect.h>



#ifndef CGAL_QT_WIDGET_ZOOMRECT_BUTTON
#define CGAL_QT_WIDGET_ZOOMRECT_BUTTON Qt::LeftButton
#endif


namespace CGAL {

class Qt_widget_zoomrect : public Qt_widget_layer
{
public:
  int first_x, first_y, x2, y2;
  bool widgetrepainted;
  bool on_first;

  Qt_widget_zoomrect() : widgetrepainted(TRUE), on_first(FALSE) {};

private:
  void draw(){
    widgetrepainted = TRUE;
  };
  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == CGAL_QT_WIDGET_ZOOMRECT_BUTTON)
    {
      if (!on_first)
      {
	      first_x = e->x();
        first_y = e->y();
        on_first = TRUE;
      }
    }
  };

  void mouseReleaseEvent(QMouseEvent *e)
  {
    if(e->button() == CGAL_QT_WIDGET_ZOOMRECT_BUTTON)		
    {
      double
	x=widget->x_real(e->x()),
	y=widget->y_real(e->y()),
	xfirst2 = widget->x_real(first_x),
	yfirst2 = widget->y_real(first_y);
			
	double	xmin, xmax, ymin, ymax;
	if(x < xfirst2) {xmin = x; xmax = xfirst2;}
	else {xmin = xfirst2; xmax = x;};
	if(y < yfirst2) {ymin = y; ymax = yfirst2;}
	else {ymin = yfirst2; ymax = y;};

	widget->set_window(xmin, xmax, ymin, ymax);
	widget->redraw();
	on_first = FALSE;
	}
  }
  void mouseMoveEvent(QMouseEvent *e)
  {
    if(on_first)
    {
      int
	      x = e->x(),
	      y = e->y();
			
      RasterOp old = widget->rasterOp();	//save the initial raster mode
      QColor old_color=widget->color();
      widget->setRasterOp(XorROP);
      widget->lock();
      widget->setColor(Qt::green);
      if(!widgetrepainted)
	      widget->painter().drawRect(first_x, first_y, x2 - first_x, y2 - first_y);
      widget->painter().drawRect(first_x, first_y, x - first_x, y - first_y);
      widget->unlock();
      widget->setColor(old_color);
      widget->setRasterOp(old);

      //save the last coordinates to redraw the screen
      x2 = x;
      y2 = y;
      widgetrepainted = FALSE;
    }
  };

  void activating()
  {
    oldcursor = widget->cursor();
    widget->setCursor(crossCursor);
    widgetrepainted = TRUE;
  };

  void deactivating()
  {
    widget->setCursor(oldcursor);
		
    RasterOp old = widget->rasterOp();	//save the initial raster mode
    widget->setRasterOp(XorROP);
    widget->lock();
    *widget << CGAL::GREEN;
    widget->unlock();
    widget->setRasterOp(old);
  };
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WIDGET_ZOOMRECT_H
