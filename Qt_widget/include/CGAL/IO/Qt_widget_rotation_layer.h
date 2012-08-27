// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

#ifndef CGAL_QT_WIDGET_ROTATION_LAYER_H
#define CGAL_QT_WIDGET_ROTATION_LAYER_H

#include <cstdio>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/pixmaps/hand.xpm>
#include <CGAL/IO/pixmaps/holddown.xpm>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <qrect.h>
#include <qcursor.h>

typedef CGAL::Simple_cartesian<double>    SCD;
typedef CGAL::Point_2<SCD>                Point;
typedef CGAL::Vector_2<SCD>               Vector;
typedef CGAL::Direction_2<SCD>            Direction;
typedef CGAL::Aff_transformation_2<SCD>   Transformation;

namespace CGAL {

class Qt_widget_rotation_layer : public Qt_widget_layer
{
public:
  Qt_widget_rotation_layer(QObject* parent = 0, const char* name = 0)
    : Qt_widget_layer(parent, name), wasrepainted(true), on_first(false){};

private:
  QCursor oldcursor;

  void draw(){
    wasrepainted = true;
  };
/*
  void timerEvent( QTimerEvent *)
  {
    if(on_first)
      widget->setCursor(QCursor( 
              QPixmap( (const char**)holddown_xpm)));
    else
      widget->setCursor(QCursor( 
              QPixmap( (const char**)hand_xpm)));
  }
*/
  bool is_pure(Qt::ButtonState s){
    if((s & Qt::ControlButton) ||
       (s & Qt::ShiftButton) ||
       (s & Qt::AltButton))
      return 0;
    else
      return 1;
  }

  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton 
       && is_pure(e->state()))
    {
      widget->setCursor(QCursor( QPixmap( (const char**)holddown_xpm)));
      if (!on_first){
	      first_x = e->x();
	      first_y = e->y();
	      on_first = true;
      }
      
      Point end_point[4] = {Point(11, 25), Point(11, 25), 
        Point(11, 25), Point(11, 25)};
      for(int i=0; i<4; i++)
        end_point[i] = (*t)(end_point[i]);
      //Point_2 first_line_end1(11, 25), first_line_end2(11, 0);
      //Point_2 second_line_end1(25, 0), second_line_end2(25, 13);
      //first_line_end1 = (*t)(first_arrow_end);
      //second_arrow_end = (*t)(second_arrow_end);
      
      //save the initial raster mode
      RasterOp old = widget->rasterOp();	
      widget->setRasterOp(XorROP);
      widget->lock();
      *widget << CGAL::WHITE;
      widget->get_painter().drawLine(11, 0, 11, 25);
      widget->get_painter().drawLine(0, 13, 25, 13);
      widget->get_painter().drawLine(11, 0, 9, 5);
      widget->get_painter().drawLine(11, 0, 13, 5);
      widget->get_painter().drawLine(25, 13, 20, 11);
      widget->get_painter().drawLine(25, 13, 20, 15);
      *widget << CGAL::RED;

      widget->unlock();
      widget->setRasterOp(old);

    }
  };

  void mouseReleaseEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton
       && is_pure(e->state()))
    {
      /*
      double x, y, xfirst2, yfirst2;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      widget->x_real(first_x, xfirst2);
      widget->y_real(first_y, yfirst2);
			
      double	xmin, xmax, ymin, ymax, distx, disty;
      if(x < xfirst2) {xmin = x; xmax = xfirst2;}
      else {xmin = xfirst2; xmax = x;};
      if(y < yfirst2) {ymin = y; ymax = yfirst2;}
      else {ymin = yfirst2; ymax = y;};
      distx = xfirst2 - x;
      disty = yfirst2 - y;
      widget->move_center(distx, disty);
      
      */
      /*
      //(*t) = tprim * (*t);
      on_first = false;


      */
      double xc = widget->x_min() + (widget->x_max() - widget->x_min())/2;
      double yc = widget->y_min() + (widget->y_max() - widget->y_min())/2;
      Transformation tprim(CGAL::ROTATION, Direction(xc+0.0001, yc+0.0001), 1, 100);
      //Transformation tprim(CGAL::TRANSLATION, Vector(-0.5, 0));
      //Transformation tprim(CGAL::SCALING, 1.2);
      (*t) = (*t) * tprim;
      on_first = false;
      widget->setCursor(QCursor( QPixmap( (const char**)hand_xpm)));
      widget->redraw();
    }
  }
  /*
  void mouseMoveEvent(QMouseEvent *e)
  {
    char tempc1[130], tempc2[40];
    double xcoord, ycoord;
    if(on_first)
    {
      int x = e->x();
      int y = e->y();
      //save the initial raster mode
      RasterOp old = widget->rasterOp();	
      widget->setRasterOp(XorROP);
      widget->lock();
        *widget << CGAL::GRAY;
      if(!wasrepainted) {
        widget->x_real(x2 - first_x, xcoord);
        widget->x_real(y2 - first_y, ycoord);
        std::sprintf(tempc1, " dx=%20.6f", xcoord);
        std::sprintf(tempc2, ", dy=%20.6f", ycoord);
        strcat(tempc1, tempc2);
        widget->get_painter().drawLine(first_x, first_y, x2, y2);
        *widget << CGAL::GREEN;
        widget->get_painter().drawText(x2, y2, tempc1, 49);
        *widget << CGAL::GRAY;
      }
      widget->x_real(x - first_x, xcoord);
      widget->x_real(y - first_y, ycoord);
      std::sprintf(tempc1, " dx=%20.6f", xcoord);
      std::sprintf(tempc2, ", dy=%20.6f", ycoord);
      strcat(tempc1, tempc2);
      widget->get_painter().drawLine(first_x, first_y, x, y);
      *widget << CGAL::GREEN;
      widget->get_painter().drawText(x, y, tempc1, 49);
      widget->unlock();
      widget->setRasterOp(old);

      //save the last coordinates to redraw the screen
      x2 = x;
      y2 = y;
      wasrepainted = false;
    }
  };
*/
  void activating()
  {
    t = widget->get_transformation();
    oldcursor = widget->cursor();
    widget->setCursor(QCursor( QPixmap( (const char**)hand_xpm)));
    wasrepainted = true;
//	  startTimer( 100 );
  };

  void deactivating()
  {
    widget->setCursor(oldcursor);
//    killTimers();
  };

  int             first_x, first_y;
  int             x2, y2;
  bool            wasrepainted;
  bool            on_first;
  Transformation  *t;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WIDGET_ROTATION_LAYER_H
