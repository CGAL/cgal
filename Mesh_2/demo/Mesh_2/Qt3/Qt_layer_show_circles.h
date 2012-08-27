// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_QT_LAYER_SHOW_CIRCLES_H
#define CGAL_QT_LAYER_SHOW_CIRCLES_H

#include <CGAL/Circle_2.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Cartesian.h>
#include <qobject.h>
#include <qcolor.h>


namespace CGAL {

// T::Geom_traits has to be a CGAL kernel!
template <class T>
class Qt_layer_show_circles : public Qt_widget_layer {
public:
  typedef typename T::Point           Point;
  typedef typename T::Segment         Segment;
  typedef typename T::Finite_faces_iterator Finite_faces_iterator;
  typedef typename T::Locate_type     Locate_type;
  typedef typename T::Face_handle     Face_handle;
  typedef typename T::Geom_traits GT;
  typedef typename GT::Circle_2 Circle;
  typedef typename GT::FT		FT;

  Qt_layer_show_circles(T* t,
			CGAL::Color c = CGAL::GRAY,
			int linewidth = 1,
			CGAL::Color fill_color = CGAL::WHITE,
			bool filled = false,
                        QObject* parent = 0, const char* name = 0) :
    Qt_widget_layer(parent, name),
    tr(t), do_erase(false), color(c), width(linewidth),
    fillcolor(fill_color), fill(filled) {};

  void draw()
  {
    Qt_widget_layer::draw();
    do_erase = false;
  };

  void mousePressEvent(QMouseEvent* e)
    {
      if (tr->dimension()<1) return;
      FT
	x=static_cast<FT>(widget->x_real(e->x())),
	y=static_cast<FT>(widget->y_real(e->y()));

      Point p(x,y);

      int li;
      Locate_type lt;
      Face_handle fh = tr->locate(p,lt,li);
      if(lt == T::FACE)
	draw_circle(fh);
    };

  void mouseMoveEvent(QMouseEvent *e)
    {
      if (tr->dimension()<1) return;
      FT
	x=static_cast<FT>(widget->x_real(e->x())),
	y=static_cast<FT>(widget->y_real(e->y()));

      Point p(x,y);

      int li;
      Locate_type lt;
      Face_handle fh = tr->locate(p,lt,li);
      if(lt == T::FACE)
	{
	  if(fh!=old_face)
	    {
	      widget->lock();

	      if(do_erase) draw_circle(old_face);
	      draw_circle(fh);
	      old_face=fh;

	      widget->unlock();
	      do_erase=true;
	    }
	}
      else
	{
	  if(do_erase)
	    draw_circle(old_face);
	  do_erase=false;
	}
    };

  void leaveEvent(QEvent* e)
    {
      Qt_widget_layer::leaveEvent(e);
      if (tr->dimension()<1) return;
      if(do_erase)
	draw_circle(old_face);
    };

private:
  void draw_circle(const Face_handle& fh) const
    {
      RasterOp oldRaster = widget->rasterOp();
      QColor oldcolor = widget->color();
      QColor oldFillColor = widget->fillColor();
      int oldwidth = widget->lineWidth();
      bool oldFilled = widget->isFilled();

      *widget << color;
      *widget << LineWidth(width) << FillColor(fillcolor);
      widget->setFilled(fill);
      widget->get_painter().setRasterOp(NotROP);

      Point v=((*fh).vertex(0))->point();
      Point c=tr->circumcenter(fh);

      *widget << Circle(c,squared_distance(v,c));
      widget->setColor(oldcolor);
      widget->setLineWidth(oldwidth);
      widget->setFillColor(oldFillColor);
      widget->setFilled(oldFilled);
      widget->setRasterOp(oldRaster);
      widget->do_paint();
    }

  T* tr;
  Face_handle  old_face;
  bool	       do_erase;

  CGAL::Color color;
  int width;
  CGAL::Color fillcolor;
  bool fill;
};//end class

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_CIRCLES_H
