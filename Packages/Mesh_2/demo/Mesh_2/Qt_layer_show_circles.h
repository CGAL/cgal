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
// file          : include/CGAL/IO/Qt_layer_show_circles.h
// package       : Qt_widget
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

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

  Qt_layer_show_circles(T &t, QLabel& l) : tr(t), label(l),
    do_erase(false) {};

  void draw() const
  {
    do_erase = false;
  };
  
  void mousePressEvent(QMouseEvent* e)
    {
      if (tr.dimension()<1) return;
      FT
	x=static_cast<FT>(widget->x_real(e->x())),
	y=static_cast<FT>(widget->y_real(e->y()));

      Point p(x,y);

      int li;
      Locate_type lt;
      Face_handle fh = tr.locate(p,lt,li);
      if(lt == T::FACE)
	draw_circle(fh);	  
    };

  void mouseMoveEvent(QMouseEvent *e)
    {
      if (tr.dimension()<1) return;
      FT
	x=static_cast<FT>(widget->x_real(e->x())),
	y=static_cast<FT>(widget->y_real(e->y()));

      Point p(x,y);

      int li;
      Locate_type lt;
      Face_handle fh = tr.locate(p,lt,li);
      if(lt == T::FACE)
	{
	  if(fh!=old_face)
	    {
	      widget->lock();
	      
	      if(do_erase) draw_circle(old_face);
	      draw_circle(fh);
	      old_face=fh;
	      
	      widget->unlock();
	      QString s;
	      s.setNum(to_double(tr.squared_minimum_sine(fh)), 'g', 5);
	      label.setText(s);
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

  void leaveEvent(QEvent *e)
    {
      if (tr.dimension()<1) return;      
      if(do_erase)
	draw_circle(old_face);
    };

private:
  void draw_circle(const Face_handle& fh) const
    {
      RasterOp oldRaster = widget->rasterOp();
      QColor oldcolor = widget->color();
      int oldwidth = widget->lineWidth();
      
      *widget << CGAL::GRAY;
      *widget << LineWidth(1);
      widget->get_painter().setRasterOp(NotROP);
      
      Point v=((*fh).vertex(0))->point();
      Point c=tr.circumcenter(fh);
      
      *widget << Circle(c,squared_distance(v,c));
      widget->setColor(oldcolor);
      widget->setLineWidth(oldwidth);
      widget->setRasterOp(oldRaster);
      widget->do_paint();
    }

  T      &tr;
  QLabel &label;
  Face_handle  old_face;
  bool	       do_erase;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_CIRCLES_H
