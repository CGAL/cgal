// ======================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium

// Copyright (c) 2002 ENS de Paris
//
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
//
// The Qt widget we provide for CGAL is distributed under the QPL,
// which is Trolltech's open source license. For more information see 
//     http://www.trolltech.com/developer/licensing/qpl.html
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// file          : include/CGAL/IO/Qt_widget_get_simple_polygon.h
// package       : Qt_widget (1.2.30)
// author(s)     : Laurent Rineau && Radu Ursu
// release       : CGAL-2.4
// release_date  : 2002, May 16
//
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

#ifndef CGAL_QT_WIDGET_GET_SIMPLE_POLYGON_H
#define CGAL_QT_WIDGET_GET_SIMPLE_POLYGON_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>  
#include <list>

#include <qcursor.h>

namespace CGAL {
template <class Polygon>
class Qt_widget_get_simple_polygon : public Qt_widget_layer
{
public:
  typedef typename Polygon::Point_2   Point_2;
  typedef typename Polygon::Segment_2 Segment_2;
  typedef typename Polygon::Edge_const_iterator  ECI;
  typedef typename Polygon::FT	      FT;

  Qt_widget_get_simple_polygon()
    : active(false), first_time(true) {}

  void draw()
  {
    if(poly.size() > 1)
    {
      ECI  it;
      widget->lock();
      RasterOp old_rasterop=widget->rasterOp();
      widget->get_painter().setRasterOp(XorROP);
      *widget << CGAL::GREEN;
      for(it = poly.edges_begin(); it != --poly.edges_end(); it++)
        *widget << *it;
      widget->setRasterOp(old_rasterop);
      widget->unlock();
    }
    return;
  };
private:

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
    if(e->button() == Qt::LeftButton && is_pure(e->state()))
    {
      FT x=static_cast<FT>(widget->x_real(e->x()));
      FT y=static_cast<FT>(widget->y_real(e->y()));
      if(!active)
      {
        active=true;
        widget->setMouseTracking(TRUE);
        last_of_poly = Point_2(x, y);
        poly.push_back(Point_2(x, y));	
      } else{
        if (last_of_poly == Point_2(x,y)) return;
        rubber_old = Point_2(x, y);
        if(is_simple()){
          poly.push_back(Point_2(x,y));	
          //show the last rubber as edge of the polygon
          widget->lock();
            RasterOp old_rasterop=widget->rasterOp();
            widget->get_painter().setRasterOp(XorROP);
	    *widget << CGAL::WHITE;
            *widget << Segment_2(rubber, last_of_poly);
	    *widget << CGAL::GREEN;
	    *widget << Segment_2(rubber, last_of_poly);
            widget->setRasterOp(old_rasterop);
          widget->unlock();
          last_of_poly = Point_2(x, y);
        }
      }
      return;
    };
    if(e->button() == Qt::RightButton && is_pure(e->state()))
    {
      if (active) {
        if(!poly.is_simple()) return;
        if(poly.is_clockwise_oriented())
          poly.reverse_orientation ();
        assert( ! poly.is_clockwise_oriented());
	  
        widget->new_object(make_object(poly));
        active = false;
        first_time = true;
        poly.erase(poly.vertices_begin(), poly.vertices_end());
        widget->redraw();
      }
    };
  };//end mousePressEvent

  // MK:: removed because it does not compile... problem with iterator
  // of polygon...
#if 0
  void keyPressEvent(QKeyEvent *e)
  {
    switch ( e->key() ) {
      case Key_Escape:			// key_escape
          if(poly.size() > 1){
            widget->lock();
              RasterOp old_rasterop=widget->rasterOp();
              widget->get_painter().setRasterOp(XorROP);
              *widget << CGAL::GREEN;
              *widget << Segment_2(*(----poly.vertices_end()), last_of_poly);
              *widget << CGAL::WHITE;
              *widget << Segment_2(rubber, last_of_poly);
              *widget << Segment_2(rubber, *(----poly.vertices_end()));
              widget->setRasterOp(old_rasterop);
            widget->unlock();
            poly.erase(--poly.vertices_end());
            last_of_poly = *(--poly.vertices_end());
          }
          break;
    }//endswitch
  }
#endif

  void mouseMoveEvent(QMouseEvent *e)
  {
    if (active)
    {
      FT x=static_cast<FT>(widget->x_real(e->x()));
      FT y=static_cast<FT>(widget->y_real(e->y()));

      rubber = Point_2(x, y);
      widget->lock();
        RasterOp old_rasterop=widget->rasterOp();
        widget->get_painter().setRasterOp(XorROP);
	// MK:: changed color...
	//        *widget << CGAL::WHITE;      	
	*widget << CGAL::ORANGE;      	
        if(!first_time)
          *widget << Segment_2(rubber_old, last_of_poly);
        *widget << Segment_2(rubber, last_of_poly);
        first_time = false;
        rubber_old = rubber;
        widget->setRasterOp(old_rasterop);
      widget->unlock();
    }
  };
  void activating()
  {	
    oldcursor = widget->cursor();
    widget->setCursor(crossCursor);
    oldpolicy = widget->focusPolicy();
    widget->setFocusPolicy(QWidget::StrongFocus);
  };
  
  void deactivating()
  {
    poly.erase(poly.vertices_begin(), poly.vertices_end());
    active = false;
    first_time = true;
    widget->setCursor(oldcursor);
    widget->setFocusPolicy(oldpolicy);
    widget->redraw();
  };
  void leaveEvent(QEvent *e)
  {
    if (active)
    {
      widget->lock();
        RasterOp old_rasterop=widget->rasterOp();
        widget->get_painter().setRasterOp(XorROP);
        *widget << CGAL::WHITE;
        *widget << Segment_2(rubber_old, last_of_poly);
        widget->setRasterOp(old_rasterop);
      widget->unlock();
      first_time = true;
    }
  }
private:
  bool is_simple()
  {
    Segment_2 rubber_segment(rubber, last_of_poly);
    if(poly.size() > 1)
    {
      ECI it;      
      for(it = poly.edges_begin(); it != ----poly.edges_end(); it++)
      {
        if(do_intersect(*it, rubber_segment))
          return false;
      }
      //if I'm out of this means that all the edges, 
      //didn't intersect the last one
      ++it;
      Object o = intersection(*it, rubber_segment);
      Point_2 p;
      if(assign(p, o))
        return true;
      else
        return false;
    }
    return true;
  }
  
protected:
  bool	active,     //true if the first point was inserted
		first_time; //true if it is the first time when 
		      //draw the rubber band
  Point_2 rubber,     //the new point of the rubber band
		  last_of_poly,	//the last point of the polygon
		  rubber_old; //the old point of the rubber band
  Polygon poly;	      //the polygon
  QWidget::FocusPolicy	oldpolicy;
  QCursor oldcursor;
};

} // namespace CGAL

#endif // CGAL_QT_WIDGET_GET_SIMPLE_POLYGON_H
