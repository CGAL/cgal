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
#include <CGAL/Cartesian.h>
#include <qobject.h>


namespace CGAL {

template <class T>
class Qt_layer_show_circles : public Qt_widget_layer {
public:
  typedef typename T::Point           Point;
  typedef typename T::Segment         Segment;
  typedef typename T::Finite_faces_iterator Finite_faces_iterator;

  Qt_layer_show_circles(T &t) : tr(t){};

  void draw()
  {  
    typedef Cartesian<double>::Circle_2 Circle_double;
    typedef Cartesian<double>::Point_2 Point_double;

    widget->lock();
    QColor oldcolor = widget->color();
    int oldwidth = widget->lineWidth();

    *widget << CGAL::GRAY;
    *widget << LineWidth(1);

    for(Finite_faces_iterator it=tr.finite_faces_begin();
	it!=tr.finite_faces_end();
	it++)
      {
	Point v=((*it).vertex(0))->point();
	double x=to_double(v.x());
	double y=to_double(v.y());
	Point_double vd(x,y);
	
	Point c=tr.circumcenter(it);
	x=to_double(c.x());
	y=to_double(c.y());
	Point_double cd(x,y);

	*widget << Circle_double(cd,squared_distance(vd,cd));
      }

    widget->setColor(oldcolor);
    widget->setLineWidth(oldwidth);
    widget->unlock();
  };
private:
  T	&tr;
  
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_CIRCLES_H
