#ifndef CGAL_QT_WIDGET_DELAUNAY_TRIANGULATION_2_H
#define CGAL_QT_WIDGET_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Delaunay_triangulation_2.h>

template < class Gt, class Tds >
Qt_widget&
operator<<(Qt_widget& w,  const Delaunay_triangulation_2<Gt,Tds> &dt)
{
  w.lock();
  dt.draw_triangulation(w);
  w.unlock();
  return w;
}

#endif
