#ifndef CGAL_QT_WIDGET_POLYGON_2_H
#define CGAL_QT_WIDGET_POLYGON_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Polygon_2.h>

namespace CGAL{

template <class Tr,class Co>
Qt_widget& operator<<(Qt_widget& w, const Polygon_2<Tr,Co>& pol)
{
  typedef Polygon_2<Tr,Co>::Vertex_const_iterator VI;
  QPointArray array;

  array.resize(pol.size());

  unsigned int n=0;
  for(VI i=pol.vertices_begin();i!=pol.vertices_end();i++)
    {
      array.setPoint(n++,w.x_pixel(to_double(i->x())),
		     w.y_pixel(to_double(i->y())));
    }
  w.get_painter().drawPolygon(array);
  w.do_paint();
  return w;
}

}//end namespace CGAL

#endif
