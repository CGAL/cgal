#ifndef CGAL_QT_LAYER_SHOW_POINTS_H
#define CGAL_QT_LAYER_SHOW_POINTS_H

#include <CGAL/IO/Qt_widget_layer.h>

namespace CGAL {

template <class C>
class Qt_layer_show_points : public Qt_widget_layer {
public:
  typedef typename C::iterator	Vertex_iterator;

  Qt_layer_show_points(C *&container, Color c=CGAL::GREEN, int pointsize=3, 
		       PointStyle pointstyle = CGAL::DISC)
    : cont(container), color(c), size(pointsize), style(pointstyle) {};

  void draw()
  {  
    *widget << color << CGAL::PointSize (size) 
	    << CGAL::PointStyle (style);

    for(Vertex_iterator it = cont->begin();
	it!=cont->end();
	++it)
      *widget << (*it);
  };
private:
  C	*&cont;
  Color color;
  int size;
  PointStyle style;
  
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_POINTS_H
