#ifndef CGAL_QT_LAYER_SHOW_POINTS_H
#define CGAL_QT_LAYER_SHOW_POINTS_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/function_objects.h>

namespace CGAL {

template <class C, class It,
  class Transform = Identity<typename It::value_type> >
class Qt_layer_show_points : public Qt_widget_layer {
public:
  typedef It Vertex_iterator;
  typedef Vertex_iterator (C::* iterator_function)() const;

  Qt_layer_show_points(C *&container,
		       iterator_function begin,
		       iterator_function end,
		       Color c = CGAL::GREEN,
		       int pointsize = 3,
		       PointStyle pointstyle = CGAL::DISC)
    : cont(container), _begin(begin), _end(end), color(c),
      size(pointsize), style(pointstyle) {};

  void draw()
  {
    *widget << color << CGAL::PointSize (size) 
	    << CGAL::PointStyle (style);

    for(Vertex_iterator it = (cont->*_begin)();
	it!=(cont->*_end)();
	++it)
      *widget << Transform()(*it);
  };
private:
  C	*&cont;
  iterator_function _begin;
  iterator_function _end;
  Color color;
  int size;
  PointStyle style;
  
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_POINTS_H
