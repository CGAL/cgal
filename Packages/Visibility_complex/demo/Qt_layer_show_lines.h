#ifndef CGAL_QT_LAYER_SHOW_LINES_H
#define CGAL_QT_LAYER_SHOW_LINES_H

#include <CGAL/IO/Qt_widget_layer.h>

namespace CGAL {

template <class C>
class Qt_layer_show_lines : public Qt_widget_layer {
public:
  typedef typename C::iterator	Vertex_iterator;

  Qt_layer_show_lines(C *&container, Color c=CGAL::GREEN, int linewidth=3)
    : cont(container), color(c), width(linewidth) {};

  void draw()
  {  
    *widget << color << CGAL::LineWidth(width);

    for(Vertex_iterator it = cont->begin();
	it!=cont->end();
	++it)
      *widget << (*it);
  };
private:
  C	*&cont;
  Color color;
  int width;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_LINES_H
