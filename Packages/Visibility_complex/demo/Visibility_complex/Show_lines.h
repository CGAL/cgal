#ifndef SHOW_LINES_H
#define SHOW_LINES_H

#include <CGAL/IO/Qt_widget_layer.h>

namespace CGAL {

template <class C, class It,
	  class Transform = Identity<typename It::value_type> >
class Show_lines : public Qt_widget_layer {
public:
  typedef It iterator;
  typedef iterator (C::* iterator_function)() const;

  Show_lines(C *&container,
	     iterator_function begin,
	     iterator_function end,
	     Color c=CGAL::GREEN,
	     int linewidth=3)
    : cont(container), _begin(begin), _end(end), color(c),
      width(linewidth) {};

  void draw()
  {  
    *widget << color << CGAL::LineWidth(width);

    for(iterator it = (cont->*_begin)();
	it!=(cont->*_end)();
	++it)
      *widget << Transform()(*it);
  };
private:
  C	*&cont;
  iterator_function _begin;
  iterator_function _end;
  Color color;
  int width;

};//end class 

} // namespace CGAL

#endif // SHOW_LINES_H
