#ifndef SHOW_LINES_H
#define SHOW_LINES_H

#include "Show_lines_base.h"

namespace CGAL {

template <class C, class It,
	  class Transform = Identity<typename It::value_type> >
class Show_lines : public Show_lines_base {
public:
  typedef Qt_widget_styled_layer::Style Style;

  typedef It iterator;
  typedef iterator (C::* iterator_function)() const;

  Show_lines(C *&container,
	     iterator_function begin,
	     iterator_function end,
	     Color c=CGAL::GREEN,
	     int linewidth=3,
	     QObject* parent = 0, const char* name = 0)
    : Show_lines_base(c, linewidth,
		      parent, name),
      cont(container), _begin(begin), _end(end) {};

  Show_lines(C *&container,
	     iterator_function begin,
	     iterator_function end,
	     Style* style,
	     QString line_color_name,
	     QString line_width_name,
	     QObject* parent = 0, const char* name = 0)
    : Show_lines_base(style, line_color_name, line_width_name,
		      parent, name),
      cont(container), _begin(begin), _end(end) {};

  void draw()
  {  
    widget->lock();
    {
      QColor old_color = widget->color();
      int old_width = widget->lineWidth();
      
      widget->setColor(style()->getColor(color));
      widget->setLineWidth(style()->getInt(width));
      
      for(iterator it = (cont->*_begin)();
	  it!=(cont->*_end)();
	  ++it)
	*widget << Transform()(*it);

      widget->setColor(old_color);
      widget->setLineWidth(old_width);
    }
    widget->unlock();
  };

private:
  C	*&cont;
  iterator_function _begin;
  iterator_function _end;
};//end class 

} // namespace CGAL

#endif // SHOW_LINES_H
