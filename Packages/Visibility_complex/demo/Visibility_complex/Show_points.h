#ifndef SHOW_POINTS_H
#define SHOW_POINTS_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/function_objects.h>

namespace CGAL {

class Show_points_base: public Qt_widget_layer {
  Q_OBJECT
public:
  Show_points_base(Color c = CGAL::GREEN,
		   int pointsize = 3,
		   PointStyle pointstyle = CGAL::DISC)
    : color(c), size(pointsize), style(pointstyle) {}

public slots:
  void setColor(QColor);
  void setPointSize(int);
  void setPointStyle(PointStyle);

protected:
  Color color;
  int size;
  PointStyle style;

}; // end Show_points_base  

template <class C, class It,
  class Transform = Identity<typename It::value_type> >
class Show_points : public Show_points_base {
public:
  typedef It iterator;
  typedef iterator (C::* iterator_function)() const;

  Show_points(C *&container,
	      iterator_function begin,
	      iterator_function end,
	      Color c = CGAL::GREEN,
	      int pointsize = 3,
	      PointStyle pointstyle = CGAL::DISC)
    : Show_points_base(c, pointstyle, pointstyle),
      cont(container), _begin(begin), _end(end) {};

  void draw()
  {
    *widget << color << CGAL::PointSize (size) 
	    << CGAL::PointStyle (style);

    for(iterator it = (cont->*_begin)();
	it!=(cont->*_end)();
	++it)
      *widget << Transform()(*it);
  };
private:
  C	*&cont;
  iterator_function _begin;
  iterator_function _end;

};//end class 

} // namespace CGAL

#endif // SHOW_POINTS_H
