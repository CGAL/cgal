#ifndef CGAL_QT_LAYER_SHOW_POINTS_H
#define CGAL_QT_LAYER_SHOW_POINTS_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/function_objects.h>
#include <qcolor.h>

namespace CGAL {

template <class C, class It,
  class Transform = Identity<typename It::value_type> >
class Qt_layer_show_points : public Qt_widget_layer {
public:
  typedef It (C::* iterator_function)() const;

  Qt_layer_show_points(C *&container,
		       iterator_function begin,
		       iterator_function end,
		       const Transform& t = Transform(),
		       Color c = CGAL::GREEN,
		       int pointsize = 3,
		       PointStyle pointstyle = CGAL::DISC)
    : cont(container), _begin(begin), _end(end), _color(c),
      size(pointsize), style(pointstyle), trans(t) {};

  Qt_layer_show_points(C *&container,
		       iterator_function begin,
		       iterator_function end,
		       Color c = CGAL::GREEN,
		       int pointsize = 3,
		       PointStyle pointstyle = CGAL::DISC)
    : cont(container), _begin(begin), _end(end), _color(c),
      size(pointsize), style(pointstyle), trans(Transform()) {};

  void draw()
  {
    widget->lock();
    QColor old_color = widget->color();
    int old_point_size = widget->pointSize();

    *widget << _color << CGAL::PointSize (size) 
	    << CGAL::PointStyle (style);

    for(It it = (cont->*_begin)();
	it!=(cont->*_end)();
	++it)
      *widget << trans(*it);
    
    widget->setPointSize(old_point_size);
    widget->setColor(old_color);
    widget->unlock();
  };

  QColor color() const { return _color; };

  int pointSize() const { return size; };

  void setColor(const QColor c)
  {
    _color = c;
    widget->redraw();
  };

  void setPointSize(const int s)
  {
    size = s;
    widget->redraw();
  };

private:
  C	*&cont;
  iterator_function _begin;
  iterator_function _end;
  Color _color;
  int size;
  PointStyle style;
  const Transform trans;
  
};//end class 

//   template <class C, class It,
// 	    class Transform = Identity<typename It::value_type> >
//   //  Qt_layer_show_points<C, It, Transform>*
//   void*
//   make_show_points_layer(C *&container,
// 			 It (C::*)() const begin,
// 			 It (C::*)() const end,
// 			 Transform& t = Transform(),
// 			 Color c = CGAL::GREEN,
// 			 int pointsize = 3,
// 			 PointStyle pointstyle = CGAL::DISC)
//   {
//     return new CGAL::Qt_layer_show_points<C, It, Transform>
//       (container, begin, end, t, c, pointsize, pointstyle);
//   };

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_POINTS_H
