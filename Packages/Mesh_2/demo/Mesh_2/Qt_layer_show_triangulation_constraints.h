#ifndef CGAL_QT_LAYER_SHOW_TRIANGULATION_CONSTRAINTS_H
#define CGAL_QT_LAYER_SHOW_TRIANGULATION_CONSTRAINTS_H

#include <CGAL/IO/Qt_widget_layer.h>

namespace CGAL {

template <class T>
class Qt_layer_show_triangulation_constraints : public Qt_widget_layer
{
public:
	
  Qt_layer_show_triangulation_constraints(T *&t,
					  CGAL::Color lc = CGAL::RED,
					  int linewidth = 1)
    : tr(t), color(lc), width(linewidth){};

  void draw()
  {
    widget->lock();

    QColor old_color = widget->color();
    int old_width = widget->lineWidth();

    widget->setColor(color);
    widget->setLineWidth(width);

    for(typename T::Edge_iterator it=tr->edges_begin();
	it!=tr->edges_end();
	it++)
      if(tr->is_constrained(*it))
	*widget << tr->segment(*it);

    widget->setLineWidth(old_width);
    widget->setColor(old_color);
    widget->unlock();
  };
	
private:
  T *&tr;
  CGAL::Color color;
  int width;
};//end class

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_TRIANGULATION_CONSTRAINTS_H
