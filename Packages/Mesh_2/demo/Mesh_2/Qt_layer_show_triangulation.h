#ifndef CGAL_QT_LAYER_SHOW_TRIANGULATION_H
#define CGAL_QT_LAYER_SHOW_TRIANGULATION_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Triangulation_2.h>

namespace CGAL {

template <class T>
class Qt_layer_show_triangulation : public Qt_widget_layer
{
public:
	
  Qt_layer_show_triangulation(T *&t,
			      CGAL::Color lc = CGAL::BLUE,
			      int linewidth = 1) 
    : tr(t), color(lc), width(linewidth) {};


  void draw()
  {
    QColor old_color = widget->color();
    int old_width = widget->lineWidth();

    widget->setColor(color);
    widget->setLineWidth(width);
      
    *widget << *tr;

    widget->setLineWidth(old_width);
    widget->setColor(old_color);
  };

private:
  T *&tr;
  CGAL::Color color;
  int width;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_TRIANGULATION_H
