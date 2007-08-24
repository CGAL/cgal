#ifndef SHOW_LINES_BASE_H
#define SHOW_LINES_BASE_H

#include "Qt_widget_styled_layer.h"

namespace CGAL {

class Show_lines_base: public Qt_widget_styled_layer {
  Q_OBJECT
public:
  typedef Qt_widget_styled_layer::Style Style;

  Show_lines_base(Color c,
		  int linewidth,
		  QObject* parent = 0, const char* name = 0);

  Show_lines_base(Style* style,
		  QString line_color_name,
		  QString line_width_name,
		  QObject* parent = 0, const char* name = 0);

public slots:
  void setColor(QColor);
  void setLineWidth(int);

protected:
  QString color;
  QString width;
}; //end Show_lines_base

} // namespace CGAL

#endif // SHOW_LINES_BASE_H
