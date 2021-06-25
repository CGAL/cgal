#include "Distribution_displayer.h"

namespace CGAL {
  class Qt_widget;
}

struct Qt_widget_displayer : public Distribution_displayer
{
  Qt_widget_displayer(CGAL::Qt_widget* widget);

  void fill_rectangle(double x1, double y1,
                      double x2, double y2,
                      CGAL::IO::Color c);

  void segment(double x1, double y1,
               double x2, double y2,
               CGAL::IO::Color c);
private:
  CGAL::Qt_widget* widget;
};
