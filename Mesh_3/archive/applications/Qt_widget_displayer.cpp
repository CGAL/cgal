#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/IO/Qt_widget.h>

#include "Qt_widget_displayer.h"

typedef CGAL::Simple_cartesian<double> Simple_kernel;
typedef Simple_kernel::Iso_rectangle_2 Rectangle_2;
typedef Simple_kernel::Segment_2 Segment_2;
typedef Simple_kernel::Point_2 Point_2;

Qt_widget_displayer::Qt_widget_displayer(CGAL::Qt_widget* w) : widget(w) {}

void Qt_widget_displayer::fill_rectangle(double x1, double y1,
                                         double x2, double y2,
                                         CGAL::IO::Color color)
{
  *widget << CGAL::FillColor(color)
          << color
          << Rectangle_2(Point_2(x1, y1), Point_2(x2, y2));
}

void Qt_widget_displayer::segment(double x1, double y1,
                                  double x2, double y2,
                                  CGAL::IO::Color color)
{
  *widget << color << Segment_2(Point_2(x1, y1),
                                Point_2(x2, y2));
}
