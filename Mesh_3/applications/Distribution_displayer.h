#ifndef DISTRIBUTION_DISPLAYER_H
#define DISTRIBUTION_DISPLAYER_H
#include <CGAL/IO/Color.h>

struct Distribution_displayer
{
  virtual void fill_rectangle(double x1, double y1,
                              double x2, double y2,
                              CGAL::Color c) = 0;

  virtual void segment(double x1, double y1,
                       double x2, double y2,
                       CGAL::Color c) = 0;

  virtual ~Distribution_displayer() {};
};

#endif
