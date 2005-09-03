#include <gd.h>
#include "Distribution_displayer.h"

struct Gd_displayer : public Distribution_displayer
{
  Gd_displayer(int sx, int sy);

  Gd_displayer(gdImagePtr gd);

  void fill_rectangle(double x1, double y1,
                      double x2, double y2,
                      CGAL::Color c);

  void segment(double x1, double y1,
               double x2, double y2,
               CGAL::Color c);

  ~Gd_displayer();

  inline 
  void set_window(const double x_min,
		  const double x_max, 
		  const double y_min, 
		  const double y_max)
  {
    set_window(x_min, x_max, y_min, y_max,
               0, 0, gdImageSX(im), gdImageSY(im));
  }

  void set_window(const double x_min,
		  const double x_max, 
		  const double y_min, 
		  const double y_max,
                  int x, int y,
                  int sx, int sy);

  int x_pixel(double x) const;
  int y_pixel(double y) const;

  bool save_png(const char* filename);

  inline gdImagePtr image()
  {
    return im;
  }

  inline 
  int gd_color(CGAL::Color c)
  {
    int i = gdImageColorExact(im, 
                              c.red(),
                              c.green(),
                              c.blue());
    if( i < 0 )
      return gdImageColorAllocate(im, 
                                  c.red(),
                                  c.green(),
                                  c.blue());
    else
      return i;
  }

private:
  double xmin, xmax, ymin, ymax; // real dimensions
  int xclip, yclip, sxclip, syclip;
  double xscal, yscal; // scaling factors

  gdImagePtr im;
};
