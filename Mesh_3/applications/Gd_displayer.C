#include "Gd_displayer.h"

Gd_displayer::Gd_displayer(int sx, int sy)
{
  im = gdImageCreate(sx, sy);
  gdImageColorAllocate(im, 255, 255, 255);
  set_window(0., 1., 0., 1.);
}

Gd_displayer::Gd_displayer(gdImagePtr gd)
{
  im = gd;
  set_window(0., 1., 0., 1.);
}

Gd_displayer::~Gd_displayer()
{
  gdImageDestroy(im);
}

void Gd_displayer::fill_rectangle(double x1, double y1,
                                  double x2, double y2,
                                  CGAL::Color c)
{
  gdImageFilledRectangle(im,
                         x_pixel(x1), y_pixel(y2),
                         x_pixel(x2), y_pixel(y1),
                         gd_color(c));
}

void Gd_displayer::segment(double x1, double y1,
                           double x2, double y2,
                           CGAL::Color c)
{
  gdImageLine(im,
              x_pixel(x1), y_pixel(y1),
              x_pixel(x2), y_pixel(y2),
              gd_color(c));
}

void Gd_displayer::set_window(const double x_min, const double x_max,
                              const double y_min, const double y_max,
                              int x, int y, int sx, int sy)
{
  xmin = x_min;
  xmax = x_max;
  ymin = y_min;
  ymax = y_max;
  xclip = x;
  yclip = y;
  sxclip = sx;
  syclip = sy;
  gdImageSetClip(im, xclip, yclip, xclip+sxclip, yclip+syclip);
  xscal = sxclip / (xmax-xmin);
  yscal = syclip / (ymax-ymin);  
}

int Gd_displayer::x_pixel(double x) const
{
  return( xclip + static_cast<int>((x-xmin)*xscal) );
}

int Gd_displayer::y_pixel(double y) const
{
  return( yclip - static_cast<int>((y-ymax)*yscal) );
}

bool Gd_displayer::save_png(const char* filename)
{
  FILE *pngout = fopen(filename, "wb");
  if (pngout != NULL)
  {
    gdImagePng(im, pngout);
    fclose(pngout);
    return true;
  }
  else 
    return false;
}
