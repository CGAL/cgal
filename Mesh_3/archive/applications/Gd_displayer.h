#include <gd.h>
#include "Distribution_displayer.h"

struct Gd_displayer : public Distribution_displayer
{
  // --- CONSTRUCTORS ---

  /** Create a Gd_displayer using with a new GD image.
      /param sx width of the image
      /param sy height of the image
  */
  Gd_displayer(int sx, int sy);

  /** Create a Gd_displayer using an existing GD Image. */
  Gd_displayer(gdImagePtr gd);

  ~Gd_displayer();


  /** \name VIRTUAL FUNCTIONS FROM Distribution_displayer */
  //@{
  void fill_rectangle(double x1, double y1,
                      double x2, double y2,
                      CGAL::IO::Color c);

  void segment(double x1, double y1,
               double x2, double y2,
               CGAL::IO::Color c);
  //@}

  /** \name FUNCTIONS SPECIFIC TO Gd_displayer */
  //@{
  /** Set world coordinates of the image.*/
  inline
  void set_window(const double x_min,
                  const double x_max,
                  const double y_min,
                  const double y_max)
  {
    set_window(x_min, x_max, y_min, y_max,
               0, 0, gdImageSX(im), gdImageSY(im));
  }

  /** Set world coordinates and a clip zone.
      \param x x coordinate (in pixel) of the upper left point of
      the clip zone
      \param y y coordinate (in pixel) of the upper left point of
      the clip zone
      \param sx width of the clip zone
      \param sy height of the clip zone
  */
  void set_window(const double x_min,
                  const double x_max,
                  const double y_min,
                  const double y_max,
                  int x, int y,
                  int sx, int sy);

  /** Returns the x coordinate, in pixel, of the world-coordinate x. */
  int x_pixel(double x) const;

  /** Returns the y coordinate, in pixel, of the world-coordinate y. */
  int y_pixel(double y) const;

  /** Save the GD image to a PNG file. */
  bool save_png(const char* filename);

  /** Returns the GD image. */
  inline gdImagePtr image()
  {
    return im;
  }

  /** Returns the index of the color c in the image palette. */
  inline
  int gd_color(CGAL::IO::Color c)
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
  //@}
  // end specific functions
private:
  /** \name real dimensions
      @{ */
  double xmin, xmax, ymin, ymax; //@}
  /** \name clip zone
      @{ */
  int xclip, yclip, sxclip, syclip; //@}
  /** \name scaling factors
      @{ */
  double xscal, yscal; //@}

  gdImagePtr im; /**< the GD image */
};
