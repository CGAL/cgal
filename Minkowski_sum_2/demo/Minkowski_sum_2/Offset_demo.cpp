
/*
 * A simple demo that computes the Minkowski sum of two polygons.
 */
#include <CGAL/basic.h>

#ifndef CGAL_USE_QT

#include <iostream>
#include <fstream>

int main ()
{
  std::cout << "Sorry, this demo needs QT..." << std::endl; 
  return (0);
}

#else

#define RWRW_STATS

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/offset_polygon_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition.h>
#include <CGAL/Polygon_convex_decomposition.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Timer.h>
#include <list>
#include <iostream>
#include <fstream>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel,
                                 Alg_kernel,Nt_traits>  Conic_traits_2;

typedef Rat_kernel::Point_2                             Point_2;
typedef CGAL::Polygon_2<Rat_kernel>                     Polygon_2;

typedef CGAL::Gps_traits_2<Conic_traits_2>              Gps_traits_2;
typedef Gps_traits_2::Polygon_2                         Offset_polygon_2;
typedef Gps_traits_2::Polygon_with_holes_2   Offset_polygon_with_holes_2;

typedef CGAL::Bbox_2                                    Bbox_2;

// Define the QT types.
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <qapplication.h>
#include <qmainwindow.h>

/*!
 * Methods of computing the offset.
 */
enum Computation_method
{
  CONVOLUTION,
  DECOMP_SSAB,
  DECOMP_OPTIMAL,
  DECOMP_HM_APPROX,
  DECOMP_GREENE,

  INVALID_METHOD = -1
};
  
/*!
 * The main window.
 */
static const int INIT_WIDTH = 600;
static const int INIT_HEIGHT = 400;

class Offset_window : public QMainWindow
{
  Q_OBJECT

private:

  CGAL::Qt_widget*             widget;
  Polygon_2                    pgn;
  Rational                     r;
  Offset_polygon_with_holes_2  offset;

public:

  Offset_window (const int& x_min, const int& y_min,
                 const int& x_max, const int& y_max,
                 const Polygon_2& polygon,
                 const Rational& radius,
                 Computation_method method, int n_reps) :
    pgn (polygon),
    r (radius)
  {
    // Define auxiliary polygon-decomposition objects.
    CGAL::Small_side_angle_bisector_decomposition<Rat_kernel> ssab_decomp;
    CGAL::Optimal_convex_decomposition<Rat_kernel>            opt_decomp;
    CGAL::Hertel_Mehlhorn_convex_decomposition<Rat_kernel>    hm_approx_decomp;
    CGAL::Greene_convex_decomposition<Rat_kernel>             greene_decomp;

    // Compute the approximated offset.
    Conic_traits_2    traits;
    CGAL::Timer       timer;
    int               i;

    timer.start();

    for (i = 0; i < n_reps; i++)
    {
      switch (method)
      {
      case CONVOLUTION:
        offset = offset_polygon_2 (traits, pgn, r);
        break;

      case DECOMP_SSAB:
	offset = offset_polygon_2 (traits, pgn, r, ssab_decomp);
	break;

      case DECOMP_OPTIMAL:
	offset = offset_polygon_2 (traits, pgn, r, opt_decomp);
	break;

      case DECOMP_HM_APPROX:
	offset = offset_polygon_2 (traits, pgn, r, hm_approx_decomp);
	break;

      case DECOMP_GREENE:
	offset = offset_polygon_2 (traits, pgn, r, greene_decomp);
	break;

      default:
	std::cerr << "Method currently not supported!." << std::endl;
	exit(1);
      }
    }
    timer.stop();
    std::cout << "Exact offset computation took "
	      << timer.time() / n_reps << " seconds." << std::endl;

    // Create the window. 
    widget = new CGAL::Qt_widget(this);
    widget->resize (INIT_WIDTH, INIT_HEIGHT);
    widget->set_window (x_min, x_max, y_min, y_max);

    connect(widget, SIGNAL(redraw_on_back()),
	    this, SLOT(redraw_win()));

    setCentralWidget(widget);
  }

private slots:

  void _draw_offset_polygon (const Offset_polygon_2& off_pgn)
  {
    typedef CGAL::Cartesian<double>                     Approx_kernel;
    typedef Approx_kernel::Point_2                      Approx_point_2;
    typedef CGAL::Polygon_2<Approx_kernel>              Approx_polygon_2;

    Offset_polygon_2::Curve_const_iterator    iter;
    Approx_polygon_2                          app_pgn;

    Nt_traits       nts;
    const double    _r = CGAL::to_double (r);
    const double    _PI = 3.14159265;
    const double    skip = 0.02; // rad ~ 1.15 deg
    double          sx, sy;
    double          tx, ty;
    double          cx, cy;
    double          alpha1, alpha2;
    double          theta;
    double          x, y;
    int             n;
    int             k;

    for (iter = off_pgn.curves_begin(); iter != off_pgn.curves_end(); ++iter)
    { 
      sx = CGAL::to_double(iter->source().x());
      sy = CGAL::to_double(iter->source().y());

      app_pgn.push_back (Approx_point_2 (sx, sy));

      if (iter->orientation() != CGAL::COLLINEAR)
      {
	tx = CGAL::to_double(iter->target().x());
	ty = CGAL::to_double(iter->target().y());
	
	// Get the center of the supporting circle and compute the start and
	// end angles of the circular arc.
        cx = CGAL::to_double (nts.convert (iter->u()) / 
                              nts.convert(-2 * iter->r()));
	cy = CGAL::to_double (nts.convert (iter->v()) / 
                              nts.convert(-2 * iter->s()));

        alpha1 = atan2 (sy - cy, sx - cx);
        alpha2 = atan2 (ty - cy, tx - cx);

	if (alpha2 < alpha1)
	  alpha2 += 2*_PI;
	
	// Create intermediate points that approximate the circular arc.
        n = static_cast<int>((alpha2 - alpha1) / skip);

	for (k = 1; k <= n; k++)
	{
          theta = alpha1 + k*skip;

	  x = cx + _r*cos(theta);
	  y = cy + _r*sin(theta);

	  app_pgn.push_back (Approx_point_2 (x, y));
	}
      }
    }

    *widget << app_pgn;
    return;
  }

  void redraw_win()
  {
    // Draw the outer boundary.
    widget->setFilled (true);
    *widget << CGAL::FillColor(CGAL::RED);
    *widget << CGAL::LineWidth(1);

    *widget << CGAL::RED;
    _draw_offset_polygon (offset.outer_boundary());

    // Draw the holes in the offset.
    Offset_polygon_with_holes_2::Hole_const_iterator   pit;

    *widget << CGAL::FillColor(CGAL::WHITE);
    for (pit = offset.holes_begin(); pit != offset.holes_end(); ++pit)
    {
      _draw_offset_polygon (*pit);
    }

    // Draw the input polygon.
    widget->setFilled (false);

    *widget << CGAL::LineWidth(2);
    *widget << CGAL::BLUE;
    *widget << pgn;

    return;
  }

};

/*!
 * Read a polygon from an input file.
 * \param filename The name of the input file.
 * \param pgn Output: The polygon.
 * \return Whether the polygon was successfuly read.
 */
bool read_polygon (const char *filename, Polygon_2& pgn)
{
  // Open the input file.
  std::ifstream          ifile (filename);

  if (! ifile.is_open())
  {
    std::cerr << "Failed to open <" << filename << ">." << std::endl;
    return (false);
  }

  // Read the polygon.
  int                    n_vertices;
  Rational               x, y;
  std::list<Point_2>     vertices;
  int                    k;

  // Read the number of polygon vertices.
  ifile >> n_vertices;

  // Read the vertices.
  for (k = 0; k < n_vertices; k++)
  {
    ifile >> x >> y;

    vertices.push_back (Point_2 (x, y));
  }
  ifile.close();

  pgn = Polygon_2 (vertices.begin(), vertices.end());

  if (! pgn.is_simple())
  {
    std::cerr << "Error - the polygon is not simple." << std::endl;
    return (false);
  }

  return (true);
}

//moc_source_file : Offset_demo.C
#include "Offset_demo.moc"

/*!
 * The main.
 */
int main (int argc, char **argv )
{
  // Read the input file.
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] 
	      << " <polygon#1> <radius> <-c|-ds|-do|-dh|-dg> [n]." 
	      << std::endl;
    return (1);
  }

  // Read the polygon from the input file.
  Polygon_2       pgn;
  Bbox_2          bbox;
  
  if (! read_polygon (argv[1], pgn))
  {
    std::cerr << "Failed to read: <" << argv[1] << ">." << std::endl;
    return (1);
  }
  bbox = pgn.bbox();

  // Read the offset radius.
  int         numer, denom;
  Rational    radius;

  if (sscanf (argv[2], "%d/%d", &numer, &denom) != 2)
  {
    std::cerr << "Invalid radius: " << argv[2] << std::endl;
    return (1);
  }
 
  radius = Rational (numer, denom);

  // Get the computation method.
  Computation_method    method = INVALID_METHOD;

  if (strcmp(argv[3], "-c") == 0)
    method = CONVOLUTION;
  else if (strcmp(argv[3], "-ds") == 0)
    method = DECOMP_SSAB;
  else if (strcmp(argv[3], "-do") == 0)
    method = DECOMP_OPTIMAL;
  else if (strcmp(argv[3], "-dh") == 0)
    method = DECOMP_HM_APPROX;
  else if (strcmp(argv[3], "-dg") == 0)
    method = DECOMP_GREENE;
  
  if (method == INVALID_METHOD)
  {
    std::cerr << "Invalid method: " << argv[3]
	      << ", plaese specify: <-c|-ds|-do|-dh|-dg>." << std::endl;
    return (1);
  }
  
  // Get the number of repetitions.
  int                   n_reps = 1;

  if (argc > 4)
  {
    if (sscanf (argv[4], "%d", &n_reps) != 1)
      n_reps = 1;
  }

  // Create the main window.
  QApplication          app (argc, argv);
  Offset_window        *w = NULL;
  const double          rad = CGAL::to_double (radius);
  const double          x_min = bbox.xmin() - rad;
  const double          x_max = bbox.xmax() + rad; 
  const double          y_min = bbox.ymin() - rad;
  const double          y_max = bbox.ymax() + rad;
  
  w = new Offset_window (static_cast<int>(x_min - 2),
			 static_cast<int>(y_min - 2),
			 static_cast<int>(x_max + 2),
			 static_cast<int>(y_max + 2),
			 pgn, radius,
			 method, n_reps);
  app.setMainWidget (w);
  w->show();
  return (app.exec());
}

#endif
