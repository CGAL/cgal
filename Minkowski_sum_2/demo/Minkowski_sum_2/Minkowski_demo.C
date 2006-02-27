/*
 * A simple demo that computes the Minkowski sum of two polygons.
 */
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
#define USE_LAZY_KERNEL

#include <CGAL/basic.h>

#ifdef USE_LAZY_KERNEL

/*
  #include <CGAL/Simple_cartesian.h>
  #include <CGAL/Lazy_kernel.h>

  typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<Rational> > Kernel;
*/
  #include <CGAL/Exact_predicates_exact_constructions_kernel.h>

  typedef CGAL::Exact_predicates_exact_constructions_kernel  Kernel;
  typedef Kernel::FT                                         Rational;

#else

  #include <CGAL/Cartesian.h>
  #include <CGAL/Gmpq.h>

  typedef CGAL::Gmpq                                         Rational;
  typedef CGAL::Cartesian<Rational>                          Kernel;

#endif

#include <CGAL/Bbox_2.h>
#include <CGAL/Timer.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition.h>
#include <CGAL/Polygon_convex_decomposition.h>

#include <list>
#include <iostream>
#include <fstream>

typedef Kernel::Point_2                             Point_2;
typedef Kernel::Segment_2                           Segment_2;
typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;
typedef CGAL::Bbox_2                                Bbox_2;

// Define the QT types.
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <qapplication.h>
#include <qmainwindow.h>

/*!
 * Methods of computing the Minkowski sum.
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

class Minkowski_sum_window : public QMainWindow
{
  Q_OBJECT

private:

  CGAL::Qt_widget*          widget;
  Polygon_2                 pgn1;
  Polygon_2                 pgn2;
  Polygon_with_holes_2      sum;
  double                    arrowhead_w;
  double                    arrowhead_h;

public:

  Minkowski_sum_window (const int& x_min, const int& y_min,
			const int& x_max, const int& y_max,
			const Polygon_2& p1, const Polygon_2& p2,
			Computation_method method, int n_reps) :
    pgn1(p1),
    pgn2(p2)
  {
    // Define auxiliary polygon-decomposition objects.
    CGAL::Small_side_angle_bisector_decomposition<Kernel>  ssab_decomp;
    CGAL::Optimal_convex_decomposition<Kernel>             opt_decomp;
    CGAL::Hertel_Mehlhorn_convex_decomposition<Kernel>     hm_approx_decomp;
    CGAL::Greene_convex_decomposition<Kernel>              greene_decomp;

    // Compute the Minkowski sum.
    CGAL::Timer       timer;
    int               i;

    timer.start();

    for (i = 0; i < n_reps; i++)
    {
      switch (method)
      {
      case CONVOLUTION:
	sum = minkowski_sum_2 (pgn1, pgn2);
	break;

      case DECOMP_SSAB:
	sum = minkowski_sum_2 (pgn1, pgn2, ssab_decomp);
	break;

      case DECOMP_OPTIMAL:
	sum = minkowski_sum_2 (pgn1, pgn2, opt_decomp);
	break;

      case DECOMP_HM_APPROX:
	sum = minkowski_sum_2 (pgn1, pgn2, hm_approx_decomp);
	break;

      case DECOMP_GREENE:
	sum = minkowski_sum_2 (pgn1, pgn2, greene_decomp);
	break;

      default:
	std::cerr << "Method currently not supported!." << std::endl;
	exit(1);
      }
    }
    timer.stop();
    std::cout << "Minkowski sum computation took "
	      << timer.time() / n_reps << " seconds." << std::endl;

    // Create the window. 
    widget = new CGAL::Qt_widget(this);
    widget->resize (INIT_WIDTH, INIT_HEIGHT);
    widget->set_window (x_min, x_max, y_min, y_max);

    connect(widget, SIGNAL(redraw_on_back()),
	    this, SLOT(redraw_win()));

    setCentralWidget(widget);

    double    scale = static_cast<double>(x_max - x_min) / INIT_WIDTH;
    if (static_cast<double>(y_max - y_min) / INIT_HEIGHT < scale)
      scale = static_cast<double>(y_max - y_min) / INIT_HEIGHT;

    arrowhead_w = 8*scale;
    arrowhead_h = 4*scale;
  }

private slots:
 
  void redraw_win()
  {
    // Draw the Minkowski sum (filled, then draw the holes).
    widget->setFilled (true);
    *widget << CGAL::FillColor(CGAL::RED);
    *widget << CGAL::LineWidth(1);

    *widget << CGAL::RED;
    *widget << sum.outer_boundary();

    Polygon_with_holes_2::Hole_const_iterator   pit;

    *widget << CGAL::FillColor(CGAL::WHITE);
    for (pit = sum.holes_begin(); pit != sum.holes_end(); ++pit)
      *widget << *pit;

    // Draw the input polygons.
    widget->setFilled (false);

    *widget << CGAL::LineWidth(2);
    *widget << CGAL::GREEN;
    *widget << pgn1;

    *widget << CGAL::LineWidth(2);
    *widget << CGAL::BLUE;
    *widget << pgn2;

    return;
  }

private:

  void _draw_arrow (const Segment_2& seg,
		    const double& arrowhead_width,
		    const double& arrowhead_height) const
  {
    const double   x1 = CGAL::to_double (seg.source().x());
    const double   y1 = CGAL::to_double (seg.source().y());
    const double   x2 = CGAL::to_double (seg.target().x());
    const double   y2 = CGAL::to_double (seg.target().y());
    const double   theta = ::atan2 (y2 - y1, x2 - x1);
    const double   alpha = 3.14159265/2 + ::atan2 (arrowhead_width, 
						   arrowhead_height);
    const double   alen = ::sqrt(arrowhead_width*arrowhead_width + 
				 arrowhead_height*arrowhead_height);
    const double   ax1 = x2 + alen * ::cos(theta + alpha);
    const double   ay1 = y2 + alen * ::sin(theta + alpha);
    const double   ax2 = x2 + alen * ::cos(theta - alpha);
    const double   ay2 = y2 + alen * ::sin(theta - alpha);
    
    widget->get_painter().moveTo (widget->x_pixel(x1), widget->y_pixel(y1));
    widget->get_painter().lineTo (widget->x_pixel(x2), widget->y_pixel(y2));
    widget->get_painter().lineTo (widget->x_pixel(ax1), widget->y_pixel(ay1));
    widget->get_painter().lineTo (widget->x_pixel(ax2), widget->y_pixel(ay2));
    widget->get_painter().lineTo (widget->x_pixel(x2), widget->y_pixel(y2));

    return;
  }
};

/*!
 * Read a polygons from an input file.
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

//moc_source_file : Minkowski_demo.C
#include "Minkowski_demo.moc"

/*!
 * The main.
 */
int main (int argc, char **argv )
{
  // Read the input file.
  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0] 
	      << " <polygon#1> <polygon#2> <-c|-ds|-do|-dh|-dg> [n]." 
	      << std::endl;
    return (1);
  }

  // Read the polygons from the input files.
  Polygon_2   pgn1, pgn2;
  Bbox_2      bbox1, bbox2;
  
  if (! read_polygon (argv[1], pgn1))
  {
    std::cerr << "Failed to read: <" << argv[1] << ">." << std::endl;
    return (1);
  }
  bbox1 = pgn1.bbox();
  
  if (! read_polygon (argv[2], pgn2))
  {
    std::cerr << "Failed to read: <" << argv[2] << ">." << std::endl;
    return (1);
  }
  bbox2 = pgn2.bbox();

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
  Minkowski_sum_window *w = NULL;
  const double          x_min = bbox1.xmin() + 
    ((bbox2.xmin() < 0) ? bbox2.xmin() : 0);
  const double          x_max = bbox1.xmax() +
    ((bbox2.xmax() > 0) ? bbox2.xmax() : 0);
  const double          y_min = bbox1.ymin() + 
    ((bbox2.ymin() < 0) ? bbox2.ymin() : 0);
  const double          y_max = bbox1.ymax() +
    ((bbox2.ymax() > 0) ? bbox2.ymax() : 0);
  
  w = new Minkowski_sum_window (static_cast<int>(x_min - 2),
				static_cast<int>(y_min - 2),
				static_cast<int>(x_max + 2),
				static_cast<int>(y_max + 2),
				pgn1, pgn2,
				method, n_reps);
  app.setMainWidget (w);
  w->show();
  return (app.exec());
}

#endif
