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

#define USE_LAZY_KERNEL

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/rotational_swept_area_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Timer.h>

#ifdef USE_LAZY_KERNEL
  typedef CGAL::Gmpq                                    Rational;
  typedef CGAL::Lazy_exact_nt<CGAL::Gmpq>               NT;
#else
  typedef CGAL::Gmpq                                    Rational;
  typedef CGAL::Gmpq                                    NT;
#endif

typedef CGAL::Cartesian<NT>                             Kernel;

#include <list>
#include <iostream>
#include <fstream>

typedef Kernel::Point_2                                 Point_2;
typedef CGAL::Polygon_2<Kernel>                         Polygon_2;

typedef CGAL::Gps_circle_segment_traits_2<Kernel>       Traits_2;
typedef Traits_2::Polygon_2                             Sv_polygon_2;

typedef CGAL::Bbox_2                                    Bbox_2;

// Define the QT types.
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <qapplication.h>
#include <qmainwindow.h>

/*!
 * The main window.
 */
static const int INIT_WIDTH = 600;
static const int INIT_HEIGHT = 400;

class RSA_window : public QMainWindow
{
  Q_OBJECT

private:

  CGAL::Qt_widget*             widget;
  Polygon_2                    pgn;
  Point_2                      pc;
  Sv_polygon_2                 sv_pgn;

public:

  RSA_window (const int& x_min, const int& y_min,
              const int& x_max, const int& y_max,
              const Polygon_2& polygon,
              const NT& sin_theta1, const NT& cos_theta1,
              const NT& sin_theta2, const NT& cos_theta2) :
    pgn (polygon)
  {
    // Locate the pair of vertices that are most distant from one another
    // and compute their midpoint as the center of rotation.
    Polygon_2::Vertex_const_iterator vi, vj;
    bool                             first_pair = true;
    Kernel                           ker;
    Kernel::Compute_squared_distance_2 dist_f = 
                                     ker.compute_squared_distance_2_object();
    NT                               dist, max_dist;
    Kernel::Construct_midpoint_2     mid_f = ker.construct_midpoint_2_object();

    for (vi = pgn.vertices_begin(); vi != pgn.vertices_end(); ++vi)
    {
      for (vj = vi, ++vj; vj != pgn.vertices_end(); ++vj)
      {
        dist = dist_f (*vi, *vj);
        if (first_pair || CGAL::compare (dist, max_dist) == CGAL::LARGER)
        {
          max_dist = dist;
          pc = mid_f (*vi, *vj);
        }
        first_pair = false;
      }
    }

    std::cout << "Rotating around (" << pc << ")." << std::endl;

    // Compute the boundary of the volume swept by the polygon when rotating
    // it from orientation theta1 to orientation theta2.
    CGAL::Timer       timer;

    timer.start();
    sv_pgn = rotational_swept_area_2 (pgn, pc,
                                      sin_theta1, cos_theta1,
                                      sin_theta2, cos_theta2);
    timer.stop();

    std::cout << "Swept-volume computation took "
	      << timer.time() << " seconds." << std::endl;

    // Create the window. 
    widget = new CGAL::Qt_widget(this);
    widget->resize (INIT_WIDTH, INIT_HEIGHT);
    widget->set_window (x_min, x_max, y_min, y_max);

    connect(widget, SIGNAL(redraw_on_back()),
	    this, SLOT(redraw_win()));

    setCentralWidget(widget);
  }

private slots:

  void redraw_win()
  {
    typedef CGAL::Cartesian<double>                     Approx_kernel;
    typedef Approx_kernel::Point_2                      Approx_point_2;
    typedef CGAL::Polygon_2<Approx_kernel>              Approx_polygon_2;

    // Draw the rotational swept volume.
    widget->setFilled (true);
    *widget << CGAL::FillColor(CGAL::RED);
    *widget << CGAL::LineWidth(1);
    *widget << CGAL::RED;

    Sv_polygon_2::Curve_const_iterator        iter;
    Approx_polygon_2                          app_pgn;

    const double    px = CGAL::to_double(pc.x());
    const double    py = CGAL::to_double(pc.y());
    const double    _PI = 3.14159265;
    const double    skip = 0.005; // rad
    double          sx, sy;
    double          tx, ty;
    double          cx, cy;
    double          alpha1, alpha2;
    double          theta;
    double          x, y;
    int             n;
    int             k;

    for (iter = sv_pgn.curves_begin(); iter != sv_pgn.curves_end(); ++iter)
    { 
      sx = CGAL::to_double(iter->source().x());
      sy = CGAL::to_double(iter->source().y());

      app_pgn.push_back (Approx_point_2 (sx, sy));

      if (iter->is_circular())
      {
        const double    _r = std::sqrt ((sx - px)*(sx - px) + 
                                        (sy - py)*(sy - py));

	tx = CGAL::to_double(iter->target().x());
	ty = CGAL::to_double(iter->target().y());
	
	// Get the center of the supporting circle and compute the start and
	// end angles of the circular arc.
	cx = CGAL::to_double(iter->supporting_circle().center().x());
	cy = CGAL::to_double(iter->supporting_circle().center().y());

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

    // Draw the input polygon.
    widget->setFilled (false);

    *widget << CGAL::LineWidth(2);
    *widget << CGAL::BLUE;
    *widget << pgn;

    // Draw the center of rotation.
    *widget << CGAL::BLACK;
    *widget << pc;

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

//moc_source_file : RSA_demo.C
#include "RSA_demo.moc"

/*!
 * The main.
 */
int main (int argc, char **argv )
{
  // Read the input file.
  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0] 
	      << " <polygon#1> <numer1>/<denom1> <numer2>/<denom2>." 
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
  int         numer1, denom1;
  int         numer2, denom2;
  NT          sin1, cos1, sin2, cos2;

  if (sscanf (argv[2], "%d/%d", &numer1, &denom1) != 2)
  {
    std::cerr << "Invalid radius: " << argv[2] << std::endl;
    return (1);
  }
  sin1 = NT (numer1) / NT (denom1);
  cos1 = NT (static_cast<int> (std::sqrt 
                               (static_cast<double> 
                                (denom1*denom1 - numer1*numer1)) + 0.001)) /
         NT (denom1);

  std::cout << "sin (theta1) = " << sin1
            << "  ,  cos(theta1) = " << cos1 << std::endl;

  if (sscanf (argv[3], "%d/%d", &numer2, &denom2) != 2)
  {
    std::cerr << "Invalid radius: " << argv[3] << std::endl;
    return (1);
  }
  sin2 = NT (numer2) / NT (denom2);
  cos2 = NT (static_cast<int> (std::sqrt 
                               (static_cast<double> 
                                (denom2*denom2 - numer2*numer2)) + 0.001)) /
         NT (denom2);
   
  std::cout << "sin (theta2) = " << sin2
            << "  ,  cos(theta2) = " << cos2 << std::endl;

  // Create the main window.
  QApplication          app (argc, argv);
  RSA_window           *w = NULL;
  const double          x_min = bbox.xmin();
  const double          x_max = bbox.xmax(); 
  const double          y_min = bbox.ymin();
  const double          y_max = bbox.ymax();
  
  w = new RSA_window (static_cast<int>(x_min - 2),
                      static_cast<int>(y_min - 2),
                      static_cast<int>(x_max + 2),
                      static_cast<int>(y_max + 2),
                      pgn,
                      sin1, cos1,
                      sin2, cos2);

  app.setMainWidget (w);
  w->show();
  return (app.exec());
}

#endif
