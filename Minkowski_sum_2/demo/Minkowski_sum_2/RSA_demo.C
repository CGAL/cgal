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

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/rotational_swept_area_2.h>
#include <CGAL/rsa_minkowski_sum_2.h>
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

typedef CGAL::Gps_circle_segment_traits_2<Rat_kernel>   Circ_seg_traits_2;
typedef Circ_seg_traits_2::Polygon_2                    Rsa_polygon_2;

typedef CGAL::Gps_traits_2<Conic_traits_2>              Gps_conic_traits_2;
typedef Gps_conic_traits_2::Polygon_2                   Sum_polygon_2;
typedef Gps_conic_traits_2::Polygon_with_holes_2        Sum_polygon_wh_2;

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
  Polygon_2                    pgn1;
  Polygon_2                    pgn2;
  Point_2                      pc;
  Rsa_polygon_2                rsa_pgn;
  Sum_polygon_wh_2             sum_pgn;

public:

  RSA_window (const int& x_min, const int& y_min,
              const int& x_max, const int& y_max,
              const Polygon_2& polygon1,
              const Polygon_2& polygon2,
              const Rational& sin_theta1, const Rational& cos_theta1,
              const Rational& sin_theta2, const Rational& cos_theta2) :
    pgn1 (polygon1),
    pgn2 (polygon2)
  {
    // Locate the pair of vertices that are most distant from one another
    // and compute their midpoint as the center of rotation.
    Polygon_2::Vertex_const_iterator vi, vj;
    bool                             first_pair = true;
    Rat_kernel                       ker;
    Rat_kernel::Compute_squared_distance_2
                              dist_f = ker.compute_squared_distance_2_object();
    Rational                         dist, max_dist;
    Rat_kernel::Construct_midpoint_2 mid_f = ker.construct_midpoint_2_object();

    for (vi = pgn2.vertices_begin(); vi != pgn2.vertices_end(); ++vi)
    {
      for (vj = vi, ++vj; vj != pgn2.vertices_end(); ++vj)
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

    // Compute the boundary of the area swept by the second polygon when
    // rotating it from orientation theta1 to orientation theta2.
    CGAL::Timer       timer;

    timer.start();
    rsa_pgn = rotational_swept_area_2 (pgn2, pc,
                                       sin_theta1, cos_theta1,
                                       sin_theta2, cos_theta2);
    timer.stop();

    std::cout << "Swept-volume computation took "
	      << timer.time() << " seconds." << std::endl;

    // Compute the Minkowski sum of this swept area with the first polygon.
    Conic_traits_2      traits;

    timer.reset();
    timer.start();
    sum_pgn = rsa_minkowski_sum_2 (traits,
                                   pgn1, rsa_pgn);
    timer.stop();

    std::cout << "Minkowski-sum computation took "
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
    widget->setFilled (false);
    //*widget << CGAL::FillColor(CGAL::RED);
    *widget << CGAL::LineWidth(3);
    *widget << CGAL::BLUE;

    Rsa_polygon_2::Curve_const_iterator       iter;
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

    for (iter = rsa_pgn.curves_begin(); iter != rsa_pgn.curves_end(); ++iter)
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

    // Draw the input polygons.
    widget->setFilled (false);

    *widget << CGAL::LineWidth(2);
    *widget << CGAL::BLACK;
    *widget << pgn1;

    *widget << CGAL::LineWidth(1);
    *widget << CGAL::BLUE;
    *widget << pgn2;

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
  if (argc < 5)
  {
    std::cerr << "Usage: " << argv[0] 
	      << " <polygon#1> <polygon#2>"
              << " <numer1>/<denom1> <numer2>/<denom2>." 
	      << std::endl;
    return (1);
  }

  // Read the polygons from the input files.
  Polygon_2       pgn1;
  Bbox_2          bbox1;
  
  if (! read_polygon (argv[1], pgn1))
  {
    std::cerr << "Failed to read: <" << argv[1] << ">." << std::endl;
    return (1);
  }
  bbox1 = pgn1.bbox();

  Polygon_2       pgn2;
  Bbox_2          bbox2;
  
  if (! read_polygon (argv[2], pgn2))
  {
    std::cerr << "Failed to read: <" << argv[2] << ">." << std::endl;
    return (1);
  }
  bbox2 = pgn2.bbox();

  // Read the offset radius.
  int         numer1, denom1;
  int         numer2, denom2;
  Rational    sin1, cos1, sin2, cos2;

  if (sscanf (argv[3], "%d/%d", &numer1, &denom1) != 2)
  {
    std::cerr << "Invalid radius: " << argv[2] << std::endl;
    return (1);
  }
  sin1 = Rational (numer1) / Rational (denom1);
  cos1 = 
    Rational (static_cast<int>
              (std::sqrt (static_cast<double> 
                          (denom1*denom1 - numer1*numer1)) + 0.001)) /
    Rational (denom1);

  std::cout << "sin (theta1) = " << sin1
            << "  ,  cos(theta1) = " << cos1 << std::endl;

  if (sscanf (argv[4], "%d/%d", &numer2, &denom2) != 2)
  {
    std::cerr << "Invalid radius: " << argv[3] << std::endl;
    return (1);
  }
  sin2 = Rational (numer2) / Rational (denom2);
  cos2 = 
    Rational (static_cast<int>
              (std::sqrt (static_cast<double> 
                          (denom2*denom2 - numer2*numer2)) + 0.001)) /
    Rational (denom2);
   
  std::cout << "sin (theta2) = " << sin2
            << "  ,  cos(theta2) = " << cos2 << std::endl;

  // Create the main window.
  QApplication          app (argc, argv);
  RSA_window           *w = NULL;
  const double          x_min = bbox1.xmin() + 
    ((bbox2.xmin() < 0) ? bbox2.xmin() : 0);
  const double          x_max = bbox1.xmax() +
    ((bbox2.xmax() > 0) ? bbox2.xmax() : 0);
  const double          y_min = bbox1.ymin() + 
    ((bbox2.ymin() < 0) ? bbox2.ymin() : 0);
  const double          y_max = bbox1.ymax() +
    ((bbox2.ymax() > 0) ? bbox2.ymax() : 0);
  
  w = new RSA_window (static_cast<int>(x_min - 2),
                      static_cast<int>(y_min - 2),
                      static_cast<int>(x_max + 2),
                      static_cast<int>(y_max + 2),
                      pgn1, pgn2,
                      sin1, cos1,
                      sin2, cos2);

  app.setMainWidget (w);
  w->show();
  return (app.exec());
}

#endif
