// ============================================================================
//
// Copyright (c) 1998, 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : rectangular_p_center_2_demo.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Demo: 2-4-Centering Axis-Parallel 2D-Rectangles
// ============================================================================


#ifdef CGAL_USE_LEDA

#include <CGAL/Cartesian.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/IO/leda_window.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/Istream_iterator.h>
#include <CGAL/rectangular_p_center_2.h>
#include <CGAL/point_generators_2.h>
//#include <CGAL/Arithmetic_filter.h>
#include <CGAL/leda_real.h>
#include <CGAL/algorithm.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cstdlib>

#ifndef _MSC_VER
using std::atoi;
using std::atof;
#endif
using std::vector;
using std::copy;
using std::back_inserter;
using std::ostream_iterator;
using std::transform;
using std::bind2nd;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using CGAL::Cartesian;
using CGAL::rectangular_p_center_2;
using CGAL::Random;
using CGAL::default_random;
using CGAL::Random_points_in_square_2;
using CGAL::Random_points_in_disc_2;
using CGAL::Istream_iterator;
using CGAL::Ostream_iterator;
using CGAL::set_pretty_mode;
using CGAL::set_ascii_mode;
using CGAL::cgalize;
using CGAL::Timer;
using CGAL::BLUE;
using CGAL::RED;
using CGAL::ORANGE;
using CGAL::GREEN;

#ifdef _MSC_VER
// that compiler cannot even distinguish between global
// and class scope, so ...
#define Base B_B_Base
#endif // _MSC_VER

#ifdef CGAL_CFG_NO_NAMESPACE
template < class P,
           class Creator =
           CGAL_STD::Creator_uniform_2< typename P::FT, P > >
class Random_p_clusters_2 : public CGAL_STD::Random_generator_base< P > {
#else
template < class P,
           class Creator =
           CGAL::Creator_uniform_2< CGAL_TYPENAME_MSVC_NULL P::FT, P > >
class Random_p_clusters_2 : public CGAL::Random_generator_base< P > {
#endif // CGAL_CFG_NO_NAMESPACE
  void generate_point() {
    typedef typename P::FT FT;
    double p = _rnd.get_double();
    Creator creator;
    if (p <= 1.0 / n)
      d_item =
        creator(FT(p0.x() + c_size * (2 * _rnd.get_double() - 1.0)),
                FT(p0.y() + c_size * (2 * _rnd.get_double() - 1.0)));
    else if (p <= 2.0 / n)
      d_item =
        creator(FT(p1.x() + c_size * (2 * _rnd.get_double() - 1.0)),
                FT(p1.y() + c_size * (2 * _rnd.get_double() - 1.0)));
    else if (p <= 3.0 / n)
      d_item =
        creator(FT(p2.x() + c_size * (2 * _rnd.get_double() - 1.0)),
                FT(p2.y() + c_size * (2 * _rnd.get_double() - 1.0)));
    else
      d_item =
        creator(FT(p3.x() + c_size * (2 * _rnd.get_double() - 1.0)),
                FT(p3.y() + c_size * (2 * _rnd.get_double() - 1.0)));
  }
public:
  typedef Random_p_clusters_2< P, Creator > This;
  typedef CGAL::Random_generator_base< P > Base;
  Random_p_clusters_2(int n_,
                      double c_size_,
                      double r = 1,
                      Random& rnd = default_random)
  : Base(r - c_size_, rnd),
    n(n_),
    c_size(c_size_),
    p0(Creator()(d_range * (2 * _rnd.get_double() - 1.0),
                 d_range * (2 * _rnd.get_double() - 1.0))),
    p1(Creator()(d_range * (2 * _rnd.get_double() - 1.0),
                 d_range * (2 * _rnd.get_double() - 1.0))),
    p2(Creator()(d_range * (2 * _rnd.get_double() - 1.0),
                 d_range * (2 * _rnd.get_double() - 1.0))),
    p3(Creator()(d_range * (2 * _rnd.get_double() - 1.0),
                 d_range * (2 * _rnd.get_double() - 1.0)))
  {
    CGAL_precondition(n >= 1 && n <= 4);
    CGAL_precondition(c_size >= 0 && c_size <= r);
    generate_point();
  }
  This& operator++() {
    generate_point();
    return *this;
  }
  This  operator++(int) {
    This tmp = *this;
    ++(*this);
    return tmp;
  }
private:
  int n;
  double c_size;
  P p0, p1, p2, p3;
};

#ifdef _MSC_VER
#undef Base
#endif // _MSC_VER

// typedefs
//typedef Filtered_exact< double, leda_real > FT;
//typedef double FT;
typedef leda_real                           FT;
typedef Cartesian< FT >                     R;
typedef CGAL::Point_2< R >                  Point;
typedef CGAL::Iso_rectangle_2< R >          Square_2;
typedef vector< Point >                     Point_cont;
typedef Point_cont::iterator                iterator;
typedef Random_points_in_square_2< Point >  Point_generator_square;
typedef Random_points_in_disc_2< Point >    Point_generator_disc;
typedef Random_p_clusters_2< Point >        Point_generator_cluster;
typedef Ostream_iterator< Point, leda_window >
  Window_stream_iterator_point;
typedef Ostream_iterator< Square_2, leda_window >
  Window_stream_iterator_square;
typedef ostream_iterator< Point >       Ostream_iterator_point;
typedef ostream_iterator< Square_2 >    Ostream_iterator_square;
typedef Istream_iterator< Point, leda_window>
  Istream_iterator_point;



// function class to construct a box
// around a point p with radius r
template < class Point, class FT, class Box >
struct Build_box
: public CGAL_STD::binary_function< Point, FT, Box >
{
  Box
  operator()(const Point& p, const FT& r) const
  {
    return Box(Point(p.x() - r, p.y() - r),
               Point(p.x() + r, p.y() + r));
  }
};

int
main(int argc, char* argv[])
{
  int number_of_points = 50;
  int number_of_clusters = 4;
  // for the cluster generators:
  double c_size = .5;     // cluster size
  Point_cont input_points;
  // indicate whether to end the interactive part
#ifndef CGAL_PCENTER_NO_SHOW
  bool done = false;
#endif // CGAL_PCENTER_NO_SHOW
  typedef Build_box< Point, FT, Square_2 >  Build_square;

#ifndef CGAL_PCENTER_NO_SHOW
  // init CGAL stuff:
  leda_window W(730, 690);
  cgalize(W);
  W.set_node_width(2);
  W.init(-1.5, 1.5, -1.2);
  int gensq_button   = W.button("Square",
                                "Generate points from the unit square.");
  int gendsk_button  = W.button("Disc",
                                "Generate points from the unit disc.");
  int gencl_button   = W.button("Cluster",
                                "Generate points from p clusters.");
  int compute_button = W.button("Compute", "Compute the p-centers.");
  int clear_button   = W.button("Clear", "Clear window.");
  int ps_button      = W.button("PS", "Generate postscript output.");
  int ascii_button   = W.button("ASCII",
                                "Generate ascii output of input points.");
  int end_button     = W.button("Quit", "Leave the program.");
  W.int_item("n", number_of_points, "Number of points.");
  W.int_item("p", number_of_clusters, "Number of clusters.");
  W.double_item("Cluster Size", c_size,
                "Size of the clusters (relevant for the cluster generator).");
  W.display();
  Window_stream_iterator_point  wout_p(W);
  Window_stream_iterator_square wout_s(W);
#endif // CGAL_PCENTER_NO_SHOW
  set_pretty_mode(cout);
  set_pretty_mode(cerr);
  set_ascii_mode(cin);
  cout.precision(10);
  cerr.precision(10);
  Ostream_iterator_point        cout_p(cerr, "\n");
  Ostream_iterator_square       cout_s(cerr, "\n");

  if (argc >= 3 && (number_of_points = atoi(argv[1])) > 0) {
    int random_seed(default_random.get_int(0, (1 << 31)));
    cerr << "3CENTER";
#ifndef CGAL_3COVER_NO_PREFILTER
    cerr << "-PREFILTER";
#endif // CGAL_3COVER_NO_PREFILTER
#ifndef CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST
    cerr << "-CHECK";
#endif // CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST
    cerr << " with " << number_of_points << " points ";

    if (argv[2][0] == 'c') {
      // random points in three clusters
      if (argc >= 4)
        c_size = atof(argv[3]);

      cerr << "in clusters of size " << c_size << endl;

      if (argc >= 5)
        // get seed from command line
        random_seed = atoi(argv[4]);
      Random my_rnd(random_seed);

      Point_generator_cluster gen(number_of_clusters, c_size, 1.0, my_rnd);
      CGAL::copy_n(gen, number_of_points, back_inserter(input_points));

    } else {
      if (argc >= 4)
        // get seed from command line
        random_seed = atoi(argv[3]);
      Random my_rnd(random_seed);

      if (argv[2][0] == 'd') {
        // random points from the unit disc
        cerr << "from the unit disc\n";
        Point_generator_disc gen(1.0, my_rnd);
        CGAL::copy_n(gen, number_of_points, back_inserter(input_points));
      } else {
        // random points from the unit square
        cerr << "from the unit square\n";
        Point_generator_square gen(1.0, my_rnd);
        CGAL::copy_n(gen, number_of_points, back_inserter(input_points));
      }
    }
    cerr << "random seed is " << random_seed << endl;
  } else if (argc >= 2 && argv[1][0] == 'i') {
#ifndef _MSC_VER
    // read parameters
    int p;
    cin >> p;
    int n;
    cin >> n;
    cout << "Computing " << p << "-centers of " << n << " points." << endl;
    
    // read input points
    for (int i = 0; i < n; ++i) {
      Point p;
      cin >> p;
      input_points.push_back(p);
    }
    number_of_points = n;
    number_of_clusters = p;
#endif
  }

#ifndef CGAL_PCENTER_NO_SHOW
  do {
    // show point set:
    W.clear();
    W << BLUE;
    copy(input_points.begin(), input_points.end(), wout_p);
#endif // CGAL_PCENTER_NO_SHOW

    FT result;
    Point_cont centers;
    if (!input_points.empty()) {
#ifndef CGAL_PCENTER_NO_SHOW
      W << GREEN
        << CGAL::bounding_box_2(input_points.begin(), input_points.end());
#endif // CGAL_PCENTER_NO_SHOW
      Timer t;
      t.start();
      rectangular_p_center_2(
        input_points.begin(),
        input_points.end(),
        back_inserter(centers),
        result,
        number_of_clusters);
      t.stop();
      cout << "[time: " << t.time() << " sec]" << endl;
      int number_of_piercing_points(centers.size());
      cerr << "Finished with diameter " << CGAL::to_double(result)
           << " and " << number_of_piercing_points
           << " points." << endl;
#ifdef CGAL_PCENTER_CHECK
      CGAL::Infinity_distance_2< R > dist;
      for (iterator i = input_points.begin(); i != input_points.end(); ++i) {
        iterator j = centers.begin();
        do {
          if (dist(*i, *j) <= result / FT(2))
            break;
          if (++j == centers.end()) {
            cerr << "!!Problemo: " << *i << endl;
#ifndef CGAL_PCENTER_NO_SHOW
            W.set_node_width(3);
            W << CGAL::VIOLET << *i;
#endif // CGAL_PCENTER_NO_SHOW
          }
        } while (j != centers.end());
      }
#endif // CGAL_PCENTER_CHECK

#ifndef CGAL_PCENTER_NO_SHOW
      // show center points
      W << RED;
      copy(centers.begin(), centers.end(), wout_p);
#endif // CGAL_PCENTER_NO_SHOW
#ifndef _MSC_VER
      copy(centers.begin(), centers.end(), cout_p);
      cerr << endl;
#endif // _MSC_VER

      // ... and the corresponding squares:
#ifndef CGAL_PCENTER_NO_SHOW
      W << ORANGE;
      transform(centers.begin(),
                centers.end(),
                wout_s,
                bind2nd(Build_square(), result / FT(2)));
#endif // CGAL_PCENTER_NO_SHOW
#ifndef _MSC_VER
      transform(centers.begin(),
                centers.end(),
                cout_s,
                bind2nd(Build_square(), result / FT(2)));
      cerr << endl;
#endif // _MSC_VER

#ifdef CGAL_PCENTER_CHECK
      // check that there is at least one square with two points
      // on opposite sides
      CGAL::Signed_x_distance_2< R > xdist;
      CGAL::Signed_y_distance_2< R > ydist;
      bool boundary = false;
      for (iterator i = centers.begin(); i != centers.end(); ++i) {
        int left = 0, right = 0, bottom = 0, top = 0;
        for (iterator j = input_points.begin(); j != input_points.end(); ++j) {
          if (xdist(*i, *j) == result / FT(2)) {
            ++left;
#ifndef CGAL_PCENTER_NO_SHOW
            W << CGAL::GREEN << *j;
#endif // CGAL_PCENTER_NO_SHOW
          }
          if (xdist(*j, *i) == result / FT(2)) {
            ++right;
#ifndef CGAL_PCENTER_NO_SHOW
            W << CGAL::GREEN << *j;
#endif // CGAL_PCENTER_NO_SHOW
          }
          if (ydist(*j, *i) == result / FT(2)) {
            ++top;
#ifndef CGAL_PCENTER_NO_SHOW
            W << CGAL::GREEN << *j;
#endif // CGAL_PCENTER_NO_SHOW
          }
          if (ydist(*i, *j) == result / FT(2)) {
            ++bottom;
#ifndef CGAL_PCENTER_NO_SHOW
            W << CGAL::GREEN << *j;
#endif // CGAL_PCENTER_NO_SHOW
          }
        }
        if (left > 0 && right > 0 || top > 0 && bottom > 0)
          boundary = true;
      }
      if (!boundary)
        cerr << "Error: No square has two points on boundary." << endl;

#endif // CGAL_PCENTER_CHECK
    } // if (!input_points.empty())

#ifndef CGAL_PCENTER_NO_SHOW
    double x, y;
    W << BLUE;
    do {
      int input = W.get_mouse(x, y);
      if (input == gensq_button) {
        // random points from the unit square
        Point_generator_square gen(1.0, default_random);
        Point_cont tmpc;
        CGAL::copy_n(gen, number_of_points, back_inserter(tmpc));
        copy(tmpc.begin(), tmpc.end(), back_inserter(input_points));
        copy(tmpc.begin(), tmpc.end(), wout_p);
      } else if (input == gendsk_button) {
        // random points from the unit disc
        Point_generator_disc gen(1.0, default_random);
        Point_cont tmpc;
        CGAL::copy_n(gen, number_of_points, back_inserter(tmpc));
        copy(tmpc.begin(), tmpc.end(), back_inserter(input_points));
        copy(tmpc.begin(), tmpc.end(), wout_p);
      } else if (input == clear_button) {
        // clear point set
        input_points.clear();
        W.clear();
      } else if (input == gencl_button) {
        // random points in three clusters
        Point_generator_cluster gen(number_of_clusters,
                                    c_size,
                                    1.0,
                                    default_random);
        Point_cont tmpc;
        CGAL::copy_n(gen, number_of_points, back_inserter(tmpc));
        copy(tmpc.begin(), tmpc.end(), back_inserter(input_points));
        copy(tmpc.begin(), tmpc.end(), wout_p);
      } else if (input == MOUSE_BUTTON(1)) {
        Point p(x, y);
        input_points.push_back(p);
        W << p;
      } else if (input == compute_button || input == MOUSE_BUTTON(2)) {
        break;
      } else if (input == ascii_button) {
        set_ascii_mode(cerr);
        cerr << number_of_clusters << endl;
        cerr << input_points.size() << endl;
        for (iterator i = input_points.begin();
             i != input_points.end();
             ++i)
          cerr << CGAL::to_double(i->x()) << " "
               << CGAL::to_double(i->y()) << "\n";
        cerr << CGAL::to_double(result) << "\n";
#ifdef _MSC_VER
        {
#endif
        for (iterator i = centers.begin(); i != centers.end(); ++i)
          cerr << CGAL::to_double(i->x()) << " "
               << CGAL::to_double(i->y()) << endl;
#ifdef _MSC_VER
        }
#endif
        set_pretty_mode(cerr);
      } else if (input == ps_button) {
        iterator xmin = std::min_element(input_points.begin(),
                                         input_points.end(),
                                         CGAL::Less_x_2< R >());
        iterator xmax = std::max_element(input_points.begin(),
                                         input_points.end(),
                                         CGAL::Less_x_2< R >());
        iterator ymin = std::min_element(input_points.begin(),
                                         input_points.end(),
                                         CGAL::Less_y_2< R >());
        iterator ymax = std::max_element(input_points.begin(),
                                         input_points.end(),
                                         CGAL::Less_y_2< R >());
        FT scale = std::max(xmax->x() - xmin->x(), ymax->y() - ymin->y());
        const int size = 500;
        const int border = 20;
        cout << "%!PS-Adobe-2.0 EPSF-1.2\n"
             << "%%Creator: rectangular_p_center_2_demo\n"
             << "%%BoundingBox: 0 0 " << size << " " << size << "\n"
             << "%%Pages: 1\n"
             << "%%EndComments\n"
             << "\n"
             << "%% global scaling factor:\n"
             << "1 1 scale\n"
             << "\n"
             << "%% Procedures:\n"
             << "\n"
             << "% square takes three arguments a b c from stack and draws\n"
             << "% a square of radius c centered at point (a,b)\n"
             << "/square\n"
             << "{ newpath 3 1 roll moveto \n"
             << "  dup 2 div neg dup rmoveto \n"
             << "  dup 0 rlineto \n"
             << "  dup 0 exch rlineto \n"
             << "  neg 0 rlineto \n"
             << "  closepath stroke} def\n"
             << "\n"
             << "% circle takes three arguments a b c from stack and draws\n"
             << "% a filled circle of radius c centered at point (a,b)\n"
             << "/circle {\n"
             << "newpath 0 360 arc closepath gsave fill grestore stroke\n"
             << "} def\n\n% input points\n0 0 0 setrgbcolor\n";
        for (iterator i = input_points.begin();
             i != input_points.end();
             ++i)
          {
            cout << int(CGAL::to_double(
                        FT(.5 + border) + (i->x() - xmin->x()) *
                        FT(size - 2.0 * border) / scale))
                 << " "
                 << int(CGAL::to_double(
                        FT(.5 + border) + (i->y() - ymin->y()) *
                        FT(size - 2.0 * border) / scale))
                 << " 2 circle\n";
          }
          cout << "\n% covering squares:\n";
#ifdef _MSC_VER
        {
#endif
        for (iterator i = centers.begin();
             i != centers.end();
             ++i)
          {
            cout << "1 0 0 setrgbcolor\n"
                 << int(CGAL::to_double(
                        FT(.5 + border) + (i->x() - xmin->x()) *
                        FT(size - 2.0 * border) / scale))
                 << " "
                 << int(CGAL::to_double(
                        FT(.5 + border) + (i->y() - ymin->y()) *
                        FT(size - 2.0 * border) / scale))
                 << " 2 circle\n";
            cout << "0 0 1 setrgbcolor\n"
                 << int(CGAL::to_double(
                        FT(.5 + border) + (i->x() - xmin->x()) *
                        FT(size - 2.0 * border) / scale))
                 << " "
                 << int(CGAL::to_double(
                        FT(.5 + border) + (i->y() - ymin->y()) *
                        FT(size - 2.0 * border) / scale))
                 << " "
                 << int(CGAL::to_double(
                        FT(.5) + result * (size - 2.0 * border) / scale))
                 << " square\n";
          }
          cout << "\nshowpage\n\n%% EOF\n" << endl;
#ifdef _MSC_VER
        }
#endif
      } else if (input == end_button || input == MOUSE_BUTTON(3))
        done = true;
    } while (!done);
  } while (!done);
#endif // CGAL_PCENTER_NO_SHOW

  return 0;
}

#else

#include <iostream>

int main()
{
  std::cerr << "This demo requires LEDA." << std::endl;
  return 0;
}

#endif

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

