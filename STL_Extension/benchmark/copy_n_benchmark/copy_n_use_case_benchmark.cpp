#include <cstdio>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>

#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/compiler_config.h>



int format_output(const char* lib, const char* generator, int n, float time) {
  return std::printf("| %s || %s || %d || %.4fM items/sec\n", lib, generator, n, time);
}

typedef CGAL::Simple_cartesian<double>         R;
typedef R::Point_2                             Point;
typedef CGAL::Creator_uniform_2<double,Point>  Creator;
typedef std::vector<Point>                     Vector;

int main(int argc, char* argv[]) {
  int n = 10000000;
  int repeats = 10;

  if(argc > 1)
    n = boost::lexical_cast<int>(argv[1]);

  if(argc > 2)
    repeats = boost::lexical_cast<int>(argv[2]);

  Vector points(n, Point());

  CGAL::Random_points_in_disc_2<Point,Creator> g( 1000.0);
  boost::timer timer;
  const char* generator = "Random_points_in_disc_2";
  float time;

  std::cout <<
    "{| \n"
    "! Library !! Generator !! #Elements !! items/sec \n"
    "|- \n";

  timer.restart();
  for (int i = 0; i < repeats; ++i) { CGAL::copy_n( g, n, points.begin()); }
  time = (double)n*repeats/timer.elapsed()/1.0E6;
  format_output("CGAL", generator , n, time);
  std::cout << "|- \n";

  timer.restart();
  for (int i = 0; i < repeats; ++i) { std::copy_n( g, n, points.begin()); }
  time = (double)n*repeats/timer.elapsed()/1.0E6;
  format_output("stdlib", generator, n, time);

  //wiki markup footer
  std::cout << "|}" << std::endl;
}
