#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Interval_skip_list.h>
#include <CGAL/Level_interval.h>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel EIK;
typedef EIK::Point_3                                   Point_3;
typedef CGAL::Projection_traits_xy_3<EIK> K;
typedef CGAL::Delaunay_triangulation_2<K>              Delaunay;
typedef Delaunay::Face_handle                          Face_handle;
typedef Delaunay::Finite_faces_iterator                Finite_faces_iterator;
typedef CGAL::Level_interval<Face_handle>              Interval;
typedef CGAL::Interval_skip_list<Interval>             Interval_skip_list;

int main()
{
  std::ifstream fin("terrain.pts"); // elevation ranges from 0 to 100
  Delaunay dt;

  dt.insert(std::istream_iterator<Point_3>(fin),
            std::istream_iterator<Point_3>());

  Interval_skip_list isl;
  for(Finite_faces_iterator fh = dt.finite_faces_begin();
      fh != dt.finite_faces_end();
      ++fh){
    isl.insert(Interval(fh));
  }
  std::list<Interval> level;
  isl.find_intervals(50, std::back_inserter(level));
  for(std::list<Interval>::iterator it = level.begin();
      it != level.end();
      ++it){
    std::cout << dt.triangle(it->face_handle()) << std::endl;
  }
  return 0;
}
