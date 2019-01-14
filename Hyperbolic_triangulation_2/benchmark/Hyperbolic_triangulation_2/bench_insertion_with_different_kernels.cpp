#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>

#include <iostream>

template<typename Point_type>
std::pair<double, double> insert_Euclidean(const std::vector<Point_type>& dpts)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Delaunay_triangulation_2< K >                 Dt;
  typedef Dt::Point_2                                         Point_2;

  CGAL::Timer timer;

  std::vector<Point_2> pts;

  for(std::size_t i=0; i<dpts.size(); ++i)
    pts.push_back(Point_2(dpts[i].x(), dpts[i].y()));

  Dt dt_end;
  timer.start();
  dt_end.insert(pts.begin(),pts.end());
  timer.stop();
  double t1 = timer.time();
  std::cout << "Euclidean   insertion, iterator input: " << t1 << std::endl;

  Dt dt_during;
  timer.reset();
  timer.start();
  for(std::size_t i=0; i<pts.size(); ++i)
    dt_during.insert(pts[i]);

  timer.stop();
  double t2 = timer.time();
  std::cout << "Euclidean   insertion, one-by-one:     " << t2 << std::endl;

  return std::pair<double, double>(t1, t2);
}

template<typename NT_point>
std::pair<double,double> insert_CK_points(const std::vector<NT_point>& dpts)
{
  typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<>   Gt;
  typedef Gt::Point_2                                             Point_2;
  typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Gt>           Dt;

  CGAL::Timer timer;
  std::vector<Point_2> pts;
  std::vector<Point_2>::iterator ip;

  for(std::size_t i=0; i<dpts.size(); ++i)
    pts.push_back(Point_2(dpts[i].x(), dpts[i].y()));

  Dt dt_end;
  timer.start();
  dt_end.insert(pts.begin(),pts.end());
  timer.stop();
  double t1 = timer.time();
  std::cout << "CK   traits insertion, iterator input: " << t1 << std::endl;

  Dt dt_during;
  timer.reset();
  timer.start();
  for(ip = pts.begin(); ip != pts.end(); ++ip)
    dt_during.insert(*ip);

  timer.stop();
  double t2 = timer.time();
  std::cout << "CK   traits insertion, one-by-one:     " << t2 << std::endl;

  return std::pair<double, double>(t1,t2);
}

template<typename NT_point>
std::pair<double,double> insert_CORE_points(const std::vector<NT_point>& dpts)
{
  typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<>      Gt;
  typedef Gt::Point_2                                             Point_2;
  typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Gt>           Dt;

  CGAL::Timer timer;
  std::vector<Point_2> pts;
  std::vector<Point_2>::iterator ip;

  for(std::size_t i=0; i<dpts.size(); ++i)
    pts.push_back(Point_2(dpts[i].x(), dpts[i].y()));

  Dt dt_end;
  timer.start();
  dt_end.insert(pts.begin(),pts.end());
  timer.stop();
  double t1 = timer.time();
  std::cout << "CORE traits insertion, iterator input: " << t1 << std::endl;

  Dt dt_during;
  timer.reset();
  timer.start();
  for(ip = pts.begin(); ip != pts.end(); ++ip)
    dt_during.insert(*ip);

  timer.stop();
  double t2 = timer.time();
  std::cout << "CORE traits insertion, on-by-one:      " << t2 << std::endl;

  return std::pair<double,double>(t1, t2);
}

int main(int argc, char** argv)
{
  typedef double                                              NT;
  typedef CGAL::Cartesian<NT>                                 NT_kernel;
  typedef NT_kernel::Point_2                                  NT_point;
  typedef CGAL::Creator_uniform_2<NT, NT_point>               Creator;

  if(argc < 2) {
    std::cout << "usage: " << argv[0] << " [number_of_points]" << std::endl;
    return -1;
  }

  int N = atoi(argv[1]);

  int iters = 1;
  if(argc >= 3)
    iters = atoi(argv[2]);

  CGAL::Random_points_in_disc_2<NT_point, Creator> in_disc(0.998);

  double timeCK1 = 0;
  double timeCK2 = 0;
  double timeCORE1 = 0;
  double timeCORE2 = 0;
  double timeEUCL1 = 0;
  double timeEUCL2 = 0;
  for(int j=0; j<iters; ++j)
  {
    std::vector<NT_point> pts;
    for(int i=0; i<N; ++i)
      pts.push_back(*(in_disc++));

    std::cout << "----------- iteration " << j << " -----------" << std::endl;
    std::pair<double, double> resEUCL = insert_Euclidean(pts);
    timeEUCL1 += resEUCL.first;
    timeEUCL2 += resEUCL.second;

    std::pair<double, double> resCK = insert_CK_points(pts);
    timeCK1 += resCK.first;
    timeCK2 += resCK.second;

    std::pair<double, double> resCORE = insert_CORE_points(pts);
    timeCORE1 += resCORE.first;
    timeCORE2 += resCORE.second;
    std::cout << std::endl;
  }

  double diters(iters);
  timeEUCL1 /= diters;
  timeCORE1 /= diters;
  timeCK1 /= diters;
  timeEUCL2 /= diters;
  timeCORE2 /= diters;
  timeCK2 /= diters;

  std::cout << "----------------------------------------------------" << std::endl;
  std::cout << "Average time for EUCL (one-at-a-time): " << timeEUCL2 << std::endl;
  std::cout << "Average time for CK   (one-at-a-time): " << timeCK2   << std::endl;
  std::cout << "Average time for CORE (one-at-a-time): " << timeCORE2 << std::endl;
  std::cout << "----------------------------------------------------" << std::endl;
  std::cout << "Average time for EUCL (iterator):      " << timeEUCL1 << std::endl;
  std::cout << "Average time for CK   (iterator):      " << timeCK1   << std::endl;
  std::cout << "Average time for CORE (iterator):      " << timeCORE1 << std::endl;
  std::cout << std::endl;

  return 0;
}
