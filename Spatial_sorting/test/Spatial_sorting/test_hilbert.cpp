#include <cassert>

#include <CGAL/hilbert_sort.h>
#include <CGAL/hilbert_sort_on_sphere.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_d.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/constructions_d.h>

#include <iostream>
#include <algorithm>
#include <vector>

#include <CGAL/Real_timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point_2;
typedef K::Point_3                                          Point_3;
typedef K::Vector_3                                         Vector_3;
typedef CGAL::Creator_uniform_2<double,Point_2>             Creator_2;
typedef CGAL::Creator_uniform_3<double,Point_3>             Creator_3;

typedef CGAL::Cartesian_d<double>                           Kd;
typedef Kd::Point_d                                         Point;
typedef CGAL::Creator_uniform_d<std::vector<double>::iterator, Point>Creator_d;

int main ()
{
  int nb_points_2 = 10000, nb_points_3 = 10000,
      nb_points_d=10000, small_nb_points_d=3;
  CGAL::Random random (42);
  CGAL::Real_timer timer;

  std::cout << "Testing Hilbert sort." << std::endl;

  {
    std::cout << "Testing 2D (median policy): Generating "<<nb_points_2<<" random points... " << std::flush;

    std::vector<Point_2> v;
    v.reserve (nb_points_2);

    CGAL::Random_points_in_square_2<Point_2> gen (1.0, random);

    for (int i = 0; i < nb_points_2 - 1; ++i)
      v.push_back (*gen++);
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point_2> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.end());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Sorting points (parallel)...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort<CGAL::Parallel_if_available_tag>(v2.begin(), v2.end());
    timer.stop();

    std::cout << "done in " << timer.time() << "seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;
    assert(v == v2);

    std::cout << "Ok" << std::endl;
  }
  {
    int size=256;               // 2^(xd)   with x=4 d=2
    double box_size = 15.0;     // 2^x -1
    std::cout << "Testing 2D (median policy): Generating "<<size<<" grid points... " << std::flush;
    std::vector<Point_2> v;
    v.reserve(size);

    CGAL::points_on_square_grid_2 (box_size, (std::size_t)size,
                                   std::back_inserter(v), Creator_2() );

    std::cout << "done." << std::endl;

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.end());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    for (int i = 0; i < size-1; ++i)
      assert(CGAL::squared_distance( v[i], v[i+1]) - 4.0 < 0.1 );

    std::cout << "OK." << std::endl;
  }
  {
    std::cout << "Testing 2D (middle policy): Generating "
              <<nb_points_2<<" random points... " << std::flush;

    std::vector<Point_2> v;
    v.reserve (nb_points_2);

    CGAL::Random_points_in_square_2<Point_2> gen (1.0, random);

    for (int i = 0; i < nb_points_2 - 1; ++i)
      v.push_back (*gen++);
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point_2> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort(v.begin(),v.end(), CGAL::Hilbert_sort_middle_policy());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    std::sort (v.begin(),  v.end(),  K().less_xy_2_object());
    std::sort (v2.begin(), v2.end(), K().less_xy_2_object());
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }
  {
    int size=256;               // 2^(xd)
    double box_size = 15.0;     // 2^x -1                with x=4 d=2
    std::cout << "Testing 2D (middle policy): Generating "
              <<size<<" grid points... " << std::flush;
    std::vector<Point_2> v;
    v.reserve(size);

    CGAL::points_on_square_grid_2 (box_size, (std::size_t)size,
                                   std::back_inserter(v), Creator_2() );

    std::cout << "done." << std::endl;

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort(v.begin(),v.end(), CGAL::Hilbert_sort_middle_policy());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    for (int i = 0; i < size-1; ++i) {
      assert(CGAL::squared_distance( v[i], v[i+1]) - 4.0 < 0.1 );
    }
    std::cout << "OK." << std::endl;
  }

  {
    std::cout << "Testing 3D (median policy): Generating "<<nb_points_3<<" random points... " << std::flush;

    std::vector<Point_3> v;
    v.reserve (nb_points_3);

    CGAL::Random_points_in_cube_3<Point_3> gen (1.0, random);

    for (int i = 0; i < nb_points_3 - 1; ++i)
      v.push_back (*gen++);
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point_3> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.end());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Sorting points (parallel)...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort<CGAL::Parallel_if_available_tag>(v2.begin(), v2.end());
    timer.stop();

    std::cout << "done in " << timer.time() << "seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }
  {
    int size=512;              // 2^(xd)   with x=3 d=3
    double box_size = 7.0;     // 2^x -1

    std::cout << "Testing 3D (median policy): Generating "<<size<<" grid points... " << std::flush;

    std::vector<Point_3> v;
    v.reserve(size);

    CGAL::points_on_cube_grid_3 (box_size, (std::size_t)size,
                                 std::back_inserter(v), Creator_3() );

    std::cout << "done." << std::endl;

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.end());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    for (int i = 0; i < size-1; ++i)
      assert(CGAL::squared_distance( v[i], v[i+1]) - 4.0 < 0.1 );

    std::cout << "OK." << std::endl;
  }

  {
    std::cout << "Testing 3D (middle policy): Generating "<<nb_points_3<<" random points... " << std::flush;

    std::vector<Point_3> v;
    v.reserve (nb_points_3);

    CGAL::Random_points_in_cube_3<Point_3> gen (1.0, random);

    for (int i = 0; i < nb_points_3 - 1; ++i)
      v.push_back (*gen++);
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point_3> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort(v.begin(),v.end(),CGAL::Hilbert_sort_middle_policy());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    std::sort (v.begin(),  v.end(),  K().less_xyz_3_object());
    std::sort (v2.begin(), v2.end(), K().less_xyz_3_object());
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }
  {
    int size=4096;              // 2^(xd)   with x=4 d=3
    double box_size = 15.0;     // 2^x -1

    std::cout << "Testing 3D (middle policy): Generating "<<size<<" grid points... " << std::flush;

    std::vector<Point_3> v;
    v.reserve(size);

    CGAL::points_on_cube_grid_3 (box_size, (std::size_t)size,
                                 std::back_inserter(v), Creator_3() );

    std::cout << "done." << std::endl;

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort(v.begin(),v.end(),CGAL::Hilbert_sort_middle_policy());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    for (int i = 0; i < size-1; ++i) {
      assert(CGAL::squared_distance( v[i], v[i+1]) - 4.0 < 0.1 );
    }

    std::cout << "OK." << std::endl;
  }
  {
    std::cout << "Testing Spherical (median policy): Generating "<<nb_points_3<<" random points... " << std::flush;

    std::vector<Point_3> v;
    v.reserve (nb_points_3);

    CGAL::Random_points_on_sphere_3<Point_3> gen (1.0, random);

    for (int i = 0; i < nb_points_3 - 1; ++i)
      v.push_back (*gen++);
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point_3> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort(v.begin(),v.end());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Sorting points (parallel)...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort<CGAL::Parallel_if_available_tag>(v2.begin(), v2.end());
    timer.stop();

    std::cout << "done in " << timer.time() << "seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }
  {
    std::cout << "Testing Spherical (median policy) + given sphere: Generating "<<nb_points_3<<" random points... " << std::flush;

    std::vector<Point_3> v;
    v.reserve (nb_points_3);

    CGAL::Random_points_on_sphere_3<Point_3> gen (2.0, random);

    for (int i = 0; i < nb_points_3 - 1; ++i)
      v.push_back (*gen++ + Vector_3(3,5,5));
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point_3> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort_on_sphere(v.begin(),v.end(), 4, CGAL::ORIGIN + Vector_3(3,5,5));
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    std::sort (v.begin(),  v.end(),  K().less_xyz_3_object());
    std::sort (v2.begin(), v2.end(), K().less_xyz_3_object());
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }
  {
    std::cout << "Testing Spherical (middle policy): Generating "<<nb_points_3<<" random points... " << std::flush;

    std::vector<Point_3> v;
    v.reserve (nb_points_3);

    CGAL::Random_points_on_sphere_3<Point_3> gen (1.0, random);

    for (int i = 0; i < nb_points_3 - 1; ++i)
      v.push_back (*gen++);
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point_3> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort_on_sphere(v.begin(),v.end(),CGAL::Hilbert_sort_middle_policy());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    std::sort (v.begin(),  v.end(),  K().less_xyz_3_object());
    std::sort (v2.begin(), v2.end(), K().less_xyz_3_object());
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }
  {
    std::cout << "Testing Spherical (middle policy) + given sphere: Generating "<<nb_points_3<<" random points... " << std::flush;

    std::vector<Point_3> v;
    v.reserve (nb_points_3);

    CGAL::Random_points_on_sphere_3<Point_3> gen (2.0, random);

    for (int i = 0; i < nb_points_3 - 1; ++i)
      v.push_back (*gen++ + Vector_3(3,5,5));
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point_3> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort_on_sphere(v.begin(),v.end(),CGAL::Hilbert_sort_middle_policy(), 4, CGAL::ORIGIN + Vector_3(3,5,5));
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    std::sort (v.begin(),  v.end(),  K().less_xyz_3_object());
    std::sort (v2.begin(), v2.end(), K().less_xyz_3_object());
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }
  {
    int dim =5;
    std::cout << "Testing "<<dim<<"D: Generating "<<nb_points_d<<" random points... " << std::flush;

    std::vector<Point> v;
    v.reserve (nb_points_d);

    CGAL::Random_points_in_cube_d<Point> gen (dim, 1.0, random);

    for (int i = 0; i < nb_points_d - 1; ++i)
      v.push_back (*gen++);
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.end(),CGAL::Hilbert_sort_median_policy());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    std::sort (v.begin(),  v.end(), Kd().less_lexicographically_d_object());
    std::sort (v2.begin(), v2.end(),Kd().less_lexicographically_d_object());
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }
  {
    int dim = 1000;
    std::cout << "Testing "<<dim<<"D: Generating "<<nb_points_d<<" random points... " << std::flush;

    std::vector<Point> v;
    v.reserve (nb_points_d);

    CGAL::Random_points_in_cube_d<Point> gen (dim, 1.0, random);

    for (int i = 0; i < nb_points_d - 1; ++i)
      v.push_back (*gen++);
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.end());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    std::sort (v.begin(),  v.end(), Kd().less_lexicographically_d_object());
    std::sort (v2.begin(), v2.end(),Kd().less_lexicographically_d_object());
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }
  {
    int dim = 10;
    std::cout << "Testing "<<dim<<"D (middle policy): Generating "<<nb_points_d<<" random points... " << std::flush;

    std::vector<Point> v;
    v.reserve (nb_points_d);

    CGAL::Random_points_in_cube_d<Point> gen (dim, 1.0, random);

    for (int i = 0; i < nb_points_d - 1; ++i)
      v.push_back (*gen++);
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point> v2 (v);

    std::cout << "            Sorting points ...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.end(),
                        CGAL::Hilbert_sort_middle_policy());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    std::sort (v.begin(),  v.end(), Kd().less_lexicographically_d_object());
    std::sort (v2.begin(), v2.end(),Kd().less_lexicographically_d_object());
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }
  {
    int dim=5;
    std::cout << "Testing "<<dim<<"D: Generating "<<small_nb_points_d<<" random points... " << std::flush;

    std::vector<Point> v;
    v.reserve (nb_points_d);

    CGAL::Random_points_in_cube_d<Point> gen (dim, 1.0, random);

    for (int i = 0; i < nb_points_d - 1; ++i)
      v.push_back (*gen++);
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.end());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    std::sort (v.begin(),  v.end(), Kd().less_lexicographically_d_object());
    std::sort (v2.begin(), v2.end(),Kd().less_lexicographically_d_object());
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }
  {
    int dim=5;
    int size=32768;            // 2^(x.dim)   with x=3
    double box_size = 7.0;     // 2^x -1

    std::cout << "Testing "<<dim<<"D: Generating "<<size<<" grid points... " << std::flush;

    std::vector<Point> v(size);

    CGAL::points_on_cube_grid_d (dim, box_size, (std::size_t)size,
                                 v.begin(), Creator_d(dim) );

    std::cout << "done." << std::endl;

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.begin()+size);
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    for (int i = 0; i < size-1; ++i)
      assert(CGAL::squared_distance( v[i], v[i+1]) - 4.0 < 0.1 );

    std::cout << "OK." << std::endl;
  }
  {
    int dim=3;
    int size=32768;            // 2^(x.dim)   with x=5
    double box_size = 31.0;     // 2^x -1

    std::cout << "Testing "<<dim<<"D (middle policy): Generating "<<size<<" grid points... " << std::flush;

    std::vector<Point> v(size);

    CGAL::points_on_cube_grid_d (dim, box_size, (std::size_t)size,
                                 v.begin(), Creator_d(dim) );

    std::cout << "done." << std::endl;

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.begin()+size,
                        CGAL::Hilbert_sort_middle_policy());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    for (int i = 0; i < size-1; ++i)
      assert(CGAL::squared_distance( v[i], v[i+1]) - 4.0 < 0.1 );

    std::cout << "OK." << std::endl;
  }

  {
    int dim=5;
    int size=32768;            // 2^(x.dim)   with x=3
    double box_size = 7.0;     // 2^x -1

    std::cout << "Testing "<<dim<<"D (middle policy): Generating "<<size<<" grid points... " << std::flush;

    std::vector<Point> v(size);

    CGAL::points_on_cube_grid_d (dim, box_size, (std::size_t)size,
                                 v.begin(), Creator_d(dim) );

    std::cout << "done." << std::endl;

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.begin()+size,
                        CGAL::Hilbert_sort_middle_policy());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    for (int i = 0; i < size-1; ++i)
      assert(CGAL::squared_distance( v[i], v[i+1]) - 4.0 < 0.1 );

    std::cout << "OK." << std::endl;
  }

  {
    int dim = 50;
    std::cout << "Testing "<<dim<<"D (median policy): Generating "<<nb_points_d<<" random points... " << std::flush;

    std::vector<Point> v;
    v.reserve (nb_points_d);

    CGAL::Random_points_in_cube_d<Point> gen (dim, 1.0, random);

    for (int i = 0; i < nb_points_d - 1; ++i)
      v.push_back (*gen++);
    v.push_back(v[0]); //insert twice the same point

    std::cout << "done." << std::endl;

    std::vector<Point> v2 (v);

    std::cout << "            Sorting points...    " << std::flush;

    timer.reset();timer.start();
    CGAL::hilbert_sort (v.begin(), v.end());
    timer.stop();

    std::cout << "done in "<<timer.time()<<"seconds." << std::endl;

    std::cout << "            Checking...          " << std::flush;

    std::sort (v.begin(),  v.end(), Kd().less_lexicographically_d_object());
    std::sort (v2.begin(), v2.end(),Kd().less_lexicographically_d_object());
    assert(v == v2);

    std::cout << "no points lost." << std::endl;
  }

  return 0;
}
