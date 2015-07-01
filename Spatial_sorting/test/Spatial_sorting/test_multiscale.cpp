#include <CGAL/basic.h>
#include <cassert>

#include <CGAL/spatial_sort.h>
#include <CGAL/spatial_sort_on_sphere.h>

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


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;

typedef CGAL::Creator_uniform_2<double,Point_2>  Creator_2;
typedef CGAL::Creator_uniform_3<double,Point_3>  Creator_3;

typedef CGAL::Cartesian_d<double>                           Kd;
typedef Kd::Point_d                                         Point;

int main ()
{
  const int nb_points_2 = 50000, nb_points_3 = 50000, nb_points_d=50000;
    CGAL::Random random (42);

    std::cout << "Testing Multiscale<Hilbert> sort." << std::endl;

    {
        std::cout << "Testing 2D: Generating points... " << std::flush;

        std::vector<Point_2> v;
        v.reserve (nb_points_2);

        CGAL::Random_points_in_square_2<Point_2,Creator_2> gen (1.0, random);

        for (int i = 0; i < nb_points_2; ++i)
            v.push_back (*gen++);

        std::cout << "done." << std::endl;

        std::vector<Point_2> v2 (v);

        std::cout << "            Sorting points...    " << std::flush;

        CGAL::spatial_sort (v.begin(), v.end());

        std::cout << "done." << std::endl;

        std::cout << "            Checking...          " << std::flush;

        std::sort (v.begin(),  v.end(),  K().less_xy_2_object());
        std::sort (v2.begin(), v2.end(), K().less_xy_2_object());
        assert(v == v2);

        std::cout << "no points lost." << std::endl;
    }

    {
        std::cout << "Testing 3D: Generating points... " << std::flush;

        std::vector<Point_3> v;
        v.reserve (nb_points_3);

        CGAL::Random_points_in_cube_3<Point_3,Creator_3> gen (1.0, random);

        for (int i = 0; i < nb_points_3; ++i)
            v.push_back (*gen++);

        std::cout << "done." << std::endl;

        std::vector<Point_3> v2 (v);

        std::cout << "            Sorting points...    " << std::flush;

        CGAL::spatial_sort (v.begin(), v.end());

        std::cout << "done." << std::endl;

        std::cout << "            Checking...          " << std::flush;

        std::sort (v.begin(),  v.end(),  K().less_xyz_3_object());
        std::sort (v2.begin(), v2.end(), K().less_xyz_3_object());
        assert(v == v2);

        std::cout << "no points lost." << std::endl;
    }

	{
        std::cout << "Testing Spherical: Generating points... " << std::flush;

        std::vector<Point_3> v;
        v.reserve (nb_points_3);

        CGAL::Random_points_on_sphere_3<Point_3,Creator_3> gen (1.0, random);

        for (int i = 0; i < nb_points_3; ++i)
            v.push_back (*gen++);

        std::cout << "done." << std::endl;

        std::vector<Point_3> v2 (v);

        std::cout << "            Sorting points...    " << std::flush;

        CGAL::spatial_sort_on_sphere (v.begin(), v.end());

        std::cout << "done." << std::endl;

        std::cout << "            Checking...          " << std::flush;

        std::sort (v.begin(),  v.end(),  K().less_xyz_3_object());
        std::sort (v2.begin(), v2.end(), K().less_xyz_3_object());
        assert(v == v2);

        std::cout << "no points lost." << std::endl;
    }

    {
        std::cout << "Testing Spherical + given sphere: Generating points... " << std::flush;

        std::vector<Point_3> v;
        v.reserve (nb_points_3);

        CGAL::Random_points_on_sphere_3<Point_3,Creator_3> gen (3.0, random);

        for (int i = 0; i < nb_points_3; ++i)
            v.push_back (*gen++ + Vector_3(10,10,3));

        std::cout << "done." << std::endl;

        std::vector<Point_3> v2 (v);

        std::cout << "            Sorting points...    " << std::flush;

        CGAL::spatial_sort_on_sphere (v.begin(), v.end(), 9, Point_3(10,10,3));

        std::cout << "done." << std::endl;

        std::cout << "            Checking...          " << std::flush;

        std::sort (v.begin(),  v.end(),  K().less_xyz_3_object());
        std::sort (v2.begin(), v2.end(), K().less_xyz_3_object());
        assert(v == v2);

        std::cout << "no points lost." << std::endl;
    }

    {
      int dim=5;
      std::cout << "Testing "<<dim<<"D: Generating points... " << std::flush;

        std::vector<Point> v;
        v.reserve (nb_points_d);

        CGAL::Random_points_in_cube_d<Point> gen (dim, 1.0, random);

        for (int i = 0; i < nb_points_d; ++i)
            v.push_back (*gen++);

        std::cout << "done." << std::endl;

        std::vector<Point> v2 (v);

        std::cout << "            Sorting points...    " << std::flush;

        CGAL::spatial_sort (v.begin(), v.end());

        std::cout << "done." << std::endl;

        std::cout << "            Checking...          " << std::flush;

        std::sort (v.begin(),  v.end(), Kd().less_lexicographically_d_object());
        std::sort (v2.begin(), v2.end(),Kd().less_lexicographically_d_object());
        assert(v == v2);

        std::cout << "no points lost." << std::endl;
    }

    return 0;
}
