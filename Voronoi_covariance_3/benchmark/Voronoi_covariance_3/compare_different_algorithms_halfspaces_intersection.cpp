/* #define CGAL_PROFILE */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/dual/halfspace_intersection_3.h>
#include <CGAL/dual/halfspace_intersection_with_constructions_3.h>

#include <CGAL/Timer.h>
#include <fstream>

/* typedef CGAL::Simple_cartesian<CGAL::Gmpq>                          K; */
/* typedef CGAL::Simple_cartesian<double>                              K; */
typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Exact_predicates_exact_constructions_kernel           EK;

typedef K::Plane_3 Plane;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;

typedef EK::Plane_3 Exact_Plane;
typedef EK::Point_3 Exact_Point;
typedef CGAL::Polyhedron_3<EK> Exact_Polyhedron_3;

#include <CGAL/point_generators_3.h>

template <typename K>
typename K::Plane_3 tangent_plane (typename K::Point_3 const& p) {
    typename K::Vector_3 v(p.x(), p.y(), p.z());
    typename K::Plane_3 plane(v.x(), v.y(), v.z(), -(p - CGAL::ORIGIN) * v);

    return plane;
}

template <class Iterator>
double test_new (Iterator first,
                Iterator beyond) {
    Polyhedron_3 poly_sphere;

    CGAL::Timer timer;
    timer.start();
    CGAL::halfspace_intersection_3(first,
                                   beyond,
                                   poly_sphere,
                                   Point(0, 0, 0));
    timer.stop();

    return timer.time();
}

template <class Iterator>
double test_old (Iterator first,
                Iterator beyond) {
    Polyhedron_3 poly_sphere;

    CGAL::Timer timer;
    timer.start();
    CGAL::halfspace_intersection_with_constructions_3(first,
                                                      beyond,
                                                      poly_sphere,
                                                      Point(0, 0, 0));
    timer.stop();

    return timer.time();
}

template <class Iterator>
double test_old_exact (Iterator first,
                       Iterator beyond) {
    Exact_Polyhedron_3 poly_sphere;

    CGAL::Timer timer;
    timer.start();
    CGAL::halfspace_intersection_with_constructions_3(first,
                                                      beyond,
                                                      poly_sphere,
                                                      Exact_Point(0, 0, 0));
    timer.stop();

    return timer.time();
}

void test (int n, double & new_time, double & old_time) {
    std::list<Plane> sphere_planes;
    CGAL::Random_points_on_sphere_3<Point> g;
    for (int i = 0; i < n; i++) {
        sphere_planes.push_back(tangent_plane<K>(*g++));
    }

    new_time = test_new(sphere_planes.begin(), sphere_planes.end());
    old_time = test_old(sphere_planes.begin(), sphere_planes.end());
}

void test_exact (int n, double & new_time, double & old_time) {
    std::list<Plane> sphere_planes;
    std::list<Exact_Plane> exact_sphere_planes;

    CGAL::Random_points_on_sphere_3<Point> g;
    CGAL::Random_points_on_sphere_3<Exact_Point> eg;

    for (int i = 0; i < n; i++) {
        sphere_planes.push_back(tangent_plane<K>(*g++));
        exact_sphere_planes.push_back(tangent_plane<EK>(*eg++));
    }

    new_time = test_new(sphere_planes.begin(), sphere_planes.end());
    old_time = test_old_exact(exact_sphere_planes.begin(), exact_sphere_planes.end());
}

int main (int argc, char *argv[]) {
    int N;
    if (argc > 1) {
        N = atoi(argv[1]);
    } else {
        N = 200;
    }

    int step;
    if (argc > 2) {
        step = atoi(argv[2]);
    } else {
        step = 2;
    }

    double times_new[N];
    double times_old[N];
    int s = 0;

    for (int i = 4; i <= N; i += step) {
        std::cout << "i = " << i << std::endl;
        double new_time, old_time;
        /* test(i, new_time, old_time); */
        test_exact(i, new_time, old_time);
        std::cout << "time for new implementation "<< new_time << std::endl;
        std::cout << "time for old implementation "<< old_time << std::endl;
        times_new[s] = new_time;
        times_old[s] = old_time;
        s++;
    }

    std::ofstream file_new("new.txt");
    for (int i = 0; i < s; i++) {
        file_new << times_new[i] << std::endl;
    }

    std::ofstream file_old("old.txt");
    for (int i = 0; i < s; i++) {
        file_old << times_old[i] << std::endl;
    }

    return 0;
}

