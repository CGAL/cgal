#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <vector>
#include <algorithm>
#include <numeric>
#include <CGAL/assertions.h>


#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/Curvatures/interpolated_curvature_measures.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <chrono>

#define N 4
namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;

int main(int argc, char* argv[])
{
    srand(time(NULL));

    std::vector<Epic::Vector_3> X(N);
    std::vector<Epic::Vector_3> U(N);

    for (int i = 0; i < N; i++) {
        X[i] = { rand() , rand() , rand() };
        U[i] = { rand() , rand() , rand() };
        U[i] = U[i] / sqrt(U[i] * U[i]);
    }
    std::cout << PMP::interpolated_mu_i_face(X, U, PMP::MU0_AREA_MEASURE) << std::endl;
    std::cout << PMP::interpolated_mu_i_face(X, U, PMP::MU1_MEAN_CURVATURE_MEASURE) << std::endl;
    std::cout << PMP::interpolated_mu_i_face(X, U, PMP::MU2_GAUSSIAN_CURVATURE_MEASURE) << std::endl;

    /*srand(time(NULL));
    Epic::Vector_3 x, y;
    Epic::FT d = 0;

    Epic::Vector_3 z;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
    {
        x * y;
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << std::endl;*/

}

