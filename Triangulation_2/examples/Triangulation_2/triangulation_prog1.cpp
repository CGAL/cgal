#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <unordered_map>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <boost/graph/graph_traits.hpp>

#include <CGAL/Polygon_mesh_processing/Curvatures/interpolated_corrected_curvature_measures.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <chrono>
#include <CGAL/Surface_mesh/IO/OFF.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Polyhedron_3<Epic> Polyhedron;
typedef CGAL::Surface_mesh<Epic::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;



int main(int argc, char* argv[])
{
    Mesh g1;
    if (!PMP::IO::read_polygon_mesh("small_bunny.obj", g1) || !CGAL::is_triangle_mesh(g1))
    {
        std::cerr << "Invalid input." << std::endl;
        return 1;
    }

    //int n1 = g1.size_of_facets();

    std::unordered_map<vertex_descriptor, Epic::Vector_3> vnm_vec;
    boost::associative_property_map< std::unordered_map<vertex_descriptor, Epic::Vector_3>> vnm(vnm_vec);

    PMP::compute_vertex_normals(g1, vnm);

    

    std::vector<Epic::FT> mu0_map, mu1_map, mu2_map;
    
    mu0_map = PMP::interpolated_corrected_measure_i(g1, PMP::MU0_AREA_MEASURE, vnm);
    mu1_map = PMP::interpolated_corrected_measure_i(g1, PMP::MU1_MEAN_CURVATURE_MEASURE, vnm);
    mu2_map = PMP::interpolated_corrected_measure_i(g1, PMP::MU2_GAUSSIAN_CURVATURE_MEASURE, vnm);

    int n = g1.faces().size();

    for (int i = 0; i < n; i++)
    {
        std::cout << mu0_map[i] << "\n";
    }

    std::cout << "\n";

    for (int i = 0; i < n; i++)
    {
        std::cout << mu1_map[i]  << "\n";
    }

    std::cout <<  "\n";

    for (int i = 0; i < n; i++)
    {
        std::cout << mu2_map[i] << "\n";
    }


    

    /*srand(time(NULL));

    CGAL::GetGeomTraits<Polyhedron> GT;

    std::vector<Epic::Vector_3> X(N);
    std::vector<Epic::Vector_3> U(N);

    for (int i = 0; i < N; i++) {
        X[i] = { rand() , rand() , rand() };
        U[i] = { rand() , rand() , rand() };
        U[i] = U[i] / sqrt(U[i] * U[i]);
    }
    

    std::cout << PMP::interpolated_mu_i_face(X, U, PMP::MU0_AREA_MEASURE) << std::endl;
    std::cout << PMP::interpolated_mu_i_face(X, U, PMP::MU1_MEAN_CURVATURE_MEASURE) << std::endl;
    std::cout << PMP::interpolated_mu_i_face(X, U, PMP::MU2_GAUSSIAN_CURVATURE_MEASURE) << std::endl;*/

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

