#include <iostream>
#include <fstream>
#include <functional>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Zp.h>
#include <CGAL/Z2.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Surface_mesh_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Filtration_lower_star.h>
#include <CGAL/HDVF/Hdvf_persistence.h>
#include <CGAL/OSM/OSM.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

//#define BUILD_TEST_DATA

//typedef int Coefficient_ring;
typedef CGAL::Z2 Coefficient_ring;
//typedef CGAL::Zp<5, char, true> Coefficient_ring;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef HDVF::Hdvf_traits_3<Kernel> Traits;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

typedef CGAL::OSM::Sub_sparse_matrix<CGAL::OSM::Sparse_chain> Sparse_matrix_struct;
typedef HDVF::Simplicial_chain_complex<Coefficient_ring, Traits, Sparse_matrix_struct> Complex;
typedef double Degree;
typedef HDVF::Filtration_lower_star<Complex, Degree> FiltrationType;
typedef HDVF::Hdvf_persistence<Complex, Degree, FiltrationType> HDVF_type;

std::function<double(const Point_3&)> sample_f = [](const Point_3& X)
{
    // Filtration : (x,y,z) -> degree x+y
    return -X[2];
} ;

int main(int argc, char **argv) {
    std::string filename;
    if (argc > 2) {
        std::cerr << "usage: test_hdvf_persistence [off_file]" << std::endl;
    }
    else if (argc == 1) filename = "data/four_triangles.obj" ;
    else filename = argv[1] ;

    // Load mesh object
    Surface_mesh sm;
    CGAL::IO::read_polygon_mesh(filename, sm);
    HDVF::Surface_mesh_io<Surface_mesh,Traits> mesh(sm) ;

    mesh.print_infos();

    // Build simplicial chain complex
    Complex complex(mesh);

    std::cout << complex;

    // -- First: Define the filtration function
    std::function<Degree(size_t)> f(HDVF::degree_function<Complex,Point_3>(complex, sample_f));


    // -- Second: build the associated lower star filtration
    FiltrationType filtration(complex, f);


    // Build empty persistent HDVF (with vtk export activated)
    HDVF_type hdvf(complex, filtration, HDVF::OPT_FULL, true);
    hdvf.write_matrices();

    // Compute a perfect HDVF
    std::cout << "--- computing" << std::endl;
    hdvf.compute_perfect_hdvf(true);
    hdvf.write_matrices();
    CGAL::IO::write_VTK(hdvf, complex, "tmp/res", true) ;

#ifdef BUILD_TEST_DATA
    // Write HDVF_persistence to a .hdvf file
    hdvf.write_hdvf_reduction("data/test_hdvf_persistence/hdvf_persist.hdvf");
#endif

    HDVF_type hdvf2(complex, filtration, HDVF::OPT_FULL, true);
    hdvf2.read_hdvf_reduction("data/test_hdvf_persistence/hdvf_persist.hdvf");

    // Compare
    bool test_rw(hdvf.compare(hdvf2));
    std::cerr << "-- Compare write/read HDVF_persistence: " << test_rw << std::endl;
    assert(test_rw);

    return 0;
}


