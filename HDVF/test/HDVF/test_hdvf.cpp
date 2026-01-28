#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Zp.h>
#include <CGAL/Z2.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/OSM/OSM.h>

#define BUILD_TEST_DATA

namespace HDVF = CGAL::Homological_discrete_vector_field;

//typedef int Coefficient_ring;
//typedef CGAL::Z2 Coefficient_ring;
typedef CGAL::Zp<5, char, true> Coefficient_ring;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef HDVF::Hdvf_traits_3<Kernel> Traits;

using Complex = HDVF::Simplicial_chain_complex<Coefficient_ring,Traits> ;
using HDVF_type = HDVF::Hdvf<Complex> ;

int main(int argc, char **argv) {
    std::string filename;
    if (argc > 2) {
        std::cerr << "usage: test_hdvf_core [off_file]" << std::endl;
    }
    else if (argc == 1) filename = "data/three_triangles.off" ;
    else filename = argv[1] ;

    std::cerr << "---- Test HDVF computation" << std::endl;

    // Load off into Mesh_object_io
    HDVF::Mesh_object_io<Traits> mesh ;
    mesh.read_off(filename);

    //    mesh.print_infos();

    // Build simplicial chain complex
    Complex complex(mesh);

    std::cout << complex;

    // Build empty HDVF
    HDVF_type hdvf(complex, HDVF::OPT_FULL) ;

    // Compute a perfect HDVF
    hdvf.compute_perfect_hdvf();
    //        hdvf.compute_rand_perfect_hdvf();
    std::cerr << std::endl;

#ifdef BUILD_TEST_DATA
    // Save HDVF to .hdvf file
    hdvf.write_hdvf_reduction("data/test_hdvf/test_hdvf.hdvf") ;
#endif

    {
        // Read HDVF from .hdvf
        HDVF_type hdvf2(complex, HDVF::OPT_FULL);
        hdvf2.read_hdvf_reduction("data/test_hdvf/test_hdvf.hdvf");

        // Compare
        bool test_true(hdvf.compare(hdvf2));
        std::cerr << "-- Test HDVF built: " << test_true << std::endl;
        assert(test_true);
    }

    // Test R
    std::cerr << "---- Test R operation" << std::endl;
    {
        bool test_valid_R(hdvf.is_valid_pair_for_R(2, 0, 1) == 1);
        std::cerr << "-- Test is_valid_pair_for_R: " << test_valid_R << std::endl;
        assert(test_valid_R);
        // Test operation
        {
            // Read HDVF from .hdvf
            HDVF_type hdvf2(complex, HDVF::OPT_FULL);
            hdvf2.read_hdvf_reduction("data/test_hdvf/test_hdvf.hdvf");
            hdvf2.R(2,0,1);
#ifdef BUILD_TEST_DATA
            // Save HDVF to .hdvf file
            hdvf2.write_hdvf_reduction("data/test_hdvf/test_hdvfR.hdvf") ;
#endif
            // Read HDVF after R from .hdvf
            HDVF_type hdvf3(complex, HDVF::OPT_FULL);
            hdvf3.read_hdvf_reduction("data/test_hdvf/test_hdvfR.hdvf");
            // Compare
            bool test_R(hdvf2.compare(hdvf3));
            std::cerr << "-- Test R operation: " << test_R << std::endl;
            assert(test_R);
        }

        bool test_invalid_R(hdvf.is_valid_pair_for_R(1, 0, 1) == 0);
        std::cerr << "-- Test incorrect is_valid_pair_for_R: " << test_invalid_R << std::endl;
        assert(test_invalid_R);
    }

    // Test M
    std::cerr << "---- Test M operation" << std::endl;
    {
        bool foundM;
        std::vector<HDVF::Cell_pair> pairsM(hdvf.find_pairs_M(1, foundM)), pairsM2;
        pairsM2.push_back(HDVF::Cell_pair({8,7,1}));
        bool test_find_pairs_M(pairsM == pairsM2);
        std::cerr << "-- Test find_pairs_M(1,foundM): " << test_find_pairs_M << std::endl;
        assert(test_find_pairs_M);

        std::vector<HDVF::Cell_pair> pairsM3(hdvf.find_pairs_M(1, foundM, 8));
        bool test_find_pairs_M3(pairsM3 == pairsM2);
        std::cerr << "-- Test find_pairs_M(1,foundM,8): " << test_find_pairs_M3 << std::endl;
        assert(test_find_pairs_M3);

        HDVF::Cell_pair cell(hdvf.find_pair_M(1, foundM)), cell2(HDVF::Cell_pair({8,7,1}));
        bool test_find_pairs_M4(cell == cell2);
        std::cerr << "-- Test find_pair_M(1,foundM): " << test_find_pairs_M4 << std::endl;
        assert(test_find_pairs_M4);

        HDVF::Cell_pair cell3(hdvf.find_pair_M(1, foundM,8));
        bool test_find_pairs_M5(cell3 == cell2);
        std::cerr << "-- Test find_pair_M(1,foundM,8): " << test_find_pairs_M5 << std::endl;
        assert(test_find_pairs_M5);

        HDVF::Cell_pair pairM(pairsM.at(0)); // Pair 8,7,1
        bool test_valid_M(hdvf.is_valid_pair_for_M(pairM.sigma, pairM.tau, pairM.dim) == 1);
        std::cerr << "-- Test is_valid_pair_for_M: " << test_valid_M << std::endl;
        assert(test_valid_M);
        // Test operation
        {
            // Read HDVF from .hdvf
            HDVF_type hdvf2(complex, HDVF::OPT_FULL);
            hdvf2.read_hdvf_reduction("data/test_hdvf/test_hdvf.hdvf");
            hdvf2.M(pairM.sigma, pairM.tau, pairM.dim);
#ifdef BUILD_TEST_DATA
            // Save HDVF to .hdvf file
            hdvf2.write_hdvf_reduction("data/test_hdvf/test_hdvfM.hdvf") ;
#endif
            // Read HDVF after R from .hdvf
            HDVF_type hdvf3(complex, HDVF::OPT_FULL);
            hdvf3.read_hdvf_reduction("data/test_hdvf/test_hdvfM.hdvf");
            // Compare
            bool test_M(hdvf2.compare(hdvf3));
            std::cerr << "-- Test M operation: " << test_M << std::endl;
            assert(test_M);
        }

        bool test_invalid_M(hdvf.is_valid_pair_for_M(8, 6, 1) == 0);
        std::cerr << "-- Test incorrect is_valid_pair_for_M: " << test_invalid_M << std::endl;
        assert(test_invalid_M);
    }

    // Test W
    std::cerr << "---- Test W operation" << std::endl;
    {
        bool foundW;
        std::vector<HDVF::Cell_pair> pairsW(hdvf.find_pairs_W(1, foundW)), pairsW2;
        pairsW2.push_back(HDVF::Cell_pair({0,7,1}));
        pairsW2.push_back(HDVF::Cell_pair({4,7,1}));
        pairsW2.push_back(HDVF::Cell_pair({3,7,1}));
        pairsW2.push_back(HDVF::Cell_pair({6,7,1}));
        bool test_find_pairs_W(pairsW == pairsW2);
        std::cerr << "-- Test find_pairs_W(1,foundW): " << test_find_pairs_W << std::endl;
        assert(test_find_pairs_W);

        std::vector<HDVF::Cell_pair> pairsW3(hdvf.find_pairs_W(1, foundW, 7));
        bool test_find_pairs_W3(pairsW3 == pairsW2);
        std::cerr << "-- Test find_pairs_W(1,foundW,7): " << test_find_pairs_W3 << std::endl;
        assert(test_find_pairs_W3);

        HDVF::Cell_pair cell(hdvf.find_pair_W(1, foundW)), cell2(HDVF::Cell_pair({0,7,1}));
        bool test_find_pairs_W4(cell == cell2);
        std::cerr << "-- Test find_pair_W(1,foundW): " << test_find_pairs_W4 << std::endl;
        assert(test_find_pairs_W4);

        HDVF::Cell_pair cell3(hdvf.find_pair_W(1, foundW,7));
        bool test_find_pairs_W5(cell3 == cell2);
        std::cerr << "-- Test find_pair_W(1,foundM,7): " << test_find_pairs_W5 << std::endl;
        assert(test_find_pairs_W5);

        HDVF::Cell_pair pairW(pairsW.at(0)); // Pair 0,7,1
        bool test_valid_W(hdvf.is_valid_pair_for_W(pairW.sigma, pairW.tau, pairW.dim) == 1);
        std::cerr << "-- Test is_valid_pair_for_W: " << test_valid_W << std::endl;
        assert(test_valid_W);
        // Test operation
        {
            // Read HDVF from .hdvf
            HDVF_type hdvf2(complex, HDVF::OPT_FULL);
            hdvf2.read_hdvf_reduction("data/test_hdvf/test_hdvf.hdvf");
            hdvf2.W(pairW.sigma, pairW.tau, pairW.dim);
#ifdef BUILD_TEST_DATA
            // Save HDVF to .hdvf file
            hdvf2.write_hdvf_reduction("data/test_hdvf/test_hdvfM.hdvf") ;
#endif
            // Read HDVF after R from .hdvf
            HDVF_type hdvf3(complex, HDVF::OPT_FULL);
            hdvf3.read_hdvf_reduction("data/test_hdvf/test_hdvfM.hdvf");
            // Compare
            bool test_W(hdvf2.compare(hdvf3));
            std::cerr << "-- Test W operation: " << test_W << std::endl;
            assert(test_W);
        }

        bool test_invalid_W(hdvf.is_valid_pair_for_W(1, 7, 1) == 0);
        std::cerr << "-- Test incorrect is_valid_pair_for_W: " << test_invalid_W << std::endl;
        assert(test_invalid_W);
    }

    //    // Test z
    CGAL::IO::write_VTK(hdvf, complex, "tmp/test_hdvf");
    //
    //    // z over a critical cell
    //    std::vector<size_t> criticals(hdvf.psc_flags(HDVF::CRITICAL,1));
    //    HDVF_type::Column_chain test_z(hdvf.z(criticals.at(0),1));
    //    CGAL::IO::write_VTK(complex, "tmp/z_critical.vtk", test_z, 1);
    //
    //    // z over a primary cell
    //    std::vector<size_t> primary(hdvf.psc_flags(HDVF::PRIMARY,1));
    //    HDVF_type::Column_chain test_z_p(hdvf.z(primary.at(0),1));
    //    Complex::chain_to_vtk(complex, "tmp/z_primary.vtk", test_z_p, 1);
    //
    //    // Test co_z
    //
    //    // co_z over a critical cell
    //    HDVF_type::Column_chain test_co_z(hdvf.co_z(criticals.at(0),1));
    //    CGAL::IO::write_VTK(complex, "tmp/co_z_critical.vtk", test_co_z, 1);
    //
    //    // co_z over a secondary cell
    //    std::vector<size_t> secondary(hdvf.psc_flags(HDVF::SECONDARY,1));
    //    HDVF_type::Column_chain test_co_z_s(hdvf.co_z(secondary.at(0),1));
    //    CGAL::IO::write_VTK(complex, "tmp/co_z_secondary.vtk", test_co_z_s, 1);

    return 0;
}


