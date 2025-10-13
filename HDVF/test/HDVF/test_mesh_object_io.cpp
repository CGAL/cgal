#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/HDVF/Mesh_object_io.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

int main() {
    // ------ Test constructor from vectors of nodes and Io_cell_types

    std::vector<HDVF::Io_node_type> nodes1, nodes2;
    std::vector<HDVF::Io_cell_type> cells1, cells2;

    // (nodes1, cells1) : "pure" mesh
    // (nodes2, cells2) : heterogeneous simplicial complex
    nodes1.push_back(HDVF::Io_node_type({0,0,0}));
    nodes1.push_back(HDVF::Io_node_type({0,1,0}));
    nodes1.push_back(HDVF::Io_node_type({0,0,1}));
    nodes2 = nodes1;
    nodes2.push_back(HDVF::Io_node_type({0,1,1}));

    cells1.push_back(HDVF::Io_cell_type({1,2,3}));
    cells1.push_back(HDVF::Io_cell_type({1,2,3}));
    cells1.push_back(HDVF::Io_cell_type({1,2,3}));
    cells2 = cells1;
    cells2.push_back(HDVF::Io_cell_type({1,4}));

    std::cerr << "Test Mesh_object_io() with d>0" << std::endl;
    std::cerr << "Construct an object with d = 2 on a 'pure' mesh" << std::endl ;
    // The constructor should raise no exception
    try {
        HDVF::Mesh_object_io mio1(2,nodes1,cells1);
        assert(mio1.nodes.size() == 3);
        assert(mio1.cells.size() == 3);
    } catch (...) {
        std::cerr << "-- check_dimension() failed" << std::endl;
        assert(false);
    }
    std::cerr << "Construct an object with d = 2 on a non homogeneous mesh" << std::endl ;
    // The constructor *should* raise an exception
    try {
        HDVF::Mesh_object_io mio1(2,nodes2,cells2);
        assert(mio1.nodes.size() == 3);
        assert(mio1.cells.size() == 3);
    } catch (...) {
        std::cerr << "-- check_dimension() properly detected the dimension changes" << std::endl;
        assert(true);
    }

    std::cerr << "Construct an object with d = -2 on a non homogeneous mesh" << std::endl ;
    // The constructor should not raise an exception
    try {
        HDVF::Mesh_object_io mio1(-2,nodes2,cells2);
        assert(mio1.nodes.size() == 4);
        assert(mio1.cells.size() == 4);
    } catch (...) {
        std::cerr << "-- check_dimension() failed" << std::endl;
        assert(false);
    }

    // Test read_off()
    {
        std::cerr << "Test read_off with data/three_triangles.off (with ordering of vertices indices)" << std::endl;
        HDVF::Mesh_object_io mio2;
        mio2.read_off("data/three_triangles.off") ;
        assert(mio2.nvertices == 6);
        assert(mio2.ncells == 3);
        assert(mio2.cells.at(0) == HDVF::Io_cell_type({0,1,2}));
        assert(mio2.cells.at(1) == HDVF::Io_cell_type({2,3,4}));
        assert(mio2.cells.at(2) == HDVF::Io_cell_type({1,3,5}));
    }

    // Test read_simp()
    {
        std::cerr << "Test read_simp with data/simple_simplicial_complex.simp (with ordering of vertices indices)" << std::endl;
        HDVF::Mesh_object_io mio2;
        mio2.read_simp("data/simple_simplicial_complex.simp") ;
        assert(mio2.nvertices == 0);
        assert(mio2.ncells == 3);
        assert(mio2.cells.at(0) == HDVF::Io_cell_type({1,2,3}));
        assert(mio2.cells.at(1) == HDVF::Io_cell_type({2,4}));
        assert(mio2.cells.at(2) == HDVF::Io_cell_type({3,4}));
    }

    return 0;
}


