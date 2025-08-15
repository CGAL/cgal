#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/HDVF/Mesh_object_io.h>

int main() {
// Test constructor for a mesh (all cells have the same dimension d -> "pure" mesh
    
    std::vector<std::vector<double> > nodes1, nodes2;
    std::vector<CGAL::HDVF::IOCellType> cells1, cells2;
    
    // (nodes1, cells1) : "pure" mesh
    // (nodes2, cells2) : heterogeneous simplicial complex
    nodes1.push_back({0,0,0});
    nodes1.push_back({0,1,0});
    nodes1.push_back({0,0,1});
    nodes2 = nodes1;
    nodes2.push_back({0,1,1});
    
    cells1.push_back(CGAL::HDVF::IOCellType({1,2,3}));
    cells1.push_back(CGAL::HDVF::IOCellType({1,2,3}));
    cells1.push_back(CGAL::HDVF::IOCellType({1,2,3}));
    cells2 = cells1;
    cells2.push_back(CGAL::HDVF::IOCellType({1,4}));
    
    std::cerr << "Test Mesh_object_io() with d>0" << std::endl;
    std::cerr << "Construct an object with d = 2 on a 'pure' mesh" << std::endl ;
    // The constructor should raise no exception
    try {
        CGAL::HDVF::Mesh_object_io mio1(2,nodes1,cells1);
        assert(mio1.nodes.size() == 3);
        assert(mio1.cells.size() == 3);
    } catch (...) {
        std::cerr << "-- check_dimension() failed" << std::endl;
        assert(false);
    }
    std::cerr << "Construct an object with d = 2 on a non homogeneous mesh" << std::endl ;
    // The constructor *should* raise an exception
    try {
        CGAL::HDVF::Mesh_object_io mio1(2,nodes2,cells2);
        assert(mio1.nodes.size() == 3);
        assert(mio1.cells.size() == 3);
    } catch (...) {
        std::cerr << "-- check_dimension() properly detected the dimension changes" << std::endl;
        assert(true);
    }
    
    std::cerr << "Construct an object with d = -2 on a non homogeneous mesh" << std::endl ;
    // The constructor should not raise an exception
    try {
        CGAL::HDVF::Mesh_object_io mio1(-2,nodes2,cells2);
        assert(mio1.nodes.size() == 4);
        assert(mio1.cells.size() == 4);
    } catch (...) {
        std::cerr << "-- check_dimension() failed" << std::endl;
        assert(false);
    }
    
    return 0;
}


