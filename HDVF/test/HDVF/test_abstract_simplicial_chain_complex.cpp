#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Abstract_simplicial_chain_complex.h>
#include <CGAL/OSM/Sparse_matrix.h>

typedef CGAL::HDVF::Z2 CoefficientType;
typedef CGAL::HDVF::Abstract_simplicial_chain_complex<CoefficientType> ChainComplexType;
typedef CGAL::OSM::Sparse_matrix<CoefficientType, CGAL::OSM::COLUMN> Column_matrix;
typedef CGAL::OSM::Sparse_chain<CoefficientType, CGAL::OSM::COLUMN> Column_chain;
typedef CGAL::OSM::Sparse_chain<CoefficientType, CGAL::OSM::ROW> Row_chain;

int main() {
    // Test constructor (from simple_simplicial_complex.simp)
    std::cerr << "Test constructor (from simple_simplicial_complex.simp)" << std::endl ;
    CGAL::HDVF::Mesh_object_io mesh;
    mesh.read_simp("data/simple_simplicial_complex.simp");
    
    ChainComplexType complex(mesh);
    
    std::cerr << "-- Check number_of_cells" << std::endl ;
    assert(complex.number_of_cells(0) == 4);
    assert(complex.number_of_cells(1) == 5);
    assert(complex.number_of_cells(2) == 1);
    
    std::vector<Column_matrix> bnds(complex.boundary_matrices());
    
//  Generate data
//    CGAL::OSM::write_matrix(bnd.at(0), "data/test_abstract_simplicial_chain_complex/bnd0.osm");
//    CGAL::OSM::write_matrix(bnd.at(1), "data/test_abstract_simplicial_chain_complex/bnd1.osm");
//    CGAL::OSM::write_matrix(bnd.at(2), "data/test_abstract_simplicial_chain_complex/bnd2.osm");
    
    std::cerr << "-- Check boundary matrices" << std::endl ;
    Column_matrix bnd;
    bool test;
    
    // Dim 0
    CGAL::OSM::read_matrix(bnd, "data/test_abstract_simplicial_chain_complex/bnd0.osm");
    test = (bnd == bnds.at(0));
    std::cerr << "--- Dimension 0: " << test << std::endl ;
    assert(test) ;
    
    // Dim 1
    CGAL::OSM::read_matrix(bnd, "data/test_abstract_simplicial_chain_complex/bnd1.osm");
    test = (bnd == bnds.at(1));
    std::cerr << "--- Dimension 1: " << test << std::endl ;
    assert(test) ;
    
    // Dim 2
    CGAL::OSM::read_matrix(bnd, "data/test_abstract_simplicial_chain_complex/bnd2.osm");
    test = (bnd == bnds.at(2));
    std::cerr << "--- Dimension 2: " << test << std::endl ;
    assert(test) ;
    
    std::cerr << "-- Check d and cod" << std::endl ;
    
    Column_chain c1, c2;
    Row_chain r1, r2;
    // Dim 0
    CGAL::OSM::read_matrix(bnd, "data/test_abstract_simplicial_chain_complex/bnd0.osm");
    c1 = CGAL::OSM::cget_column(bnd, 0);
    c2 = complex.d(0,0);
    test = (c1 == c2);
    std::cerr << "--- Dimension 0 (d(0)): " << test << std::endl ;
    assert(test) ;
    r1 = CGAL::OSM::get_row(bnd, 0);
    r2 = complex.cod(0,0);
    test = (c1 == c2);
    std::cerr << "--- Dimension 0 (cod(0))): " << test << std::endl ;
    assert(test) ;
    
    // Dim 0
    CGAL::OSM::read_matrix(bnd, "data/test_abstract_simplicial_chain_complex/bnd0.osm");
    c1 = CGAL::OSM::cget_column(bnd, 0);
    c2 = complex.d(0,0);
    test = (c1 == c2);
    std::cerr << "--- Dimension 0 (d(0)): " << test << std::endl ;
    assert(test) ;
    r1 = CGAL::OSM::get_row(bnd, 0);
    r2 = complex.cod(0,0);
    test = (c1 == c2);
    std::cerr << "--- Dimension 0 (cod(0))): " << test << std::endl ;
    assert(test) ;
    
    // Dim 1
    CGAL::OSM::read_matrix(bnd, "data/test_abstract_simplicial_chain_complex/bnd1.osm");
    c1 = CGAL::OSM::cget_column(bnd, 0);
    c2 = complex.d(0,1);
    test = (c1 == c2);
    std::cerr << "--- Dimension 1 (d(0)): " << test << std::endl ;
    assert(test) ;
    r1 = CGAL::OSM::get_row(bnd, 0);
    r2 = complex.cod(0,1);
    test = (c1 == c2);
    std::cerr << "--- Dimension 1 (cod(0))): " << test << std::endl ;
    assert(test) ;
    
    // Dim 2
    CGAL::OSM::read_matrix(bnd, "data/test_abstract_simplicial_chain_complex/bnd2.osm");
    c1 = CGAL::OSM::cget_column(bnd, 0);
    c2 = complex.d(0,2);
    test = (c1 == c2);
    std::cerr << "--- Dimension 2 (d(0)): " << test << std::endl ;
    assert(test) ;
    r1 = CGAL::OSM::get_row(bnd, 0);
    r2 = complex.cod(0,2);
    test = (c1 == c2);
    std::cerr << "--- Dimension 2 (cod(0))): " << test << std::endl ;
    assert(test) ;
    
    return 0;
}


