#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/HDVF/Triangulation_3_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

typedef int Coefficient_ring;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
 
typedef CGAL::Triangulation_3<K>      Triangulation;
 
typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;

typedef HDVF::Hdvf_traits_3<K> Traits;
typedef HDVF::Triangulation_3_io<Triangulation,Traits> Triangulation_3_io;
typedef HDVF::Simplicial_chain_complex<Coefficient_ring, Traits> Complex;

int main() {
    // construction from a list of points :
    std::list<Point> L;
    L.push_front(Point(0,0,0));
    L.push_front(Point(1,0,0));
    L.push_front(Point(0,1,0));
    L.push_front(Point(0,0,1));
    L.push_front(Point(1,1,1));
    L.push_front(Point(2,2,2));
   
    Triangulation T(L.begin(), L.end());

    Triangulation_3_io tri3_io(T);
    Complex complex(tri3_io);
    
    std::cout << "Complex built: " << complex << std::endl;
    Complex::chain_complex_to_vtk(complex, "tmp/complex.vtk");

    return 0;
}


