#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/HDVF/Surface_mesh_io.h>
#include <CGAL/Zp.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> SurfaceMesh;
typedef HDVF::Hdvf_traits_3<K> Traits;
typedef CGAL::Homological_discrete_vector_field::Surface_mesh_io<SurfaceMesh,Traits> SurfaceMeshIO;
typedef CGAL::Zp<5, char, true> Coefficient_ring;
typedef HDVF::Simplicial_chain_complex<Coefficient_ring,Traits> Complex;

namespace HDVF = CGAL::Homological_discrete_vector_field;

int main() {
    SurfaceMesh mesh;
    std::string filename("data/three_triangles.off");

    CGAL::IO::read_polygon_mesh(filename, mesh);

    SurfaceMeshIO mesh_io(mesh);

    // Build the associated Simplicial_chain_complex;

    Complex complex(mesh_io);

    // Test Surface_mesh_io IO chains

    SurfaceMesh::face_index fh(*(mesh.faces().begin()));
    SurfaceMesh::halfedge_index he(mesh.halfedge(fh));
    size_t i(CGAL::Homological_discrete_vector_field::halfedge_handle_to_cell_index(he, 2, mesh_io, complex));
    SurfaceMesh::halfedge_index he2(CGAL::Homological_discrete_vector_field::cell_index_to_halfedge_handle(i, 2, mesh_io, complex));

    std::cout << "he == he2: " << (he==he2) << std::endl;

    CGAL::OSM::Sparse_chain<Coefficient_ring, CGAL::OSM::COLUMN> chain(complex.number_of_cells(1));
    chain.set_coefficient(0, Coefficient_ring(1.));
    chain.set_coefficient(2, Coefficient_ring(1.));
    SurfaceMeshIO::Simplex_chain<Coefficient_ring> simplex_chain(HDVF::chain_to_simplex_chain(chain, 1, mesh_io, complex));

    return 0;
}


