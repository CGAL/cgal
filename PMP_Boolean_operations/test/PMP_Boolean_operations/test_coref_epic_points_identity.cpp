#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;

int main(int argc, char** argv)
{
  for(int i=0; i< (argc-1)/2;++i)
  {
    Surface_mesh sm1, sm2;
    std::ifstream input(argv[2*i+1]);
    assert(input);
    input >> sm1;
    input.close();
    input.open(argv[2*(i+1)]);
    assert(input);
    input >> sm2;
    input.close();

    CGAL::Polygon_mesh_processing::corefine(sm1, sm2);

    assert(sm1.is_valid());
    assert(sm2.is_valid());

    std::size_t nbv1=sm1.number_of_vertices(), nbv2=sm2.number_of_vertices();
    std::size_t nbe1=sm1.number_of_edges(), nbe2=sm2.number_of_edges();
    std::size_t nbf1=sm1.number_of_faces(), nbf2=sm2.number_of_faces();

    CGAL::Polygon_mesh_processing::corefine(sm1, sm2);

    assert( nbv1 == sm1.number_of_vertices() );
    assert( nbv2 == sm2.number_of_vertices() );
    assert( nbe1 == sm1.number_of_edges() );
    assert( nbe2 == sm2.number_of_edges() );
    assert( nbf1 == sm1.number_of_faces() );
    assert( nbf2 == sm2.number_of_faces() );
  }
}
