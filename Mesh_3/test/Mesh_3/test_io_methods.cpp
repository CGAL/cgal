/*
  This test file tests <CGAL/Mesh_3/IO.h>, in ascii and binary mode.
*/

#include <CGAL/Volume_mesher_default_triangulation_3.h>

// c2t3
#include <CGAL/Complex_2_in_triangulation_3.h>

#include <iostream>  // std::cerr
#include <algorithm> // std::distance

#include <fstream>
#include <string>

#include <boost/array.hpp>

// default triangulation
typedef CGAL::Volume_mesher_default_triangulation_3 Tr;

// Kernel
typedef CGAL::Kernel_traits<Tr::Point>::Kernel K;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

enum Mode { ASCII, BINARY};

const std::string filename = "test_io_methods.out";

bool test(Mode mode)
{
  Tr tr;
  C2t3 c2t3(tr);

  boost::array<K::Point_3, 4> points;

  points[0]=K::Point_3(0, 0, 0);
  points[1]=K::Point_3(1, 0, 0);
  points[2]=K::Point_3(0, 1, 0);
  points[3]=K::Point_3(0, 0, 1);

  std::vector<Tr::Vertex_handle> v(points.size());
  typedef std::vector<Tr::Vertex_handle>::size_type size_type;

  for(size_type i = 0; i < points.size(); ++i)
    v[i] = tr.insert(points[i]);

  // Set 2 triangles in the complex
  // So it will be 2 facets, and 5 edges in the complex
  for(int i = 0; i < 2; ++i)
  {
    Tr::Cell_handle cell;
    int i1, i2, i3;
    tr.is_facet(v[i+1], v[(i+2)&3], v[(i+3)&3],
                cell,
                i1, i2, i3);
    c2t3.set_in_complex(Tr::Facet(cell, 6 - i1 - i2 - i3));
  }

  std::string filename2 = filename;
  if(mode == BINARY) 
    filename2+=".binary";

  std::ofstream out(filename2.c_str());

  if(mode == BINARY)
    CGAL::set_binary_mode(out);

  CGAL::Mesh_3::output_mesh(out, c2t3);

  out.close();

  std::ifstream in(filename2.c_str());
  CGAL::Mesh_3::input_mesh(in, c2t3,
                           true, // "true" means debugging on
                           &std::cerr);

  Tr::size_type number_of_edges = std::distance(c2t3.edges_begin(),
                                                c2t3.edges_end());

  Tr::size_type number_of_facets = std::distance(c2t3.facets_begin(),
                                                 c2t3.facets_end());

  Tr::size_type number_of_boundary_edges = 
    std::distance(c2t3.boundary_edges_begin(),
                  c2t3.boundary_edges_end());

  std::cerr << "Number of edges: " << number_of_edges
            << "\nNumber of facets: " << number_of_facets
            << "\nNumber of boundary edges: " << number_of_boundary_edges
            << "\n";

  std::cerr << "checking points coordinates... ";
  bool check = true;
  for(size_type i = 0; i < points.size(); ++i)
  {
    Tr::Vertex_handle dummy_v;
    check = check && tr.is_vertex(points[i], dummy_v);
  }
  if(check)
    std::cerr << "ok\n";
  else
    std::cerr << "error\n";

  // Excepted results:
  return 
    check &&
    ( number_of_boundary_edges == 4 ) &&
    ( number_of_edges == 5 ) &&
    ( number_of_facets == 2 );  
} // end test

int main()
{
  std::cerr << "testing ASCII mode...\n";
  bool result = test(ASCII);
  std::cerr << "\ntesting BINARY mode...\n";
  result = result && test(BINARY);

  return (result ? EXIT_SUCCESS : EXIT_FAILURE);
}
