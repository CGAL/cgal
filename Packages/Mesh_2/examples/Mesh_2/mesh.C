#include <CGAL/basic.h>
#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Mesh_2.h>
#include <CGAL/Mesh_face_base_2.h>
#include <CGAL/Mesh_default_traits_2.h>

typedef CGAL::Simple_cartesian<double> K1;
typedef CGAL::Filtered_kernel<K1> K2;
struct K : public K2 {};

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Mesh_default_traits_2<K> Meshtraits;
typedef CGAL::Constrained_Delaunay_triangulation_2<Meshtraits, Tds,
  CGAL::Exact_predicates_tag> Tr;

typedef CGAL::Mesh_2<Tr> Mesh;

typedef K::Point_2 Point;

Mesh mesh;

int main(int argc, char** argv)
{
  if(argc<2 || argc> 3)
    {
      std::cerr << "Usage: " << std::endl
		<< argv[0] << " input.poly [output.poly]" <<
	std::endl;
      return 1;
    };
  std::ifstream input(argv[1]);
  if(input)
    {
      mesh.read_poly(input);
      mesh.refine();
    }
  else
    {
      std::cerr << "Bad file: " << argv[1] << std::endl;
      return 1;
    }
  
  if(argc==2)
    mesh.write_poly(std::cout);
  else
    {
      std::ofstream output(argv[2]);
      mesh.write_poly(output);
    }
};
