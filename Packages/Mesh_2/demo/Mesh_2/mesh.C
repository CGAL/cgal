#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Delaunay_mesh_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_area_traits_2.h>

typedef CGAL::Simple_cartesian<double> K1;
typedef CGAL::Filtered_kernel<K1> K2;
struct K : public K2 {};

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Delaunay_mesh_area_traits_2<K> Meshtraits;
typedef CGAL::Constrained_Delaunay_triangulation_2<Meshtraits, Tds,
  CGAL::Exact_predicates_tag> Tr;

typedef K::Point_2 Point;

void usage(char** argv)
{
  std::cerr << "Usage: " << std::endl
	    << argv[0] << " [-Q] [-a area] input.poly [output.poly]" << std::endl;
}

int main(int argc, char** argv)
{
  int arg_count = 1;
  bool terminal_output = true;
  Meshtraits traits;
  Tr t;

  if(argc < 2)
    {
      usage(argv);
      return 1;
    }

  while(argv[arg_count][0] == '-' && argv[arg_count] != "--")
    {
      if(std::string(argv[arg_count]) == "-Q")
	terminal_output = false;
      else if(std::string(argv[arg_count]) == "-a")
	{
	  double area_bound;
	    if( (argc > arg_count+1) && 
		std::istringstream(argv[arg_count+1]) >> area_bound )
	      {
		traits.set_area_bound(area_bound);
		++arg_count;
	      }
	    else
	      {
		std::cerr << "The area " << argv[arg_count+1]
			  << " is not a double." << std::endl;
		usage(argv);
		return 1;
	      }
	}
      else
	{
	  std::cerr << "Unknown option " << argv[arg_count] << std::endl;
	  usage(argv);
	  return 1;
	}
      ++arg_count;
    }
  if(argv[arg_count] == "--")
    ++arg_count;

  if(argc < arg_count+1 || argc > arg_count+2)
    {
      usage(argv);
      return 1;
    };

  std::ifstream input(argv[arg_count]);
  if(input)
    {
      CGAL::read_poly(t, input);
      CGAL::refine_Delaunay_mesh_2(t, traits);

      if(argc==arg_count+1)
	{
	  if(terminal_output)
	    CGAL::write_poly(t, std::cout);
	}
      else
	{
	  std::ofstream output(argv[arg_count+1]);
	  CGAL::write_poly(t, output);
	}
      if(terminal_output)
	std::cerr 
	  << "Mesh points: " << t.number_of_vertices() << std::endl
	  << "Mesh triangles: " << t.number_of_faces () << std::endl;
      
    }
  else
    {
      std::cerr << "Bad file: " << argv[arg_count] << std::endl;
      usage(argv);
      return 1;
    }
  
  return 0;
};
