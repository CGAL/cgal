#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Conform_2.h>
#include <CGAL/Read_write.h>
#include <CGAL/Mesh_default_traits_2.h>

typedef CGAL::Simple_cartesian<double> K1;
typedef CGAL::Filtered_kernel<K1> K2;
struct K : public K2 {};

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Mesh_default_traits_2<K> Meshtraits;
typedef CGAL::Constrained_Delaunay_triangulation_2<Meshtraits, Tds,
  CGAL::Exact_predicates_tag> Tr;

typedef CGAL::Conform_triangulation_2<Tr> Conform;

typedef K::Point_2 Point;

Conform conform;

void usage(char** argv)
{
  std::cerr << "Usage: " << std::endl
	    << argv[0] << " [-Q] input.poly [output.poly]" << std::endl;
}

int main(int argc, char** argv)
{
  int arg_count = 1;
  bool terminal_output = true;

  if(argc < 2)
    {
      usage(argv);
      return 1;
    }

  while(argv[arg_count][0] == '-' && argv[arg_count] != "--")
    {
      if(std::string(argv[arg_count]) == "-Q")
	terminal_output = false;
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
      read_poly(conform, input);

      conform.init(CGAL::Gabriel_conform_policy_2());
      
      typedef Conform::cluster_vertices_iterator Vh_iterator;
      for(Vh_iterator it = conform.clusters_vertices_begin();
	  it != conform.clusters_vertices_end();
	  ++it)
	{
	  typedef std::pair<Conform::vertices_in_cluster_iterator,
	    Conform::vertices_in_cluster_iterator> It_pair;

	  int n = conform.number_of_clusters_at_vertex(*it);

	  std::cout << "Point(" << (*it)->point() << ")" << std::endl
		    << "  " << n
		    << " cluster(s)." << std::endl;
	  
	  for(int i = 0; i<n; ++i)
	    {
	      std::cout << "  cluster " << i << ":" << std::endl;
	      const It_pair& it_pair =
		conform.vertices_in_cluster_sequence(*it, i);
	      for(Conform::vertices_in_cluster_iterator
		    vic_it = it_pair.first;
		  vic_it != it_pair.second;
		  ++vic_it)
		std::cout << "    " << (*vic_it)->point() << std::endl;
	    }
	};

      conform.gabriel_conform();
    }
  else
    {
      std::cerr << "Bad file: " << argv[arg_count] << std::endl;
      usage(argv);
      return 1;
    }
  
  if(argc==arg_count+1)
    {
      if(terminal_output)
	write_poly(conform, std::cout);
    }
  else
    {
      std::ofstream output(argv[arg_count+1]);
      write_poly(conform, output);
    }
  if(terminal_output)
    std::cerr 
      << "Mesh points: " << conform.number_of_vertices() << std::endl
      << "Mesh triangles: " << conform.number_of_faces () << std::endl;

  return 0;
};
