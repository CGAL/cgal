#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <CGAL/Timer.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_area_criteria_2.h>

#include <CGAL/IO/File_poly.h>

typedef CGAL::Simple_cartesian<double> K1;
typedef CGAL::Filtered_kernel<K1> K2;

#ifdef USE_FILTRED
struct K: public K2 {};
#else
struct K : public K1 {};
#endif

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> Tr;
typedef CGAL::Delaunay_mesh_area_criteria_2<Tr> Mesh_criteria;

typedef CGAL::Delaunay_mesher_2<Tr, Mesh_criteria> Mesher;

typedef K::Point_2 Point;

void usage(char** argv)
{
  std::cerr << "Usage: " << std::endl
	    << argv[0] << " [-Q] [-a area] input.poly" << std::endl;
}

int main(int argc, char** argv)
{
  int arg_count = 1;
  bool terminal_output = true;

  Tr tr;

  Mesher mesher(tr);

  Mesh_criteria criteria;

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
		criteria.set_area_bound(area_bound);
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
      CGAL::Timer t;
      t.start();
      CGAL::read_triangle_poly_file(tr, input);
      t.stop();
      if(terminal_output)
	std::cout << "Delaunay time: " << t.time() << std::endl;
      
      mesher.set_criteria(criteria);
      
      t.reset(); t.start();
      mesher.refine_mesh();
      t.stop();
      if(terminal_output)
	std::cout << "Meshing time: " << t.time() << std::endl;
    }
  else
    {
      std::cerr << "Bad file: " << argv[arg_count] << std::endl;
      usage(argv);
      return 1;
    }
  
  if(argc==arg_count+1)
    {
    }
  else
    {
      std::ofstream output(argv[arg_count+1]);
      CGAL::write_triangle_poly_file(tr, output);
    }
  if(terminal_output)
    std::cout 
      << "Mesh points: " << tr.number_of_vertices() << std::endl
      << "Mesh triangles: " << tr.number_of_faces () << std::endl;

  return 0;
};
