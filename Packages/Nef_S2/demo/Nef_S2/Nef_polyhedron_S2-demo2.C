#include <LOCAL/CGALH3.h>
#include <CGAL/random_selection.h>
#include <CGAL/point_generators_3.h>
#include "Nef_polyhedron_S2.h"
#include "Nef_S2/Nef_polyhedron_S2_OGLUT_stream.h"
#include <LEDA/param_handler.h>

//#define TRACEV(t) std::cerr << #t << " = " << t <<std::endl

typedef CGAL::Nef_polyhedron_S2<HKernel> Nef_polyhedron_S2;
typedef Nef_polyhedron_S2::Sphere_point   Sphere_point;
typedef Nef_polyhedron_S2::Sphere_segment Sphere_segment;
typedef Nef_polyhedron_S2::Sphere_circle  Sphere_circle;
typedef Nef_polyhedron_S2::Explorer Explorer;

typedef CGAL::Creator_uniform_3<RT,Point_3>  Creator;
typedef CGAL::Random_points_in_cube_3<Point_3,Creator> Point_source;


int main(int argc, char **argv)
{
  CGAL::set_pretty_mode ( std::cerr );
  SETDTHREAD(911);
  // Sphere_geometry 11 
  // Sphere_geometry_OGL 13
  // Segment_overlay 23
  // SM_overlayer 53
  Point_3 p(0,0,0);

  leda_string input_file;
  leda_param_handler H(argc,argv,".nd",false);
  H.add_parameter("file_of_circles:-i:string:");
  leda_param_handler::init_all();
  H.get_parameter("-i",input_file);

  std::list<Sphere_circle> L;
  if ( input_file == "" ) { // create random input:
    L.push_back( Sphere_circle(1,0,0) );
    L.push_back( Sphere_circle(0,1,0) );
    L.push_back( Sphere_circle(0,0,1) );
    L.push_back( Sphere_circle(1,1,1) );
    L.push_back( Sphere_circle(-1,1,1) );
    L.push_back( Sphere_circle(1,-1,1) );
    L.push_back( Sphere_circle(1,1,-1) );
  } else { // read input from file:
    std::ifstream input(input_file);
    CGAL_assertion_msg(input,"no input log.");
    Sphere_circle c;
    while ( input >> c ) L.push_back(c);
  }

  // output log:
  std::ofstream output("nef2.log");
  std::list<Sphere_circle>::iterator it;
  CGAL_forall_iterators(it,L) output << *it << ' ';
  output << std::endl;
  output.close();

  // partition input into two lists
  Nef_polyhedron_S2 Ni, N;
  bool first(false);
  CGAL_forall_iterators(it,L) {
    if ( first ) {
      N = Nef_polyhedron_S2(*it);
      first = false;
    } else {
      Ni = Nef_polyhedron_S2(*it);
      N = N ^ Ni;
    }
  }

  //std::cerr << Ni << N std::endl;
  CGAL::ogl << N; 
  CGAL::ogl << "Symmetric Difference"; 
  CGAL::ogl.display();
  return 0;

}


