#include <LOCAL/CGALH3.h>
#include <CGAL/copy_n.h>
#include <CGAL/random_selection.h>
#include <CGAL/point_generators_3.h>
#include <LEDA/param_handler.h>
#include <CGAL/Nef_S2/Sphere_map.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_io_parser.h>
#include "SM_overlayer.h"

typedef CGAL::Sphere_geometry<HKernel> SKernel;
typedef CGAL::Sphere_map<SKernel> Sphere_map;

typedef Sphere_map::Vertex_handle   Vertex_handle;
typedef Sphere_map::Halfedge_handle Halfedge_handle;
typedef Sphere_map::Face_handle     Face_handle;
typedef CGAL::SM_decorator<Sphere_map,SKernel> SM_decorator;
typedef CGAL::SM_overlayer<SM_decorator>  SM_overlayer;

typedef CGAL::Creator_uniform_3<RT,Point_3>  Creator;
typedef CGAL::Random_points_in_cube_3<Point_3,Creator> Point_source;
typedef SKernel::Sphere_point   SPoint;
typedef SKernel::Sphere_segment SSegment;

struct OR {
bool operator()(bool b1, bool b2) const
{ return b1||b2; }
};

int main(int argc, char **argv)
{
  CGAL::set_pretty_mode ( std::cerr );
  SETDTHREAD(911);
  // Sphere_geometry 11 
  // Sphere_geometry_OGL 13
  // Segment_overlay 23
  // SM_overlayer 53
  Point_3 p(0,0,0);
  Sphere_map E1,E2,E3;
  SM_decorator D1(E1),D2(E2),D3(E3);
  SM_overlayer O1(D1),O2(D2),O3(D3);

  int n;
  leda_string input_file;
  leda_param_handler H(argc,argv,".sg",false);
  H.add_parameter("number_of_lines:-n:int:10");
  H.add_parameter("file_of_segments:-i:string:");
  leda_param_handler::init_all();
  H.get_parameter("-n",n);
  H.get_parameter("-i",input_file);

  CGAL_assertion_msg(n>0,"-n value must be postive.");

  std::list<SSegment> L;
  if ( input_file == "" ) {
    Point_source S(5);
    Point_3 p1,p2,ph;
    Point_3 o(0,0,0);
    while ( n-- > 0 ) {
      do { ph = *S++; } while ( ph == o );
      Plane_3 h(o,(ph-CGAL::ORIGIN).direction()); 
      do { p1 = *S++; } 
      while ( p1 == o || h.projection(p1) == o );
      do { p2 = *S++; } 
      while ( p2 == o || h.projection(p2) == o );
      SPoint p3(h.projection(p1)),p4(h.projection(p2));
      int which = CGAL::default_random.get_int(0,3);
      if ( p3 == p4 ) which = 3;
      if ( p3 == p4.opposite() ) which = 2;
      switch ( which ) {
        case 0: // short
          L.push_back( SSegment(p3,p4,true) ); break;
        case 1: // long
          L.push_back( SSegment(p3,p4,false) ); break;
        case 2: // halfcircle
          L.push_back( SSegment(p3,p3.opposite(),h) ); break;
        case 3: // trivial
          L.push_back( SSegment(p3,p3,h) ); break;
      }
    }
  } else {
    std::ifstream input(input_file);
    CGAL_assertion_msg(input,"no input log.");
    SSegment s;
    while ( input >> s ) L.push_back(s);
  }
  std::ofstream output("smo-demo.log");
  std::list<SSegment>::iterator it;
  forall_iterators(it,L) output << *it;
  output << std::endl;
  output.close();
  std::list<SSegment> L1,L2;
  int b=0;
  forall_iterators(it,L) {
    if ( b == 0 ) L1.push_back(*it);
    else          L2.push_back(*it);
    b = 1-b;
  }
  O1.create_from_segments(L1.begin(),L1.end());
  O1.simplify(); // O1.dump(std::cerr);
  O2.create_from_segments(L2.begin(),L2.end());
  O2.simplify(); // O2.dump(std::cerr);
  O3.subdivide(E1,E2); 
  O3.select(OR());
  O3.simplify(); // O3.dump(std::cerr);
  return 0;
}


