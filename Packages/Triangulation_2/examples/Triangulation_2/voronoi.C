#include <CGAL/basic.h>
#include <fstream>


// Define shorter names to please linker (g++/egcs)
#define Cartesian Cart

#include <CGAL/Cartesian.h>

#include <CGAL/squared_distance_2.h>   // to avoid a g++ problem
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Segment_2.h>

#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

class Options {
public:
    Options()
        :  file_input(false)

    {}

    char program[100];
    char fname[100];
    bool file_input;
};



void usage(char* program)
{
  std::cerr << "\nNAME\n     "
       << program << " - Triangulation of a point set\n\n";
  std::cerr << "SYNOPSIS\n     "
       << program << " [-file fname]\n";

  std::cerr << "\nDESCRIPTION\n"
       << "     Triangulates a point set that comes from a file or stdin.\n";
  std::cerr << "\nOPTIONS\n"
       << "     All options can be abbreviated by their first character\n\n";
}


bool
parse(int argc, char* argv[], Options &opt)
{
    strcpy(opt.program, argv[0]);
    --argc;
    argv++;

    while ((argc > 0) && (argv[0][0] == '-')){
        if ((!strcmp(argv[0], "-f")) || (!strcmp(argv[0], "-file"))) {
          strcpy(opt.fname, argv[1]);
          opt.file_input = true;
          argv += 2;
          argc -= 2;
      }
      else if ((!strcmp(argv[0], "-?")) ||
               (!strcmp(argv[0], "-h")) ||
               (!strcmp(argv[0], "-help"))) {
          usage(opt.program);
          return false;
      }
      else {
          std::cerr << "Unrecognized option " << argv[0] << std::endl;
          usage(opt.program);
          return false;
      }
    }
  if(argc > 0){
      std::cerr << "Unrecognized option " << argv[0] << std::endl;
      usage(opt.program);
      return false;
  }
  return true;
}


typedef double coord_type;
typedef CGAL::Cartesian<coord_type>  Rpst;

typedef CGAL::Triangulation_euclidean_traits_2<Rpst> Gt;
typedef Gt::Point_2       Point;
typedef Gt::Segment_2     Segment;
typedef Gt::Ray_2         Ray;
typedef Gt::Line_2        Line;
typedef Gt::Triangle_2    Triangle;

typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds>  Triangulation;

typedef Triangulation::Face  Face;
typedef Triangulation::Vertex Vertex;
typedef Triangulation::Face_handle  Face_handle;
typedef Triangulation::Vertex_handle Vertex_handle;

typedef Triangulation::Face_circulator  Face_circulator;
typedef Triangulation::Vertex_circulator  Vertex_circulator;

typedef Triangulation::Locate_type Locate_type;

typedef Triangulation::Face_iterator  Face_iterator;
typedef Triangulation::Vertex_iterator  Vertex_iterator;
typedef Triangulation::Edge_iterator  Edge_iterator;
typedef Triangulation::Line_face_circulator  Line_face_circulator;




void input_from_file(Triangulation &T,
                     const Options& opt)
{
    if(! opt.file_input){
        return;
    }

    std::ifstream is(opt.fname);
    CGAL::set_ascii_mode(is);

    int n;
    is >> n;
    std::cout << "Reading " << n << " points" << std::endl;

    std::istream_iterator<Point> begin(is);
    std::istream_iterator<Point> end;
    T.insert(begin, end);
}

int
main(int argc, char* argv[])
{
    Options opt;
    parse(argc, argv, opt);

    Triangulation T;

    input_from_file(T, opt);
    
    std::cout << "Vertices :" << std::endl << "=============" << std::endl;
    Face_iterator fit, fbegin=T.faces_begin(), fend=T.faces_end(); 
    for (fit=fbegin; fit != fend; ++fit)
	std::cout << T.dual(fit) << std::endl;

    std::cout << "Segment :" << std::endl << "=============" << std::endl;
      {
	Edge_iterator eit, ebegin=T.edges_begin(), eend=T.edges_end(); 
	for (eit=ebegin; eit != eend; ++eit)
	  {
	    CGAL::Object o = T.dual(eit);
	    Triangulation::Segment s;
	    if (CGAL::assign(s,o)) std::cout << s << std::endl;
	  }
      }

    std::cout << "Rays :" << std::endl << "=============" << std::endl;
      {
	Edge_iterator eit, ebegin=T.edges_begin(), eend=T.edges_end(); 
	for (eit=ebegin; eit != eend; ++eit)
	  {
	    CGAL::Object o = T.dual(eit);
	    Ray r;
	    if (CGAL::assign(r,o)) std::cout << r << std::endl;
	  }
      }

    return 0;
}
