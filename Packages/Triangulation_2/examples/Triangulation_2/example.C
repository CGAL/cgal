
#include <CGAL/basic.h>
#include <cassert>
#include <iostream>
#include <fstream>


// Define shorter names to please linker (g++/egcs)
#define Cartesian Ca
#define Homogeneous Ho

#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>
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
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Triangulation_2<Gt,Tds>  Triangulation;
//typedef CGAL::Delaunay_triangulation_2<Gt,Tds>  Delaunay_triangulation_2;

typedef Gt::Point_2          Point;
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
//     Point p;
//     for(int i=0; i<n; i++) {
//       is >> p; 
//       T.insert(p);
//     }
}

void input_from_range(Triangulation &T)
{
    std::list<Point> L;
    L.push_front(Point(0,0));
    L.push_front(Point(1,0));
    L.push_front(Point(1,1));

    int n =T.insert(L.begin(), L.end());
    std::cout << n << " points inserted from a list." << std::endl;

    std::vector<Point> V(3);
    V[0] = Point(0, 0);
    V[1] = Point(0.4, 0.4);
    V[2] = Point(0.3, 0.3);

    n = T.insert(V.begin(), V.end());
    std::cout << n << " points inserted from a vector." << std::endl;
}


void faces_along_line(Triangulation &T)
{
    Point p(0.2, 0.6), q(0.7, 0.4);

    Face_handle f = T.locate(p);
    Line_face_circulator lfc = T.line_walk(p, q, f),
                         done(lfc);
    if(lfc == (CGAL_NULL_TYPE) NULL){
        std::cout << "Line does not intersect convex hull" << std::endl;
    } else {
        int count = 0;
        do{
            if(! T.is_infinite(lfc)){
                count++;
            }
        }while(++lfc != done);
        std::cout << "The line intersects " 
		  << count << " finite faces" << std::endl;
    }
}

void convex_hull(Triangulation &T)
{
    Point p, q;

    Vertex_circulator chc = T.infinite_vertex()->incident_vertices(),
                      done(chc);
    if(chc == (CGAL_NULL_TYPE)NULL) {
        std::cout << "convex hull is empty" << std::endl;
    } else {
        p = chc->point();
        do {
            --chc;
            q = chc->point();
            p = q;
        } while(chc != done);
        
    }
}

void fileIO(Triangulation &T,
            const Options& opt)
{
    std::cout << "The triangulation will be written to a file and read again\n";
    {
        std::ofstream out("tr");
        CGAL::set_ascii_mode(out);
        out << T << std::endl;
    }
    Triangulation T2;

    std::ifstream in("tr");
    CGAL::set_ascii_mode(in);
    in >> T2;
    assert( T2.is_valid() );
}

int
main(int argc, char* argv[])
{
    Options opt;
    parse(argc, argv, opt);

    Triangulation T;

    input_from_range(T);
    input_from_file(T, opt);
 
    std::list<Point> L;
    L.push_front(Point(0,0));
    L.push_front(Point(1,0));
    L.push_front(Point(1,1));

    int n = T.insert(L.begin(), L.end());
    std::cout << n << " points inserted from a list." << std::endl;
    

    T.insert(Point(0,0));
    T.is_valid();
    
    faces_along_line(T);

    convex_hull(T);

    fileIO(T, opt);

    return 0;
}

