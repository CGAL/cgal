#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <list>

#include <CGAL/Simple_cartesian.h>

#include <CEP/Visibility_complex/Shortest_path_segment_traits.h>
#include <CEP/Visibility_complex/Shortest_path_2.h>
//
// Number type used for input
typedef int                                         NT;  
// Number type used when computing distances
typedef double                                      SNT; 

typedef CGAL::Simple_cartesian<NT>                  SC;
struct Rep : public SC {};

typedef CGAL::Shortest_path_segment_traits<Rep,SNT> Gt;
typedef Gt::Disk                                    Segment;
typedef Gt::Point_2                                 Point;
typedef CGAL::Visibility_complex_2<Gt>              VC;
typedef VC::Vertex                                  Vertex;

int main()
{
    // Reading segments from file
    std::list<Segment> O;
    std::ifstream ifs("segments");
    std::istream_iterator<Segment> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,std::back_inserter(O));

    // Two corners of the window
    Point p(0,0);
    Point q(1000,800);

    shortest_path_2(O.begin(),O.end(), 
		    p,q, 
		    std::ostream_iterator<Segment>(std::cout,"\n"),
		    Gt());

    return 0;
}
