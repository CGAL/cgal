#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <list>

#include <CGAL/Simple_cartesian.h>

#include <CEP/Visibility_complex/Shortest_path_segment_traits.h>
#include <CEP/Visibility_complex/Shortest_path_2.h>

#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Ostream_iterator.h>
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

typedef CGAL::Window_stream                         Window_stream;

int main()
{
    // Reading segments from file
    std::list<Segment> O;
    std::ifstream ifs("segments");
    std::istream_iterator<Segment> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,std::back_inserter(O));

    // Open a LEDA window
    Window_stream W(500,400); // physical window size
    W.init(0, 800., 0.);   // logical window size
    CGAL::cgalize(W);
    W.display();

    // Two corners of the window
    Point p(0,0);
    Point q(1000,800);

    // Otput the segments of the shortest path from p to q in reverse order.
    W << CGAL::RED;
    W.set_line_width(3);
    shortest_path_2(O.begin(),O.end(), 
		    p,q, 
		    CGAL::Ostream_iterator<Segment,Window_stream>(W),
		    Gt());

    // Output the segments in BLACK
    W << CGAL::BLACK;
    W.set_line_width(1);
    std::copy(O.begin(),O.end(),
	      CGAL::Ostream_iterator<Segment,Window_stream>(W));

    // Show window
    std::cout << "Click in window to quit" << std::endl;
    W.read_mouse();
    return 0;
}
