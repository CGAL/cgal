#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <list>

#include <CGAL/Simple_cartesian.h>

#include <CEP/Visibility_complex/Visibility_complex_segment_traits.h>
#include <CEP/Visibility_complex/Visibility_complex_2.h>

#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Ostream_iterator.h>

typedef int                                          FT;
typedef CGAL::Simple_cartesian <FT>                  Rep;
typedef CGAL::Visibility_complex_segment_traits<Rep> Gt;
typedef CGAL::Visibility_complex_2<Gt>               Visibility_complex;
typedef Visibility_complex::Vertex                   Vertex;
typedef Gt::Disk                                     Segment;

typedef CGAL::Window_stream                          Window_stream;

int main()
{
    // Reading segments from file
    std::list<Segment> O;
    std::ifstream ifs("segments");
    std::istream_iterator<Segment> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,std::back_inserter(O));

    // Computing visibility graph
    Visibility_complex V(O.begin(),O.end());

    // Open a LEDA window
    Window_stream W(500,400); // physical window size
    W.init(0, 800., 0.);   // logical window size
    CGAL::cgalize(W);
    W.display();

    // Output the bitangents in RED
    W << CGAL::RED;
    std::copy(V.vertices_begin(),V.vertices_end(),
	      CGAL::Ostream_iterator<Vertex,Window_stream>(W));

    // Output the segments in BLACK
    W << CGAL::BLACK;
    std::copy(O.begin(),O.end(),
	      CGAL::Ostream_iterator<Segment,Window_stream>(W));

    // Show window
    std::cout << "Click in window to quit" << std::endl;
    W.read_mouse();
    return 0;
}
