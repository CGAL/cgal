#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <list>

#include <CGAL/Simple_cartesian.h>

#include <CEP/Visibility_complex/Visibility_complex_point_traits.h>
#include <CEP/Visibility_complex/Visibility_complex_2.h>

#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Ostream_iterator.h>

typedef int                                          FT;
typedef CGAL::Simple_cartesian <FT>                  Rep;
typedef CGAL::Visibility_complex_point_traits<Rep> Gt;
typedef CGAL::Visibility_complex_2<Gt>               Visibility_complex;
typedef Visibility_complex::Vertex                   Vertex;
typedef Visibility_complex::Antichain                Antichain;
typedef Gt::Disk                                     Point;

typedef CGAL::Window_stream                          Window_stream;

int main()
{
    // Reading segments from file
    std::list<Point> O;
    std::ifstream ifs("segments");
    std::istream_iterator<Point> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,std::back_inserter(O));

    // Computing visibility graph
    Visibility_complex V(O.begin(),O.end());

    // Computing the first antichain
    std::list<Vertex> empty;
    Antichain A(O.begin(),O.end(), empty.begin(), empty.end());

    // Open a LEDA window
    Window_stream W(500,400); // physical window size
    W.init(0, 800., 0.);   // logical window size
    CGAL::cgalize(W);
    W.display();

    // Output the bitangents in RED
    W << CGAL::RED;
    std::copy(A.vertices_begin(),A.vertices_end(),
	      CGAL::Ostream_iterator<Vertex,Window_stream>(W));

    // Output the segments in BLACK
    W << CGAL::BLACK;
    std::copy(O.begin(),O.end(),
	      CGAL::Ostream_iterator<Point,Window_stream>(W));

    // Show window
    std::cout << "Click in window to quit" << std::endl;
    W.read_mouse();
    return 0;
}
