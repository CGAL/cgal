#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/IO/Istream_iterator.h>
#include <CGAL/IO/Window_stream.h>
#include <iostream>
#include <algorithm>

typedef CGAL::Point_2< CGAL::Cartesian<double> >           Point;
typedef CGAL::Istream_iterator<Point, CGAL::Window_stream> Iterator;

void init_window( leda_window& W) {
    CGAL::cgalize( W);
    W.set_fg_color( leda_green);
    W.display();
    W.init(-1.0, 1.0, -1.0);
}

int main () {
    CGAL::Window_stream window( 512, 512);
    init_window(window);
    std::copy( Iterator(window), Iterator(),
               std::ostream_iterator<Point>(std::cout,"\n"));
    return 0;
}
