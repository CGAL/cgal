#ifndef CGAL_USE_LEDA
#include <iostream>
int main(){ std::cout << "This demo needs LEDA" << std::endl; return 0;}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/IO/Istream_iterator.h>
#include <CGAL/IO/Window_stream.h>
#include <iostream>
#include <algorithm>

typedef CGAL::Cartesian<double>::Point_2                   Point;
typedef CGAL::Istream_iterator<Point, CGAL::Window_stream> Iterator;

#ifdef CGAL_USE_CGAL_WINDOW
#define leda_window CGAL::window
#define leda_green  CGAL::green
#endif

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

#endif
