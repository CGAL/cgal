// file: demo/Min_circle_2/Min_circle_2_simple_demo.C

#include <CGAL/Cartesian.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Cartesian<double>              K;
typedef CGAL::Min_circle_2_traits_2<K>       Traits;
typedef CGAL::Min_circle_2<Traits>           Min_circle;

int main ()
{
    std::istream_iterator<K::Point_2>  in_start(std::cin);
    std::istream_iterator<K::Point_2>  in_end;

    Min_circle  mc (in_start, in_end);   // smallest enclosing circle

    CGAL::Window_stream  W;
    W.init(-256.0, 255.0, -256.0);
    W.display();    
    W << CGAL::BLACK << mc;              // output circle + points
    W << CGAL::GREEN << mc.circle();     // output circle

    W.read_mouse();                      // wait for mouse click in window
    return 0;
}

