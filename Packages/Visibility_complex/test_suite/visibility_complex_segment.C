#include <CGAL/basic.h>

#include <fstream>
#include <list>

#include <CGAL/Simple_cartesian.h>

#include <CEP/Visibility_complex/Visibility_complex_segment_traits.h>
#include <CEP/Visibility_complex/Visibility_complex_2.h>

typedef int                                         FT;
typedef CGAL::Simple_cartesian <FT>                 Rep;
typedef CGAL::Visibility_complex_segment_traits<Rep> Gt;
typedef CGAL::Visibility_complex_2<Gt>              Visibility_complex;
typedef Visibility_complex::Vertex                  Vertex;
typedef Gt::Disk                                    Segment;

int main()
{
    std::list<Segment> O;
    std::ifstream ifs("segments");
    std::istream_iterator<Segment> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,std::back_inserter(O));

    // Computing visibility graph
    Visibility_complex V(O.begin(),O.end());
    if (V.size() != 872) exit(1);
    exit(0);
}
