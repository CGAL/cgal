#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <list>

#include <CGAL/Simple_cartesian.h>

#include <CEP/Visibility_complex/Visibility_complex_polygon_traits.h>
#include <CEP/Visibility_complex/Visibility_complex_2.h>

typedef int                                         FT;
typedef CGAL::Simple_cartesian <FT>                  Rep;
typedef CGAL::Visibility_complex_polygon_traits<Rep> Gt;
typedef CGAL::Visibility_complex_antichain<Gt>      Antichain;
typedef CGAL::Visibility_complex_2<Gt>              Visibility_complex;
typedef Visibility_complex::Vertex                  Vertex;
typedef Gt::Disk                                    Disk;
typedef Gt::Segment_2                               Segment;

int main()
{
    std::list<Disk> O;
    std::ifstream ifs("polygon");
    std::istream_iterator<Disk> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,std::back_inserter(O));

    // Computing visibility graph
    Visibility_complex V(O.begin(),O.end());
    if (V.size() != 1596) exit(1);
    exit(0);
}
