#include <CGAL/basic.h>

#include <fstream>
#include <list>

#include <CGAL/Simple_cartesian.h>

#include <CEP/Visibility_complex/Shortest_path_segment_traits.h>
#include <CEP/Visibility_complex/shortest_path_2.h>

// Number type used for input
typedef int                                         NT;  
// Number type used when computing distances
typedef double                                      SNT; 

typedef CGAL::Simple_cartesian<NT>                  R;

typedef CGAL::Shortest_path_segment_traits<R,SNT>   Gt;
typedef Gt::Disk                                    Segment;
typedef Gt::Point_2                                 Point;

int main()
{
    // Reading segments from file
    std::list<Segment> O;
    std::ifstream ifs("segments");
    std::istream_iterator<Segment> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,std::back_inserter(O));

    // We want to ccompute the shortest path between p and q
    Point p(0,0);
    Point q(1000,800);

    // Output the bitangents of the shortest path to cout
    CGAL::shortest_path_2(O.begin(),O.end(), 
			  p,q, 
			  std::ostream_iterator<Segment>(cout , "\n"),
			  Gt());
    return 0;
}
