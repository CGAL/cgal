#include <list>
#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Visibility_complex_2/Segment_traits.h>
#include <CGAL/Shortest_path_2.h>

typedef CGAL::Simple_cartesian<int>  K;
typedef CGAL::Visibility_complex_2_segment_traits<K>  Gt;
typedef Gt::Disk                                     Segment;
typedef CGAL::Point_2<K> Point_2;
typedef CGAL::Shortest_path_2<Gt> Shp;
typedef Shp::Visibility_complex_2 VC;

int main()
{
    std::list<Segment> D;
    Gt::Make_convex_from_point mcfp;
    D.push_back(mcfp(Point_2(3,6)));
    D.push_back(mcfp(Point_2(6,-20)));
    D.push_back(Segment(Point_2(0,0),Point_2(10,10)));
    D.push_back(Segment(Point_2(5,2),Point_2(11,5)));
    D.push_back(Segment(Point_2(-4,-10),Point_2(6,-3)));
    VC vc(D.begin(),D.end());
    Shp shp(vc);
    shp.compute_shortest_paths(*vc.disks_begin());
    shp.get_path_bitangents(vc.disks_begin()[1],
      std::ostream_iterator<VC::Bitangent_2>(std::cout,"\n"));
    return 0;
}
