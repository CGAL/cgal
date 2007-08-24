#include <iostream>
#include <fstream>
#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Visibility_complex_2.h>
#include <CGAL/Visibility_complex_2/Segment_traits.h>
#include <CGAL/Shortest_path_2.h>

typedef CGAL::Simple_cartesian<int>                  K;
typedef CGAL::Visibility_complex_2_segment_traits<K>  Gt;

int main()
{
  std::ifstream di("data/segment.d");
  std::istream_iterator<Gt::Disk> disk_it(di),disk_end;
  std::list<Gt::Disk> disks;
  std::copy(disk_it,disk_end,std::back_inserter(disks));

  CGAL::Compute_free_bitangents_2<Gt>()(disks.begin(),disks.end(),
    std::ostream_iterator<Gt::Bitangent_2>(std::cout,"\n"));
  return 0;
}
