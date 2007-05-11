#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesher/Surface_mesher_edges_level.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/iterator.h>

#include <iostream>  // std::cerr
#include <ext/algorithm> // std::copy_n
#include <cstdlib>   // macros EXIT_*

#include <boost/format.hpp>

// default triangulation
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// generator
typedef CGAL::Random_points_in_cube_3<Tr::Point> Point_generator;

int main(int , char**)
{
  Tr tr;
  Point_generator pointgen(1);

  __gnu_cxx::copy_n(pointgen, 100, CGAL::inserter(tr));

  int result = EXIT_SUCCESS;

  Tr::size_type counter = 0;
  for(Tr::Edge_iterator eit = tr.edges_begin();
      eit != tr.edges_end();
      ++eit)
  {
    ++counter;
    Tr::Edge canon = CGAL::Surface_mesher::canonical_edge(tr, *eit);
    if(canon != *eit )
    {
      std::cerr <<
	::boost::format("edge (%1%, %2%, %3%)\n"
			"  canonical edge (%4%, %5%, %6%)\n")
	% &*eit->first % eit->second % eit->third
	% &*canon.first % canon.second % canon.third;
      result = EXIT_FAILURE;
    }
  }
  std::cerr << counter << " edges tested.\n";
  return result;
}
