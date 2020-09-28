#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Random.h>

#define CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE

// Define the kernel.
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point_3;

// Define the structure.
typedef CGAL::Delaunay_triangulation_3<Kernel>                               DT;
typedef DT::Segment_simplex_iterator                   Segment_simplex_iterator;
typedef CGAL::Triangulation_simplex_3<DT::Triangulation_data_structure> Simplex;

int main(int argc, char* argv[])
{
  const char* fname = (argc>1) ? argv[1] : "data/blobby.xyz";

  std::vector<Point_3> points;
  std::ifstream stream(fname);
  if (!stream ||
      !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // Construct the Delaunay triangulation.
  DT dt( points.begin(), points.end() );
  assert( dt.is_valid() );

  CGAL::Random rng;
  std::cout << "Random seed is " << CGAL::get_default_random().get_seed() << std::endl;

  unsigned int nb_cells = 0, nb_facets = 0, nb_edges = 0, nb_vertex = 0;

  ////////////////////////////////////////////////////////////
  // Construct a traverser and use begin/end iterators
  ////////////////////////////////////////////////////////////
  Point_3 p1(rng.get_double(-0.48, 0.31),
             rng.get_double(-0.22, 0.22),
             rng.get_double(-0.19, 0.19));
  Point_3 p2(rng.get_double(-0.48, 0.31),
             rng.get_double(-0.22, 0.22),
             rng.get_double(-0.19, 0.19));

  Segment_simplex_iterator st = dt.segment_traverser_simplices_begin(p1, p2);
  Segment_simplex_iterator stend = dt.segment_traverser_simplices_end();

  for (; st != stend; ++st)
  {
    if (st->dimension() == 3)          ++nb_cells;
    else if (st->dimension() == 2)     ++nb_facets;
    else if (st->dimension() == 1)     ++nb_edges;
    else if (st->dimension() == 0)     ++nb_vertex;
  }

#ifdef CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE
  std::cout << "While traversing from " << p1 << " to " << p2 << std::endl;
  std::cout << "\tcells    : " << nb_cells << std::endl;
  std::cout << "\tfacets   : " << nb_facets << std::endl;
  std::cout << "\tedges    : " << nb_edges << std::endl;
  std::cout << "\tvertices : " << nb_vertex << std::endl;
  std::cout << std::endl << std::endl;
#endif

  ////////////////////////////////////////////////////////////
  // Construct a traverser and use range-iterator
  ////////////////////////////////////////////////////////////

  nb_cells = 0, nb_facets = 0, nb_edges = 0, nb_vertex = 0;
  for (const Simplex& s : dt.segment_traverser_simplices(p1, p2))
  {
    if (s.dimension() == 3)          ++nb_cells;
    else if (s.dimension() == 2)     ++nb_facets;
    else if (s.dimension() == 1)     ++nb_edges;
    else if (s.dimension() == 0)     ++nb_vertex;
  }

#ifdef CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE
  std::cout << "While traversing from " << p1 << " to " << p2 << std::endl;
  std::cout << "\tcells    : " << nb_cells << std::endl;
  std::cout << "\tfacets   : " << nb_facets << std::endl;
  std::cout << "\tedges    : " << nb_edges << std::endl;
  std::cout << "\tvertices : " << nb_vertex << std::endl;
  std::cout << std::endl << std::endl;
#endif

  return 0;
}
